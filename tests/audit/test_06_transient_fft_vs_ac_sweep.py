"""
Audit case 6: transient with sinusoidal bias (FFT) vs ac_sweep at the
same frequency.

The AC small-signal admittance Y(omega) is by construction the
fundamental frequency component of the current response to a small
sinusoidal voltage perturbation. This test runs the transient runner
with V(t) = V_DC + dV * sin(omega t) materialized as a `voltage_t`
table (M16.7) on the anode contact, takes the FFT of the terminal
current at omega, and compares to ac_sweep at that frequency.

Windowing
---------
The transient runs for an integer number of cycles, so a rectangular
window already gives clean spectral peaks at the fundamental. We apply
a Hann window before the FFT anyway: the off-fundamental leakage from
solver noise is small but non-zero, and Hann's flat top keeps the
fundamental amplitude estimate stable across slight phase mismatches
between V(t) and the sample grid. The complex transfer function
`Y_transient_fft = FFT(I)[bin_F] / FFT(V)[bin_F]` is window-invariant
because the same window multiplies both numerator and denominator.

Tuning
------
F_HZ = 1 MHz, dV = 1 mV (small-signal), V_DC = 0.4 V on the
`rc_ac_sweep` benchmark (1D pn diode, 800 cells, 20 um). N_periods = 4
cycles at 200 samples per cycle gives 800 BDF1 steps with dt = 5 ns.
The audit gate is 5%; over-sampling further is wasted CI time.

The acceptance assertion is `|Y_transient_fft - Y_ac| / |Y_ac| < 0.05`
per IMPROVEMENT_GUIDE § M16.7 acceptance test 1.
"""
from __future__ import annotations

import copy

import numpy as np
import pytest

from ._helpers import (
    load_benchmark,
    require_dolfinx,
    write_csv,
    write_markdown,
)

CASE = "06_transient_fft_vs_ac_sweep"
V_DC = 0.4
DV = 1.0e-3
F_HZ = 1.0e6

# CI-budget tuning. 4 cycles at 200 samples/cycle = 800 BDF1 steps.
N_PERIODS = 4
SAMPLES_PER_PERIOD = 200
N_SAMPLES = N_PERIODS * SAMPLES_PER_PERIOD
DT = 1.0 / (F_HZ * SAMPLES_PER_PERIOD)
T_END = N_PERIODS / F_HZ


@pytest.mark.audit
def test_transient_fft_vs_ac_sweep():
    require_dolfinx()

    from semi.runners.ac_sweep import run_ac_sweep
    from semi.runners.transient import run_transient

    cfg = load_benchmark("rc_ac_sweep/rc_ac_sweep.json")

    # AC reference value at f = F_HZ.
    cfg_ac = copy.deepcopy(cfg)
    cfg_ac["contacts"][0]["voltage"] = V_DC
    cfg_ac["solver"]["dc_bias"] = {"contact": "anode", "voltage": V_DC}
    cfg_ac["solver"]["ac"]["frequencies"] = {"type": "list", "values": [F_HZ]}
    ac_res = run_ac_sweep(cfg_ac)
    Y_ac = complex(ac_res.Y[0])

    # Build the transient configuration. Sample V(t) on a uniform grid
    # at SAMPLES_PER_PERIOD per cycle for N_PERIODS cycles. The grid
    # ends one dt before T_END so the last BDF step lands exactly on
    # `times[-1]` and the runner returns at that t.
    times = np.arange(N_SAMPLES) * DT
    values = V_DC + DV * np.sin(2.0 * np.pi * F_HZ * times)

    cfg_tr = copy.deepcopy(cfg)
    cfg_tr["contacts"][0]["voltage"] = V_DC  # sets ramp target & is the
                                              # static fallback
    cfg_tr["contacts"][0]["voltage_t"] = {
        "type": "table",
        "times": times.tolist(),
        "values": values.tolist(),
    }
    cfg_tr["solver"] = {
        "type": "transient",
        "t_end": float(T_END),
        "dt": float(DT),
        "order": 1,  # BDF1 is enough for small-signal current; matches
                     # the AC linearization point well.
        "max_steps": 4 * N_SAMPLES,
        "bc_ramp_steps": 20,
        "output_every": N_SAMPLES + 1,
    }

    tr_res = run_transient(cfg_tr)

    # Pull I(t) on the anode contact at every recorded step. The IV
    # rows are in chronological order; one row per (t, contact).
    iv_rows = [r for r in tr_res.iv if r["contact"] == "anode"]
    iv_rows.sort(key=lambda r: r["t"])
    t_arr = np.array([r["t"] for r in iv_rows])
    v_arr = np.array([r["V"] for r in iv_rows])
    j_arr = np.array([r["J"] for r in iv_rows])

    # Drop the t=0 entry so the sample count matches the table grid
    # (the runner records an initial entry before the first BDF step).
    if len(t_arr) == N_SAMPLES + 1:
        t_arr = t_arr[1:]
        v_arr = v_arr[1:]
        j_arr = j_arr[1:]

    n = len(t_arr)
    assert n == N_SAMPLES, (
        f"expected {N_SAMPLES} transient samples, got {n}"
    )

    # Demean and apply a Hann window before the FFT. The window is
    # applied to both V and J so it cancels in the admittance ratio.
    window = np.hanning(n)
    v_ac = (v_arr - V_DC) * window
    j_ac = (j_arr - np.mean(j_arr)) * window
    V_fft = np.fft.rfft(v_ac)
    J_fft = np.fft.rfft(j_ac)

    # The fundamental sits at bin = N_PERIODS (integer cycles in the
    # window). Pick the strongest bin near it as a robustness guard.
    expected_bin = N_PERIODS
    search = slice(max(0, expected_bin - 1), expected_bin + 2)
    bin_F = expected_bin + (
        int(np.argmax(np.abs(V_fft[search]))) - (expected_bin - search.start)
    )

    Y_transient_fft = J_fft[bin_F] / V_fft[bin_F]
    rel_err = abs(Y_transient_fft - Y_ac) / abs(Y_ac)

    write_csv(
        CASE,
        [
            "V_DC", "f_Hz", "Y_ac_re", "Y_ac_im",
            "Y_transient_fft_re", "Y_transient_fft_im",
            "rel_err", "status",
        ],
        [[
            V_DC, F_HZ, Y_ac.real, Y_ac.imag,
            Y_transient_fft.real, Y_transient_fft.imag,
            float(rel_err), "passed" if rel_err < 0.05 else "failed",
        ]],
    )
    # As in Cases 02 and 05: per-run floats live in the CSV; the
    # markdown is a stable summary so the CI sync check (`git diff
    # --exit-code docs/PHYSICS_AUDIT.md`) is reproducible. The MUMPS
    # LU pivot ordering on the indefinite forward-bias AC block
    # contributes the same ~1-2% Re(Y) noise floor as Case 05; the
    # transient FFT inherits it via Y_ac in the rel_err numerator.
    write_markdown(
        CASE,
        "Case 06 - transient (FFT) vs ac_sweep",
        f"At V_DC = {V_DC} V, f = 1 MHz, dV = 1 mV on the "
        f"`rc_ac_sweep` benchmark, the FFT of I(t) under "
        f"V(t) = V_DC + dV*sin(omega t) (sampled into a "
        f"`voltage_t` table; M16.7) recovers the AC small-signal "
        f"admittance Y(omega) within the 5% audit gate. Per-run "
        f"values (Y_ac, Y_transient_fft, rel_err) are written to "
        f"`/tmp/audit/{CASE}.csv`. The transient runner uses BDF1 "
        f"over 4 cycles at 200 samples per cycle (800 timesteps, "
        f"dt = 5 ns) with `bc_ramp_steps = 20` ramping to V_DC "
        f"before the time loop. A Hann window applied to both V(t) "
        f"and I(t) before the rfft cancels in the admittance ratio "
        f"Y(omega) = J_fft[bin] / V_fft[bin].\n\n"
        f"CSV: `/tmp/audit/{CASE}.csv`",
    )

    assert rel_err < 0.05, (
        f"|Y_transient_fft - Y_ac| / |Y_ac| = {rel_err:.3e} exceeds "
        f"5% gate; Y_ac = {Y_ac}, Y_transient_fft = {Y_transient_fft}"
    )
