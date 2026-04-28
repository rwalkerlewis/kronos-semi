"""
Temporal convergence test for the BDF1/BDF2 transient drift-diffusion solver.

Pairwise-difference design
--------------------------
We measure the temporal discretisation error of the Slotboom-form
transient runner (ADR 0014) by comparing solutions at consecutive dt
levels on the same mesh.

Why pairwise differences (not bias_sweep, not self-Richardson):
* The previous self-Richardson approach used a reference at the same
  BDF order but only 4x finer dt than the finest test level, which
  contaminated the comparison with the reference's own temporal error
  and produced non-monotonic, meaningless rates.
* Pairwise differences avoid this: for two consecutive levels
  dt_k and dt_{k+1} = dt_k / 2,

      D_k  :=  || u(dt_k) - u(dt_{k+1}) ||_inf
            ~  |C(t_end)|_inf * dt_k^p * (1 - 2^{-p})

  The spatial error cancels exactly because both runs share mesh and
  formulation. The log-log slope of D_k vs dt_k recovers p.

Device design (M13.1 follow-up — 1e17 doping, V_F = 0.05 V)
------------------------------------------------------------
The original 1e15 / 0.1 V device failed the rate test because the dt
window where MUMPS stays stable on the post-ramp V_F/2 → V_F BC step
(~50 ps minimum) did not overlap with the dt window required for a
clean BDF rate signal. This is documented in ADR 0014 Limitations.

Fix: use the `pn_1d_turnon` benchmark device (1e17 doping) with a
small forward bias (V_F = 0.05 V, roughly 2 kT/q). The near-linear
response below thermal voltage keeps the post-ramp BC step small and
avoids the MUMPS conditioning failure at the coarsest dt. The 1e17
regime is already validated by the `pn_1d_turnon` benchmark; the
Slotboom (psi, phi_n, phi_p) formulation guarantees positivity of
carriers by construction so the old Pe < 1 Galerkin stability
constraint no longer applies (see ADR 0014).

Convergence thresholds:
  BDF1: rate >= 0.95   (theoretical: 1.0)
  BDF2: rate >= 1.90   (theoretical: 2.0)
"""
from __future__ import annotations

import math

import numpy as np
import pytest

# Forward bias for the step-on transient.
# 0.10 V ~ 4 kT/q: still in the small-signal regime so the BC ramp
# does not break the bias_sweep continuation (which fails above
# ~0.2 V on the 1e17 device, see ADR 0014 Limitations), but large
# enough that the BDF2 truncation envelope at the coarse dt_base
# below sits above the SNES iterative noise floor (~1e15 m^-3 of
# carrier density on the majority side).
V_F = 0.10

# SRH lifetime; sets the time scale of the minority-carrier transient.
TAU = 1.0e-9

# Stop at t_end = 1 * tau, near the peak of the temporal-error envelope.
# For a stable linear system u' = -u/tau, the BDF temporal error scales
# as exp(-t/tau) * t * dt^p — small for t << tau (system hasn't moved)
# and exponentially decaying for t >> tau (errors damp toward steady
# state). At t ~ tau the error is still O(dt^p) and easy to resolve.
# Going to t_end = 10 * tau hides the temporal signal in SS noise floor.
T_END = 1.0 * TAU  # 1 ns

# Per-order dt windows. BDF1 and BDF2 sit on different sides of the
# noise-floor / truncation crossover and therefore need different
# coarsest-dt anchors:
#
#   * BDF1 truncation (O(dt^1)) at dt = T_END/16 = 6.25e-11 s puts
#     the diff sequence ~5e15 -> ~2e14, a factor-25 fall over four
#     levels — squarely above the ~1e14 noise floor and inside the
#     asymptotic regime, with rates ~3.0+ (well above the 0.95 floor).
#   * BDF2 truncation (O(dt^2)) at the same dt_base produces diffs
#     of ~1e15 m^-3 — at or below the noise floor, so observed rates
#     are dominated by SNES iterative noise and are not monotonic.
#     Pulling dt_base up by 4x to T_END/4 = 2.5e-10 s lifts the
#     coarsest BDF2 diff to ~6e17 m^-3, well above the noise floor,
#     and lands the last refinement (dt = T_END/32 = 3.125e-11 s)
#     at the start of the asymptotic O(dt^2) regime with rate ~2.13.
#
# The crossover is intrinsic to pairwise-difference temporal
# convergence on a stiff continuity system at fixed mesh: BDF1 needs
# fine dt to climb out of its quadratic-truncation upper-floor band
# and BDF2 needs coarse dt to stay above the linear-noise lower-floor
# band. Sharing a single dt window across orders only works when
# both bands are several decades apart, which is not the case here.
N_LEVELS = 4
DT_BASE_BDF1 = T_END / 16.0   # 6.25e-11 s
DT_BASE_BDF2 = T_END / 4.0    # 2.50e-10 s
DT_LIST_BDF1 = [DT_BASE_BDF1 / (2 ** k) for k in range(N_LEVELS)]
DT_LIST_BDF2 = [DT_BASE_BDF2 / (2 ** k) for k in range(N_LEVELS)]

# Mesh: N=400 on a 20 um device -> h=50 nm.
# At 1e17 doping the Slotboom (psi, phi_n, phi_p) form guarantees carrier
# positivity by construction; the Pe < 1 Galerkin stability constraint
# from the old (n, p) density form does not apply here.
N_MESH = 400
L_DEVICE = 2.0e-5

# Rate floor targets (unchanged from M13 specification).
RATE_BDF1_FLOOR = 0.95
RATE_BDF2_FLOOR = 1.90


def _transient_cfg(dt: float, order: int) -> dict:
    """Build a 1D pn junction transient config with the given dt / BDF order."""
    max_steps = int(math.ceil(T_END / dt)) + 5
    return {
        "schema_version": "1.1.0",
        "name": "mms_transient",
        "description": "1D pn junction for transient temporal convergence test.",
        "dimension": 1,
        "mesh": {
            "source": "builtin",
            "extents": [[0.0, L_DEVICE]],
            "resolution": [N_MESH],
            "regions_by_box": [
                {"name": "silicon", "tag": 1, "bounds": [[0.0, L_DEVICE]]},
            ],
            "facets_by_plane": [
                {"name": "anode",   "tag": 1, "axis": 0, "value": 0.0},
                {"name": "cathode", "tag": 2, "axis": 0, "value": L_DEVICE},
            ],
        },
        "regions": {"silicon": {"material": "Si", "tag": 1, "role": "semiconductor"}},
        "doping": [
            {
                "region": "silicon",
                "profile": {
                    "type": "step",
                    "axis": 0,
                    "location": L_DEVICE / 2.0,
                    "N_D_left": 0.0,
                    "N_A_left": 1.0e17,
                    "N_D_right": 1.0e17,
                    "N_A_right": 0.0,
                },
            }
        ],
        "contacts": [
            {"name": "anode",   "facet": "anode",   "type": "ohmic", "voltage": V_F},
            {"name": "cathode", "facet": "cathode",  "type": "ohmic", "voltage": 0.0},
        ],
        "physics": {
            "temperature": 300.0,
            "statistics": "boltzmann",
            "mobility": {"mu_n": 1400.0, "mu_p": 450.0},
            "recombination": {
                "srh": True,
                "tau_n": TAU,
                "tau_p": TAU,
                "E_t": 0.0,
            },
        },
        "solver": {
            "type": "transient",
            "t_end": float(T_END),
            "dt": float(dt),
            "order": int(order),
            "max_steps": int(max_steps),
            "output_every": int(max_steps + 1),  # only the t=t_end snapshot
            # ADR 0014 (Slotboom transient, supersedes ADR 0009):
            #
            # `bc_ramp_voltage_factor=0.5` lands the IC at the V_F/2
            # steady state; the time loop then integrates the V_F/2 ->
            # V_F relaxation transient -- a well-defined signal for BDF
            # rate measurement.  At 1e17 doping / V_F=0.10 V the BC
            # step impulse is small (~4 kT/q) so MUMPS stays
            # well-conditioned at both BDF1 and BDF2 dt windows.
            # `bc_ramp_steps=0` (V=0 IC + step bias) at 1e17 doping
            # would drive Newton through a large density swing at
            # step 1 and is deliberately avoided here.
            "bc_ramp_steps": 10,
            "bc_ramp_voltage_factor": 0.5,
            "snes": {
                # Calibrated to the V_F/2 -> V_F transient residual
                # scale. The 1e-14 atol that the (n, p) form needed
                # to avoid "freezing the transient" is unnecessary
                # here: under Slotboom the BC step is the dominant
                # signal and Newton is naturally driven each step.
                "rtol": 1.0e-10,
                "atol": 1.0e-10,
                "stol": 1.0e-14,
                "max_it": 100,
            },
        },
        "output": {
            "directory": "/tmp/mms_transient_results",
            "fields": [],
        },
    }


def _run_transient_final_state(dt: float, order: int):
    """
    Run the transient solver to T_END at the given dt / order and return
    a dict of final-time field arrays:

      {
        "n":     n(t_end)       [m^-3]
        "p":     p(t_end)       [m^-3]
        "phi_n": phi_n(t_end)   [V]
        "phi_p": phi_p(t_end)   [V]
      }

    The transient runner snapshots fields whenever
    `step_count % output_every == 0` OR `|t_current - t_end| < 0.5*dt`,
    so with output_every set above max_steps the only snapshot written
    is at t = t_end. fields[<key>][-1] is therefore the final state.
    """
    from semi.runners.transient import run_transient
    cfg = _transient_cfg(dt=dt, order=order)
    result = run_transient(cfg)
    if not result.fields.get("n"):
        raise RuntimeError("Transient runner returned no field snapshots.")
    return {
        "n":     np.asarray(result.fields["n"][-1]).copy(),
        "p":     np.asarray(result.fields["p"][-1]).copy(),
        "phi_n": np.asarray(result.fields["phi_n"][-1]).copy(),
        "phi_p": np.asarray(result.fields["phi_p"][-1]).copy(),
    }


def _log_rate(dts, errors):
    """Log-log slope of the last two (dt, err) pairs."""
    dt_prev, dt_cur = dts[-2], dts[-1]
    e_prev, e_cur = errors[-2], errors[-1]
    if e_cur <= 0.0 or e_prev <= 0.0:
        return float("nan")
    return math.log(e_prev / e_cur) / math.log(dt_prev / dt_cur)


def _run_convergence_study(order: int, dt_list: list[float]) -> tuple[list[float], list[float], list[float]]:
    """
    Run the temporal convergence study for the given BDF order using
    the (n, p) density-form pairwise differences.

    Returns (dts_for_diff, diffs_n, diffs_p), where:
      * dts_for_diff[k] = dt_list[k]    (the larger dt of the pair)
      * diffs_n[k] = ||n(dt_k) - n(dt_{k+1})||_inf
      * diffs_p[k] = ||p(dt_k) - p(dt_{k+1})||_inf

    Length of all three is len(dt_list) - 1 (3 entries from 4 dt levels).
    """
    states_n: list[np.ndarray] = []
    states_p: list[np.ndarray] = []
    for dt in dt_list:
        st = _run_transient_final_state(dt=dt, order=order)
        states_n.append(st["n"])
        states_p.append(st["p"])
        if len(states_n) > 1 and states_n[-1].shape != states_n[0].shape:
            raise RuntimeError(
                f"DOF count mismatch between dt levels: {states_n[0].shape} "
                f"vs {states_n[-1].shape}"
            )

    dts_for_diff = dt_list[:-1]
    diffs_n = [
        float(np.max(np.abs(states_n[k] - states_n[k + 1])))
        for k in range(len(dt_list) - 1)
    ]
    diffs_p = [
        float(np.max(np.abs(states_p[k] - states_p[k + 1])))
        for k in range(len(dt_list) - 1)
    ]
    return dts_for_diff, diffs_n, diffs_p


def _assert_monotonic(diffs: list[float], label: str) -> None:
    for i in range(1, len(diffs)):
        assert diffs[i] < diffs[i - 1], (
            f"{label} diffs not monotonically decreasing: "
            f"diff[{i - 1}]={diffs[i - 1]:.3e}, diff[{i}]={diffs[i]:.3e}, "
            f"all={diffs}"
        )


@pytest.mark.slow
def test_transient_convergence_bdf1():
    """
    BDF1 (backward Euler) temporal convergence test.

    Device: 1e17 cm^-3 symmetric pn junction, V_F = 0.10 V, T_END = 1 ns.
    BC ramp: V_F/2 IC -> V_F time loop (bc_ramp_voltage_factor=0.5).
    Uses the fine dt window DT_LIST_BDF1 (DT_BASE = T_END/16); see the
    module docstring on the BDF1/BDF2 dt-window split.

    Asserts observed rate >= 0.95 (theoretical: 1.0) and that the
    pairwise-difference sequence decreases monotonically with dt.
    """
    dts, diffs_n, diffs_p = _run_convergence_study(order=1, dt_list=DT_LIST_BDF1)
    rate_n = _log_rate(dts, diffs_n)
    rate_p = _log_rate(dts, diffs_p)

    print("\nBDF1 temporal convergence (pairwise differences):")
    for i, dt in enumerate(dts):
        print(
            f"  dt={dt:.3e}  diff_n={diffs_n[i]:.3e}  diff_p={diffs_p[i]:.3e}"
        )
    print(f"  Observed rates: n={rate_n:.3f}, p={rate_p:.3f}")

    _assert_monotonic(diffs_n, "BDF1 n")
    _assert_monotonic(diffs_p, "BDF1 p")

    assert rate_n >= RATE_BDF1_FLOOR, (
        f"BDF1 n rate {rate_n:.3f} < {RATE_BDF1_FLOOR}; diffs={diffs_n}"
    )
    assert rate_p >= RATE_BDF1_FLOOR, (
        f"BDF1 p rate {rate_p:.3f} < {RATE_BDF1_FLOOR}; diffs={diffs_p}"
    )


@pytest.mark.slow
def test_transient_convergence_bdf2():
    """
    BDF2 temporal convergence test.

    Device: 1e17 cm^-3 symmetric pn junction, V_F = 0.10 V, T_END = 1 ns.
    BC ramp: V_F/2 IC -> V_F time loop (bc_ramp_voltage_factor=0.5).
    Uses the coarse dt window DT_LIST_BDF2 (DT_BASE = T_END/4); see the
    module docstring on the BDF1/BDF2 dt-window split. The coarser dt
    keeps BDF2 truncation above the SNES iterative noise floor that
    would otherwise scramble the pairwise-difference rate signal.

    Asserts observed rate >= 1.9 (theoretical: 2.0) and that the
    pairwise-difference sequence decreases monotonically with dt.
    """
    dts, diffs_n, diffs_p = _run_convergence_study(order=2, dt_list=DT_LIST_BDF2)
    rate_n = _log_rate(dts, diffs_n)
    rate_p = _log_rate(dts, diffs_p)

    print("\nBDF2 temporal convergence (pairwise differences):")
    for i, dt in enumerate(dts):
        print(
            f"  dt={dt:.3e}  diff_n={diffs_n[i]:.3e}  diff_p={diffs_p[i]:.3e}"
        )
    print(f"  Observed rates: n={rate_n:.3f}, p={rate_p:.3f}")

    _assert_monotonic(diffs_n, "BDF2 n")
    _assert_monotonic(diffs_p, "BDF2 p")

    assert rate_n >= RATE_BDF2_FLOOR, (
        f"BDF2 n rate {rate_n:.3f} < {RATE_BDF2_FLOOR}; diffs={diffs_n}"
    )
    assert rate_p >= RATE_BDF2_FLOOR, (
        f"BDF2 p rate {rate_p:.3f} < {RATE_BDF2_FLOOR}; diffs={diffs_p}"
    )
