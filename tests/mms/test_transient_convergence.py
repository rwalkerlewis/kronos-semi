"""
MMS temporal convergence test for the transient drift-diffusion solver.

Reference strategy
------------------
The transient solver starts from equilibrium at V=0 and immediately
applies the bias V=0.1V at t=0.  As t -> infinity the solution must
converge to the steady-state drift-diffusion solution at V=0.1V.  We
obtain that steady-state reference independently via run_bias_sweep on
the same mesh.  This reference is independent of dt and BDF order.

For large enough t_end (t_end >> tau_p), the transient solution is deep
in steady state and the dominant error is accumulated time-stepping
error, which scales as dt^p where p is the BDF order.

Convergence study
-----------------
- Run transient to t_end = 50 * tau_p at four dt refinement levels.
- Compute err_n = ||n_transient(t_end; dt_k) - n_bias_sweep||_inf
- Assert observed convergence rates:
    * BDF1: rate >= 0.95  (expected 1.0)
    * BDF2: rate >= 1.9   (expected 2.0)

Design parameters (Pe < 1)
--------------------------
- 1D pn junction, 20 um, N_A = N_D = 1e15 cm^-3 (low doping -> Pe < 1)
- N_MESH = 400 elements
- V_bias = 0.1 V (well below 0.6 V to keep injection moderate)
- tau_n = tau_p = 1e-9 s (short lifetime -> fast relaxation to SS)
- t_end = 50 * tau_p = 5.0e-8 s
- DT_BASE = t_end / 200 = 2.5e-10 s; four levels 200/400/800/1600 steps

This test runs with the real transient runner on a minimal 1D problem.
It is placed in tests/mms/ and is collected by pytest.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

# Minority-carrier lifetime; sets the relaxation timescale
_TAU_P = 1.0e-9   # s

# Run deep into steady state so accumulated time-stepping error dominates
T_END = 50.0 * _TAU_P   # 5e-8 s

# Four dt refinement levels: DT_BASE / 2^k for k in [0,1,2,3]
# Step counts: 200, 400, 800, 1600
DT_BASE = T_END / 200.0   # 2.5e-10 s
N_LEVELS = 4
DT_LIST = [DT_BASE / (2 ** k) for k in range(N_LEVELS)]

# Mesh: fine enough that spatial error << temporal error at the coarsest dt
N_MESH = 400
L_DEVICE = 2.0e-5   # 20 um
V_BIAS = 0.1        # V, applied forward bias

# Rate floor targets
RATE_BDF1_FLOOR = 0.95
RATE_BDF2_FLOOR = 1.90


def _base_cfg() -> dict:
    """
    Build a minimal 1D pn junction config (Pe < 1: low doping, fine mesh).
    The solver block is omitted; callers fill it in.
    """
    return {
        "schema_version": "1.1.0",
        "name": "mms_transient",
        "description": "Minimal 1D pn junction for transient MMS convergence.",
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
                    "N_A_left": 1.0e15,
                    "N_D_right": 1.0e15,
                    "N_A_right": 0.0,
                },
            }
        ],
        "contacts": [
            {"name": "anode",   "facet": "anode",   "type": "ohmic", "voltage": V_BIAS},
            {"name": "cathode", "facet": "cathode",  "type": "ohmic", "voltage": 0.0},
        ],
        "physics": {
            "temperature": 300.0,
            "statistics": "boltzmann",
            "mobility": {"mu_n": 1400.0, "mu_p": 450.0},
            "recombination": {
                "srh": True,
                "tau_n": _TAU_P,
                "tau_p": _TAU_P,
                "E_t": 0.0,
            },
        },
        "output": {
            "directory": "/tmp/mms_transient_results",
            "fields": [],
        },
    }


def _minimal_transient_cfg(dt: float, t_end: float, order: int) -> dict:
    """
    Build a 1D transient config for the pn junction with the given
    timestep and BDF order.
    """
    max_steps = int(math.ceil(t_end / dt)) + 5
    cfg = _base_cfg()
    cfg["solver"] = {
        "type": "transient",
        "t_end": float(t_end),
        "dt": float(dt),
        "order": int(order),
        "max_steps": int(max_steps),
        "output_every": int(max_steps + 1),  # only snapshot at t_end
        "snes": {
            "rtol": 1.0e-10,
            "atol": 1.0e-7,
            "stol": 1.0e-14,
            "max_it": 50,
        },
    }
    return cfg


def _compute_ss_reference() -> tuple[np.ndarray, np.ndarray]:
    """
    Compute the steady-state reference solution at V=V_BIAS via
    run_bias_sweep on the same mesh.  This reference is independent of
    dt and BDF order.

    Returns (n_ref, p_ref) in physical units (cm^-3 equivalent in SI).
    """
    from semi.runners.bias_sweep import run_bias_sweep

    ss_cfg = _base_cfg()
    ss_cfg["solver"] = {
        "type": "bias_sweep",
        "snes": {"rtol": 1.0e-14, "atol": 1.0e-14, "stol": 1.0e-14, "max_it": 60},
        "continuation": {
            "min_step": 1.0e-3,
            "max_halvings": 6,
            "max_step": 0.1,
            "easy_iter_threshold": 4,
            "grow_factor": 2.0,
        },
    }
    ss_cfg["contacts"][0]["voltage_sweep"] = {
        "start": 0.0,
        "stop": V_BIAS,
        "step": 0.05,
    }
    ss_result = run_bias_sweep(ss_cfg)
    return ss_result.n_phys, ss_result.p_phys


def _run_and_get_final_state(dt: float, order: int) -> tuple[np.ndarray, np.ndarray]:
    """
    Run the transient solver to T_END and return the final (n, p) arrays
    in physical units.
    """
    from semi.runners.transient import run_transient

    cfg = _minimal_transient_cfg(dt=dt, t_end=T_END, order=order)
    result = run_transient(cfg)

    if result.fields.get("n") and result.fields["n"]:
        n_final = result.fields["n"][-1]
        p_final = result.fields["p"][-1]
    else:
        raise RuntimeError("No field snapshots were written by the transient runner.")
    return n_final, p_final


def _log_rate(dts, errors):
    """Log-log slope of last two (dt, err) pairs."""
    dt_prev, dt_cur = dts[-2], dts[-1]
    e_prev, e_cur = errors[-2], errors[-1]
    if e_cur <= 0.0 or e_prev <= 0.0:
        return float("nan")
    return math.log(e_prev / e_cur) / math.log(dt_prev / dt_cur)


def _run_convergence_study(order: int) -> tuple[list[float], list[float]]:
    """
    Run temporal convergence study for the given BDF order.
    Returns (errors_n_inf, errors_p_inf) lists (one per dt level).

    Reference: bias_sweep steady-state at V=V_BIAS (independent of dt).
    Error metric: ||n_transient(t_end; dt_k) - n_bias_sweep||_inf
    """
    n_ref, p_ref = _compute_ss_reference()

    errors_n = []
    errors_p = []
    for dt in DT_LIST:
        n_k, p_k = _run_and_get_final_state(dt=dt, order=order)
        err_n = float(np.max(np.abs(n_k - n_ref)))
        err_p = float(np.max(np.abs(p_k - p_ref)))
        errors_n.append(err_n)
        errors_p.append(err_p)

    return errors_n, errors_p


@pytest.mark.slow
def test_transient_convergence_bdf1():
    """
    BDF1 (backward Euler) temporal convergence test.
    Assert observed rate >= 0.95 (theoretical: 1.0).
    """
    errors_n, errors_p = _run_convergence_study(order=1)
    rate_n = _log_rate(DT_LIST, errors_n)
    rate_p = _log_rate(DT_LIST, errors_p)

    print("\nBDF1 temporal convergence:")
    for i, dt in enumerate(DT_LIST):
        print(f"  dt={dt:.2e}  err_n={errors_n[i]:.3e}  err_p={errors_p[i]:.3e}")
    print(f"  Observed rates: n={rate_n:.3f}, p={rate_p:.3f}")

    assert rate_n >= RATE_BDF1_FLOOR, (
        f"BDF1 n rate {rate_n:.3f} < {RATE_BDF1_FLOOR}; "
        f"errors={errors_n}"
    )
    assert rate_p >= RATE_BDF1_FLOOR, (
        f"BDF1 p rate {rate_p:.3f} < {RATE_BDF1_FLOOR}; "
        f"errors={errors_p}"
    )


@pytest.mark.slow
def test_transient_convergence_bdf2():
    """
    BDF2 temporal convergence test.
    Assert observed rate >= 1.9 (theoretical: 2.0).
    """
    errors_n, errors_p = _run_convergence_study(order=2)
    rate_n = _log_rate(DT_LIST, errors_n)
    rate_p = _log_rate(DT_LIST, errors_p)

    print("\nBDF2 temporal convergence:")
    for i, dt in enumerate(DT_LIST):
        print(f"  dt={dt:.2e}  err_n={errors_n[i]:.3e}  err_p={errors_p[i]:.3e}")
    print(f"  Observed rates: n={rate_n:.3f}, p={rate_p:.3f}")

    assert rate_n >= RATE_BDF2_FLOOR, (
        f"BDF2 n rate {rate_n:.3f} < {RATE_BDF2_FLOOR}; "
        f"errors={errors_n}"
    )
    assert rate_p >= RATE_BDF2_FLOOR, (
        f"BDF2 p rate {rate_p:.3f} < {RATE_BDF2_FLOOR}; "
        f"errors={errors_p}"
    )
