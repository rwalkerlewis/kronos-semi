"""
MMS temporal convergence test for the transient drift-diffusion solver.

Reference strategy (BDF2 self-Richardson)
------------------------------------------
The transient solver is run at four dt refinement levels (BDF1 or BDF2)
and compared against a BDF2 reference computed at dt_ref = DT_LIST[-1]/8.
Because both reference and test use the same density-form discretisation
(run_transient), the spatial mismatch between the two solutions is
identically zero.  The only source of difference is temporal integration
error, which scales as dt^p where p is the BDF order.

Key design choice: long lifetime + large dt to isolate slow mode
-----------------------------------------------------------------
The pn junction has two timescales:
  * Fast: tau_transport ~ 9e-11 s (depletion-edge dynamics)
  * Slow: tau_p (minority carrier recombination)

The slow mode gives clean dt^p convergence; the fast mode does not, because
its spatial profile (exp(±psi/VT) at the depletion edge) is non-smooth and
BDF errors there do not follow asymptotic scaling at practical dt values.

By choosing tau_p = 1e-7 s and DT_LIST = [1e-8, ...]:
  * T_END = tau_p = 1e-7 s >> tau_transport * 1111  (fast mode fully damped)
  * DT_LIST[0] = 1e-8 s < tau_p  (slow mode still active => errors measurable)
  * All DT_LIST values >> tau_transport  (fast mode stiff; A-stable BDF damps it)

With fast-mode errors vanishing exponentially at T_END (exp(-1111) ~ 0), the
global temporal error is dominated entirely by the slow minority-carrier mode,
which scales as dt^p for all four refinement levels.

Why not bias_sweep reference?
------------------------------
A bias_sweep (Slotboom-form) reference was tried but fails for the pn
junction: the Slotboom and density-form discretisations converge to *different*
discrete steady states (mismatch ~1e19 m^-3 at the depletion edge, versus
temporal errors ~1e12 m^-3).  No choice of t_end makes the temporal signal
visible above this floor.

Convergence study
-----------------
- T_END = 1 * tau_p (transient is active; temporal errors accumulating).
- Four dt refinement levels: DT_BASE / 2^k for k in [0,1,2,3].
  Step counts: 10 / 20 / 40 / 80 (fast wall-clock time per run).
- Reference: BDF2 at dt_ref = DT_LIST[-1] / 8  ->  640 steps.
  Computed once and cached across both test functions.
- Error metric: ||u(t_end; dt_k) - u_ref(t_end; dt_ref)||_inf
- Assert observed convergence rates:
    * BDF1: rate >= 0.95  (expected ~1.0)
    * BDF2: rate >= 1.9   (expected ~2.0)

Design parameters (Pe < 1)
--------------------------
- 1D pn junction, 20 um, N_A = N_D = 1e15 cm^-3 (Pe < 1 in depletion)
- N_MESH = 400 elements
- V_bias = 0.1 V (well below built-in potential)
- tau_n = tau_p = 1e-7 s (long lifetime so tau_p >> tau_transport)
- T_END = 1 * tau_p = 1e-7 s
- DT_BASE = T_END / 10 = 1e-8 s

This test runs with the real transient runner on a minimal 1D problem.
It is placed in tests/mms/ and is collected by pytest.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

# Long minority-carrier lifetime: ensures tau_p >> tau_transport ~ 9e-11 s
# so the fast depletion-edge mode is fully decayed (exp(-1111) ~ 0) by T_END.
_TAU_P = 1.0e-7   # s

# Simulation window: 1 time constant keeps the transient active so
# temporal errors from the slow mode are large enough to measure cleanly.
T_END = 1.0 * _TAU_P   # 1e-7 s

# Four dt refinement levels: DT_BASE / 2^k for k in [0,1,2,3]
# All DT values >> tau_transport so fast depletion mode is stiff/damped.
# Step counts: 10, 20, 40, 80 (fast per run)
DT_BASE = T_END / 10.0   # 1e-8 s
N_LEVELS = 4
DT_LIST = [DT_BASE / (2 ** k) for k in range(N_LEVELS)]

# Reference: BDF2 at 8x the finest test dt.
# BDF2 ref error (O(dt_ref^2)) << BDF1 test error (O(dt)):
#   contamination = (dt_ref/DT_LIST[-1])^2 = 1/64 ~ 1.6%
DT_REF = DT_LIST[-1] / 8.0   # 1.5625e-10 s  ->  640 steps

# Mesh: fine enough for Pe < 1 in the depletion region at 1e15 doping
N_MESH = 400
L_DEVICE = 2.0e-5   # 20 um
V_BIAS = 0.1        # V, applied forward bias

# Rate floor targets
RATE_BDF1_FLOOR = 0.95
RATE_BDF2_FLOOR = 1.90

# Module-level cache: compute the BDF2 reference once for both tests
_cached_reference: tuple[np.ndarray, np.ndarray] | None = None


def _base_cfg() -> dict:
    """
    Build a minimal 1D pn junction config (Pe < 1: 1e15 doping, N=400).
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
    timestep, end time, and BDF order.
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


def _run_and_get_final_state(
    dt: float, order: int, t_end: float = T_END
) -> tuple[np.ndarray, np.ndarray]:
    """
    Run the transient solver to *t_end* and return the final (n, p) arrays
    in physical units.
    """
    from semi.runners.transient import run_transient

    cfg = _minimal_transient_cfg(dt=dt, t_end=t_end, order=order)
    result = run_transient(cfg)

    if result.fields.get("n") and result.fields["n"]:
        n_final = result.fields["n"][-1]
        p_final = result.fields["p"][-1]
    else:
        raise RuntimeError("No field snapshots were written by the transient runner.")
    return n_final, p_final


def _compute_transient_reference() -> tuple[np.ndarray, np.ndarray]:
    """
    Compute the BDF2 reference solution at dt_ref = DT_LIST[-1] / 8.

    Using the same density-form transient runner eliminates any spatial
    mismatch between reference and test.  The BDF2 reference error
    O(dt_ref^2) is negligible relative to BDF1 test errors O(dt) and
    causes only ~1.6% contamination of BDF2 test errors at the finest level.

    Result is cached at module level so both test functions share it.
    """
    global _cached_reference
    if _cached_reference is None:
        _cached_reference = _run_and_get_final_state(
            dt=DT_REF, order=2, t_end=T_END
        )
    return _cached_reference


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

    Reference: BDF2 transient at DT_REF = DT_LIST[-1]/8, same density-form
    discretisation as the test runs.  Spatial mismatch is identically zero.
    Error metric: ||u(t_end; dt_k) - u_ref(t_end; dt_ref)||_inf
    """
    n_ref, p_ref = _compute_transient_reference()

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
