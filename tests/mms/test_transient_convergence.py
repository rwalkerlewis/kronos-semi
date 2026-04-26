"""
Temporal convergence test for the BDF1/BDF2 transient drift-diffusion solver.

Pairwise-difference design (M13 final)
--------------------------------------
We measure the temporal discretisation error of the (psi, n, p) form
transient runner by comparing solutions at consecutive dt levels on the
same mesh.

Why pairwise differences (not bias_sweep, not self-Richardson):
* The transient runner discretises the (psi, n, p) form. The bias_sweep
  runner discretises the Slotboom (psi, phi_n, phi_p) form. On the same
  mesh, the two formulations give *different* discrete steady-state
  solutions because their Galerkin residuals are different (the test
  function for the carrier-density variable differs). The difference
  ~3e18 m^-3 in n at the depletion edge does not vanish as dt -> 0, so
  bias_sweep cannot serve as a reference for the n-form temporal error.
* The previous self-Richardson approach used a reference at the same
  BDF order but only 4x finer dt than the finest test level, which
  contaminated the comparison with the reference's own temporal error
  and produced non-monotonic, meaningless rates.

Pairwise differences avoid both pitfalls:

  Let u(dt; t_end) be the n-form transient solution at the same mesh,
  same BDF order, computed with timestep dt. Its error vs the unknown
  true solution u_exact(t_end) of the *continuous* (n, p) PDE is

      u(dt; t_end) - u_exact(t_end)  =  C(t_end) * dt^p  +  e_spatial

  where p is the BDF order and e_spatial is the (dt-independent)
  spatial discretisation error. Therefore for two consecutive levels
  dt_k and dt_{k+1} = dt_k / 2,

      D_k  :=  || u(dt_k) - u(dt_{k+1}) ||_inf
            =  || C(t_end) * (dt_k^p - dt_{k+1}^p) ||_inf
            ~  |C(t_end)|_inf * dt_k^p * (1 - 2^{-p})

  The spatial error cancels exactly because both runs share mesh and
  formulation. The log-log slope of D_k vs dt_k recovers p.

Pe < 1 design (preserved from M13 round 3)
------------------------------------------
Device parameters keep the element Peclet number below 1 so the
standard Galerkin discretisation maintains positive carrier densities
through Newton iteration:

  N_A = N_D = 1e15 cm^-3 (Si, 300 K) -> V_bi ~ 0.577 V, W ~ 1.23 um
  ψ_bi ~ 22.2,  h_max = 2*W/ψ_bi ~ 111 nm
  L = 20 um, N = 400 -> h = 50 nm < 111 nm  -> Pe ~ 0.45

Convergence thresholds (unchanged from M13 specification):
  BDF1: rate >= 0.95
  BDF2: rate >= 1.90
"""
from __future__ import annotations

import math

import numpy as np
import pytest

# Forward bias for the step-on transient.
V_F = 0.1

# SRH lifetime; sets the time scale of the minority-carrier transient.
TAU = 1.0e-9

# Stop at t_end = 1 * tau, near the peak of the temporal-error envelope.
# For a stable linear system u' = -u/tau, the BDF temporal error scales
# as exp(-t/tau) * t * dt^p — small for t << tau (system hasn't moved)
# and exponentially decaying for t >> tau (errors damp toward steady
# state). At t ~ tau the error is still O(dt^p) and easy to resolve.
# Going to t_end = 10 * tau hides the temporal signal in SS noise floor.
T_END = 1.0 * TAU  # 1 ns

# DT_LIST is chosen so the coarsest level still has enough steps that
# BDF2 (which uses one BDF1 step to seed its history) is not dominated
# by that startup step. With T_END / 16 = 6.25e-11 s as DT_BASE, the
# coarsest run does 16 steps (15 of them BDF2 in BDF2 mode) and the
# finest does 128 steps.
DT_BASE = T_END / 16.0  # 6.25e-11 s
N_LEVELS = 4
DT_LIST = [DT_BASE / (2 ** k) for k in range(N_LEVELS)]

# Mesh: N=400 on a 20 um device -> h=50 nm -> Pe ~ 0.45 < 1 for 1e15 doping.
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
                    "N_A_left": 1.0e15,
                    "N_D_right": 1.0e15,
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
            "snes": {
                # Tight atol because the transient residual at this device
                # (1e15 doping, V=0.1 V) is small in scaled units; the M13
                # default atol=1e-7 lets SNES "converge" at iteration 0
                # without updating the solution, freezing the transient.
                # 1e-14 forces real Newton work each step.
                "rtol": 1.0e-12,
                "atol": 1.0e-14,
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
    (n_final, p_final) as numpy arrays in physical units (m^-3).

    The transient runner snapshots fields whenever
    `step_count % output_every == 0` OR `|t_current - t_end| < 0.5*dt`,
    so with output_every set above max_steps the only snapshot written
    is at t = t_end. fields["n"][-1] is therefore the final state.
    """
    from semi.runners.transient import run_transient
    cfg = _transient_cfg(dt=dt, order=order)
    result = run_transient(cfg)
    if not result.fields.get("n"):
        raise RuntimeError("Transient runner returned no field snapshots.")
    return (
        np.asarray(result.fields["n"][-1]).copy(),
        np.asarray(result.fields["p"][-1]).copy(),
    )


def _log_rate(dts, errors):
    """Log-log slope of the last two (dt, err) pairs."""
    dt_prev, dt_cur = dts[-2], dts[-1]
    e_prev, e_cur = errors[-2], errors[-1]
    if e_cur <= 0.0 or e_prev <= 0.0:
        return float("nan")
    return math.log(e_prev / e_cur) / math.log(dt_prev / dt_cur)


def _run_convergence_study(order: int) -> tuple[list[float], list[float], list[float]]:
    """
    Run the temporal convergence study for the given BDF order.

    Returns (dts_for_diff, diffs_n, diffs_p), where:
      * dts_for_diff[k] = DT_LIST[k]    (the larger dt of the pair)
      * diffs_n[k] = ||n(dt_k) - n(dt_{k+1})||_inf
      * diffs_p[k] = ||p(dt_k) - p(dt_{k+1})||_inf

    Length of all three is N_LEVELS - 1 (3 entries from 4 dt levels).
    """
    states_n: list[np.ndarray] = []
    states_p: list[np.ndarray] = []
    for dt in DT_LIST:
        n_k, p_k = _run_transient_final_state(dt=dt, order=order)
        states_n.append(n_k)
        states_p.append(p_k)
        if len(states_n) > 1 and states_n[-1].shape != states_n[0].shape:
            raise RuntimeError(
                f"DOF count mismatch between dt levels: {states_n[0].shape} "
                f"vs {states_n[-1].shape}"
            )

    dts_for_diff = DT_LIST[:-1]
    diffs_n = [
        float(np.max(np.abs(states_n[k] - states_n[k + 1])))
        for k in range(N_LEVELS - 1)
    ]
    diffs_p = [
        float(np.max(np.abs(states_p[k] - states_p[k + 1])))
        for k in range(N_LEVELS - 1)
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
    Asserts observed rate >= 0.95 (theoretical: 1.0) and that the
    pairwise-difference sequence decreases monotonically with dt.
    """
    dts, diffs_n, diffs_p = _run_convergence_study(order=1)
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
    Asserts observed rate >= 1.9 (theoretical: 2.0) and that the
    pairwise-difference sequence decreases monotonically with dt.
    """
    dts, diffs_n, diffs_p = _run_convergence_study(order=2)
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
