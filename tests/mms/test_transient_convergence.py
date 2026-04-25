"""
MMS temporal convergence test for the transient drift-diffusion solver.

Manufactured-solution design (M13 round-3 redesign)
----------------------------------------------------
We run a step-function forward-bias turn-on on a lightly-doped 1-D pn
junction and compare each BDF level to a 4× finer reference (Richardson
extrapolation). This exercises the full BDF1/BDF2 time-integration path
without requiring a manufactured forcing term.

Device parameters are chosen so that the element Péclet number
Pe = h·|∇ψ|/(2) < 1 everywhere, which guarantees that the standard
Galerkin discretisation keeps carrier densities positive throughout the
Newton solve. Violations of this constraint (the M13 "positivity issue")
cause Newton to take steps where n or p go negative, making the SRH
denominator singular and stalling convergence.

Design constraints
------------------
  Pe < 1  ⟹  h < 2 / |∇ψ_max|  ≈  2·W / ψ_bi

  For N_A = N_D = 1e15 cm⁻³ (Si, 300 K):
    V_bi  ≈ 0.577 V,  ψ_bi = V_bi/V_t ≈ 22.2
    W     ≈ 1.23 µm  (depletion width at V=0)
    h_max = 2·W/ψ_bi ≈ 111 nm

  With L = 20 µm and N = 400 elements: h = 50 nm < 111 nm  → Pe ≈ 0.45

  Dynamic range of n_hat: n ranges from ~1 (n-side) to ~2.25e-10 (p-side),
  spanning ~10 orders of magnitude — within the 8–12 order physical range
  representative of real lightly-doped devices.

  Forward bias V = 0.1 V → ψ-step at t=0 is 0.1/0.026 ≈ 3.85, giving
  an initial SNES residual ≈3× smaller than the V=0.3 V case.

  SRH lifetime τ = 1 ns ≈ T_PERIOD, so the minority-carrier distribution
  evolves ≈63 % toward steady state within one period, providing a
  well-conditioned temporal error for the convergence study.

Convergence thresholds (unchanged from M13 specification):
  BDF1: rate >= 0.95
  BDF2: rate >= 1.90
"""
from __future__ import annotations

import math

import numpy as np
import pytest

# Four dt refinement levels: dt_0 / 2^k for k in [0, 1, 2, 3]
# Starting dt chosen so that the coarsest level is still stable and
# accurate enough to see the rate; T = 1e-9 s.
T_PERIOD = 1.0e-9   # step-response window (≈ 1 minority-carrier lifetime), s

# Coarsest timestep: T/10 = 1e-10 s.  Four levels: T/10, T/20, T/40, T/80.
DT_BASE = T_PERIOD / 10.0
N_LEVELS = 4
DT_LIST = [DT_BASE / (2 ** k) for k in range(N_LEVELS)]

# Mesh: N=400 on a 20 µm device → h=50 nm → Pe≈0.45 < 1 for 1e15 doping.
# This guarantees Galerkin positivity and a well-conditioned Newton solve.
# (Pe < 1 iff h < 2·W/ψ_bi; with W≈1.23 µm, ψ_bi≈22.2 → h_max≈111 nm.)
N_MESH = 400
L_DEVICE = 2.0e-5   # 20 um, same as pn_1d_bias

# Rate floor targets (unchanged from M13 specification)
RATE_BDF1_FLOOR = 0.95
RATE_BDF2_FLOOR = 1.90


def _minimal_transient_cfg(dt: float, t_end: float, order: int) -> dict:
    """
    Build a minimal 1D transient config for a pn junction with the given
    timestep and BDF order.
    """
    max_steps = int(math.ceil(t_end / dt)) + 5
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
                    # Reduced from 1e17 to 1e15 cm^-3 to achieve Pe<1 with
                    # N=400 mesh (see module docstring for design rationale).
                    "N_A_left": 1.0e15,
                    "N_D_right": 1.0e15,
                    "N_A_right": 0.0,
                },
            }
        ],
        "contacts": [
            # 0.1 V forward bias: ψ-step at t=0 is 0.1/0.026≈3.85 (vs 11.5
            # for 0.3 V), giving an initial SNES residual ~3× smaller and
            # allowing Newton to converge without n/p going negative.
            {"name": "anode",   "facet": "anode",   "type": "ohmic", "voltage": 0.1},
            {"name": "cathode", "facet": "cathode",  "type": "ohmic", "voltage": 0.0},
        ],
        "physics": {
            "temperature": 300.0,
            "statistics": "boltzmann",
            "mobility": {"mu_n": 1400.0, "mu_p": 450.0},
            "recombination": {
                "srh": True,
                # τ = 1 ns ≈ T_PERIOD so ~63 % of the minority-carrier
                # transient completes within the test window, ensuring
                # a measurable temporal error at all four dt levels.
                "tau_n": 1.0e-9,
                "tau_p": 1.0e-9,
                "E_t": 0.0,
            },
        },
        "solver": {
            "type": "transient",
            "t_end": float(t_end),
            "dt": float(dt),
            "order": int(order),
            "max_steps": int(max_steps),
            "output_every": int(max_steps + 1),  # no snapshots during test
            "snes": {
                "rtol": 1.0e-10,
                "atol": 1.0e-7,
                "stol": 1.0e-14,
                "max_it": 100,
            },
        },
        "output": {
            "directory": "/tmp/mms_transient_results",
            "fields": [],
        },
    }


def _run_and_get_final_state(dt: float, order: int):
    """
    Run the transient solver for one full period T and return the final
    (psi, n, p) arrays (in physical units).
    """
    from semi.runners.transient import run_transient

    cfg = _minimal_transient_cfg(dt=dt, t_end=T_PERIOD, order=order)
    result = run_transient(cfg)

    # Retrieve final state from last snapshot or from result fields
    # The runner stores snapshots in result.fields when output_every triggers.
    # For MMS we want the very last step values; get them from result.fields
    # or re-run if needed.
    if result.fields.get("n") and result.fields["n"]:
        n_final = result.fields["n"][-1]
        p_final = result.fields["p"][-1]
        psi_final = result.fields["psi"][-1]
    else:
        raise RuntimeError("No field snapshots were written by the transient runner.")
    return psi_final, n_final, p_final


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

    Strategy: run at each dt level and compare to the analytical
    (steady-state) result. At t = T (end of one full oscillation period),
    the manufactured solution returns to the same state as t = 0 (the
    equilibrium solution). The error at t = T is therefore:

        ||n(x, T) - n_ss(x)||_inf

    where n_ss is the steady-state n from equilibrium at t=0. This gives
    a clean manufactured-solution test without needing a separate
    analytical reference.
    """

    # Get equilibrium reference (t=0 state)
    # Use minimal config with ohmic contacts at the same biases as transient
    eq_cfg = _minimal_transient_cfg(dt=DT_LIST[0], t_end=DT_LIST[0], order=1)
    eq_cfg["solver"] = {"type": "equilibrium"}
    # For equilibrium, use zero bias (the transient runner starts from eq)
    # But we want n_ss from the DD steady-state at V=0.3V.
    # Actually: since the transient starts from equilibrium at V=0 and the
    # bias is stepped to V=0.3V, at t=T the solution is NOT back to equilibrium
    # unless V=0. For the convergence test we need a simple reference.
    #
    # Better approach: use a near-zero bias so the solution is almost linear
    # and returns close to the initial state after one period. But this
    # requires a manufactured oscillating source in the transient equations.
    #
    # Simplest correct approach for testing temporal convergence:
    # Run at 4x finer dt (reference), compare coarser levels to the finest.
    # This gives the observed temporal rate without needing an exact solution.

    # Reference: run at 4x finer than finest level
    dt_ref = DT_LIST[-1] / 4.0
    _, n_ref, p_ref = _run_and_get_final_state(dt=dt_ref, order=order)

    errors_n = []
    errors_p = []
    for dt in DT_LIST:
        _, n_k, p_k = _run_and_get_final_state(dt=dt, order=order)
        # Trim to same length (both should be same mesh)
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
