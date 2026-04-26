"""
Steady-state limit test for the transient solver.

Run a 1D pn junction under the transient solver for long enough that
the solution converges to steady state, then compare the final IV to
the result from run_bias_sweep on the same configuration.

Acceptance criterion: relative error in J at each contact < 1e-4.

Physical basis
--------------
For t >> tau_p (minority-carrier lifetime), the transient solution must
converge to the steady-state drift-diffusion solution at the same applied
bias. This test verifies that the transient formulation in (n, p) form
gives the same terminal currents as the Slotboom-form steady-state solver.

Configuration used
------------------
- 1D symmetric pn junction, 2 um, N_A = N_D = 1e17 cm^-3
- V_anode = 0.3 V (well below 0.6 V to keep moderate injection)
- tau_n = tau_p = 1e-9 s (short lifetime so t_end doesn't need to be long)
- dt = 5e-11 s, t_end = 3e-8 s (30 * tau)
"""
from __future__ import annotations

import pytest

_L = 2.0e-6    # device length, m
_TAU = 1.0e-9  # short lifetime for fast convergence to steady state, s
_V_F = 0.3     # forward bias, V


def _make_cfg(solver_block: dict) -> dict:
    """Build a minimal 1D pn junction config with the given solver block."""
    return {
        "schema_version": "1.1.0",
        "name": "test_transient_ss_limit",
        "description": "Steady-state limit test for transient solver.",
        "dimension": 1,
        "mesh": {
            "source": "builtin",
            "extents": [[0.0, _L]],
            "resolution": [100],
            "regions_by_box": [
                {"name": "silicon", "tag": 1, "bounds": [[0.0, _L]]},
            ],
            "facets_by_plane": [
                {"name": "anode",   "tag": 1, "axis": 0, "value": 0.0},
                {"name": "cathode", "tag": 2, "axis": 0, "value": _L},
            ],
        },
        "regions": {"silicon": {"material": "Si", "tag": 1, "role": "semiconductor"}},
        "doping": [
            {
                "region": "silicon",
                "profile": {
                    "type": "step",
                    "axis": 0,
                    "location": _L / 2.0,
                    "N_D_left": 0.0,
                    "N_A_left": 1.0e17,
                    "N_D_right": 1.0e17,
                    "N_A_right": 0.0,
                },
            }
        ],
        "contacts": [
            {"name": "anode",   "facet": "anode",   "type": "ohmic", "voltage": _V_F},
            {"name": "cathode", "facet": "cathode",  "type": "ohmic", "voltage": 0.0},
        ],
        "physics": {
            "temperature": 300.0,
            "statistics": "boltzmann",
            "mobility": {"mu_n": 1400.0, "mu_p": 450.0},
            "recombination": {
                "srh": True,
                "tau_n": _TAU,
                "tau_p": _TAU,
                "E_t": 0.0,
            },
        },
        "solver": solver_block,
        "output": {
            "directory": "/tmp/test_transient_ss_limit",
            "fields": [],
        },
    }


@pytest.mark.xfail(
    reason=(
        "Known M13.1 issue: the (n,p) Galerkin transient does not "
        "converge to the same discrete steady state as the Slotboom "
        "bias_sweep, even with tight SNES tolerance and refined mesh. "
        "Tracked in GH issue #34 (Scharfetter-Gummel "
        "edge-flux discretization). See ADR 0009 'Known limitation "
        "(M13.1)'."
    ),
    strict=False,
)
def test_transient_steady_state_limit():
    """
    Run the transient solver until steady state and compare IV to
    run_bias_sweep on the same configuration.

    Asserts: max relative error in J at anode < 1e-4.
    """
    from semi.runners.bias_sweep import run_bias_sweep
    from semi.runners.transient import run_transient

    # Steady-state reference: bias_sweep to V_F in one step
    ss_cfg = _make_cfg({
        "type": "bias_sweep",
        "snes": {"rtol": 1.0e-14, "atol": 1.0e-14, "stol": 1.0e-14, "max_it": 60},
        "continuation": {
            "min_step": 1.0e-3, "max_halvings": 6,
            "max_step": 0.1, "easy_iter_threshold": 4, "grow_factor": 2.0,
        },
    })
    # Add the voltage_sweep so bias_sweep walks to V_F
    ss_cfg["contacts"][0]["voltage_sweep"] = {
        "start": 0.0, "stop": _V_F, "step": 0.05,
    }
    ss_result = run_bias_sweep(ss_cfg)

    # Extract J at anode at V = V_F from steady-state IV
    j_ss_rows = [r for r in ss_result.iv if abs(r.get("V", 0.0) - _V_F) < 1.0e-4]
    assert j_ss_rows, (
        f"No SS IV row found near V={_V_F} V; available: {[r.get('V') for r in ss_result.iv]}"
    )
    J_ss = float(j_ss_rows[-1]["J"])

    # Transient run: run for 30 * tau (time for SS limit)
    dt = 5.0e-11  # 50 ps
    t_end = 30.0 * _TAU  # 30 ns >> 1 lifetime
    max_steps = int(t_end / dt) + 10
    tr_cfg = _make_cfg({
        "type": "transient",
        "t_end": t_end,
        "dt": dt,
        "order": 2,
        "max_steps": max_steps,
        "output_every": max_steps + 1,  # snapshot only at end
        "snes": {"rtol": 1.0e-10, "atol": 1.0e-12, "stol": 1.0e-14, "max_it": 100},
    })
    tr_result = run_transient(tr_cfg)

    # Extract J at anode from the last transient IV entry
    final_iv = [r for r in tr_result.iv if r.get("contact") == "anode"]
    assert final_iv, "No IV data for anode in transient result."
    J_tr = float(final_iv[-1]["J"])

    # Compare
    if abs(J_ss) < 1.0e-30:
        abs_err = abs(J_tr - J_ss)
        assert abs_err < 1.0e-10, (
            f"J_ss ~ 0: absolute error {abs_err:.3e} A/m^2 (J_tr={J_tr:.3e})"
        )
    else:
        rel_err = abs(J_tr - J_ss) / abs(J_ss)
        assert rel_err < 1.0e-4, (
            f"Transient steady-state J={J_tr:.6e} A/m^2 differs from "
            f"bias_sweep J={J_ss:.6e} A/m^2 by {rel_err:.3e} (threshold 1e-4)"
        )
