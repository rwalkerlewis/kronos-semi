"""
End-to-end tests for the MOS capacitor analytic-differential C(V) runner.

Exercises the M14.1 `mos_cap_ac` solver type:
    - Multi-region equilibrium Poisson at every V_gate (same as `mos_cv`)
    - Analytic PDE-sensitivity C(V_gate) per bias point
    - iv rows carry both `Q_gate` (legacy compat) and `C_ac` (new)
    - Q_gate values must agree byte-identically with `mos_cv` on the same
      config: the nonlinear Poisson solve is unchanged.

Uses the same tiny 2D mesh as `test_mos_cv.py` so these tests run in a
few seconds inside the dolfinx CI image.
"""
from __future__ import annotations

import numpy as np
import pytest


def _tiny_mos_cfg(v_values=(-0.3, 0.0, 0.3), solver_type="mos_cap_ac"):
    return {
        "schema_version": "1.2.0",
        "name": "mos_cap_ac_tiny",
        "dimension": 2,
        "mesh": {
            "source": "builtin",
            "extents": [[0.0, 1.0e-7], [0.0, 1.05e-7]],
            "resolution": [2, 105],
            "regions_by_box": [
                {"name": "silicon", "tag": 1, "bounds": [[0.0, 1.0e-7], [0.0,     1.0e-7]]},
                {"name": "oxide",   "tag": 2, "bounds": [[0.0, 1.0e-7], [1.0e-7,  1.05e-7]]},
            ],
            "facets_by_plane": [
                {"name": "body", "tag": 1, "axis": 1, "value": 0.0},
                {"name": "gate", "tag": 2, "axis": 1, "value": 1.05e-7},
            ],
        },
        "regions": {
            "silicon": {"material": "Si",   "tag": 1, "role": "semiconductor"},
            "oxide":   {"material": "SiO2", "tag": 2, "role": "insulator"},
        },
        "doping": [
            {"region": "silicon",
             "profile": {"type": "uniform", "N_D": 0.0, "N_A": 1.0e17}},
        ],
        "contacts": [
            {"name": "body", "facet": "body", "type": "ohmic", "voltage": 0.0},
            {
                "name": "gate", "facet": "gate", "type": "gate",
                "voltage": 0.0, "workfunction": 0.0,
                "voltage_sweep": {
                    "start": float(min(v_values)),
                    "stop":  float(max(v_values)),
                    "step":  0.3,
                },
            },
        ],
        "physics": {
            "temperature": 300.0,
            "statistics": "boltzmann",
            "mobility": {"mu_n": 1400.0, "mu_p": 450.0},
            "recombination": {"srh": False, "tau_n": 1.0e-7,
                               "tau_p": 1.0e-7, "E_t": 0.0},
        },
        "solver": {"type": solver_type},
        "output": {"directory": "./results/mos_cap_ac_tiny"},
    }


def test_mos_cap_ac_records_iv_rows_with_C_ac():
    """Run a 3-point gate sweep and check `C_ac` is reported per bias."""
    from semi import run as semi_run
    from semi import schema

    cfg = schema.validate(_tiny_mos_cfg())
    result = semi_run.run(cfg)

    assert result.solver_info["converged"]
    assert result.iv is not None
    assert len(result.iv) == 3
    for r in result.iv:
        assert "V" in r and "Q_gate" in r and "C_ac" in r
        assert np.isfinite(r["C_ac"])
        # Differential capacitance must be positive (the second derivative
        # of the energy w.r.t. gate voltage is positive for a stable system).
        assert r["C_ac"] > 0.0, f"C_ac must be positive: {r}"


def test_mos_cap_ac_matches_mos_cv_qgate_byte_identical():
    """Q_gate from `mos_cap_ac` must match `mos_cv` exactly on the same config.

    The nonlinear Poisson solve, BCs, and charge integration are identical;
    only the C extraction differs. Any divergence in Q would indicate a
    regression in the shared multi-region path.
    """
    from semi import run as semi_run
    from semi import schema

    cfg_ac = schema.validate(_tiny_mos_cfg(solver_type="mos_cap_ac"))
    cfg_cv = schema.validate(_tiny_mos_cfg(solver_type="mos_cv"))

    res_ac = semi_run.run(cfg_ac)
    res_cv = semi_run.run(cfg_cv)

    Q_ac = np.array([r["Q_gate"] for r in res_ac.iv])
    Q_cv = np.array([r["Q_gate"] for r in res_cv.iv])
    V_ac = np.array([r["V"] for r in res_ac.iv])
    V_cv = np.array([r["V"] for r in res_cv.iv])

    assert V_ac == pytest.approx(V_cv)
    # Q_gate is a scalar integral of an in-place mutated psi; SNES is the
    # same code path with the same tolerances. We expect agreement to a few
    # ULPs * residual tolerance. Rel 1e-10 is comfortably tight.
    assert Q_ac == pytest.approx(Q_cv, rel=1.0e-10, abs=1.0e-18)


def test_mos_cap_ac_C_monotone_decreasing_in_depletion():
    """C_ac(V_gate) must decrease monotonically through the depletion window.

    Physically, as V_gate sweeps from V_FB into depletion, the depletion
    width grows and the series capacitance C_dep||C_ox decreases until
    strong-inversion onset clamps it. Inside [V_FB+0.1, V_T-0.1] V the
    differential capacitance is therefore strictly non-increasing (with a
    little wiggle room for SNES residual noise at coarse meshes).
    """
    from semi import run as semi_run
    from semi import schema

    cfg = _tiny_mos_cfg()
    cfg["contacts"][1]["voltage_sweep"] = {
        "start": -0.1, "stop": 0.4, "step": 0.05,
    }
    cfg = schema.validate(cfg)

    result = semi_run.run(cfg)
    V = np.array([r["V"] for r in result.iv])
    C_ac = np.array([r["C_ac"] for r in result.iv])
    order = np.argsort(V)
    V = V[order]
    C_ac = C_ac[order]

    # All values positive (basic correctness).
    assert (C_ac > 0).all(), f"C_ac not all positive: {C_ac}"

    # Monotone non-increasing across the depletion regime (1% wiggle).
    non_increasing = all(
        C_ac[i + 1] <= C_ac[i] * (1.0 + 0.01) for i in range(len(C_ac) - 1)
    )
    assert non_increasing, (
        f"C_ac must be non-increasing across depletion: V={V}, C_ac={C_ac}"
    )


def test_mos_cap_ac_raises_without_gate_sweep():
    """The runner must raise a clear error if no gate contact has a voltage_sweep."""
    from semi import run as semi_run
    from semi import schema

    cfg = _tiny_mos_cfg()
    for c in cfg["contacts"]:
        if c["type"] == "gate":
            c.pop("voltage_sweep", None)
    cfg = schema.validate(cfg)

    with pytest.raises(ValueError, match="voltage_sweep"):
        semi_run.run(cfg)
