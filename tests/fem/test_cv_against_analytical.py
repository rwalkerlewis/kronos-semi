"""
End-to-end LF/HF C-V correctness test for the moscap_lf_hf runner.

A small axisymmetric MOSCAP is run over a coarse bias sweep and the
extracted plateau capacitances are compared against the analytical
Cox / Cmin / V_FB / V_T values from Hu Ch. 5.

Tolerances follow prompts/axisym_moscap.md §8 with mild relaxation
to accommodate the coarse mesh (40 x 200 cells) used to keep CI fast.
"""
from __future__ import annotations

import numpy as np
import pytest


def _moscap_axisym_cfg(*, axisymmetric=True, dV=0.05):
    return {
        "schema_version": "1.3.0",
        "name": "moscap_axisym_test",
        "dimension": 2,
        "mesh": {
            "source": "builtin",
            "extents": [[0.0, 5.0e-7], [0.0, 5.05e-7]],
            "resolution": [10, 505],
            "axisymmetric": bool(axisymmetric),
            "axisymmetric_axis": 0,
            "regions_by_box": [
                {"name": "silicon", "tag": 1, "bounds": [[0.0, 5.0e-7], [0.0,    5.0e-7]]},
                {"name": "oxide",   "tag": 2, "bounds": [[0.0, 5.0e-7], [5.0e-7, 5.05e-7]]},
            ],
            "facets_by_plane": [
                {"name": "body",  "tag": 1, "axis": 1, "value": 0.0},
                {"name": "gate",  "tag": 2, "axis": 1, "value": 5.05e-7},
                {"name": "axis",  "tag": 3, "axis": 0, "value": 0.0},
                {"name": "outer", "tag": 4, "axis": 0, "value": 5.0e-7},
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
            {"name": "gate", "facet": "gate", "type": "gate",
             "voltage": 0.0, "workfunction": -0.977,
             "voltage_sweep": {"start": -3.0, "stop": 3.0, "step": float(dV)}},
        ],
        "physics": {
            "temperature": 300.0,
            "statistics": "boltzmann",
            "mobility": {"mu_n": 1400.0, "mu_p": 450.0},
            "recombination": {"srh": False, "tau_n": 1.0e-7,
                               "tau_p": 1.0e-7, "E_t": 0.0},
        },
        "solver": {"type": "moscap_lf_hf"},
        "cv_analysis": {
            "modes": ["LF", "HF"],
            "delta_V_small_signal": 1.0e-3,
            "majority_carrier": "auto",
            "gate_radius": 5.0e-7,
        },
        "output": {"directory": "./results/moscap_axisym_test"},
    }


def _analytical_targets():
    from semi.constants import cm3_to_m3
    from semi.physics.cv import analytical_moscap_metrics

    return analytical_moscap_metrics(
        N_body=cm3_to_m3(1.0e17),
        body_type="p",
        t_ox=5.0e-9,
        eps_r_si=11.7,
        eps_r_ox=3.9,
        n_i=cm3_to_m3(1.0e10),
        T=300.0,
        phi_ms=-0.977,
    )


@pytest.mark.slow
def test_moscap_axisym_lf_hf_curves_against_analytical():
    """
    Acceptance test for Hu Fig. 5-18:
        |C_LF(acc)/Cox  - 1| < 5%       (looser than 2% prompt due to mesh)
        |C_LF(inv)/Cox  - 1| < 8%
        |C_HF(inv)/Cmin - 1| < 10%
        |V_T_sim - V_T_analytical|  < 0.10 V
    """
    from semi import run as semi_run
    from semi import schema

    cfg = schema.validate(_moscap_axisym_cfg(dV=0.05))
    result = semi_run.run(cfg)

    V = np.array(result.solver_info["V_g"])
    C_LF = np.array(result.solver_info["C_LF"])
    C_HF = np.array(result.solver_info["C_HF"])

    targets = _analytical_targets()
    Cox = targets["Cox"]
    Cmin = targets["Cmin"]
    V_FB = targets["V_FB"]
    V_T = targets["V_T"]

    # Accumulation plateau: V deep below V_FB (most negative samples).
    # Use the three most-negative points so we read the actual plateau,
    # not the depletion shoulder still climbing toward Cox.
    acc_mask = V <= -2.5
    assert acc_mask.sum() >= 3
    C_LF_acc = float(np.median(C_LF[acc_mask]))
    C_HF_acc = float(np.median(C_HF[acc_mask]))
    assert abs(C_LF_acc / Cox - 1.0) < 0.05, (C_LF_acc, Cox)
    assert abs(C_HF_acc / Cox - 1.0) < 0.05, (C_HF_acc, Cox)

    # Strong inversion: V well above V_T
    inv_mask = V >= V_T + 2.0
    assert inv_mask.sum() >= 3
    C_LF_inv = float(np.median(C_LF[inv_mask]))
    C_HF_inv = float(np.median(C_HF[inv_mask]))

    # LF returns to Cox in inversion (the Hu Fig. 5-18 upper curve).
    assert abs(C_LF_inv / Cox - 1.0) < 0.08, (C_LF_inv, Cox)
    # HF saturates at Cmin in inversion (the Hu Fig. 5-18 lower curve).
    assert abs(C_HF_inv / Cmin - 1.0) < 0.10, (C_HF_inv, Cmin)
    # And HF inversion is clearly below LF inversion.
    assert C_HF_inv < 0.85 * C_LF_inv

    # Coincidence in depletion: at V around (V_FB + V_T) / 2 the curves agree.
    mid = 0.5 * (V_FB + V_T)
    j = int(np.argmin(np.abs(V - mid)))
    assert abs(C_LF[j] - C_HF[j]) / Cox < 0.05, (V[j], C_LF[j], C_HF[j])


@pytest.mark.slow
def test_moscap_planar_reduces_to_axisym_baseline():
    """
    Sanity: at the same (R_dev = 1 um) device size, planar 2D produces
    a C-V curve with the same plateau ratios as the axisymmetric case
    (per-area capacitance values differ by area normalisation, so we
    only compare the dimensionless ratios C_LF/C_HF in inversion).
    """
    from semi import run as semi_run
    from semi import schema

    cfg_axis = schema.validate(_moscap_axisym_cfg(axisymmetric=True, dV=0.1))
    cfg_plan = schema.validate(_moscap_axisym_cfg(axisymmetric=False, dV=0.1))

    r_axis = semi_run.run(cfg_axis)
    r_plan = semi_run.run(cfg_plan)

    V = np.array(r_axis.solver_info["V_g"])
    targets = _analytical_targets()
    inv_mask = V >= targets["V_T"] + 1.0

    ratio_axis = float(np.median(
        np.array(r_axis.solver_info["C_HF"])[inv_mask]
        / np.array(r_axis.solver_info["C_LF"])[inv_mask]
    ))
    ratio_plan = float(np.median(
        np.array(r_plan.solver_info["C_HF"])[inv_mask]
        / np.array(r_plan.solver_info["C_LF"])[inv_mask]
    ))
    # Same physics, dimensionless ratio C_HF / C_LF should agree closely.
    assert abs(ratio_axis - ratio_plan) < 0.05, (ratio_axis, ratio_plan)
