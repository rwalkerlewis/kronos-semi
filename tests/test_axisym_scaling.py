"""
Pure-Python checks for analytical MOSCAP reference values used by the
M14.2 LF/HF C-V tests.

Validates:
    * `analytical_moscap_metrics` returns sensible Cox / Cmin for a
      P-body 1e17 / 5 nm SiO2 device.
    * `infer_majority_carrier` picks holes for P-body and electrons for
      N-body doping.
    * Scaling.lambda2 stays O(1e-6 .. 1e-2) for the axisymmetric mesh,
      i.e. the radial weight does not break nondimensionalisation.
"""
from __future__ import annotations

import math

import pytest


def test_analytical_metrics_p_body():
    from semi.constants import cm3_to_m3
    from semi.physics.cv import analytical_moscap_metrics

    out = analytical_moscap_metrics(
        N_body=cm3_to_m3(1.0e17),
        body_type="p",
        t_ox=5.0e-9,
        eps_r_si=11.7,
        eps_r_ox=3.9,
        n_i=cm3_to_m3(1.0e10),
        T=300.0,
        phi_ms=-0.977,
    )
    Cox = out["Cox"]
    Cmin = out["Cmin"]
    # Cox = eps0 * 3.9 / 5e-9 ~ 6.9 mF/m^2
    assert 6.5e-3 < Cox < 7.5e-3, Cox
    # Cmin should be much less than Cox for NA=1e17 (~Cox/8 for this stack)
    assert 0.0 < Cmin < Cox
    assert 5.0e-4 < Cmin < 2.0e-3, Cmin
    # phi_F ~ kT/q * ln(1e17/1e10) ~ 0.0259 * 16.1 ~ 0.417 V
    assert 0.40 < out["phi_F"] < 0.43, out["phi_F"]
    # V_T - V_FB > 2 phi_F (depletion charge term contributes)
    assert out["V_T"] - out["V_FB"] > 2.0 * out["phi_F"]


def test_infer_majority_carrier():
    from semi.physics.cv import infer_majority_carrier

    cfg_p = {"doping": [
        {"region": "si", "profile": {"type": "uniform", "N_D": 0.0, "N_A": 1.0e17}}
    ]}
    cfg_n = {"doping": [
        {"region": "si", "profile": {"type": "uniform", "N_D": 1.0e17, "N_A": 0.0}}
    ]}
    assert infer_majority_carrier(cfg_p) == "holes"
    assert infer_majority_carrier(cfg_n) == "electrons"


def test_scaling_axisym_lambda2_well_conditioned():
    """The axisymmetric weight is geometric and must not change scaling."""
    from semi.materials import get_material
    from semi.scaling import make_scaling_from_config

    cfg = {
        "physics": {"temperature": 300.0},
        "mesh": {
            "source": "builtin",
            "extents": [[0.0, 1.0e-6], [0.0, 1.005e-6]],
            "resolution": [40, 201],
        },
        "doping": [
            {"region": "silicon",
             "profile": {"type": "uniform", "N_D": 0.0, "N_A": 1.0e17}}
        ],
    }
    sc = make_scaling_from_config(cfg, get_material("Si"))
    # lambda2 should be small (Debye << device) but not pathological.
    assert 1.0e-9 < sc.lambda2 < 1.0e-1, sc.lambda2
    # L_D2 = lambda2 * L0^2 must be non-degenerate.
    L_D2 = sc.lambda2 * sc.L0 ** 2
    assert math.isfinite(L_D2) and L_D2 > 0
