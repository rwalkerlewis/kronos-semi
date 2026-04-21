"""Tests for semi.scaling — no dolfinx required."""
import math

import pytest

from semi import materials, scaling
from semi.constants import EPS0, Q, thermal_voltage


def test_thermal_voltage_300k():
    assert thermal_voltage(300.0) == pytest.approx(0.025852, abs=1e-5)


def test_scaling_basic():
    sc = scaling.Scaling(L0=1.0e-6, C0=1.0e23, T=300.0, mu0=0.14, n_i=1.0e16)
    assert sc.V0 == pytest.approx(0.025852, abs=1e-5)
    assert sc.D0 == pytest.approx(sc.V0 * 0.14)
    assert sc.t0 == pytest.approx(sc.L0**2 / sc.D0)
    assert sc.J0 == pytest.approx(Q * sc.D0 * sc.C0 / sc.L0)


def test_lambda_squared_matches_debye_formula():
    """lambda^2 * eps_r should equal (Debye/L)^2 for reference density C0."""
    eps_r = 11.7
    C0 = 1.0e23
    L0 = 1.0e-6
    sc = scaling.Scaling(L0=L0, C0=C0, T=300.0, mu0=0.14, n_i=1.0e16)
    L_D = math.sqrt(eps_r * EPS0 * sc.V0 / (Q * C0))
    assert (sc.lambda2 * eps_r) == pytest.approx((L_D / L0)**2, rel=1e-10)


def test_debye_length_reasonable_for_si_1e17():
    """Textbook: L_D(Si, 1e17 cm^-3, 300K) ~ 13 nm."""
    sc = scaling.Scaling(L0=2.0e-6, C0=1.0e23, T=300.0, mu0=0.14, n_i=1.0e16)
    # Debye with eps_r accounted separately:
    import math
    L_D_si = math.sqrt(11.7 * EPS0 * sc.V0 / (Q * 1e23))
    assert 10e-9 < L_D_si < 16e-9


def test_make_scaling_from_config():
    si = materials.get_material("Si")
    cfg = {
        "physics": {"temperature": 300.0},
        "mesh": {"source": "builtin", "extents": [[0.0, 2.0e-6]], "resolution": [100]},
        "doping": [{
            "region": "si",
            "profile": {"type": "step", "axis": 0, "location": 1e-6,
                        "N_D_left": 0.0, "N_A_left": 1e17,
                        "N_D_right": 1e17, "N_A_right": 0.0},
        }],
    }
    sc = scaling.make_scaling_from_config(cfg, si)
    assert sc.L0 == pytest.approx(2.0e-6)
    assert sc.C0 == pytest.approx(1.0e23)  # peak in m^-3
    assert sc.T == 300.0


def test_make_scaling_floor():
    """C0 should have a 1e16 cm^-3 floor so tiny doping doesn't produce huge lambda^2."""
    si = materials.get_material("Si")
    cfg = {
        "physics": {"temperature": 300.0},
        "mesh": {"source": "builtin", "extents": [[0.0, 1e-6]], "resolution": [50]},
        "doping": [{
            "region": "si",
            "profile": {"type": "uniform", "N_D": 1e10, "N_A": 0.0},  # way below floor
        }],
    }
    sc = scaling.make_scaling_from_config(cfg, si)
    # floor is 1e16 cm^-3 = 1e22 m^-3
    assert sc.C0 == pytest.approx(1.0e22)
