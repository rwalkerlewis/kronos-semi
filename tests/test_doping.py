"""Tests for semi.doping — no dolfinx required."""
import numpy as np
import pytest

from semi import doping
from semi.constants import cm3_to_m3


def test_uniform():
    spec = [{"region": "x", "profile": {"type": "uniform", "N_D": 1e17, "N_A": 0.0}}]
    f = doping.build_profile(spec)
    x = np.array([[0.0, 0.5, 1.0]])
    result = f(x)
    assert np.allclose(result, cm3_to_m3(1e17))


def test_uniform_with_acceptor():
    spec = [{"region": "x", "profile": {"type": "uniform", "N_D": 0.0, "N_A": 1e17}}]
    f = doping.build_profile(spec)
    x = np.array([[0.5]])
    # N_D - N_A = -1e17 cm^-3 = -1e23 m^-3
    assert np.allclose(f(x), cm3_to_m3(-1e17))


def test_step():
    spec = [{
        "region": "x",
        "profile": {
            "type": "step", "axis": 0, "location": 0.5,
            "N_D_left": 0.0, "N_A_left": 1e17,
            "N_D_right": 1e17, "N_A_right": 0.0,
        },
    }]
    f = doping.build_profile(spec)
    x = np.array([[0.25, 0.75]])  # left side, right side
    result = f(x)
    assert result[0] == pytest.approx(cm3_to_m3(-1e17))  # p-side
    assert result[1] == pytest.approx(cm3_to_m3(+1e17))  # n-side


def test_gaussian_donor():
    spec = [{
        "region": "x",
        "profile": {
            "type": "gaussian",
            "center": [0.5], "sigma": [0.1],
            "peak": 1e18, "dopant": "donor",
        },
    }]
    f = doping.build_profile(spec)
    x = np.array([[0.5, 0.6, 1.0]])  # peak, 1 sigma, far
    result = f(x)
    # At peak: should match peak value (no background)
    assert result[0] == pytest.approx(cm3_to_m3(1e18))
    # At 1 sigma: exp(-0.5) = 0.6065
    assert result[1] == pytest.approx(cm3_to_m3(1e18) * np.exp(-0.5), rel=1e-6)
    # Far away: essentially zero
    assert abs(result[2]) < cm3_to_m3(1e18) * 1e-5


def test_gaussian_acceptor():
    spec = [{
        "region": "x",
        "profile": {
            "type": "gaussian",
            "center": [0.0], "sigma": [0.1],
            "peak": 1e18, "dopant": "acceptor",
        },
    }]
    f = doping.build_profile(spec)
    x = np.array([[0.0]])
    # Acceptor -> contributes NEGATIVE to N_D - N_A
    assert f(x)[0] == pytest.approx(-cm3_to_m3(1e18))


def test_gaussian_2d():
    spec = [{
        "region": "x",
        "profile": {
            "type": "gaussian",
            "center": [0.0, 0.0], "sigma": [0.1, 0.2],
            "peak": 1e18, "dopant": "donor",
        },
    }]
    f = doping.build_profile(spec)
    x = np.array([[0.1, 0.0], [0.0, 0.2]])  # shape (2, 2): two 2D points
    result = f(x)
    # Both points are at exactly 1 sigma in one direction
    expected = cm3_to_m3(1e18) * np.exp(-0.5)
    assert np.allclose(result, expected, rtol=1e-6)


def test_multiple_profiles_superpose():
    """Background uniform + Gaussian implant should sum."""
    spec = [
        {"region": "x", "profile": {"type": "uniform", "N_D": 0.0, "N_A": 1e16}},
        {"region": "x", "profile": {
            "type": "gaussian",
            "center": [0.0], "sigma": [0.1],
            "peak": 1e20, "dopant": "donor",
        }},
    ]
    f = doping.build_profile(spec)
    # At origin: 1e20 donor - 1e16 acceptor = ~1e20 (donor dominates)
    r = f(np.array([[0.0]]))[0]
    assert r == pytest.approx(cm3_to_m3(1e20 - 1e16), rel=1e-4)
    # Far away: just the background
    r_far = f(np.array([[10.0]]))[0]
    assert r_far == pytest.approx(-cm3_to_m3(1e16), rel=1e-4)


def test_input_shape_1d_accepted():
    """Should handle x of shape (N,) by treating it as 1D."""
    spec = [{"region": "x", "profile": {"type": "uniform", "N_D": 1e17, "N_A": 0.0}}]
    f = doping.build_profile(spec)
    result = f(np.array([0.0, 0.5, 1.0]))
    assert result.shape == (3,)
    assert np.allclose(result, cm3_to_m3(1e17))
