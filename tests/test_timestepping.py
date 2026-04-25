"""
Pure-Python tests for BDFCoefficients.

Tests the coefficient class without requiring dolfinx.
"""
from __future__ import annotations

import numpy as np
import pytest

from semi.timestepping import BDFCoefficients


def test_bdf1_coeffs():
    bdf = BDFCoefficients(1)
    assert bdf.order == 1
    assert bdf.coeffs == (1.0, -1.0)


def test_bdf2_coeffs():
    bdf = BDFCoefficients(2)
    assert bdf.order == 2
    assert bdf.coeffs == (1.5, -2.0, 0.5)


def test_unsupported_order_raises():
    with pytest.raises(ValueError, match="not supported"):
        BDFCoefficients(3)


def test_bdf1_apply_constant():
    """BDF1 time derivative of constant should be zero."""
    bdf = BDFCoefficients(1)
    u0 = np.array([5.0, 5.0])
    u1 = np.array([5.0, 5.0])
    result = bdf.apply([u0, u1], dt=0.1)
    np.testing.assert_allclose(result, 0.0)


def test_bdf1_apply_linear():
    """BDF1 time derivative of linear function is the slope."""
    bdf = BDFCoefficients(1)
    dt = 0.01
    slope = 3.0
    u0 = np.array([0.0])
    u1 = np.array([slope * dt])
    result = bdf.apply([u0, u1], dt=dt)
    np.testing.assert_allclose(result, [slope], rtol=1.0e-12)


def test_bdf2_apply_quadratic():
    """BDF2 is exact for quadratic polynomials u(t) = t^2.

    At t_n+1 = 2*dt: du/dt = 2*t = 4*dt.
    BDF2 uses values at t_n+1=2dt, t_n=dt, t_n-1=0.
    """
    bdf = BDFCoefficients(2)
    dt = 0.5
    u0 = np.array([0.0])    # t=0
    u1 = np.array([dt**2])  # t=dt
    u2 = np.array([(2*dt)**2])  # t=2*dt
    result = bdf.apply([u0, u1, u2], dt=dt)
    # du/dt at t=2*dt for u=t^2 is 2*2*dt = 4*dt = 2.0
    expected = np.array([2.0 * 2 * dt])
    np.testing.assert_allclose(result, expected, rtol=1.0e-12)


def test_apply_insufficient_history_raises():
    bdf = BDFCoefficients(2)
    with pytest.raises(ValueError, match="requires at least"):
        bdf.apply([np.array([1.0])], dt=0.1)  # need 3 values for BDF2
