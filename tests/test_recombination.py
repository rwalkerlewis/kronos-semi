"""Tests for the SRH recombination kernel (pure-Python path)."""
from __future__ import annotations

import numpy as np
import pytest

from semi.physics import recombination

N_I = 1.0e16              # m^-3, Si at 300 K (order of magnitude)
TAU_N = 1.0e-7            # s
TAU_P = 1.0e-7            # s


def _np_rate(n, p, E_t_over_Vt=0.0):
    return recombination.srh_rate_np(n, p, N_I, TAU_N, TAU_P, E_t_over_Vt)


def test_srh_zero_at_equilibrium_intrinsic():
    """R = 0 when n p = n_i^2 at the intrinsic point."""
    n = N_I
    p = N_I
    assert _np_rate(n, p) == pytest.approx(0.0, abs=1e-30)


def test_srh_zero_at_equilibrium_extrinsic():
    """R = 0 for any (n, p) obeying mass action n p = n_i^2."""
    n = 1.0e22
    p = N_I ** 2 / n
    assert _np_rate(n, p) == pytest.approx(0.0, rel=0, abs=1.0e-10)


def test_srh_positive_under_forward_bias():
    """Above equilibrium (n p > n_i^2) recombination is strictly positive."""
    n = 1.0e18
    p = 1.0e16
    assert n * p > N_I ** 2
    assert _np_rate(n, p) > 0.0


def test_srh_negative_under_reverse_bias():
    """Below equilibrium (n p < n_i^2) the kernel models net generation."""
    n = 1.0e6
    p = 1.0e6
    assert n * p < N_I ** 2
    assert _np_rate(n, p) < 0.0


def test_srh_low_level_limit_n_type():
    """
    Low-level injection in n-type (n0 >> p, delta_p << n0): R ~ delta_p / tau_p.
    """
    n0 = 1.0e22
    p0 = N_I ** 2 / n0
    delta_p = 1.0e12              # << n0
    n = n0                        # delta_n ~ delta_p << n0 so n stays
    p = p0 + delta_p
    R = _np_rate(n, p)
    assert R == pytest.approx(delta_p / TAU_P, rel=1e-3)


def test_srh_low_level_limit_p_type():
    """Low-level injection in p-type: R ~ delta_n / tau_n."""
    p0 = 1.0e22
    n0 = N_I ** 2 / p0
    delta_n = 1.0e12
    n = n0 + delta_n
    p = p0
    R = _np_rate(n, p)
    assert R == pytest.approx(delta_n / TAU_N, rel=1e-3)


def test_srh_high_level_saturation():
    """
    High-level injection (delta_n = delta_p >> doping): R ~ delta_n / (tau_n + tau_p).
    """
    delta = 1.0e22
    n = delta
    p = delta
    R = _np_rate(n, p)
    # At delta >> n_i, num ~ delta^2 and den ~ (tau_p + tau_n) delta.
    assert R == pytest.approx(delta / (TAU_N + TAU_P), rel=1e-6)


def test_srh_deep_trap_suppresses_rate():
    """Traps far from mid-gap make recombination less efficient."""
    n = 1.0e18
    p = 1.0e16
    R_mid = _np_rate(n, p, E_t_over_Vt=0.0)
    R_deep = _np_rate(n, p, E_t_over_Vt=10.0)
    assert 0.0 < R_deep < R_mid


def test_srh_asymmetric_lifetimes_limit_on_slower_tau():
    """
    With tau_n >> tau_p in strongly n-type material (n >> p), R is set by
    the minority lifetime tau_p, not tau_n.
    """
    n0 = 1.0e22
    delta_p = 1.0e14
    p = N_I ** 2 / n0 + delta_p
    R = recombination.srh_rate_np(n0, p, N_I, tau_n=1.0e-3, tau_p=1.0e-7)
    assert R == pytest.approx(delta_p / 1.0e-7, rel=2e-2)


def test_srh_vector_input_matches_scalar():
    """Vectorized evaluation equals element-wise scalar evaluation."""
    rng = np.random.default_rng(1)
    n = rng.uniform(1.0e16, 1.0e20, size=20)
    p = rng.uniform(1.0e10, 1.0e16, size=20)
    R_vec = _np_rate(n, p)
    R_loop = np.array([_np_rate(float(nn), float(pp)) for nn, pp in zip(n, p, strict=True)])
    np.testing.assert_allclose(R_vec, R_loop, rtol=1e-12)


def test_scaled_tau_conversion():
    t0 = 1.0e-9
    assert recombination.scaled_tau(1.0e-7, t0) == pytest.approx(100.0)
