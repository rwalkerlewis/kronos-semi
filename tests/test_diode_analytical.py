"""
Tests for the pure-Python diode analytical helpers in
semi.diode_analytical. No dolfinx dependency.
"""
from __future__ import annotations

import numpy as np
import pytest

from semi.constants import Q
from semi.diode_analytical import (
    depletion_width,
    shockley_long_diode_saturation,
    sns_total_reference,
    srh_generation_reference,
)

V_T = 0.025852
N_A = 1.0e23
N_D = 1.0e23
N_I = 1.0e16
EPS = 11.7 * 8.8541878128e-12
MU_N = 0.14
MU_P = 0.045
TAU_N = 1.0e-8
TAU_P = 1.0e-8


def test_shockley_saturation_matches_hand_calc():
    J_s, L_n, L_p = shockley_long_diode_saturation(
        N_A, N_D, N_I, MU_N, MU_P, TAU_N, TAU_P, V_T,
    )
    D_n = V_T * MU_N
    D_p = V_T * MU_P
    assert L_n == pytest.approx((D_n * TAU_N) ** 0.5)
    assert L_p == pytest.approx((D_p * TAU_P) ** 0.5)
    expected = Q * N_I ** 2 * (D_n / (L_n * N_A) + D_p / (L_p * N_D))
    assert J_s == pytest.approx(expected)


def test_depletion_width_vanishes_near_vbi():
    V_bi = V_T * np.log(N_A * N_D / N_I ** 2)
    W_near = depletion_width(N_A, N_D, N_I, EPS, V_T, V_bi - 1.0e-8)
    W_mid = depletion_width(N_A, N_D, N_I, EPS, V_T, 0.3)
    W_zero = depletion_width(N_A, N_D, N_I, EPS, V_T, 0.0)
    assert float(W_near) < float(W_mid) < float(W_zero)


def test_depletion_width_grows_under_reverse_bias():
    W_zero = float(depletion_width(N_A, N_D, N_I, EPS, V_T, 0.0))
    W_rev = float(depletion_width(N_A, N_D, N_I, EPS, V_T, -2.0))
    # W(V) scales as sqrt(V_bi - V); for V=-2 vs V=0 with V_bi~0.83,
    # the ratio is sqrt((0.83+2)/0.83) ~ 1.85.
    assert 1.6 < W_rev / W_zero < 2.1


def test_sns_reference_zero_at_v_zero():
    J_total, J_diff, J_rec, J_s = sns_total_reference(
        N_A, N_D, N_I, EPS, MU_N, MU_P, TAU_N, TAU_P, V_T, V=0.0,
    )
    assert float(J_total) == pytest.approx(0.0, abs=1e-12)
    assert float(J_diff) == pytest.approx(0.0, abs=1e-12)
    assert float(J_rec) == pytest.approx(0.0, abs=1e-12)
    assert J_s > 0.0


def test_sns_reference_forward_dominated_by_rec_at_low_bias():
    """At V=0.15 V with tau=1e-8 s the SNS recombination branch is
    many orders of magnitude above pure Shockley diffusion, which is
    the whole reason the M3 verifier needed this term."""
    V = 0.15
    J_total, J_diff, J_rec, _J_s = sns_total_reference(
        N_A, N_D, N_I, EPS, MU_N, MU_P, TAU_N, TAU_P, V_T, V=V,
    )
    assert float(J_rec) > 100.0 * float(J_diff)
    assert float(J_total) == pytest.approx(float(J_diff) + float(J_rec))


def test_sns_reference_vectorised_matches_scalar():
    Vs = np.array([0.1, 0.2, 0.3, 0.5])
    J_vec, _, _, _ = sns_total_reference(
        N_A, N_D, N_I, EPS, MU_N, MU_P, TAU_N, TAU_P, V_T, V=Vs,
    )
    for V, J_v in zip(Vs, J_vec, strict=True):
        J_s, _, _, _ = sns_total_reference(
            N_A, N_D, N_I, EPS, MU_N, MU_P, TAU_N, TAU_P, V_T, V=float(V),
        )
        assert float(J_s) == pytest.approx(float(J_v), rel=1e-12)


def test_srh_generation_zero_at_v_zero():
    J_net, W0, V_bi = srh_generation_reference(
        N_A, N_D, N_I, EPS, TAU_N, TAU_P, V_T, V=0.0,
    )
    assert float(J_net) == pytest.approx(0.0, abs=1e-12)
    assert W0 > 0.0
    assert V_bi == pytest.approx(V_T * np.log(N_A * N_D / N_I ** 2))


def test_srh_generation_grows_with_reverse_bias():
    Vs = np.array([-0.5, -1.0, -1.5, -2.0])
    J_net, _W0, _V_bi = srh_generation_reference(
        N_A, N_D, N_I, EPS, TAU_N, TAU_P, V_T, V=Vs,
    )
    # The helper returns the signed (q n_i / 2 tau_eff)(W(V) - W(0))
    # which is positive for V < 0 because W(V) > W(0); the verifier
    # compares |J_sim| against |J_net|.
    assert np.all(J_net > 0.0)
    # |J_net| must be non-decreasing as V gets more negative.
    assert np.all(np.diff(np.abs(J_net)) > 0.0)


def test_srh_generation_order_of_magnitude_deep_reverse():
    """For this device (tau=1e-8, N=1e17) the deep-reverse generation
    current is dominated by SRH, not Shockley diffusion, and is on the
    order of 1e-2 A/m^2. The M2 reverse-bias spec (J saturating to
    J_s ~ 1e-7 A/m^2) does not hold; the M3 verifier compares to
    this quantity instead."""
    J_net, _W0, _V_bi = srh_generation_reference(
        N_A, N_D, N_I, EPS, TAU_N, TAU_P, V_T, V=-2.0,
    )
    assert 5.0e-3 < abs(float(J_net)) < 5.0e-2
