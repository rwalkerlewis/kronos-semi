"""
Unit tests for the Scharfetter-Gummel Bernoulli + 1D edge-flux primitives
(M13.1, ADR 0012).

Pure-Python tests: no dolfinx import. Tests run in the pure-python CI
matrix on Python 3.10 / 3.11 / 3.12.

Test coverage
-------------
1. Bernoulli regimes:
    - B(0) = 1 exactly.
    - Taylor branch (|x| < 1e-3) matches mpmath to 1e-12 relative.
    - Mid-range branch (1e-3 <= |x| <= 30) matches mpmath to 1e-12 relative.
    - Large positive (x > 30) matches mpmath to 1e-12 relative.
    - Large negative (x < -30) matches mpmath asymptote `B(x) -> -x` to 1e-12.
    - 1000-point log-uniform-in-|x| sample on [-50, 50], mixed signs:
      relative error vs mpmath stays below 1e-12 throughout.

2. Algebraic identities:
    - `B(x) - B(-x) = -x` (verified in scaled units across the sample).
      This is the corrected identity; the prompt-draft form `B(x) + B(-x) = -x`
      is wrong (the sum is `x coth(x/2)`).

3. Vectorised `bernoulli_array` agrees with scalar `bernoulli` to bit-identity.

4. Edge flux sign / scaling guards (ADR 0012 test plan #2 and #3):
    - Electron sign+scaling at small Peclet: SG flux agrees with the
      midpoint-Galerkin reference to within 5% relative.
    - Electron two-direction sign: swap (i, j), flux flips sign with
      same magnitude.
    - **Hole guard 2b**: device-equation symmetry under
      `(n <-> p, psi -> -psi)`. Hole flux equals the negative of the
      mirrored electron flux to 1e-10 relative. The hole sign is the
      most error-prone part of SG; this guard is non-negotiable.
    - Hole two-direction sign: swap (i, j), flux flips sign with
      same magnitude.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

mpmath = pytest.importorskip("mpmath")

from semi.fem.scharfetter_gummel import (  # noqa: E402
    bernoulli,
    bernoulli_array,
    midpoint_galerkin_flux_n,
    sg_edge_flux_n,
    sg_edge_flux_n_array,
    sg_edge_flux_p,
    sg_edge_flux_p_array,
)

# --------------------------------------------------------------------------- #
# mpmath reference                                                            #
# --------------------------------------------------------------------------- #

def _bernoulli_mp(x: float, dps: int = 50) -> float:
    """High-precision reference: B(x) = x / (exp(x) - 1) at `dps` decimals."""
    if x == 0.0:
        return 1.0
    with mpmath.workdps(dps):
        xm = mpmath.mpf(repr(x))  # repr to avoid float -> mpf rounding surprises
        return float(xm / (mpmath.exp(xm) - 1))


# --------------------------------------------------------------------------- #
# Bernoulli regime tests                                                      #
# --------------------------------------------------------------------------- #

def test_bernoulli_at_zero_is_exactly_one():
    assert bernoulli(0.0) == 1.0


def test_bernoulli_at_zero_array_is_exactly_one():
    out = bernoulli_array(np.array([0.0, 0.0, 0.0]))
    assert (out == 1.0).all()


@pytest.mark.parametrize("x", [
    1.0e-12, 1.0e-9, 1.0e-6, 1.0e-4, 5.0e-4,
    -1.0e-12, -1.0e-9, -1.0e-6, -1.0e-4, -5.0e-4,
])
def test_bernoulli_taylor_regime_matches_mpmath(x):
    rel = abs(bernoulli(x) - _bernoulli_mp(x)) / abs(_bernoulli_mp(x))
    assert rel < 1.0e-12, f"x={x}: rel={rel}"


@pytest.mark.parametrize("x", [
    1.0e-3, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 29.99,
    -1.0e-3, -0.1, -0.5, -1.0, -2.0, -5.0, -10.0, -20.0, -29.99,
])
def test_bernoulli_mid_regime_matches_mpmath(x):
    rel = abs(bernoulli(x) - _bernoulli_mp(x)) / abs(_bernoulli_mp(x))
    assert rel < 1.0e-12, f"x={x}: rel={rel}"


@pytest.mark.parametrize("x", [30.001, 35.0, 40.0, 50.0, 100.0])
def test_bernoulli_large_positive_matches_mpmath(x):
    # B(x) -> 0 as x -> +inf; mpmath is the reference.
    ref = _bernoulli_mp(x)
    got = bernoulli(x)
    if abs(ref) > 0:
        rel = abs(got - ref) / abs(ref)
        assert rel < 1.0e-12, f"x={x}: rel={rel}"
    else:
        assert got == 0.0 or abs(got) < 1.0e-300


@pytest.mark.parametrize("x", [-30.001, -35.0, -40.0, -50.0, -100.0])
def test_bernoulli_large_negative_matches_mpmath(x):
    ref = _bernoulli_mp(x)
    rel = abs(bernoulli(x) - ref) / abs(ref)
    assert rel < 1.0e-12, f"x={x}: rel={rel}"


def test_bernoulli_1000_random_points_logspace_match_mpmath():
    """
    1000-point log-uniform-in-|x| sample on [-50, 50], mixed signs.
    Relative error vs mpmath stays below 1e-12.
    """
    rng = np.random.default_rng(seed=0xB1B2)
    abs_x = 10.0 ** rng.uniform(-12.0, math.log10(50.0), size=1000)
    signs = rng.choice([-1.0, 1.0], size=1000)
    xs = signs * abs_x

    worst = 0.0
    for x in xs:
        ref = _bernoulli_mp(float(x))
        got = bernoulli(float(x))
        if abs(ref) == 0.0:
            assert abs(got) < 1.0e-300, f"x={x}: got={got} but ref=0"
            continue
        rel = abs(got - ref) / abs(ref)
        worst = max(worst, rel)
    assert worst < 1.0e-12, f"worst relative error = {worst}"


# --------------------------------------------------------------------------- #
# Identity tests                                                              #
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("x", [
    -25.0, -10.0, -1.0, -0.1, -1.0e-4, 1.0e-4, 0.1, 1.0, 10.0, 25.0,
])
def test_bernoulli_minus_minus_x_equals_minus_x(x):
    """
    Algebraic identity B(x) - B(-x) = -x.

    Numerator e^x + e^{-x} - 2 = 2 cosh(x) - 2.
    Denominator (e^x - 1)(e^{-x} - 1) = -(2 cosh(x) - 2).
    Ratio: -1, so B(x) - B(-x) = -x.

    The earlier prompt draft stated B(x) + B(-x) = -x, which is wrong;
    the sum is x coth(x/2).
    """
    diff = bernoulli(x) - bernoulli(-x)
    assert diff == pytest.approx(-x, rel=1.0e-12, abs=1.0e-14)


def test_bernoulli_sum_is_x_coth_half_x_not_minus_x():
    """The 'wrong' identity from the earlier prompt draft, formalised
    so future readers don't reintroduce it."""
    x = 1.0
    correct_sum = x * (1.0 / math.tanh(x / 2.0))
    assert bernoulli(x) + bernoulli(-x) == pytest.approx(correct_sum, rel=1.0e-12)
    # And it is NOT -x.
    assert abs(bernoulli(x) + bernoulli(-x) - (-x)) > 0.5


# --------------------------------------------------------------------------- #
# Vectorisation tests                                                         #
# --------------------------------------------------------------------------- #

def test_bernoulli_array_matches_scalar_bitwise():
    xs = np.array([
        -50.0, -30.001, -29.999, -1.0, -1.0e-3, -1.0e-9,
        0.0,
        1.0e-9, 1.0e-3, 1.0, 29.999, 30.001, 50.0,
    ])
    expected = np.array([bernoulli(float(x)) for x in xs])
    got = bernoulli_array(xs)
    np.testing.assert_array_equal(got, expected)


# --------------------------------------------------------------------------- #
# Edge-flux sign and scaling guards (ADR 0012 test plan #2 and #3)            #
# --------------------------------------------------------------------------- #

# Test edges:
#
#   * SMALL_PECLET edge: dpsi = 0.005 V at room T gives dpsi_scaled = 0.193
#     and cell Peclet |dpsi_scaled|/2 = 0.097 << 1. SG and midpoint-Galerkin
#     agree at this Peclet because both reduce to centred-difference
#     diffusion plus a near-linear drift correction.
#   * USER_SPEC edge: dpsi = 0.1 V (the user's prompt value). dpsi_scaled
#     = 3.87 puts the cell at Peclet 1.93 which is firmly in the upwind
#     regime where SG and Galerkin diverge by design. We use this edge
#     for the two-direction sign-flip tests, which are robust at any
#     Peclet because they exercise antisymmetry of the SG flux under
#     (i <-> j, dpsi -> -dpsi).
#
# **Deviation from prompt spec**: the prompt asserted that Peclet at
# dpsi=0.1V was "small enough" for the Galerkin agreement test. It is
# not (Pe = 1.93). For the agreement guard we use dpsi = 0.005 V instead;
# the sign-flip and hole-symmetry tests retain dpsi = 0.1 V per the spec.
N_I = 1.0e16
N_J = 1.0e17
H_IJ = 1.0e-7  # 100 nm
MU_N = 1400.0e-4  # m^2/Vs (Si electron, room-T)
MU_P = 450.0e-4   # m^2/Vs (Si hole, room-T)
V_T = 0.025852  # k_B T / q at 300 K

# Small-Peclet edge for the agreement test.
DPSI_VOLTS_SMALL = 0.005
DPSI_SCALED_SMALL = DPSI_VOLTS_SMALL / V_T  # ~ 0.193, Peclet ~ 0.097

# User-spec edge for the sign-flip and hole-symmetry tests.
DPSI_VOLTS_USER = 0.1
DPSI_SCALED_USER = DPSI_VOLTS_USER / V_T  # ~ 3.87, Peclet ~ 1.93


def test_sg_edge_flux_n_agrees_with_midpoint_galerkin_at_small_peclet():
    """
    Independent reference: midpoint-Galerkin flux. Different derivation
    (centred difference) so this catches sign-convention errors and
    factor-of-2 errors. Slotboom 2-node closed form is NOT used because
    it derives from the same B(x).

    Cell Peclet at dpsi=0.005V: |dpsi_scaled|/2 = 0.097, well below 1.
    SG-vs-Galerkin disagreement scales as O(Pe^2) ~ 1%. Tolerance: 5%.
    Sign error: orders of magnitude. Factor-of-2 error: 2x.
    """
    sg = sg_edge_flux_n(N_I, N_J, DPSI_SCALED_SMALL, H_IJ, MU_N, V_T)
    ref = midpoint_galerkin_flux_n(N_I, N_J, DPSI_SCALED_SMALL, H_IJ, MU_N, V_T)
    rel = abs(sg - ref) / abs(ref)
    assert rel < 0.05, (
        f"SG vs midpoint-Galerkin disagreement {rel*100:.3f}% > 5%; "
        f"sg={sg:.6e}, ref={ref:.6e}, dpsi_scaled={DPSI_SCALED_SMALL:.4f}"
    )


def test_sg_edge_flux_n_two_direction_sign():
    """
    Swap (i, j): flux flips sign with same magnitude. The SG flux is
    naturally antisymmetric in the (i <-> j, dpsi -> -dpsi) swap. This
    holds at any Peclet number; we use the user-spec dpsi = 0.1 V.
    """
    forward = sg_edge_flux_n(N_I, N_J, DPSI_SCALED_USER, H_IJ, MU_N, V_T)
    reverse = sg_edge_flux_n(N_J, N_I, -DPSI_SCALED_USER, H_IJ, MU_N, V_T)
    assert reverse == pytest.approx(-forward, rel=1.0e-10, abs=0.0)


def test_sg_edge_flux_p_two_direction_sign():
    """Same antisymmetry for holes at the user-spec edge."""
    p_i, p_j = N_J, N_I  # hole densities mirror electrons here
    forward = sg_edge_flux_p(p_i, p_j, DPSI_SCALED_USER, H_IJ, MU_P, V_T)
    reverse = sg_edge_flux_p(p_j, p_i, -DPSI_SCALED_USER, H_IJ, MU_P, V_T)
    assert reverse == pytest.approx(-forward, rel=1.0e-10, abs=0.0)


def test_hole_flux_satisfies_device_equation_symmetry_guard_2b():
    """
    Guard 2b (the most important sign test).

    The device equations are symmetric under
    `(n <-> p, psi -> -psi, anode <-> cathode)`. Take an electron-flux
    configuration and mirror it under that transformation; the
    resulting hole-flux configuration must give `F_p = -F_n` to
    machine precision (1e-10 relative).

    Mirror map for our test edge:
        electrons : (n_i, n_j, dpsi)              -> F_n
        holes (mirror): (p_i = n_i, p_j = n_j, -dpsi) -> F_p
    Then F_p must equal -F_n exactly (using mu_p = mu_n here so the
    comparison is purely sign-driven).

    The hole sign is the place SG implementations most often go wrong.
    This guard is non-negotiable.
    """
    mu_shared = MU_N  # use the same mobility so the comparison is purely sign-driven
    f_n = sg_edge_flux_n(N_I, N_J, DPSI_SCALED_USER, H_IJ, mu_shared, V_T)
    f_p = sg_edge_flux_p(N_I, N_J, -DPSI_SCALED_USER, H_IJ, mu_shared, V_T)
    assert f_p == pytest.approx(-f_n, rel=1.0e-10, abs=0.0), (
        f"Hole-flux device-equation symmetry violated: F_n={f_n:.6e}, "
        f"F_p_mirror={f_p:.6e}, expected F_p = -F_n. Hole sign is wrong."
    )


def test_sg_edge_flux_n_array_matches_scalar():
    """Vectorised electron flux must agree with scalar evaluation
    point-by-point."""
    rng = np.random.default_rng(seed=0xE001)
    n_i = rng.uniform(1e14, 1e18, size=10)
    n_j = rng.uniform(1e14, 1e18, size=10)
    dpsi = rng.uniform(-5.0, 5.0, size=10)
    h = rng.uniform(0.5e-7, 5e-7, size=10)

    vec = sg_edge_flux_n_array(n_i, n_j, dpsi, h, MU_N, V_T)
    scalar = np.array([
        sg_edge_flux_n(float(n_i[k]), float(n_j[k]), float(dpsi[k]),
                       float(h[k]), MU_N, V_T)
        for k in range(10)
    ])
    np.testing.assert_allclose(vec, scalar, rtol=1.0e-14, atol=0.0)


def test_sg_edge_flux_p_array_matches_scalar():
    """Same for holes."""
    rng = np.random.default_rng(seed=0xE002)
    p_i = rng.uniform(1e14, 1e18, size=10)
    p_j = rng.uniform(1e14, 1e18, size=10)
    dpsi = rng.uniform(-5.0, 5.0, size=10)
    h = rng.uniform(0.5e-7, 5e-7, size=10)

    vec = sg_edge_flux_p_array(p_i, p_j, dpsi, h, MU_P, V_T)
    scalar = np.array([
        sg_edge_flux_p(float(p_i[k]), float(p_j[k]), float(dpsi[k]),
                       float(h[k]), MU_P, V_T)
        for k in range(10)
    ])
    np.testing.assert_allclose(vec, scalar, rtol=1.0e-14, atol=0.0)


# --------------------------------------------------------------------------- #
# High-Peclet sanity (regime where Galerkin fails and SG must hold up)        #
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("dpsi_scaled", [10.0, 50.0, -10.0, -50.0])
def test_sg_edge_flux_n_finite_at_high_peclet(dpsi_scaled):
    """The SG flux must be finite (no overflow / NaN) at large |dpsi|.
    This is the regime where the M13 Galerkin discretisation drove
    n into 1e24 m^-3 territory; SG should not."""
    f = sg_edge_flux_n(1.0e16, 1.0e17, dpsi_scaled, 1.0e-7, MU_N, V_T)
    assert math.isfinite(f), f"SG flux not finite at dpsi={dpsi_scaled}"
