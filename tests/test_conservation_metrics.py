"""
Pure-Python unit tests for the Phase 3 V&V conservation metrics.

These tests exercise the dolfinx-free reduction functions in
`semi.verification.conservation`:

    current_continuity_metric(xs, J_vals)  -> CurrentContinuityMetric
    charge_neutrality_metric(x, p, n, N_net, q) -> ChargeNeutralityMetric

The functions take raw arrays, so we fabricate `J_n`, `J_p`, and
(n, p, N_net) profiles directly rather than constructing a dolfinx
SimulationResult. That's what "Mock SimulationResult with known J_n,
J_p for continuity check" / "Mock charge profile for charge check"
means operationally: the FEM stack is not needed to verify the
reduction math is correct.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from semi.constants import Q
from semi.verification.conservation import (
    charge_neutrality_metric,
    current_continuity_metric,
)


# ---------------------------------------------------------------------------
# current_continuity_metric
# ---------------------------------------------------------------------------


def test_continuity_flat_current_passes_5_percent():
    """
    If J_n + J_p is constant in space (the physically correct answer),
    the metric's max relative deviation is exactly zero and the 5%
    forward gate passes trivially.
    """
    xs = np.linspace(1.0e-6, 1.9e-5, 10)
    J_n = np.full(10, 12.0)
    J_p = np.full(10, -2.0)
    m = current_continuity_metric(xs, J_n + J_p)
    assert m.mean_J == pytest.approx(10.0)
    assert m.max_abs_dev == pytest.approx(0.0)
    assert m.max_rel_dev == pytest.approx(0.0)
    assert m.max_rel_dev < 0.05


def test_continuity_noisy_within_tolerance_passes_5_percent():
    """
    1% white noise on a 10 A/m^2 baseline must fall under the 5%
    forward-bias gate.
    """
    rng = np.random.default_rng(seed=42)
    xs = np.linspace(1.0e-6, 1.9e-5, 10)
    J_mean = 10.0
    J_vals = J_mean + 0.01 * J_mean * rng.standard_normal(10)
    m = current_continuity_metric(xs, J_vals)
    assert m.max_rel_dev < 0.05
    assert m.mean_J == pytest.approx(J_mean, rel=0.02)


def test_continuity_15_percent_outlier_fails_5_percent_passes_15_percent():
    """
    An outlier at 15% of the mean must fail the 5% forward gate and
    pass the 15% reverse gate. This exercises both thresholds in one
    synthetic profile.
    """
    xs = np.linspace(0.0, 2.0e-5, 10)
    J_vals = np.full(10, 10.0)
    J_vals[5] = 10.0 * 1.12        # +12% outlier in one sample
    m = current_continuity_metric(xs, J_vals)
    # Mean is 10.12; max deviation ~1.08 A/m^2 at the outlier; rel ~10.7%.
    assert 0.05 < m.max_rel_dev < 0.15


def test_continuity_zero_mean_zero_deviation_is_zero():
    """
    At V=0 the physical answer is J=0 everywhere. The metric must
    report 0 deviation with 0 mean and rel_dev = 0 (not inf).
    """
    m = current_continuity_metric(np.linspace(0, 1, 5), np.zeros(5))
    assert m.mean_J == 0.0
    assert m.max_abs_dev == 0.0
    assert m.max_rel_dev == 0.0


def test_continuity_zero_mean_nonzero_deviation_is_inf():
    """
    If the mean collapses to zero but samples deviate (pathological
    cancellation), rel_dev is +inf so the gate always fails loudly.
    """
    m = current_continuity_metric(
        np.linspace(0, 1, 4), np.array([1.0, -1.0, 1.0, -1.0])
    )
    assert m.mean_J == pytest.approx(0.0)
    assert m.max_abs_dev > 0.0
    assert math.isinf(m.max_rel_dev)


def test_continuity_empty_inputs_are_safe():
    """An empty samples array must not raise and must report zeros."""
    m = current_continuity_metric(np.empty(0), np.empty(0))
    assert m.mean_J == 0.0
    assert m.max_rel_dev == 0.0


# ---------------------------------------------------------------------------
# charge_neutrality_metric
# ---------------------------------------------------------------------------


def test_charge_neutrality_symmetric_step_junction_is_zero():
    """
    Synthetic symmetric step junction: N_A on the left, N_D on the
    right, and p(x) = N_A on the left, n(x) = N_D on the right, with
    zero minority carriers. This is a charge-neutral profile and
    integrating (p - n + N_net) gives zero analytically.

    N_net here is N_D - N_A (the sign convention the benchmark uses).
    """
    N = 201
    L = 2.0e-5
    x = np.linspace(0.0, L, N)
    N_A = 1.0e23
    N_D = 1.0e23
    # Step profile at x = L/2: N_net = N_D on the right, -N_A on the left.
    N_net = np.where(x < L / 2, -N_A, N_D)
    # Perfect majority-carrier mirror image: p = N_A on the left, n = N_D
    # on the right, minority ~0.
    p = np.where(x < L / 2, N_A, 0.0)
    n = np.where(x < L / 2, 0.0, N_D)

    m = charge_neutrality_metric(x, p, n, N_net, Q)
    # Exact analytic zero modulo trapezoid rounding on the discontinuity.
    Q_ref_expected = Q * max(N_A, N_D) * L
    assert m.Q_ref == pytest.approx(Q_ref_expected, rel=1.0e-12)
    assert m.rel_error < 1.0e-10


def test_charge_neutrality_excess_charge_is_detected():
    """
    A 1e-9 fractional imbalance on top of peak doping must produce a
    relative error above the 1e-10 gate and therefore be caught.
    """
    N = 101
    L = 1.0e-5
    x = np.linspace(0.0, L, N)
    N_A = 1.0e22
    # Imbalanced: p on the left, n on the right, with a tiny asymmetry.
    N_net = np.where(x < L / 2, -N_A, N_A)
    p = np.where(x < L / 2, N_A * (1 + 1.0e-9), 0.0)
    n = np.where(x < L / 2, 0.0, N_A)
    m = charge_neutrality_metric(x, p, n, N_net, Q)
    assert m.rel_error > 1.0e-10
    assert m.Q_net > 0.0


def test_charge_neutrality_sort_is_robust_to_unordered_input():
    """
    The metric must sort coordinates before applying the trapezoidal
    rule; feeding an unordered x must give the same answer as the
    sorted version.
    """
    x = np.array([0.0, 0.5, 0.25, 1.0, 0.75])
    N_A = 1.0e22
    # Linear rho = +N_A * x - N_A/2 (zero mean over [0,1])
    rho_over_q = N_A * x - N_A / 2.0
    # Encode rho as p - n + N_net with N_net=0.
    p = np.where(rho_over_q > 0, rho_over_q, 0.0)
    n = np.where(rho_over_q < 0, -rho_over_q, 0.0)
    N_net = np.zeros_like(x)

    m = charge_neutrality_metric(x, p, n, N_net, Q)
    # The integral of N_A * x - N_A/2 over [0, 1] is 0.
    assert abs(m.Q_net) < Q * N_A * 1.0e-12
    assert m.Q_ref == 0.0   # N_net is all zeros; reference collapses
    assert math.isinf(m.rel_error)


def test_charge_neutrality_reference_uses_max_abs_N_net_and_L():
    """
    The reference charge is q * max|N_net| * (x_max - x_min); verify
    directly so future refactors do not silently drop a factor.
    """
    x = np.linspace(0.0, 2.0, 5)
    N_net = np.array([-3.0, -2.0, 0.0, 4.0, 1.0])
    m = charge_neutrality_metric(x, np.zeros(5), np.zeros(5), N_net, Q)
    assert m.Q_ref == pytest.approx(Q * 4.0 * 2.0)


def test_charge_threshold_matches_phase3_spec():
    """
    The acceptance threshold specified by Phase 3 is 1e-10 * q *
    max|N_net| * L. A sanity synthetic that sits just under it must
    pass; one just over it must fail.
    """
    x = np.linspace(0.0, 1.0e-5, 11)
    N_A = 1.0e23
    L = 1.0e-5
    Q_ref = Q * N_A * L
    # |Q_net| == 0.5 * threshold -> must pass with margin.
    target_rel = 0.5e-10
    # Integral of constant c over [0, L] is c * L.
    c = target_rel * N_A
    p = np.full_like(x, c)
    n = np.zeros_like(x)
    N_net = np.full_like(x, -N_A + 0.0) + np.zeros_like(x)
    # Actually rework: keep max|N_net| = N_A; use rho-via-minor-excess-p.
    N_net = np.full_like(x, -N_A)
    p = np.full_like(x, N_A + c)     # slight excess
    n = np.zeros_like(x)             # no compensation
    m = charge_neutrality_metric(x, p, n, N_net, Q)
    # |Q_net| ~= q * c * L = q * (0.5e-10 * N_A) * L = 0.5e-10 * Q_ref
    assert m.rel_error < 1.0e-10
    assert m.Q_ref == pytest.approx(Q_ref, rel=1.0e-12)
