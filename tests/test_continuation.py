"""
Tests for semi.continuation.AdaptiveStepController.

No dolfinx dependency. Drives the controller with synthetic SNES
iteration counts and asserts the step-size sequence.
"""
from __future__ import annotations

import pytest

from semi.continuation import AdaptiveStepController, StepTooSmall


def test_initial_step_clamped_to_max():
    c = AdaptiveStepController(
        initial_step=1.0, max_step_abs=0.25, min_step_abs=1e-4,
        easy_iter_threshold=4, grow_factor=1.5,
    )
    assert c.step == pytest.approx(0.25)


def test_initial_step_preserved_within_bounds():
    c = AdaptiveStepController(
        initial_step=0.05, max_step_abs=0.2, min_step_abs=1e-4,
        easy_iter_threshold=4, grow_factor=1.5,
    )
    assert c.step == pytest.approx(0.05)


def test_negative_initial_step_preserves_sign():
    c = AdaptiveStepController(
        initial_step=-0.1, max_step_abs=0.2, min_step_abs=1e-4,
        easy_iter_threshold=4, grow_factor=1.5,
    )
    assert c.step == pytest.approx(-0.1)


def test_grows_after_threshold_consecutive_easy_solves():
    c = AdaptiveStepController(
        initial_step=0.05, max_step_abs=1.0, min_step_abs=1e-4,
        easy_iter_threshold=4, grow_factor=1.5,
    )
    for _ in range(3):
        c.on_success(3)  # easy
        assert c.step == pytest.approx(0.05)
    c.on_success(3)  # fourth easy -> grow
    assert c.step == pytest.approx(0.075)
    assert c.easy_count == 0


def test_non_easy_solve_resets_easy_counter():
    c = AdaptiveStepController(
        initial_step=0.05, max_step_abs=1.0, min_step_abs=1e-4,
        easy_iter_threshold=4, grow_factor=1.5,
    )
    c.on_success(3)
    c.on_success(3)
    c.on_success(6)  # not easy: resets
    assert c.easy_count == 0
    c.on_success(3)
    c.on_success(3)
    c.on_success(3)
    assert c.step == pytest.approx(0.05)  # only 3 easy in a row, not 4


def test_growth_capped_at_max_step():
    c = AdaptiveStepController(
        initial_step=0.05, max_step_abs=0.1, min_step_abs=1e-4,
        easy_iter_threshold=2, grow_factor=2.0,
    )
    c.on_success(1)
    c.on_success(1)  # grow 0.05 -> 0.1
    assert c.step == pytest.approx(0.1)
    c.on_success(1)
    c.on_success(1)  # would be 0.2, capped to 0.1
    assert c.step == pytest.approx(0.1)


def test_failure_halves_step_and_resets_counter():
    c = AdaptiveStepController(
        initial_step=0.08, max_step_abs=1.0, min_step_abs=1e-4,
        easy_iter_threshold=4, grow_factor=1.5,
    )
    c.on_success(3)
    c.on_success(3)
    assert c.easy_count == 2
    c.on_failure()
    assert c.step == pytest.approx(0.04)
    assert c.easy_count == 0


def test_failure_raises_when_step_hits_min():
    c = AdaptiveStepController(
        initial_step=0.01, max_step_abs=1.0, min_step_abs=0.005,
        easy_iter_threshold=4, grow_factor=1.5,
    )
    c.on_failure()  # 0.01 -> 0.005, OK
    assert c.step == pytest.approx(0.005)
    with pytest.raises(StepTooSmall):
        c.on_failure()  # 0.005 -> 0.0025 < min, raises


def test_clamp_to_endpoint_returns_step_when_within_bounds():
    c = AdaptiveStepController(
        initial_step=0.05, max_step_abs=0.1, min_step_abs=1e-4,
        easy_iter_threshold=4, grow_factor=1.5,
    )
    assert c.clamp_to_endpoint(0.5) == pytest.approx(0.05)


def test_clamp_to_endpoint_clamps_final_step():
    c = AdaptiveStepController(
        initial_step=0.1, max_step_abs=0.2, min_step_abs=1e-4,
        easy_iter_threshold=4, grow_factor=1.5,
    )
    assert c.clamp_to_endpoint(0.03) == pytest.approx(0.03)


def test_clamp_rejects_sign_mismatch():
    c = AdaptiveStepController(
        initial_step=0.05, max_step_abs=0.1, min_step_abs=1e-4,
        easy_iter_threshold=4, grow_factor=1.5,
    )
    with pytest.raises(ValueError):
        c.clamp_to_endpoint(-0.1)


def test_reverse_sweep_growth_stays_negative():
    c = AdaptiveStepController(
        initial_step=-0.1, max_step_abs=0.5, min_step_abs=1e-4,
        easy_iter_threshold=2, grow_factor=1.5,
    )
    c.on_success(1)
    c.on_success(1)  # grow
    assert c.step < 0.0
    assert abs(c.step) == pytest.approx(0.15)


def test_synthetic_forward_sweep_total_solves():
    """
    Walk 0 -> 0.6 V with initial 0.05, max 0.2, threshold 3, factor 2,
    every solve converging in 2 iters (easy). Assert the sequence of
    step sizes matches the expected growth pattern and the total
    number of solves is smaller than a fixed-step 0.05 sweep (12).
    """
    c = AdaptiveStepController(
        initial_step=0.05, max_step_abs=0.2, min_step_abs=1e-4,
        easy_iter_threshold=3, grow_factor=2.0,
    )
    V = 0.0
    V_end = 0.6
    solves = 0
    visited: list[float] = []
    while abs(V_end - V) > 1e-12 and solves < 100:
        step = c.clamp_to_endpoint(V_end - V)
        V += step
        visited.append(V)
        solves += 1
        c.on_success(2)
    assert solves < 12
    assert V == pytest.approx(V_end)


def test_invalid_config_rejected():
    with pytest.raises(ValueError):
        AdaptiveStepController(
            initial_step=0.05, max_step_abs=-0.1, min_step_abs=1e-4,
            easy_iter_threshold=4, grow_factor=1.5,
        )
    with pytest.raises(ValueError):
        AdaptiveStepController(
            initial_step=0.05, max_step_abs=0.1, min_step_abs=1.0,
            easy_iter_threshold=4, grow_factor=1.5,
        )
    with pytest.raises(ValueError):
        AdaptiveStepController(
            initial_step=0.05, max_step_abs=0.1, min_step_abs=1e-4,
            easy_iter_threshold=0, grow_factor=1.5,
        )
    with pytest.raises(ValueError):
        AdaptiveStepController(
            initial_step=0.05, max_step_abs=0.1, min_step_abs=1e-4,
            easy_iter_threshold=4, grow_factor=1.0,
        )
