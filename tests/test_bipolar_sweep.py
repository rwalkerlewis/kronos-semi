"""
Pure-Python tests for the bipolar-sweep leg computation.

The M7 resistor benchmark introduced symmetric sweeps that cross
zero. `compute_bipolar_legs` converts such a sweep into the two
single-direction endpoints the `AdaptiveStepController` can walk
(V = 0 -> most-negative -> most-positive). Unipolar sweeps (the
pn-junction and MOS benchmarks) return an empty list so the original
single-endpoint ramp path runs unchanged.

This module lives in the pure-Python tier: no dolfinx, no mesh.
"""
from __future__ import annotations

import pytest


def test_bipolar_legs_populated_when_sweep_crosses_zero():
    from semi.runners.bias_sweep import compute_bipolar_legs

    v_sweep_list = [0.0, -0.010, -0.005, 0.005, 0.010]
    legs = compute_bipolar_legs(v_sweep_list)

    assert legs, "expected bipolar_legs to be populated for a sign-spanning sweep"
    assert len(legs) == 2
    assert legs[0] == pytest.approx(-0.010)
    assert legs[1] == pytest.approx(0.010)


def test_bipolar_legs_empty_for_unipolar_positive_sweep():
    from semi.runners.bias_sweep import compute_bipolar_legs

    assert compute_bipolar_legs([0.0, 0.1, 0.2, 0.3, 0.6]) == []


def test_bipolar_legs_empty_for_unipolar_negative_sweep():
    from semi.runners.bias_sweep import compute_bipolar_legs

    assert compute_bipolar_legs([0.0, -0.5, -1.0, -2.0]) == []


def test_bipolar_legs_ignores_numerical_zero():
    from semi.runners.bias_sweep import compute_bipolar_legs

    assert compute_bipolar_legs([0.0, 1.0e-15, 0.1]) == []
    assert compute_bipolar_legs([0.0, -1.0e-15, -0.1]) == []
