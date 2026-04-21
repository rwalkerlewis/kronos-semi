"""Tests for semi.constants — no dolfinx required."""
import pytest

from semi import constants


def test_constants_match_2019_si():
    assert constants.Q == pytest.approx(1.602176634e-19, rel=1e-12)
    assert constants.KB == pytest.approx(1.380649e-23, rel=1e-12)


def test_thermal_voltage():
    # 300 K -> ~25.85 mV
    assert constants.thermal_voltage(300.0) == pytest.approx(0.025852, abs=1e-5)
    # Scales linearly
    assert constants.thermal_voltage(600.0) == pytest.approx(2 * 0.025852, abs=2e-5)


def test_cm3_roundtrip():
    assert constants.m3_to_cm3(constants.cm3_to_m3(1.0e17)) == pytest.approx(1.0e17)


def test_cm2_roundtrip():
    assert constants.m2_to_cm2(constants.cm2_to_m2(1400.0)) == pytest.approx(1400.0)


def test_cm3_to_m3_factor():
    # 1 cm^-3 = 1e6 m^-3
    assert constants.cm3_to_m3(1.0) == pytest.approx(1.0e6)


def test_cm2_to_m2_factor():
    # 1 cm^2/Vs = 1e-4 m^2/Vs
    assert constants.cm2_to_m2(1.0) == pytest.approx(1.0e-4)
