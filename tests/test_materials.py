"""Tests for semi.materials — no dolfinx required."""
import pytest

from semi import materials
from semi.constants import EPS0, cm2_to_m2, cm3_to_m3


def test_get_material_known():
    si = materials.get_material("Si")
    assert si.name == "Si"
    assert si.is_semiconductor()
    assert si.epsilon_r == pytest.approx(11.7)


def test_get_material_unknown():
    with pytest.raises(KeyError, match="Unknown material"):
        materials.get_material("Unobtainium")


def test_list_materials_includes_core():
    names = materials.list_materials()
    for required in ("Si", "Ge", "GaAs", "SiO2", "HfO2", "Si3N4"):
        assert required in names


def test_si_parameters_match_modern_tcad():
    """Values chosen to align with Sentaurus defaults. Don't let them drift."""
    si = materials.get_material("Si")
    assert si.Eg == pytest.approx(1.12)
    assert si.n_i == pytest.approx(cm3_to_m3(1.0e10))  # Altermatt consensus
    assert si.mu_n == pytest.approx(cm2_to_m2(1400.0))
    assert si.mu_p == pytest.approx(cm2_to_m2(450.0))


def test_insulators_have_zero_carrier_params():
    for name in ("SiO2", "HfO2", "Si3N4"):
        mat = materials.get_material(name)
        assert mat.is_insulator()
        assert not mat.is_semiconductor()
        assert mat.n_i == 0.0
        assert mat.mu_n == 0.0
        assert mat.Eg == 0.0


def test_epsilon_absolute():
    si = materials.get_material("Si")
    assert si.epsilon == pytest.approx(11.7 * EPS0)


def test_sio2_matches_nominal():
    sio2 = materials.get_material("SiO2")
    assert sio2.epsilon_r == pytest.approx(3.9)


def test_gaas_reasonable():
    gaas = materials.get_material("GaAs")
    assert 1.4 < gaas.Eg < 1.45
    # Mobility should be much higher than Si
    assert gaas.mu_n > materials.get_material("Si").mu_n
