"""Tests for semi.schema — no dolfinx required."""
import json
from pathlib import Path

import pytest

from semi import schema


@pytest.fixture
def minimal_cfg():
    return {
        "schema_version": "1.0.0",
        "name": "test",
        "dimension": 1,
        "mesh": {
            "source": "builtin",
            "extents": [[0.0, 1.0e-6]],
            "resolution": [100],
            "facets_by_plane": [
                {"name": "L", "tag": 1, "axis": 0, "value": 0.0},
                {"name": "R", "tag": 2, "axis": 0, "value": 1.0e-6},
            ],
        },
        "regions": {"si": {"material": "Si", "tag": 1, "role": "semiconductor"}},
        "doping": [
            {"region": "si", "profile": {"type": "uniform", "N_D": 1e17, "N_A": 0.0}},
        ],
        "contacts": [
            {"name": "L", "facet": "L", "type": "ohmic", "voltage": 0.0},
            {"name": "R", "facet": "R", "type": "ohmic", "voltage": 0.0},
        ],
    }


def test_minimal_validates(minimal_cfg):
    result = schema.validate(minimal_cfg)
    # Defaults should now be present
    assert result["physics"]["temperature"] == 300.0
    assert result["solver"]["atol"] == pytest.approx(1.0e-10)


def test_missing_required(minimal_cfg):
    del minimal_cfg["dimension"]
    with pytest.raises(schema.SchemaError, match=r"dimension"):
        schema.validate(minimal_cfg)


def test_bad_dimension(minimal_cfg):
    minimal_cfg["dimension"] = 5
    with pytest.raises(schema.SchemaError):
        schema.validate(minimal_cfg)


def test_bad_mesh_source(minimal_cfg):
    minimal_cfg["mesh"]["source"] = "cookie"
    with pytest.raises(schema.SchemaError):
        schema.validate(minimal_cfg)


def test_defaults_do_not_override(minimal_cfg):
    """Explicit temperature should survive default-filling."""
    minimal_cfg["physics"] = {"temperature": 77.0}
    result = schema.validate(minimal_cfg)
    assert result["physics"]["temperature"] == 77.0


def test_load_benchmark_file():
    """The shipped benchmark JSON must validate."""
    path = Path(__file__).parent.parent / "benchmarks" / "pn_1d" / "pn_junction.json"
    cfg = schema.load(path)
    assert cfg["name"] == "pn_junction_1d"
    assert cfg["dimension"] == 1
    assert "_source_dir" in cfg


def test_dumps_strips_internal():
    cfg = {"name": "x", "dimension": 1, "_source_path": "/tmp/x.json", "_internal": "no"}
    s = schema.dumps(cfg)
    parsed = json.loads(s)
    assert "_source_path" not in parsed
    assert "_internal" not in parsed
    assert parsed["name"] == "x"


def test_recombination_E_t_default(minimal_cfg):
    result = schema.validate(minimal_cfg)
    rec = result["physics"]["recombination"]
    assert rec["E_t"] == 0.0
    assert rec["tau_n"] > 0.0
    assert rec["tau_p"] > 0.0


def test_voltage_sweep_accepted(minimal_cfg):
    minimal_cfg["contacts"][1]["voltage_sweep"] = {
        "start": 0.0, "stop": 0.6, "step": 0.05,
    }
    result = schema.validate(minimal_cfg)
    sweep = result["contacts"][1]["voltage_sweep"]
    assert sweep["start"] == 0.0
    assert sweep["stop"] == pytest.approx(0.6)
    assert sweep["step"] == pytest.approx(0.05)


def test_voltage_sweep_rejects_nonpositive_step(minimal_cfg):
    minimal_cfg["contacts"][1]["voltage_sweep"] = {
        "start": 0.0, "stop": 0.5, "step": 0.0,
    }
    with pytest.raises(schema.SchemaError):
        schema.validate(minimal_cfg)


def test_solver_type_drift_diffusion(minimal_cfg):
    minimal_cfg["solver"] = {"type": "drift_diffusion"}
    result = schema.validate(minimal_cfg)
    assert result["solver"]["type"] == "drift_diffusion"


def test_solver_type_bias_sweep(minimal_cfg):
    minimal_cfg["solver"] = {"type": "bias_sweep"}
    result = schema.validate(minimal_cfg)
    assert result["solver"]["type"] == "bias_sweep"


def test_continuation_defaults(minimal_cfg):
    result = schema.validate(minimal_cfg)
    cont = result["solver"]["continuation"]
    assert cont["max_halvings"] >= 0
    assert cont["min_step"] > 0.0
    assert cont["easy_iter_threshold"] >= 1
    assert cont["grow_factor"] > 1.0


def test_continuation_adaptive_overrides(minimal_cfg):
    minimal_cfg["solver"] = {
        "type": "bias_sweep",
        "continuation": {
            "max_step": 0.2,
            "easy_iter_threshold": 3,
            "grow_factor": 2.0,
        },
    }
    result = schema.validate(minimal_cfg)
    cont = result["solver"]["continuation"]
    assert cont["max_step"] == pytest.approx(0.2)
    assert cont["easy_iter_threshold"] == 3
    assert cont["grow_factor"] == pytest.approx(2.0)


def test_continuation_rejects_bad_grow_factor(minimal_cfg):
    minimal_cfg["solver"] = {
        "type": "bias_sweep",
        "continuation": {"grow_factor": 1.0},
    }
    with pytest.raises(schema.SchemaError):
        schema.validate(minimal_cfg)


def test_continuation_rejects_nonpositive_max_step(minimal_cfg):
    minimal_cfg["solver"] = {
        "type": "bias_sweep",
        "continuation": {"max_step": 0.0},
    }
    with pytest.raises(schema.SchemaError):
        schema.validate(minimal_cfg)


def test_doping_gaussian_schema():
    cfg = {
        "schema_version": "1.0.0",
        "name": "g", "dimension": 2,
        "mesh": {"source": "builtin", "extents": [[0, 1], [0, 1]], "resolution": [10, 10]},
        "regions": {"a": {"material": "Si"}},
        "doping": [{
            "region": "a",
            "profile": {
                "type": "gaussian",
                "center": [0.5, 0.5], "sigma": [0.1, 0.1],
                "peak": 1e18, "dopant": "donor",
            },
        }],
        "contacts": [{"name": "c", "facet": 1, "type": "ohmic"}],
    }
    schema.validate(cfg)  # should not raise


# --- M15: axisymmetric / gate-extras / cv_sweep ----------------------------

def _moscap_axi_cfg():
    """Minimal axisymmetric MOSCAP config for cross-field invariant tests."""
    return {
        "schema_version": "1.3.0",
        "name": "moscap_axi",
        "dimension": 2,
        "coordinate_system": "axisymmetric",
        "mesh": {
            "source": "builtin",
            "extents": [[0.0, 15.0e-6], [0.0, 2.005e-6]],
            "resolution": [40, 40],
            "regions_by_box": [
                {"name": "silicon", "tag": 1,
                 "bounds": [[0.0, 15.0e-6], [0.0, 2.0e-6]]},
                {"name": "oxide", "tag": 2,
                 "bounds": [[0.0, 15.0e-6], [2.0e-6, 2.005e-6]]},
            ],
            "facets_by_plane": [
                {"name": "body",  "tag": 1, "axis": 1, "value": 0.0},
                {"name": "gate",  "tag": 2, "axis": 1, "value": 2.005e-6},
                {"name": "outer", "tag": 3, "axis": 0, "value": 15.0e-6},
            ],
        },
        "regions": {
            "silicon": {"material": "Si", "tag": 1, "role": "semiconductor"},
            "oxide":   {"material": "SiO2", "tag": 2, "role": "insulator"},
        },
        "doping": [{
            "region": "silicon",
            "profile": {"type": "uniform", "N_D": 0.0, "N_A": 1.0e17},
        }],
        "contacts": [
            {"name": "body", "facet": "body", "type": "ohmic", "voltage": 0.0},
            {"name": "gate", "facet": "gate", "type": "gate",
             "voltage": 0.0, "oxide_thickness": 5.0e-9,
             "gate_workfunction_eV": 4.05, "oxide_region": "oxide"},
        ],
        "cv_sweep": {
            "Vg_min": -2.0, "Vg_max": 2.0, "n_points": 161,
            "frequency_mode": "both",
        },
    }


def test_axisymmetric_validates():
    cfg = schema.validate(_moscap_axi_cfg())
    assert cfg["coordinate_system"] == "axisymmetric"
    assert cfg["cv_sweep"]["frequency_mode"] == "both"


def test_coordinate_system_default_is_cartesian(minimal_cfg):
    out = schema.validate(minimal_cfg)
    assert out["coordinate_system"] == "cartesian"


def test_axisymmetric_rejects_dim_not_2():
    cfg = _moscap_axi_cfg()
    cfg["dimension"] = 1
    cfg["mesh"]["extents"] = [[0.0, 1.0e-6]]
    cfg["mesh"]["resolution"] = [40]
    with pytest.raises(schema.SchemaError, match=r"dimension == 2"):
        schema.validate(cfg)


def test_axisymmetric_rejects_nonzero_r_origin():
    cfg = _moscap_axi_cfg()
    cfg["mesh"]["extents"][0][0] = 1.0e-7  # r must start at exactly 0
    with pytest.raises(schema.SchemaError, match=r"extents\[0\]\[0\] == 0\.0"):
        schema.validate(cfg)


def test_gate_contact_requires_oxide_thickness():
    cfg = _moscap_axi_cfg()
    del cfg["contacts"][1]["oxide_thickness"]
    with pytest.raises(schema.SchemaError, match=r"oxide_thickness"):
        schema.validate(cfg)


def test_gate_workfunction_eV_optional():
    cfg = _moscap_axi_cfg()
    del cfg["contacts"][1]["gate_workfunction_eV"]
    schema.validate(cfg)  # should not raise


def test_cv_sweep_requires_exactly_one_gate():
    cfg = _moscap_axi_cfg()
    cfg["contacts"][1]["type"] = "ohmic"  # remove the only gate
    del cfg["contacts"][1]["oxide_thickness"]
    del cfg["contacts"][1]["gate_workfunction_eV"]
    with pytest.raises(schema.SchemaError, match=r"exactly one contact"):
        schema.validate(cfg)


def test_cv_sweep_rejects_inverted_range():
    cfg = _moscap_axi_cfg()
    cfg["cv_sweep"]["Vg_min"] = 2.0
    cfg["cv_sweep"]["Vg_max"] = -2.0
    with pytest.raises(schema.SchemaError, match=r"Vg_min.*Vg_max"):
        schema.validate(cfg)


def test_cv_sweep_rejects_too_few_points():
    cfg = _moscap_axi_cfg()
    cfg["cv_sweep"]["n_points"] = 3
    with pytest.raises(schema.SchemaError):
        schema.validate(cfg)


def test_cv_sweep_frequency_mode_default_is_both():
    cfg = _moscap_axi_cfg()
    del cfg["cv_sweep"]["frequency_mode"]
    out = schema.validate(cfg)
    assert out["cv_sweep"]["frequency_mode"] == "both"


def test_cv_sweep_rejects_bad_frequency_mode():
    cfg = _moscap_axi_cfg()
    cfg["cv_sweep"]["frequency_mode"] = "ultra"
    with pytest.raises(schema.SchemaError):
        schema.validate(cfg)
