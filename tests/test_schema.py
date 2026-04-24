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
