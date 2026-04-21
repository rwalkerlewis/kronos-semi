"""Tests for semi.schema — no dolfinx required."""
import json
from pathlib import Path

import pytest

from semi import schema


@pytest.fixture
def minimal_cfg():
    return {
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


def test_doping_gaussian_schema():
    cfg = {
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
