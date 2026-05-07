"""
Tests for the M16.7 schema additions: contacts[].voltage_t for the
transient runner. Pure-Python; no dolfinx required.

Mirrors the M16.5 / M16.6 schema-test pattern: every existing
benchmark JSON keeps validating, the new field is optional, and the
loader rejects mutually exclusive or solver-mismatched configs at
validate time so users see the failure before the FEM path.
"""
from __future__ import annotations

import copy
import json
import math
from pathlib import Path

import pytest

from semi import schema


@pytest.fixture
def transient_cfg():
    return {
        "schema_version": "2.7.0",
        "name": "voltage_t_test",
        "dimension": 1,
        "mesh": {
            "source": "builtin",
            "extents": [[0.0, 1.0e-6]],
            "resolution": [10],
            "regions_by_box": [
                {"name": "silicon", "tag": 1, "bounds": [[0.0, 1.0e-6]]},
            ],
            "facets_by_plane": [
                {"name": "anode",   "tag": 1, "axis": 0, "value": 0.0},
                {"name": "cathode", "tag": 2, "axis": 0, "value": 1.0e-6},
            ],
        },
        "regions": {
            "silicon": {"material": "Si", "tag": 1, "role": "semiconductor"},
        },
        "doping": [
            {
                "region": "silicon",
                "profile": {"type": "uniform", "N_D": 1.0e17, "N_A": 0.0},
            }
        ],
        "contacts": [
            {"name": "anode",   "facet": "anode",   "type": "ohmic",
             "voltage": 0.0},
            {"name": "cathode", "facet": "cathode", "type": "ohmic",
             "voltage": 0.0},
        ],
        "solver": {"type": "transient", "t_end": 1.0e-9, "dt": 1.0e-10},
    }


def test_voltage_t_table_validates(transient_cfg):
    transient_cfg["contacts"][0]["voltage_t"] = {
        "type": "table",
        "times": [0.0, 1.0e-10, 2.0e-10],
        "values": [0.0, 0.3, 0.6],
    }
    result = schema.validate(transient_cfg)
    vt = result["contacts"][0]["voltage_t"]
    assert vt["type"] == "table"
    assert vt["times"] == [0.0, 1.0e-10, 2.0e-10]
    assert vt["values"] == [0.0, 0.3, 0.6]


def test_voltage_t_step_validates(transient_cfg):
    transient_cfg["contacts"][0]["voltage_t"] = {
        "type": "step", "t0": 5.0e-10, "v0": 0.0, "v1": 0.6,
    }
    result = schema.validate(transient_cfg)
    vt = result["contacts"][0]["voltage_t"]
    assert vt["type"] == "step"
    assert vt["t0"] == pytest.approx(5.0e-10)
    assert vt["v0"] == pytest.approx(0.0)
    assert vt["v1"] == pytest.approx(0.6)


def test_voltage_t_table_non_monotonic_rejected(transient_cfg):
    transient_cfg["contacts"][0]["voltage_t"] = {
        "type": "table",
        "times": [0.0, 2.0e-10, 1.0e-10],
        "values": [0.0, 0.3, 0.6],
    }
    with pytest.raises(schema.SchemaError, match="monotonically"):
        schema.validate(transient_cfg)


def test_voltage_t_table_length_mismatch_rejected(transient_cfg):
    transient_cfg["contacts"][0]["voltage_t"] = {
        "type": "table",
        "times": [0.0, 1.0e-10, 2.0e-10],
        "values": [0.0, 0.3],
    }
    with pytest.raises(schema.SchemaError, match="equal length"):
        schema.validate(transient_cfg)


def test_voltage_t_with_voltage_sweep_rejected(transient_cfg):
    transient_cfg["contacts"][0]["voltage_t"] = {
        "type": "step", "t0": 5.0e-10, "v0": 0.0, "v1": 0.6,
    }
    transient_cfg["contacts"][0]["voltage_sweep"] = {
        "start": 0.0, "stop": 0.6, "step": 0.1,
    }
    with pytest.raises(schema.SchemaError, match="mutually exclusive"):
        schema.validate(transient_cfg)


def test_voltage_t_on_bias_sweep_rejected(transient_cfg):
    transient_cfg["solver"] = {"type": "bias_sweep"}
    transient_cfg["contacts"][0]["voltage_t"] = {
        "type": "step", "t0": 5.0e-10, "v0": 0.0, "v1": 0.6,
    }
    with pytest.raises(schema.SchemaError, match="transient"):
        schema.validate(transient_cfg)


def test_voltage_t_on_equilibrium_rejected(transient_cfg):
    transient_cfg["solver"] = {"type": "equilibrium"}
    transient_cfg["contacts"][0]["voltage_t"] = {
        "type": "step", "t0": 5.0e-10, "v0": 0.0, "v1": 0.6,
    }
    with pytest.raises(schema.SchemaError, match="transient"):
        schema.validate(transient_cfg)


def test_voltage_t_on_ac_sweep_rejected(transient_cfg):
    transient_cfg["solver"] = {"type": "ac_sweep"}
    transient_cfg["contacts"][0]["voltage_t"] = {
        "type": "step", "t0": 5.0e-10, "v0": 0.0, "v1": 0.6,
    }
    with pytest.raises(schema.SchemaError, match="transient"):
        schema.validate(transient_cfg)


def test_voltage_t_unknown_variant_rejected(transient_cfg):
    transient_cfg["contacts"][0]["voltage_t"] = {
        "type": "sine", "amplitude": 0.05, "frequency": 1.0e6,
    }
    with pytest.raises(schema.SchemaError):
        schema.validate(transient_cfg)


def test_voltage_t_step_missing_required_field_rejected(transient_cfg):
    transient_cfg["contacts"][0]["voltage_t"] = {
        "type": "step", "t0": 5.0e-10, "v0": 0.0,
    }
    with pytest.raises(schema.SchemaError, match="v1"):
        schema.validate(transient_cfg)


def test_voltage_t_table_too_few_points_rejected(transient_cfg):
    transient_cfg["contacts"][0]["voltage_t"] = {
        "type": "table",
        "times": [0.0],
        "values": [0.3],
    }
    with pytest.raises(schema.SchemaError):
        schema.validate(transient_cfg)


def test_voltage_t_unknown_property_rejected(transient_cfg):
    transient_cfg["contacts"][0]["voltage_t"] = {
        "type": "step", "t0": 5.0e-10, "v0": 0.0, "v1": 0.6,
        "frequency": 1.0e6,
    }
    with pytest.raises(schema.SchemaError):
        schema.validate(transient_cfg)


def test_default_fill_leaves_voltage_t_unset(transient_cfg):
    """Contacts without voltage_t do not gain the field by default,
    preserving v0.22.0 byte-identity on every existing benchmark."""
    result = schema.validate(transient_cfg)
    assert "voltage_t" not in result["contacts"][0]
    assert "voltage_t" not in result["contacts"][1]


def test_voltage_t_two_contacts_independent(transient_cfg):
    transient_cfg["contacts"][0]["voltage_t"] = {
        "type": "step", "t0": 5.0e-10, "v0": 0.0, "v1": 0.6,
    }
    transient_cfg["contacts"][1]["voltage_t"] = {
        "type": "table",
        "times": [0.0, 1.0e-10, 2.0e-10],
        "values": [0.0, 0.0, 0.0],
    }
    result = schema.validate(transient_cfg)
    assert result["contacts"][0]["voltage_t"]["type"] == "step"
    assert result["contacts"][1]["voltage_t"]["type"] == "table"


def test_voltage_t_sampled_sinusoid_validates(transient_cfg):
    """The audit-test waveform shape (sinusoid sampled into a table)
    validates and round-trips through the loader."""
    n = 64
    f = 1.0e6
    times = [i / (n * f) for i in range(n)]
    values = [0.4 + 1.0e-3 * math.sin(2.0 * math.pi * f * t) for t in times]
    transient_cfg["contacts"][0]["voltage_t"] = {
        "type": "table", "times": times, "values": values,
    }
    result = schema.validate(transient_cfg)
    vt = result["contacts"][0]["voltage_t"]
    assert len(vt["times"]) == n
    assert len(vt["values"]) == n


def test_existing_benchmarks_validate_under_v2_7_0():
    """Every shipped benchmark JSON must continue to validate against
    the v2.7.0 strict schema (no benchmark uses voltage_t today, so
    every prior config is bit-identical to v0.22.0)."""
    benchmarks_dir = Path(__file__).parent.parent / "benchmarks"
    for json_path in sorted(benchmarks_dir.rglob("*.json")):
        if json_path.parent.name == "results":
            continue
        with json_path.open() as f:
            cfg = json.load(f)
        if not str(cfg.get("schema_version", "")).startswith("2."):
            continue
        schema.validate(copy.deepcopy(cfg))


def test_supported_minor_advanced_to_7():
    assert schema.SCHEMA_SUPPORTED_MINOR >= 7
