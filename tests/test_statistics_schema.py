"""
Schema-side tests for the M16.4 Fermi-Dirac statistics dispatch.

The schema additive minor bump v2.3.0 -> v2.4.0 widens the
`physics.statistics` enum from `["boltzmann"]` to
`["boltzmann", "fermi_dirac"]`. Default stays `"boltzmann"`. v2.0.0,
v2.1.0, v2.2.0, and v2.3.0 inputs continue to validate (additive).

These tests live in their own module (mirroring
`tests/test_mobility_schema.py` and `tests/test_recombination.py`)
so the M16.x schema-side surface stays browseable per slice.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from semi import schema


@pytest.fixture
def minimal_cfg_v2():
    return {
        "schema_version": "2.0.0",
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


def test_supported_minor_is_at_least_4():
    """The engine must accept the v2.4.0 minor introduced by M16.4."""
    assert schema.SCHEMA_SUPPORTED_MINOR >= 4


def test_default_statistics_is_boltzmann(minimal_cfg_v2):
    """An input without physics.statistics defaults to boltzmann."""
    result = schema.validate(minimal_cfg_v2)
    assert result["physics"]["statistics"] == "boltzmann"


def test_explicit_boltzmann_statistics(minimal_cfg_v2):
    minimal_cfg_v2["physics"] = {"statistics": "boltzmann"}
    result = schema.validate(minimal_cfg_v2)
    assert result["physics"]["statistics"] == "boltzmann"


def test_fermi_dirac_statistics_validates(minimal_cfg_v2):
    """The fermi_dirac branch added in v2.4.0 must validate."""
    minimal_cfg_v2["schema_version"] = "2.4.0"
    minimal_cfg_v2["physics"] = {"statistics": "fermi_dirac"}
    result = schema.validate(minimal_cfg_v2)
    assert result["physics"]["statistics"] == "fermi_dirac"


def test_unknown_statistics_rejected(minimal_cfg_v2):
    minimal_cfg_v2["physics"] = {"statistics": "foo"}
    with pytest.raises(schema.SchemaError):
        schema.validate(minimal_cfg_v2)


def test_v2_3_0_input_still_validates(minimal_cfg_v2):
    """v2.3.0 inputs must continue to validate against the v2.4.0 schema."""
    minimal_cfg_v2["schema_version"] = "2.3.0"
    schema.validate(minimal_cfg_v2)


def test_v2_2_0_input_still_validates(minimal_cfg_v2):
    minimal_cfg_v2["schema_version"] = "2.2.0"
    schema.validate(minimal_cfg_v2)


def test_v2_1_0_input_still_validates(minimal_cfg_v2):
    minimal_cfg_v2["schema_version"] = "2.1.0"
    schema.validate(minimal_cfg_v2)


def test_v2_0_0_input_still_validates(minimal_cfg_v2):
    minimal_cfg_v2["schema_version"] = "2.0.0"
    schema.validate(minimal_cfg_v2)


def test_schema_v2_4_0_examples_present():
    """The schema_version examples list must advertise v2.4.0."""
    schema_v2 = schema.get_schema(2)
    examples = schema_v2["properties"]["schema_version"]["examples"]
    assert "2.4.0" in examples
    # Backwards-compat advertisement of older minors must remain.
    assert "2.3.0" in examples
    assert "2.0.0" in examples


def test_schema_statistics_enum_widened():
    """The physics.statistics enum carries both branches."""
    schema_v2 = schema.get_schema(2)
    enum = schema_v2["properties"]["physics"]["properties"]["statistics"]["enum"]
    assert "boltzmann" in enum
    assert "fermi_dirac" in enum


def test_every_benchmark_validates_unchanged():
    """Every shipped benchmark JSON must validate against v2.4.0 without edits."""
    bench_dir = Path(__file__).parent.parent / "benchmarks"
    json_paths = sorted(bench_dir.glob("*/*.json"))
    assert json_paths, "no benchmark JSONs found"
    for path in json_paths:
        with path.open() as f:
            cfg = json.load(f)
        # Only schema-validating files; some bench dirs may carry
        # mesh-spec JSONs without schema_version. Skip those.
        if "schema_version" not in cfg:
            continue
        schema.validate(cfg)
