"""
Tests for the M14.3 strict-mode schema (v2.0.0).

Coverage:
    - schemas/input.v2.json is a valid Draft-07 schema.
    - Every benchmark JSON in benchmarks/ validates as strict v2.
    - An intentional typo (`voltag` instead of `voltage`) in a contact
      entry is rejected by v2 with a message naming the offending field.
    - A v1 input still loads but emits a DeprecationWarning.
    - v3+ majors are rejected (sanity check on the supported-majors gate).
"""
from __future__ import annotations

import glob
import json
import warnings
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SCHEMAS_DIR = REPO_ROOT / "schemas"
V1_PATH = SCHEMAS_DIR / "input.v1.json"
V2_PATH = SCHEMAS_DIR / "input.v2.json"
BENCHMARKS_DIR = REPO_ROOT / "benchmarks"


def _v2_validator():
    import jsonschema

    schema = json.loads(V2_PATH.read_text())
    return jsonschema.Draft7Validator(schema)


def _base_v2_cfg() -> dict:
    return {
        "schema_version": "2.0.0",
        "name": "strict_test",
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


def test_v2_schema_file_exists_and_is_valid_draft7():
    import jsonschema

    assert V2_PATH.exists(), "schemas/input.v2.json must ship with M14.3"
    schema = json.loads(V2_PATH.read_text())
    jsonschema.Draft7Validator.check_schema(schema)
    # v2 must lock its schema_version to 2.x.y so v1 inputs cannot
    # accidentally be accepted under the strict gate.
    pat = schema["properties"]["schema_version"]["pattern"]
    assert pat.startswith("^2\\.")


def test_every_benchmark_validates_strict_v2():
    v = _v2_validator()
    paths = sorted(glob.glob(str(BENCHMARKS_DIR / "*" / "*.json")))
    assert paths, "no benchmark JSONs found under benchmarks/*/*.json"
    for p in paths:
        cfg = json.loads(Path(p).read_text())
        errors = sorted(v.iter_errors(cfg), key=lambda e: list(e.path))
        assert not errors, (
            f"{p} fails strict v2 validation:\n"
            + "\n".join(
                f"  at {'.'.join(str(x) for x in e.absolute_path) or '<root>'}: {e.message}"
                for e in errors[:5]
            )
        )


def test_strict_v2_rejects_contact_typo():
    """Acceptance test 4: `voltag` (typo for `voltage`) is rejected by v2.

    The v1 schema silently ignores the typo (additionalProperties is
    not constrained), which is the bug v2 fixes. The v2 error must
    name the offending field so a UI form-builder can highlight it.
    """
    v = _v2_validator()
    cfg = _base_v2_cfg()
    cfg["contacts"][0]["voltag"] = 0.0  # intentional typo

    errors = sorted(v.iter_errors(cfg), key=lambda e: list(e.path))
    assert errors, "v2 strict mode must reject extra properties"
    messages = " | ".join(e.message for e in errors)
    assert "voltag" in messages, (
        f"strict-v2 error should name the offending field 'voltag'; got: {messages!r}"
    )


def test_strict_v2_rejects_extra_top_level_key():
    """Sanity check: additionalProperties:false also bites at the root."""
    v = _v2_validator()
    cfg = _base_v2_cfg()
    cfg["nonsense_field"] = 42

    errors = sorted(v.iter_errors(cfg), key=lambda e: list(e.path))
    assert errors
    messages = " | ".join(e.message for e in errors)
    assert "nonsense_field" in messages


def test_v1_input_loads_with_deprecation_warning():
    """Acceptance test 4 part 2: v1 inputs still load, but log a
    DeprecationWarning so callers can plan their migration.
    """
    from semi import schema as schema_mod

    cfg = _base_v2_cfg()
    cfg["schema_version"] = "1.4.0"

    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        schema_mod.validate(cfg)
    deprecations = [w for w in captured if issubclass(w.category, DeprecationWarning)]
    assert deprecations, "v1 schema validation must emit a DeprecationWarning"
    assert any("v1" in str(w.message) for w in deprecations)


def test_v2_input_loads_without_warning():
    from semi import schema as schema_mod

    cfg = _base_v2_cfg()

    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        schema_mod.validate(cfg)
    deprecations = [w for w in captured if issubclass(w.category, DeprecationWarning)]
    assert not deprecations, (
        f"v2 (strict) inputs must not emit DeprecationWarning; got {deprecations}"
    )


def test_validate_rejects_v3_or_higher():
    from semi import schema as schema_mod

    cfg = _base_v2_cfg()
    cfg["schema_version"] = "3.0.0"
    with pytest.raises(schema_mod.SchemaError, match=r"major"):
        schema_mod.validate(cfg)


def test_get_schema_returns_correct_major():
    """The new `get_schema(major)` accessor returns each shipped schema."""
    from semi import schema as schema_mod

    s1 = schema_mod.get_schema(1)
    s2 = schema_mod.get_schema(2)
    assert s1["title"].endswith("v1")
    assert s2["title"].endswith("v2")
    # Strict-mode v2 must have additionalProperties:false at the root.
    assert s2.get("additionalProperties") is False
