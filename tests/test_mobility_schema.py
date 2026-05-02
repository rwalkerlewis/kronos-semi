"""
Tests for the M16.1 schema 2.1.0 mobility dispatch.

Coverage:
    - schemas/input.v2.json declares the new caughey_thomas enum value
      and the four new vsat_n / vsat_p / beta_n / beta_p parameters.
    - Every existing benchmark JSON validates unchanged against v2.1.0
      (additive minor; no benchmark migration required).
    - A caughey_thomas input with default vsat_* / beta_* parameters
      validates.
    - A caughey_thomas input with an explicit vsat_n override validates.
    - An unknown mobility model string is rejected with an enum error.
    - An extra mobility property (strict-v2 additionalProperties:false)
      is rejected.
    - A v2.0.0 input (no vsat_* / beta_*, model "constant") still
      validates against the v2.1.0 schema (additive bump invariant).

Acceptance test 3 in `docs/IMPROVEMENT_GUIDE.md` § M16.1 is covered by
the benchmark-validation loop.
"""
from __future__ import annotations

import glob
import json
from copy import deepcopy
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SCHEMAS_DIR = REPO_ROOT / "schemas"
V2_PATH = SCHEMAS_DIR / "input.v2.json"
BENCHMARKS_DIR = REPO_ROOT / "benchmarks"


def _v2_validator():
    import jsonschema

    schema = json.loads(V2_PATH.read_text())
    return jsonschema.Draft7Validator(schema)


def _base_v2_cfg() -> dict:
    return {
        "schema_version": "2.1.0",
        "name": "mobility_test",
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


def test_v2_schema_declares_caughey_thomas_enum():
    """The v2.1.0 schema must list caughey_thomas as a valid model value."""
    schema = json.loads(V2_PATH.read_text())
    mob = schema["properties"]["physics"]["properties"]["mobility"]
    enum = mob["properties"]["model"]["enum"]
    assert "constant" in enum
    assert "caughey_thomas" in enum
    for key in ("vsat_n", "vsat_p", "beta_n", "beta_p"):
        assert key in mob["properties"], f"v2.1.0 schema must declare {key}"


def test_every_existing_benchmark_validates_unchanged():
    """Acceptance test 3: every shipped benchmark JSON validates against
    v2.1.0 without modification (additive minor bump)."""
    v = _v2_validator()
    paths = sorted(glob.glob(str(BENCHMARKS_DIR / "*" / "*.json")))
    assert paths, "no benchmark JSONs found under benchmarks/*/*.json"
    for p in paths:
        cfg = json.loads(Path(p).read_text())
        errors = sorted(v.iter_errors(cfg), key=lambda e: list(e.path))
        assert not errors, (
            f"{p} fails v2.1.0 validation:\n"
            + "\n".join(
                f"  at {'.'.join(str(x) for x in e.absolute_path) or '<root>'}: {e.message}"
                for e in errors[:5]
            )
        )


def test_caughey_thomas_default_parameters_validate():
    v = _v2_validator()
    cfg = _base_v2_cfg()
    cfg["physics"] = {"mobility": {"model": "caughey_thomas"}}
    errors = list(v.iter_errors(cfg))
    assert not errors, errors


def test_caughey_thomas_explicit_vsat_validates():
    v = _v2_validator()
    cfg = _base_v2_cfg()
    cfg["physics"] = {
        "mobility": {
            "model": "caughey_thomas",
            "vsat_n": 1.2e7,
            "vsat_p": 8.5e6,
            "beta_n": 2.0,
            "beta_p": 1.0,
        }
    }
    errors = list(v.iter_errors(cfg))
    assert not errors, errors


def test_unknown_mobility_model_rejected():
    v = _v2_validator()
    cfg = _base_v2_cfg()
    cfg["physics"] = {"mobility": {"model": "foo"}}
    errors = list(v.iter_errors(cfg))
    assert errors, "unknown mobility model must be rejected by the enum"
    messages = " | ".join(e.message for e in errors)
    assert "foo" in messages or "enum" in messages.lower()


def test_extra_mobility_property_rejected():
    """Strict v2 additionalProperties:false applies to physics.mobility."""
    v = _v2_validator()
    cfg = _base_v2_cfg()
    cfg["physics"] = {"mobility": {"model": "constant", "foo": 1.0}}
    errors = list(v.iter_errors(cfg))
    assert errors, "extra mobility key must be rejected under strict v2"
    messages = " | ".join(e.message for e in errors)
    assert "foo" in messages


def test_v2_0_0_input_still_validates_against_v2_1_0_schema():
    """Additive minor bump invariant: a v2.0.0 input (model=constant
    without vsat_* / beta_*) continues to validate."""
    v = _v2_validator()
    cfg = _base_v2_cfg()
    cfg["schema_version"] = "2.0.0"
    cfg["physics"] = {"mobility": {"model": "constant", "mu_n": 1400.0, "mu_p": 450.0}}
    errors = list(v.iter_errors(cfg))
    assert not errors, errors


def test_validate_accepts_v2_1_0():
    """The engine validator dispatches v2.1.0 to input.v2.json without
    a deprecation warning."""
    import warnings

    from semi import schema as schema_mod

    cfg = _base_v2_cfg()
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        schema_mod.validate(deepcopy(cfg))
    deprecations = [w for w in captured if issubclass(w.category, DeprecationWarning)]
    assert not deprecations, (
        f"v2.1.0 inputs must not emit DeprecationWarning; got {deprecations}"
    )


def test_validate_fills_caughey_thomas_defaults():
    """_fill_defaults preserves the caughey_thomas model when set, and
    keeps the constant defaults intact otherwise."""
    from semi import schema as schema_mod

    cfg = _base_v2_cfg()
    cfg["physics"] = {"mobility": {"model": "caughey_thomas"}}
    out = schema_mod.validate(cfg)
    mob = out["physics"]["mobility"]
    assert mob["model"] == "caughey_thomas"
    assert mob["mu_n"] == pytest.approx(1400.0)
    assert mob["mu_p"] == pytest.approx(450.0)
