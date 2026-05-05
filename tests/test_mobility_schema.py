"""
Tests for the M16.1 schema 2.1.0 and M16.2 schema 2.2.0 mobility
dispatch.

M16.1 coverage:
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

M16.2 coverage:
    - schemas/input.v2.json declares the new lombardi enum value, the
      bulk_model selector, the lombardi sub-object (B_*, C_*,
      lambda_*, delta_*), and the interface_facet_tag (nullable).
    - Every existing benchmark JSON validates unchanged against v2.2.0.
    - A lombardi input with a non-null interface_facet_tag and default
      lombardi parameters validates.
    - A lombardi input with a null interface_facet_tag is rejected by
      the schema loader (conditional check in semi/schema.py, not
      pure JSON Schema).
    - An extra property under physics.mobility.lombardi is rejected
      (additionalProperties:false).
    - A bulk_model value not in the enum is rejected.
    - A v2.0.0 / v2.1.0 input continues to validate against v2.2.0.

Acceptance test 3 in `docs/IMPROVEMENT_GUIDE.md` § M16.1 / § M16.2 is
covered by the benchmark-validation loop.
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


# ---------------------------------------------------------------------------
# M16.2 schema 2.2.0 lombardi dispatch
# ---------------------------------------------------------------------------


def test_v2_schema_declares_lombardi_enum():
    """The v2.2.0 schema must list lombardi in the model enum and
    declare the bulk_model selector, the lombardi sub-object, and the
    interface_facet_tag."""
    schema = json.loads(V2_PATH.read_text())
    mob = schema["properties"]["physics"]["properties"]["mobility"]
    enum = mob["properties"]["model"]["enum"]
    assert enum == ["constant", "caughey_thomas", "lombardi"]
    bulk = mob["properties"]["bulk_model"]
    assert bulk["enum"] == ["constant", "caughey_thomas"]
    assert bulk["default"] == "constant"
    tag = mob["properties"]["interface_facet_tag"]
    assert tag["type"] == ["integer", "null"]
    assert tag["default"] is None
    lombardi = mob["properties"]["lombardi"]
    assert lombardi["additionalProperties"] is False
    for key in (
        "B_n",
        "B_p",
        "C_n",
        "C_p",
        "lambda_n",
        "lambda_p",
        "delta_n",
        "delta_p",
    ):
        assert key in lombardi["properties"]


def test_every_existing_benchmark_validates_against_v2_2_0():
    """Acceptance: every shipped benchmark JSON validates against
    v2.2.0 without modification (additive minor bump)."""
    v = _v2_validator()
    paths = sorted(glob.glob(str(BENCHMARKS_DIR / "*" / "*.json")))
    assert paths
    for p in paths:
        cfg = json.loads(Path(p).read_text())
        errors = sorted(v.iter_errors(cfg), key=lambda e: list(e.path))
        assert not errors, (
            f"{p} fails v2.2.0 validation:\n"
            + "\n".join(
                f"  at {'.'.join(str(x) for x in e.absolute_path) or '<root>'}: {e.message}"
                for e in errors[:5]
            )
        )


def _lombardi_v2_cfg(*, with_facet_tag: bool = True) -> dict:
    cfg = _base_v2_cfg()
    cfg["schema_version"] = "2.2.0"
    mob: dict = {"model": "lombardi", "bulk_model": "caughey_thomas"}
    if with_facet_tag:
        mob["interface_facet_tag"] = 1
    else:
        mob["interface_facet_tag"] = None
    cfg["physics"] = {"mobility": mob}
    return cfg


def test_lombardi_with_facet_tag_validates_in_pure_json_schema():
    """JSON Schema layer accepts a lombardi block with an integer
    facet tag; the conditional check on null is enforced by the
    loader in semi/schema.py, not by pure JSON Schema."""
    v = _v2_validator()
    cfg = _lombardi_v2_cfg(with_facet_tag=True)
    errors = list(v.iter_errors(cfg))
    assert not errors, errors


def test_lombardi_null_facet_tag_passes_pure_json_schema():
    """The pure JSON Schema accepts null because the conditional
    requirement lives in the loader (project convention: no
    if/then/else in the JSON Schema; the loader carries the rule)."""
    v = _v2_validator()
    cfg = _lombardi_v2_cfg(with_facet_tag=False)
    errors = list(v.iter_errors(cfg))
    assert not errors, errors


def test_lombardi_null_facet_tag_rejected_by_loader():
    """semi.schema.validate must raise SchemaError when
    physics.mobility.model='lombardi' and interface_facet_tag is None."""
    from semi import schema as schema_mod

    cfg = _lombardi_v2_cfg(with_facet_tag=False)
    with pytest.raises(schema_mod.SchemaError, match="interface_facet_tag"):
        schema_mod.validate(deepcopy(cfg))


def test_lombardi_extra_property_rejected():
    """Strict v2 additionalProperties:false applies to the lombardi
    sub-object."""
    v = _v2_validator()
    cfg = _lombardi_v2_cfg(with_facet_tag=True)
    cfg["physics"]["mobility"]["lombardi"] = {"foo": 1.0}
    errors = list(v.iter_errors(cfg))
    assert errors
    messages = " | ".join(e.message for e in errors)
    assert "foo" in messages


def test_lombardi_unknown_bulk_model_rejected():
    v = _v2_validator()
    cfg = _lombardi_v2_cfg(with_facet_tag=True)
    cfg["physics"]["mobility"]["bulk_model"] = "lombardi"
    errors = list(v.iter_errors(cfg))
    assert errors, "bulk_model='lombardi' must be rejected by enum"
    messages = " | ".join(e.message for e in errors).lower()
    assert "enum" in messages or "lombardi" in messages


def test_v2_1_0_input_still_validates_against_v2_2_0_schema():
    """Additive minor bump invariant: a v2.1.0 caughey_thomas input
    (no bulk_model / lombardi sub-object) continues to validate
    under v2.2.0."""
    v = _v2_validator()
    cfg = _base_v2_cfg()
    cfg["schema_version"] = "2.1.0"
    cfg["physics"] = {"mobility": {"model": "caughey_thomas"}}
    errors = list(v.iter_errors(cfg))
    assert not errors, errors


def test_validate_fills_lombardi_defaults():
    """_fill_defaults populates the lombardi sub-object with the
    Lombardi 1988 / Sentaurus defaults when model='lombardi' and
    leaves other branches untouched."""
    from semi import schema as schema_mod

    cfg = _lombardi_v2_cfg(with_facet_tag=True)
    out = schema_mod.validate(deepcopy(cfg))
    mob = out["physics"]["mobility"]
    assert mob["model"] == "lombardi"
    assert mob["bulk_model"] == "caughey_thomas"
    lombardi = mob["lombardi"]
    assert lombardi["B_n"] == pytest.approx(4.75e7)
    assert lombardi["B_p"] == pytest.approx(9.93e6)
    assert lombardi["C_n"] == pytest.approx(1.74e5)
    assert lombardi["C_p"] == pytest.approx(8.84e5)
    assert lombardi["lambda_n"] == pytest.approx(0.125)
    assert lombardi["lambda_p"] == pytest.approx(0.0317)
    assert lombardi["delta_n"] == pytest.approx(5.82e14)
    assert lombardi["delta_p"] == pytest.approx(2.0546e14)


def test_validate_does_not_inject_lombardi_block_for_other_models():
    """When model != 'lombardi', _fill_defaults must not inject a
    lombardi sub-object (keeps constant / caughey_thomas defaults
    bit-equivalent to v0.17.0)."""
    from semi import schema as schema_mod

    cfg = _base_v2_cfg()
    cfg["schema_version"] = "2.2.0"
    cfg["physics"] = {"mobility": {"model": "caughey_thomas"}}
    out = schema_mod.validate(deepcopy(cfg))
    mob = out["physics"]["mobility"]
    assert "lombardi" not in mob
    assert "bulk_model" not in mob

    cfg = _base_v2_cfg()
    cfg["schema_version"] = "2.2.0"
    out = schema_mod.validate(deepcopy(cfg))
    mob = out["physics"]["mobility"]
    assert mob["model"] == "constant"
    assert "lombardi" not in mob
    assert "bulk_model" not in mob
