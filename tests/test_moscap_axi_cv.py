"""
Schema and unit tests for the axisymmetric MOSCAP C-V feature.

Tests:
1. Schema validation of the benchmark config
2. Schema rejects invalid configs (bad dimension, bad contact type)
3. Mesh building (requires dolfinx)
4. C-V smoke test (requires dolfinx)
"""
from __future__ import annotations

import json
import pathlib
import pytest

from semi.schema import validate

BENCH_DIR = pathlib.Path(__file__).parent.parent / "benchmarks" / "moscap_axi"
MOSCAP_JSON = BENCH_DIR / "moscap_2d_axi.json"


def load_cfg():
    with open(MOSCAP_JSON) as f:
        return json.load(f)


# ------------------------------------------------------------------ #
# Schema tests (pure Python, no dolfinx)
# ------------------------------------------------------------------ #


def test_schema_accepts_axisymmetric_2d():
    """The benchmark config should pass schema validation."""
    cfg = load_cfg()
    validate(cfg)  # must not raise


def test_schema_accepts_mos_gate_contact():
    """Schema accepts 'mos_gate' contact type."""
    cfg = load_cfg()
    # Find the gate contact and verify it is mos_gate
    gate_contacts = [c for c in cfg["contacts"] if c["type"] == "mos_gate"]
    assert len(gate_contacts) >= 1, "Expected at least one mos_gate contact"
    validate(cfg)


def test_schema_accepts_cv_lf_solver():
    """Schema accepts solver.type == 'cv_lf'."""
    cfg = load_cfg()
    cfg["solver"] = {"type": "cv_lf"}
    validate(cfg)


def test_schema_accepts_cv_hf_solver():
    """Schema accepts solver.type == 'cv_hf'."""
    cfg = load_cfg()
    cfg["solver"] = {"type": "cv_hf"}
    validate(cfg)


def test_schema_rejects_invalid_dimension():
    """Schema should reject unknown dimension value."""
    from semi.schema import SchemaError
    cfg = load_cfg()
    cfg["dimension"] = "invalid_dim"
    with pytest.raises(SchemaError):
        validate(cfg)


def test_schema_rejects_invalid_contact_type():
    """Schema should reject unknown contact type."""
    from semi.schema import SchemaError
    cfg = load_cfg()
    cfg["contacts"][1]["type"] = "bad_type"
    with pytest.raises(SchemaError):
        validate(cfg)


def test_schema_accepts_integer_dimensions():
    """Schema still accepts legacy integer dimensions 1, 2, 3."""
    cfg = load_cfg()
    for dim in [1, 2, 3]:
        cfg_copy = dict(cfg)
        cfg_copy["dimension"] = dim
        validate(cfg_copy)


# ------------------------------------------------------------------ #
# Mesh test (requires dolfinx)
# ------------------------------------------------------------------ #


def test_axisymmetric_mesh_builds():
    """Build the axisymmetric mesh and check basic properties."""
    try:
        import dolfinx.mesh  # noqa: F401
        from semi.mesh import build_mesh
    except ImportError:
        pytest.skip("dolfinx not available")

    cfg = load_cfg()
    msh, cell_tags, facet_tags = build_mesh(cfg)

    assert msh is not None, "Mesh should not be None"
    assert cell_tags is not None, "Cell tags should not be None"
    assert facet_tags is not None, "Facet tags should not be None"

    # Check cell tags cover both regions
    import numpy as np
    unique_cell_tags = set(int(t) for t in np.unique(cell_tags.values))
    assert 1 in unique_cell_tags, "Silicon cells (tag=1) missing"
    assert 2 in unique_cell_tags, "Oxide cells (tag=2) missing"

    # Check facet tags
    unique_facet_tags = set(int(t) for t in np.unique(facet_tags.values))
    assert 12 in unique_facet_tags, "Body facet (tag=12) missing"
    assert 13 in unique_facet_tags, "Gate facet (tag=13) missing"


def test_bcs_resolve_with_axisymmetric_config():
    """BCs should resolve correctly for the axisymmetric config."""
    try:
        from semi.bcs import resolve_contacts
    except ImportError:
        pytest.skip("dolfinx not available")

    cfg = load_cfg()
    # resolve without mesh (pure-Python path)
    contacts = resolve_contacts(cfg, facet_tags=None, voltages={"gate": 0.5})

    names = {c.name for c in contacts}
    kinds = {c.kind for c in contacts}
    assert "body" in names
    assert "gate" in names
    assert "ohmic" in kinds
    assert "mos_gate" in kinds

    gate_bc = next(c for c in contacts if c.name == "gate")
    assert abs(gate_bc.V_applied - 0.5) < 1e-12
    assert gate_bc.work_function is not None
    assert abs(gate_bc.work_function - 4.05) < 1e-9
