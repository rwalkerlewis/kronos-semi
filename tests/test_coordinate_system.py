"""
Schema-level tests for the axisymmetric MOSCAP coordinate_system field.

These tests live alongside the existing pure-Python schema tests and
do not require dolfinx; they verify the cross-field validation rules
introduced in schema 1.3.0:

  - axisymmetric requires dimension == 2,
  - axisymmetric requires r_min >= 0 in builtin meshes,
  - no Dirichlet (ohmic, gate) contact may pin the symmetry axis r = 0,
  - cartesian (default) is unchanged.
"""
import copy

import pytest

from semi.schema import SchemaError, validate


def _base_axisym_cfg() -> dict:
    return {
        "schema_version": "1.3.0",
        "name": "moscap_axisym_test",
        "dimension": 2,
        "coordinate_system": "axisymmetric",
        "mesh": {
            "source": "builtin",
            "extents": [[0.0, 200.0e-6], [-5.0e-6, 1.0e-8]],
            "resolution": [40, 60],
            "regions_by_box": [
                {"name": "silicon", "tag": 1, "bounds": [[0.0, 200.0e-6], [-5.0e-6, 0.0]]},
                {"name": "oxide", "tag": 2, "bounds": [[0.0, 50.0e-6], [0.0, 1.0e-8]]},
            ],
            "facets_by_plane": [
                {"name": "axis", "tag": 1, "axis": 0, "value": 0.0},
                {"name": "outer", "tag": 2, "axis": 0, "value": 200.0e-6},
                {"name": "body", "tag": 3, "axis": 1, "value": -5.0e-6},
                {"name": "gate", "tag": 4, "axis": 1, "value": 1.0e-8},
            ],
        },
        "regions": {
            "silicon": {"material": "Si", "tag": 1, "role": "semiconductor"},
            "oxide": {"material": "SiO2", "tag": 2, "role": "insulator"},
        },
        "doping": [
            {"region": "silicon", "profile": {"type": "uniform", "N_D": 0.0, "N_A": 5.0e16}},
        ],
        "contacts": [
            {"name": "body", "facet": "body", "type": "ohmic", "voltage": 0.0},
            {"name": "gate", "facet": "gate", "type": "gate", "voltage": 0.0},
        ],
    }


def test_axisymmetric_valid_minimal():
    cfg = _base_axisym_cfg()
    out = validate(cfg)
    assert out["coordinate_system"] == "axisymmetric"


def test_default_coordinate_system_cartesian():
    cfg = _base_axisym_cfg()
    cfg.pop("coordinate_system")
    out = validate(cfg)
    assert out["coordinate_system"] == "cartesian"


def test_axisymmetric_rejects_dimension_3():
    cfg = _base_axisym_cfg()
    cfg["dimension"] = 3
    cfg["mesh"]["extents"].append([0.0, 1.0e-6])
    cfg["mesh"]["resolution"].append(10)
    with pytest.raises(SchemaError, match="requires dimension=2"):
        validate(cfg)


def test_axisymmetric_rejects_negative_radial_extent():
    cfg = _base_axisym_cfg()
    cfg["mesh"]["extents"][0][0] = -1.0e-6
    with pytest.raises(SchemaError, match="r_min >= 0"):
        validate(cfg)


def test_axisymmetric_rejects_dirichlet_on_axis():
    cfg = _base_axisym_cfg()
    # Replace the body contact (planar bulk) with one pinned at r = 0.
    cfg["contacts"] = [
        {"name": "body", "facet": "body", "type": "ohmic", "voltage": 0.0},
        {"name": "axis_pin", "facet": "axis", "type": "ohmic", "voltage": 0.0},
        {"name": "gate", "facet": "gate", "type": "gate", "voltage": 0.0},
    ]
    with pytest.raises(SchemaError, match="symmetry axis"):
        validate(cfg)


def test_axisymmetric_unknown_value_rejected():
    cfg = _base_axisym_cfg()
    cfg["coordinate_system"] = "spherical"
    with pytest.raises(SchemaError):
        validate(cfg)


def test_cartesian_3d_still_works():
    """Adding the field must not break unrelated cartesian inputs."""
    cfg = _base_axisym_cfg()
    cfg["coordinate_system"] = "cartesian"
    cfg["dimension"] = 3
    cfg["mesh"]["extents"].append([0.0, 1.0e-6])
    cfg["mesh"]["resolution"].append(10)
    out = validate(copy.deepcopy(cfg))
    assert out["coordinate_system"] == "cartesian"
