"""
Schema validation for the M14.2 axisymmetric MOSCAP additions.

Pure-Python (no dolfinx). Exercises:
    * mesh.axisymmetric flag is accepted on dimension == 2 builtin meshes
    * cv_analysis block validates with both LF and HF modes
    * solver.type "moscap_lf_hf" is accepted
"""
from __future__ import annotations

import pytest

from semi import schema


def _base_axisym_cfg():
    return {
        "schema_version": "1.3.0",
        "name": "moscap_axisym_schema",
        "dimension": 2,
        "mesh": {
            "source": "builtin",
            "extents": [[0.0, 1.0e-6], [0.0, 1.005e-6]],
            "resolution": [40, 201],
            "axisymmetric": True,
            "axisymmetric_axis": 0,
            "regions_by_box": [
                {"name": "silicon", "tag": 1, "bounds": [[0.0, 1.0e-6], [0.0,    1.0e-6]]},
                {"name": "oxide",   "tag": 2, "bounds": [[0.0, 1.0e-6], [1.0e-6, 1.005e-6]]},
            ],
            "facets_by_plane": [
                {"name": "body",  "tag": 1, "axis": 1, "value": 0.0},
                {"name": "gate",  "tag": 2, "axis": 1, "value": 1.005e-6},
                {"name": "axis",  "tag": 3, "axis": 0, "value": 0.0},
                {"name": "outer", "tag": 4, "axis": 0, "value": 1.0e-6},
            ],
        },
        "regions": {
            "silicon": {"material": "Si",   "tag": 1, "role": "semiconductor"},
            "oxide":   {"material": "SiO2", "tag": 2, "role": "insulator"},
        },
        "doping": [
            {"region": "silicon",
             "profile": {"type": "uniform", "N_D": 0.0, "N_A": 1.0e17}},
        ],
        "contacts": [
            {"name": "body", "facet": "body", "type": "ohmic", "voltage": 0.0},
            {"name": "gate", "facet": "gate", "type": "gate",
             "voltage": 0.0, "workfunction": -0.977,
             "voltage_sweep": {"start": -2.0, "stop": 2.5, "step": 0.05}},
        ],
        "solver": {"type": "moscap_lf_hf"},
        "cv_analysis": {
            "modes": ["LF", "HF"],
            "delta_V_small_signal": 1.0e-3,
            "majority_carrier": "auto",
            "gate_radius": 1.0e-6,
        },
    }


def test_axisym_cfg_validates():
    cfg = schema.validate(_base_axisym_cfg())
    assert cfg["mesh"]["axisymmetric"] is True
    assert cfg["mesh"]["axisymmetric_axis"] == 0
    assert cfg["solver"]["type"] == "moscap_lf_hf"
    assert cfg["cv_analysis"]["modes"] == ["LF", "HF"]


def test_cv_analysis_majority_carrier_validates_enum():
    cfg = _base_axisym_cfg()
    cfg["cv_analysis"]["majority_carrier"] = "not-a-valid-option"
    with pytest.raises(schema.SchemaError):
        schema.validate(cfg)


def test_cv_analysis_modes_only_LF_HF():
    cfg = _base_axisym_cfg()
    cfg["cv_analysis"]["modes"] = ["LF", "BAD"]
    with pytest.raises(schema.SchemaError):
        schema.validate(cfg)


def test_axisymmetric_axis_only_0_or_1():
    cfg = _base_axisym_cfg()
    cfg["mesh"]["axisymmetric_axis"] = 2
    with pytest.raises(schema.SchemaError):
        schema.validate(cfg)


def test_cv_analysis_optional():
    """cv_analysis is optional; runner falls back to defaults."""
    cfg = _base_axisym_cfg()
    cfg.pop("cv_analysis")
    schema.validate(cfg)


def test_axisymmetric_flag_optional_default_false():
    cfg = _base_axisym_cfg()
    cfg["mesh"].pop("axisymmetric")
    cfg["mesh"].pop("axisymmetric_axis")
    cfg = schema.validate(cfg)
    # Schema does not auto-fill; absence implies planar.
    assert "axisymmetric" not in cfg["mesh"] or cfg["mesh"]["axisymmetric"] is False
