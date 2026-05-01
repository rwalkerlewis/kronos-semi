"""Tests for schema 1.4.0: solver.backend and solver.compute (M15 Phase B).

Pure-Python; no dolfinx required. Verifies:
  - existing benchmark JSONs still validate against 1.4.0
  - cross-field rules between solver.backend and solver.compute.device fire
  - default-fill is byte-equivalent to today (no compute block, backend
    defaults to cpu-mumps)
  - schema-version major-gate still accepts 1.4.0 and rejects 2.0.0
"""
from __future__ import annotations

import copy
import glob
import json
from pathlib import Path

import pytest

from semi import schema

BENCHMARKS_DIR = Path(__file__).resolve().parents[1] / "benchmarks"


@pytest.fixture
def minimal_cfg():
    return {
        "schema_version": "1.4.0",
        "name": "compute_test",
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
            {"region": "si", "profile": {"type": "uniform", "N_D": 1.0e17, "N_A": 0.0}},
        ],
        "contacts": [
            {"name": "L", "facet": "L", "type": "ohmic", "voltage": 0.0},
            {"name": "R", "facet": "R", "type": "ohmic", "voltage": 0.0},
        ],
    }


def test_supported_minor_bumped_to_4():
    assert schema.SCHEMA_SUPPORTED_MINOR == 4


def test_schema_version_140_accepted(minimal_cfg):
    minimal_cfg["schema_version"] = "1.4.0"
    result = schema.validate(minimal_cfg)
    assert result["solver"]["backend"] == "cpu-mumps"


def test_schema_version_major_two_rejected(minimal_cfg):
    minimal_cfg["schema_version"] = "2.0.0"
    with pytest.raises(schema.SchemaError, match=r"major"):
        schema.validate(minimal_cfg)


def test_default_fill_no_compute_block(minimal_cfg):
    """Default-fill: no backend, no compute -> backend=cpu-mumps, no compute."""
    result = schema.validate(minimal_cfg)
    assert result["solver"]["backend"] == "cpu-mumps"
    assert "compute" not in result["solver"], (
        "default-fill must NOT add a compute block; today's CPU runs stay "
        "byte-equivalent only when the compute object is absent"
    )
    # Existing linear_solver defaults still apply.
    assert result["solver"]["linear_solver"]["pc_type"] == "lu"
    assert result["solver"]["linear_solver"]["factor_mat_solver_type"] == "mumps"


def test_default_fill_round_trips(minimal_cfg):
    """Validating an already-validated config is idempotent."""
    once = schema.validate(copy.deepcopy(minimal_cfg))
    twice = schema.validate(copy.deepcopy(once))
    assert once["solver"]["backend"] == twice["solver"]["backend"]
    assert ("compute" in once["solver"]) == ("compute" in twice["solver"])


def test_explicit_cpu_mumps_with_cpu_device_ok(minimal_cfg):
    minimal_cfg["solver"] = {
        "backend": "cpu-mumps",
        "compute": {"device": "cpu"},
    }
    result = schema.validate(minimal_cfg)
    assert result["solver"]["backend"] == "cpu-mumps"
    assert result["solver"]["compute"]["device"] == "cpu"


def test_explicit_cpu_mumps_with_auto_device_ok(minimal_cfg):
    minimal_cfg["solver"] = {
        "backend": "cpu-mumps",
        "compute": {"device": "auto"},
    }
    schema.validate(minimal_cfg)


def test_cpu_mumps_rejects_cuda_device(minimal_cfg):
    minimal_cfg["solver"] = {
        "backend": "cpu-mumps",
        "compute": {"device": "cuda"},
    }
    with pytest.raises(schema.SchemaError, match=r"cpu-mumps"):
        schema.validate(minimal_cfg)


def test_cpu_mumps_rejects_hip_device(minimal_cfg):
    minimal_cfg["solver"] = {
        "backend": "cpu-mumps",
        "compute": {"device": "hip"},
    }
    with pytest.raises(schema.SchemaError, match=r"cpu-mumps"):
        schema.validate(minimal_cfg)


def test_gpu_amgx_with_cuda_ok(minimal_cfg):
    minimal_cfg["solver"] = {
        "backend": "gpu-amgx",
        "compute": {"device": "cuda"},
    }
    schema.validate(minimal_cfg)


def test_gpu_amgx_with_auto_ok(minimal_cfg):
    minimal_cfg["solver"] = {
        "backend": "gpu-amgx",
        "compute": {"device": "auto"},
    }
    schema.validate(minimal_cfg)


def test_gpu_amgx_rejects_cpu_device(minimal_cfg):
    minimal_cfg["solver"] = {
        "backend": "gpu-amgx",
        "compute": {"device": "cpu"},
    }
    with pytest.raises(schema.SchemaError, match=r"gpu-amgx"):
        schema.validate(minimal_cfg)


def test_gpu_hypre_rejects_cpu_device(minimal_cfg):
    minimal_cfg["solver"] = {
        "backend": "gpu-hypre",
        "compute": {"device": "cpu"},
    }
    with pytest.raises(schema.SchemaError, match=r"gpu-hypre"):
        schema.validate(minimal_cfg)


def test_gpu_hypre_with_hip_ok(minimal_cfg):
    minimal_cfg["solver"] = {
        "backend": "gpu-hypre",
        "compute": {"device": "hip"},
    }
    schema.validate(minimal_cfg)


def test_backend_auto_validation_time_passthrough(minimal_cfg):
    """`auto` is resolved at solve time (Phase C); validation accepts any
    legal device (including device left unspecified)."""
    minimal_cfg["solver"] = {"backend": "auto"}
    schema.validate(minimal_cfg)
    minimal_cfg["solver"] = {"backend": "auto", "compute": {"device": "cuda"}}
    schema.validate(minimal_cfg)
    minimal_cfg["solver"] = {"backend": "auto", "compute": {"device": "cpu"}}
    schema.validate(minimal_cfg)


def test_unknown_backend_rejected(minimal_cfg):
    minimal_cfg["solver"] = {"backend": "tpu-special"}
    with pytest.raises(schema.SchemaError):
        schema.validate(minimal_cfg)


def test_compute_precision_only_float64(minimal_cfg):
    minimal_cfg["solver"] = {"compute": {"precision": "float32"}}
    with pytest.raises(schema.SchemaError):
        schema.validate(minimal_cfg)


def test_compute_additional_properties_rejected(minimal_cfg):
    minimal_cfg["solver"] = {"compute": {"device": "cpu", "extra_key": 42}}
    with pytest.raises(schema.SchemaError):
        schema.validate(minimal_cfg)


def test_compute_default_values_present_in_schema():
    s = schema.SCHEMA["properties"]["solver"]["properties"]
    assert s["backend"]["default"] == "cpu-mumps"
    assert "cpu-mumps" in s["backend"]["enum"]
    assert "gpu-amgx" in s["backend"]["enum"]
    assert "gpu-hypre" in s["backend"]["enum"]
    assert "auto" in s["backend"]["enum"]
    c = s["compute"]["properties"]
    assert c["device"]["default"] == "cpu"
    assert c["precision"]["default"] == "float64"
    assert c["preconditioner"]["default"] == "auto"
    assert c["linear_solver"]["default"] == "auto"


def test_all_existing_benchmarks_still_validate():
    """Every benchmark JSON validates unchanged against schema 1.4.0."""
    paths = sorted(glob.glob(str(BENCHMARKS_DIR / "*" / "*.json")))
    assert paths, "no benchmark JSONs found under benchmarks/*/*.json"
    for p in paths:
        cfg = json.loads(Path(p).read_text())
        result = schema.validate(cfg)
        assert result["solver"]["backend"] == "cpu-mumps", (
            f"{p}: default-fill must keep cpu-mumps (no behavior change)"
        )
        # Compute object must not be auto-added; today's CPU path is
        # byte-equivalent only when no compute block is present.
        assert "compute" not in result["solver"], (
            f"{p}: compute block was injected; this would not be a "
            f"byte-equivalent change"
        )
