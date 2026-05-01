"""Tests for backend -> PETSc options translation (M15 Phase D).

Pure-Python; no dolfinx required. Verifies the data-only mapping
function in semi.compute that runners use to build the petsc_options
dict for solve_nonlinear[_block].
"""
from __future__ import annotations

import pytest

from semi import compute


def test_cpu_mumps_returns_empty_overrides():
    """cpu-mumps must add no overrides, so the existing
    DEFAULT_PETSC_OPTIONS direct-LU stays bit-identical to pre-M15."""
    opts = compute.petsc_options_for_backend("cpu-mumps", device="cpu")
    assert opts == {}


def test_auto_input_rejected():
    with pytest.raises(compute.ConfigError, match=r"resolve_backend"):
        compute.petsc_options_for_backend("auto")


def test_unknown_backend_rejected():
    with pytest.raises(compute.ConfigError, match=r"unknown backend"):
        compute.petsc_options_for_backend("tpu-special")


def test_gpu_amgx_cuda_defaults():
    opts = compute.petsc_options_for_backend("gpu-amgx", device="cuda")
    assert opts["mat_type"] == "aijcusparse"
    assert opts["vec_type"] == "cuda"
    assert opts["ksp_type"] == "gmres"
    assert opts["pc_type"] == "amgx"
    # MUMPS factor solver default must be cancelled, otherwise PETSc
    # would try to apply it on top of the GPU PC.
    assert opts["pc_factor_mat_solver_type"] is None


def test_gpu_amgx_hip_uses_hipsparse_and_hip_vec():
    opts = compute.petsc_options_for_backend("gpu-amgx", device="hip")
    assert opts["mat_type"] == "aijhipsparse"
    assert opts["vec_type"] == "hip"


def test_gpu_amgx_rejects_cpu_device():
    with pytest.raises(compute.ConfigError, match=r"device"):
        compute.petsc_options_for_backend("gpu-amgx", device="cpu")


def test_gpu_amgx_rejects_invalid_pc():
    with pytest.raises(compute.ConfigError, match=r"preconditioner"):
        compute.petsc_options_for_backend(
            "gpu-amgx", device="cuda", preconditioner="hypre-boomeramg",
        )


def test_gpu_amgx_user_can_pick_gamg():
    opts = compute.petsc_options_for_backend(
        "gpu-amgx", device="cuda", preconditioner="gamg",
    )
    assert opts["pc_type"] == "gamg"


def test_gpu_amgx_user_can_override_ksp():
    opts = compute.petsc_options_for_backend(
        "gpu-amgx", device="cuda", linear_solver="bcgs",
    )
    assert opts["ksp_type"] == "bcgs"


def test_gpu_hypre_cuda_defaults():
    opts = compute.petsc_options_for_backend("gpu-hypre", device="cuda")
    assert opts["mat_type"] == "aijcusparse"
    assert opts["vec_type"] == "cuda"
    assert opts["ksp_type"] == "bcgs"
    assert opts["pc_type"] == "hypre"
    assert opts["pc_hypre_type"] == "boomeramg"
    assert opts["pc_hypre_boomeramg_relax_type_all"] == "Chebyshev"
    assert opts["pc_hypre_use_gpu"] == "true"
    assert opts["pc_factor_mat_solver_type"] is None


def test_gpu_hypre_rejects_invalid_pc():
    with pytest.raises(compute.ConfigError, match=r"preconditioner"):
        compute.petsc_options_for_backend(
            "gpu-hypre", device="cuda", preconditioner="amgx",
        )


def test_auto_device_defaults_to_cuda():
    """When the runner cannot determine the device, default to CUDA;
    HIP must be opted into explicitly."""
    opts = compute.petsc_options_for_backend("gpu-amgx", device="auto")
    assert opts["mat_type"] == "aijcusparse"
    assert opts["vec_type"] == "cuda"


def _minimal_cfg_with_solver(solver_block: dict) -> dict:
    return {
        "schema_version": "1.4.0",
        "name": "backend_settings_test",
        "dimension": 1,
        "mesh": {"source": "builtin", "extents": [[0.0, 1.0e-6]],
                 "resolution": [10],
                 "facets_by_plane": [
                     {"name": "L", "tag": 1, "axis": 0, "value": 0.0},
                     {"name": "R", "tag": 2, "axis": 0, "value": 1.0e-6},
                 ]},
        "regions": {"si": {"material": "Si", "tag": 1, "role": "semiconductor"}},
        "doping": [{"region": "si", "profile": {"type": "uniform", "N_D": 1e17, "N_A": 0.0}}],
        "contacts": [
            {"name": "L", "facet": "L", "type": "ohmic", "voltage": 0.0},
            {"name": "R", "facet": "R", "type": "ohmic", "voltage": 0.0},
        ],
        "solver": solver_block,
    }


def test_backend_settings_default_cpu_path():
    """The default-fill path produces a cpu-mumps backend with empty
    petsc_options overrides."""
    from semi import schema as schema_mod
    cfg = schema_mod.validate(_minimal_cfg_with_solver({}))
    info = compute.backend_settings_from_cfg(cfg)
    assert info["requested"] == "cpu-mumps"
    assert info["resolved"] == "cpu-mumps"
    assert info["device"] == "cpu"
    assert info["petsc_options"] == {}


def test_backend_settings_explicit_gpu_amgx_unavailable_raises():
    """Acceptance test A3: requesting an unavailable GPU backend must
    fail fast at backend resolution, not silently fall back."""
    from semi import schema as schema_mod
    cfg = schema_mod.validate(_minimal_cfg_with_solver({
        "backend": "gpu-amgx",
        "compute": {"device": "cuda"},
    }))
    # Force the available list to be CPU-only to simulate a CPU-only
    # CI host even if PETSc is somehow GPU-built.
    import unittest.mock as _mock
    with _mock.patch.object(
        compute, "available_backends", return_value=["cpu-mumps"],
    ):
        with pytest.raises(compute.ConfigError, match=r"gpu-amgx"):
            compute.backend_settings_from_cfg(cfg)


def test_backend_settings_auto_falls_back_to_cpu_when_no_gpu():
    """auto on a CPU-only build resolves to cpu-mumps (this is a legal
    user request that means 'use the best available'). Distinct from
    A3, where the user explicitly named a GPU backend."""
    from semi import schema as schema_mod
    cfg = schema_mod.validate(_minimal_cfg_with_solver({
        "backend": "auto",
    }))
    import unittest.mock as _mock
    with _mock.patch.object(
        compute, "available_backends", return_value=["cpu-mumps"],
    ):
        info = compute.backend_settings_from_cfg(cfg)
    assert info["requested"] == "auto"
    assert info["resolved"] == "cpu-mumps"
    assert info["device"] == "cpu"
    assert info["petsc_options"] == {}


def test_backend_settings_gpu_amgx_when_available():
    from semi import schema as schema_mod
    cfg = schema_mod.validate(_minimal_cfg_with_solver({
        "backend": "gpu-amgx",
        "compute": {"device": "cuda", "linear_solver": "gmres"},
    }))
    import unittest.mock as _mock
    with _mock.patch.object(
        compute, "available_backends",
        return_value=["cpu-mumps", "gpu-amgx"],
    ):
        info = compute.backend_settings_from_cfg(cfg)
    assert info["resolved"] == "gpu-amgx"
    assert info["device"] == "cuda"
    assert info["petsc_options"]["mat_type"] == "aijcusparse"
    assert info["petsc_options"]["pc_type"] == "amgx"
    assert info["petsc_options"]["ksp_type"] == "gmres"
