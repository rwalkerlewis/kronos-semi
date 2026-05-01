"""GPU backend acceptance tests (M15).

These tests exercise the full GPU linear-solver path. They skip
gracefully on any host where PETSc was built without the requested
GPU backend, so the suite remains green on CPU-only CI.

A1 (correctness): solve the pn_1d benchmark on cpu-mumps and on every
available GPU backend. The psi solution must agree to within 1e-8
relative L2.
"""
from __future__ import annotations

import json
import pathlib

import numpy as np
import pytest

dolfinx = pytest.importorskip("dolfinx")  # noqa: F841

from semi import compute  # noqa: E402
from semi.run import run as run_simulation  # noqa: E402

_REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
_PN1D_INPUT = _REPO_ROOT / "benchmarks" / "pn_1d" / "input.json"


def _gpu_backends_available() -> list[str]:
    return [b for b in compute.available_backends() if b.startswith("gpu-")]


def _load_pn1d_cfg() -> dict:
    with _PN1D_INPUT.open("r", encoding="utf-8") as f:
        return json.load(f)


@pytest.mark.parametrize("gpu_backend", _gpu_backends_available() or ["__none__"])
def test_pn1d_gpu_matches_cpu_mumps(gpu_backend):
    """Acceptance test A1: GPU result matches CPU-MUMPS reference."""
    if gpu_backend == "__none__":
        pytest.skip(
            f"no GPU backend in PETSc build: available={compute.available_backends()}"
        )

    cfg_cpu = _load_pn1d_cfg()
    cfg_cpu.setdefault("solver", {})["backend"] = "cpu-mumps"
    res_cpu = run_simulation(cfg_cpu)

    cfg_gpu = _load_pn1d_cfg()
    cfg_gpu.setdefault("solver", {})["backend"] = gpu_backend
    cfg_gpu["solver"].setdefault("compute", {})
    res_gpu = run_simulation(cfg_gpu)

    psi_cpu = np.asarray(res_cpu.psi.x.array)
    psi_gpu = np.asarray(res_gpu.psi.x.array)
    rel = np.linalg.norm(psi_gpu - psi_cpu) / np.linalg.norm(psi_cpu)
    assert rel < 1.0e-8, (
        f"GPU psi disagrees with CPU-MUMPS reference: "
        f"backend={gpu_backend} rel L2 = {rel:.3e}"
    )

    # Backend metadata must propagate through the solver layer to the
    # SimulationResult so the manifest is correctly populated.
    info = res_gpu.solver_info
    assert info["backend_resolved"] == gpu_backend
    assert info["device"] in ("cuda", "hip")
