"""
Gated regression test for the axisymmetric MOSCAP benchmark.

This test is skipped when dolfinx is not available, but its
schema-validation + form-builder smoke tests run anywhere dolfinx
is installed (the CI image, the kronos-semi Docker image).

The full LF-vs-HF C-V curve check (against
``benchmarks/moscap_axisym_2d/reference_cv.csv``) is exercised end-to-end
by the notebook and by the benchmark runner. Here we only verify:

  - the JSON validates under schema 1.3.0,
  - the gmsh mesh, when present, can be loaded,
  - the r-weighted form builders accept that mesh and produce a
    well-typed UFL form (no silent dimension or measure mismatches).
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

dolfinx = pytest.importorskip("dolfinx")
gmshio = pytest.importorskip("dolfinx.io.gmshio")
ufl = pytest.importorskip("ufl")

from mpi4py import MPI

from semi.schema import validate
from semi.physics.axisymmetric import (
    build_equilibrium_poisson_form_axisym_mr,
)
from semi.scaling import ScalingContext


BENCH_DIR = Path(__file__).resolve().parent.parent / "benchmarks" / "moscap_axisym_2d"
CFG_PATH = BENCH_DIR / "moscap_axisym.json"
MSH_PATH = BENCH_DIR / "moscap_axisym.msh"


def test_axisym_moscap_json_validates():
    cfg = json.loads(CFG_PATH.read_text())
    out = validate(cfg)
    assert out["coordinate_system"] == "axisymmetric"
    assert out["dimension"] == 2


@pytest.mark.skipif(not MSH_PATH.exists(), reason="gmsh mesh not generated; run gmsh -2 first")
def test_axisym_form_builder_accepts_gmsh_mesh():
    msh, cell_tags, _ = gmshio.read_from_msh(
        str(MSH_PATH), MPI.COMM_WORLD, gdim=2,
    )
    V = dolfinx.fem.functionspace(msh, ("Lagrange", 1))
    psi = dolfinx.fem.Function(V)
    N_hat = dolfinx.fem.Function(V)

    sc = ScalingContext.default(temperature=300.0)
    eps_r_fn = dolfinx.fem.Function(dolfinx.fem.functionspace(msh, ("DG", 0)))
    eps_r_fn.x.array[:] = 11.7  # placeholder; multi-region map is filled by the runner

    F = build_equilibrium_poisson_form_axisym_mr(
        V, psi, N_hat, sc, eps_r_fn, cell_tags, semi_tag=1,
    )
    # The form must be a UFL Form over the right mesh and have the
    # right rank; we don't assemble it here (the runner does).
    assert F.ufl_domain().ufl_cargo() is msh
    assert F.arguments()  # has a TestFunction
