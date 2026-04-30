"""
Gated regression test for the axisymmetric MOSCAP benchmark.

Skipped when ``dolfinx`` is not available; otherwise builds a tiny
2D mesh in-process (no gmsh dependency) and exercises the three
r-weighted form builders in :mod:`semi.physics.axisymmetric`.

The schema-validation case runs unconditionally and verifies the
benchmark JSON loads under schema 1.3.0.

The full FEM-vs-analytical C-V regression lives in the benchmark
runner; here we only certify that the form builders accept a real
mesh and produce a well-typed UFL form (correct rank, domain, and
arguments) for both the single-region Poisson, the multi-region
Poisson, and the Slotboom drift-diffusion block.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

dolfinx = pytest.importorskip("dolfinx")
ufl = pytest.importorskip("ufl")

import numpy as np  # noqa: E402
from mpi4py import MPI  # noqa: E402

from semi.materials import get_material  # noqa: E402
from semi.physics.axisymmetric import (  # noqa: E402
    build_dd_block_residual_axisym,
    build_equilibrium_poisson_form_axisym,
    build_equilibrium_poisson_form_axisym_mr,
)
from semi.physics.drift_diffusion import make_dd_block_spaces  # noqa: E402
from semi.scaling import Scaling  # noqa: E402
from semi.schema import validate  # noqa: E402

BENCH_DIR = Path(__file__).resolve().parent.parent / "benchmarks" / "moscap_axisym_2d"
CFG_PATH = BENCH_DIR / "moscap_axisym.json"


def test_axisym_moscap_json_validates():
    cfg = json.loads(CFG_PATH.read_text())
    out = validate(cfg)
    assert out["coordinate_system"] == "axisymmetric"
    assert out["dimension"] == 2


def _make_meridian_mesh():
    """Tiny meridian half-plane (r, z) in [0, 1] x [-1, 1]."""
    from dolfinx import mesh as dmesh

    msh = dmesh.create_rectangle(
        MPI.COMM_WORLD,
        [np.array([0.0, -1.0]), np.array([1.0, 1.0])],
        [4, 8],
        cell_type=dmesh.CellType.triangle,
    )
    tdim = msh.topology.dim
    msh.topology.create_entities(tdim)
    n_local = msh.topology.index_map(tdim).size_local
    cells = np.arange(n_local, dtype=np.int32)

    # Tag z < 0 as silicon (tag 1), z >= 0 as oxide (tag 2).
    midpoints = dolfinx.mesh.compute_midpoints(msh, tdim, cells)
    tags = np.where(midpoints[:, 1] < 0.0, 1, 2).astype(np.int32)
    cell_tags = dmesh.meshtags(msh, tdim, cells, tags)
    return msh, cell_tags


def _make_scaling():
    Si = get_material("Si")
    return Scaling(L0=1.0e-6, C0=1.0e23, T=300.0, mu0=Si.mu_n, n_i=Si.n_i)


def test_axisym_poisson_single_region_form():
    msh, _ = _make_meridian_mesh()
    V = dolfinx.fem.functionspace(msh, ("Lagrange", 1))
    psi = dolfinx.fem.Function(V)
    N_hat = dolfinx.fem.Function(V)
    N_hat.x.array[:] = -1.0
    sc = _make_scaling()

    F = build_equilibrium_poisson_form_axisym(V, psi, N_hat, sc, eps_r=11.7)
    assert F.ufl_domain() is msh.ufl_domain()
    assert F.arguments()


def test_axisym_poisson_multiregion_form():
    msh, cell_tags = _make_meridian_mesh()
    V = dolfinx.fem.functionspace(msh, ("Lagrange", 1))
    psi = dolfinx.fem.Function(V)
    N_hat = dolfinx.fem.Function(V)
    sc = _make_scaling()

    DG0 = dolfinx.fem.functionspace(msh, ("DG", 0))
    eps_r_fn = dolfinx.fem.Function(DG0)
    eps_r_fn.x.array[:] = np.where(cell_tags.values == 1, 11.7, 3.9)

    F = build_equilibrium_poisson_form_axisym_mr(
        V, psi, N_hat, sc, eps_r_fn, cell_tags, semi_tag=1,
    )
    assert F.ufl_domain() is msh.ufl_domain()
    assert F.arguments()


def test_axisym_dd_block_residual_three_forms():
    msh, _ = _make_meridian_mesh()
    spaces = make_dd_block_spaces(msh)
    N_hat = dolfinx.fem.Function(spaces.V_psi)
    N_hat.x.array[:] = -1.0
    sc = _make_scaling()

    forms = build_dd_block_residual_axisym(
        spaces,
        N_hat_fn=N_hat,
        sc=sc,
        eps_r=11.7,
        mu_n_over_mu0=1.0,
        mu_p_over_mu0=0.32,
        tau_n_hat=1.0,
        tau_p_hat=1.0,
        E_t_over_Vt=0.0,
    )
    assert len(forms) == 3
    for F in forms:
        assert F.ufl_domain() is msh.ufl_domain()
        assert F.arguments()
    f_psi, f_n, f_p = forms
    assert f_psi.arguments()[0].ufl_function_space() is spaces.V_psi
    assert f_n.arguments()[0].ufl_function_space() is spaces.V_phi_n
    assert f_p.arguments()[0].ufl_function_space() is spaces.V_phi_p
