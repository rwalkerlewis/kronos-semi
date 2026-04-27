"""
Unit tests for the per-edge SG residual assembly (M13.1 Path B1).

Verifies that `assemble_sg_residual_1d` correctly:

1. Computes the per-edge SG flux matching the closed-form
   `sg_edge_flux_n` / `sg_edge_flux_p` (algebraic agreement at every
   Peclet, not just zero field).
2. Scatters the per-edge contributions into per-vertex residual entries
   with the correct sign convention: `+L0_sq * F` at the left vertex
   of each edge, `-L0_sq * F` at the right vertex (divergence form of
   `dn/dt = -div(F_n)` after IBP).
3. Returns zero net residual at interior vertices when all edges are
   in equilibrium (uniform doping, balanced n*p = n_i^2 case).
4. Has correct boundary contribution: residual at boundary vertices
   accumulates the single-edge flux from the unique adjacent cell.
"""
from __future__ import annotations

import numpy as np
import pytest

dolfinx = pytest.importorskip("dolfinx")

import dolfinx.mesh  # noqa: E402,F811
from mpi4py import MPI  # noqa: E402

from semi.fem.scharfetter_gummel import sg_edge_flux_n, sg_edge_flux_p  # noqa: E402
from semi.fem.sg_assembly import (  # noqa: E402
    assemble_sg_jacobian_correction_1d,
    assemble_sg_residual_1d,
    compute_edge_topology_1d,
    vertex_to_dof_map,
)

MU_N_HAT = 1.0  # scaled mobilities (test uses raw scaling factor 1)
MU_P_HAT = 0.32  # roughly mu_p / mu_n for Si
L0_SQ = 1.0  # neutralised so test compares directly


def _build_test_mesh_and_spaces(n_cells: int, L: float = 1.0):
    """Create a 1D unit-interval mesh with `n_cells` cells, plus three P1 spaces."""
    msh = dolfinx.mesh.create_interval(MPI.COMM_WORLD, n_cells, [0.0, L])
    from dolfinx.fem import functionspace
    V = functionspace(msh, ("Lagrange", 1))
    return msh, V


def test_compute_edge_topology_1d_returns_left_right_ordered():
    """Per-cell vertex pairs must be sorted by x: cell_vertices[k][0] is left."""
    msh, V = _build_test_mesh_and_spaces(10)
    cell_vertices, h_cells, cell_x = compute_edge_topology_1d(msh)
    assert cell_vertices.shape == (10, 2)
    assert h_cells.shape == (10,)
    # Each cell has length L/n_cells = 0.1
    np.testing.assert_allclose(h_cells, 0.1, rtol=1e-10)
    # Left x < right x for every cell
    assert (cell_x[:, 0] < cell_x[:, 1]).all()


def test_assemble_sg_residual_zero_at_zero_carriers_zero_field():
    """Trivial case: n = p = 0, psi = 0 -> residual = 0 everywhere."""
    msh, V = _build_test_mesh_and_spaces(10)
    cell_v, h, _ = compute_edge_topology_1d(msh)
    v2d = vertex_to_dof_map(V)

    n_dofs = V.dofmap.index_map.size_local
    psi_arr = np.zeros(n_dofs)
    n_arr = np.zeros(n_dofs)
    p_arr = np.zeros(n_dofs)

    R_n, R_p = assemble_sg_residual_1d(
        psi_arr, n_arr, p_arr,
        cell_v, h, v2d, v2d, v2d,
        MU_N_HAT, MU_P_HAT, L0_SQ,
    )
    np.testing.assert_array_equal(R_n, 0.0)
    np.testing.assert_array_equal(R_p, 0.0)


def test_assemble_sg_residual_uniform_n_zero_field():
    """
    Uniform n with zero psi gradient: per-edge flux is zero (no drift,
    no diffusion), so residual is zero at every vertex.
    """
    msh, V = _build_test_mesh_and_spaces(10)
    cell_v, h, _ = compute_edge_topology_1d(msh)
    v2d = vertex_to_dof_map(V)

    n_dofs = V.dofmap.index_map.size_local
    psi_arr = np.zeros(n_dofs)
    n_arr = np.full(n_dofs, 1.0e16)  # uniform n
    p_arr = np.full(n_dofs, 1.0e10)  # uniform p

    R_n, R_p = assemble_sg_residual_1d(
        psi_arr, n_arr, p_arr,
        cell_v, h, v2d, v2d, v2d,
        MU_N_HAT, MU_P_HAT, L0_SQ,
    )
    np.testing.assert_allclose(R_n, 0.0, atol=1e-6)
    np.testing.assert_allclose(R_p, 0.0, atol=1e-6)


def test_assemble_sg_residual_single_edge_matches_closed_form():
    """
    On a 2-cell, 3-vertex mesh with non-trivial state, the per-vertex
    residuals must match per-edge scalar flux scattered with sign
    convention {+L0_sq*F at left vertex, -L0_sq*F at right vertex}.
    """
    msh, V = _build_test_mesh_and_spaces(2)
    cell_v, h, cell_x = compute_edge_topology_1d(msh)
    v2d = vertex_to_dof_map(V)

    # Two cells, three vertices. Set state.
    n_dofs = V.dofmap.index_map.size_local
    psi_arr = np.zeros(n_dofs)
    n_arr = np.zeros(n_dofs)
    p_arr = np.zeros(n_dofs)
    # vertex 0: left, vertex 1: middle, vertex 2: right (assuming dolfinx ordering)
    # Use vertex coordinates to identify
    x_to_v = {round(float(msh.geometry.x[v, 0]), 6): v for v in range(n_dofs)}
    v_left = x_to_v[0.0]
    v_mid = x_to_v[0.5]
    v_right = x_to_v[1.0]
    psi_arr[v2d[v_left]] = 0.0
    psi_arr[v2d[v_mid]] = 0.5
    psi_arr[v2d[v_right]] = 1.5
    n_arr[v2d[v_left]] = 1.0e16
    n_arr[v2d[v_mid]] = 1.0e17
    n_arr[v2d[v_right]] = 1.0e18
    p_arr[v2d[v_left]] = 1.0e15
    p_arr[v2d[v_mid]] = 1.0e14
    p_arr[v2d[v_right]] = 1.0e13

    R_n, R_p = assemble_sg_residual_1d(
        psi_arr, n_arr, p_arr,
        cell_v, h, v2d, v2d, v2d,
        MU_N_HAT, MU_P_HAT, L0_SQ,
    )

    # Compute reference via the scalar `sg_edge_flux_n`/`p` per cell.
    R_n_ref = np.zeros(n_dofs)
    R_p_ref = np.zeros(n_dofs)
    for k in range(2):
        i, j = cell_v[k]
        h_k = h[k]
        dpsi = psi_arr[v2d[j]] - psi_arr[v2d[i]]
        F_n = sg_edge_flux_n(
            n_arr[v2d[i]], n_arr[v2d[j]], dpsi, h_k, MU_N_HAT, 1.0,
        )
        F_p = sg_edge_flux_p(
            p_arr[v2d[i]], p_arr[v2d[j]], dpsi, h_k, MU_P_HAT, 1.0,
        )
        # Note: sg_edge_flux_n uses V_t in its arg; we pass V_t=1 so the
        # raw mu_n*V_t/h matches our (mu/h) convention here.
        R_n_ref[v2d[i]] += L0_SQ * F_n
        R_n_ref[v2d[j]] -= L0_SQ * F_n
        R_p_ref[v2d[i]] += L0_SQ * F_p
        R_p_ref[v2d[j]] -= L0_SQ * F_p

    np.testing.assert_allclose(R_n, R_n_ref, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(R_p, R_p_ref, rtol=1e-12, atol=1e-12)


def test_assemble_sg_residual_interior_vertex_is_div_of_flux():
    """
    On a 3-cell mesh with non-trivial state, the residual at the
    interior vertex `v` is `+F(left edge) - F(right edge)`, i.e., the
    discrete divergence of the flux.
    """
    msh, V = _build_test_mesh_and_spaces(3)
    cell_v, h, cell_x = compute_edge_topology_1d(msh)
    v2d = vertex_to_dof_map(V)

    n_dofs = V.dofmap.index_map.size_local
    psi_arr = np.linspace(0.0, 1.0, n_dofs)  # linearly increasing psi
    # Re-order so it's actually monotone in x:
    coord_order = np.argsort(msh.geometry.x[:, 0])
    psi_dof_in_x = np.linspace(0.0, 1.0, n_dofs)
    psi_arr_reordered = np.zeros(n_dofs)
    for k, v in enumerate(coord_order):
        psi_arr_reordered[v2d[v]] = psi_dof_in_x[k]
    psi_arr = psi_arr_reordered
    n_arr = np.full(n_dofs, 1.0e16)
    p_arr = np.full(n_dofs, 1.0e10)

    R_n, R_p = assemble_sg_residual_1d(
        psi_arr, n_arr, p_arr,
        cell_v, h, v2d, v2d, v2d,
        MU_N_HAT, MU_P_HAT, L0_SQ,
    )

    # All cells should have IDENTICAL flux (uniform n, linearly varying psi
    # so dpsi is the same on every cell). Then the interior vertex residuals
    # are F - F = 0 (div = 0 of constant flux). Boundary vertices get +-F.
    interior_verts = [v for k, v in enumerate(coord_order) if 0 < k < n_dofs - 1]
    # Relative tolerance: flux magnitude is ~mu*n/h ~ 1e16; allow 1e-12 relative.
    flux_scale = MU_N_HAT * 1.0e16 * 3.0  # rough order of magnitude
    rel_tol = 1.0e-12
    abs_tol = flux_scale * rel_tol
    for v in interior_verts:
        dof = v2d[v]
        assert abs(R_n[dof]) < abs_tol, (
            f"interior vertex {v}: R_n = {R_n[dof]} (tol {abs_tol:.3e})"
        )
        assert abs(R_p[dof]) < abs_tol, (
            f"interior vertex {v}: R_p = {R_p[dof]} (tol {abs_tol:.3e})"
        )


def _galerkin_convdiff_residual_1d(
    psi_arr, n_arr, p_arr,
    cell_vertices, h_cells, v2d_psi, v2d_n, v2d_p,
    mu_n_hat, mu_p_hat, L0_sq,
):
    """
    Reference Galerkin convection-diffusion residual matching the M13 form
    `inner(grad(*), grad(v)) -+ *grad(psi).grad(v)` integrated cell by cell.

    Used inside the FD-verification regression test. Independent of UFL
    so the regression test does not need to assemble a dolfinx form.
    """
    R_gal_n = np.zeros_like(n_arr)
    R_gal_p = np.zeros_like(p_arr)
    for k in range(len(h_cells)):
        i, j = int(cell_vertices[k, 0]), int(cell_vertices[k, 1])
        h = float(h_cells[k])
        psi_i = psi_arr[v2d_psi[i]]
        psi_j = psi_arr[v2d_psi[j]]
        n_i = n_arr[v2d_n[i]]
        n_j = n_arr[v2d_n[j]]
        p_i = p_arr[v2d_p[i]]
        p_j = p_arr[v2d_p[j]]
        dpsi = psi_j - psi_i
        n_avg = 0.5 * (n_i + n_j)
        p_avg = 0.5 * (p_i + p_j)
        gal_n_i = -L0_sq * mu_n_hat * ((n_j - n_i) / h - n_avg * dpsi / h)
        gal_p_i = -L0_sq * mu_p_hat * ((p_j - p_i) / h + p_avg * dpsi / h)
        R_gal_n[v2d_n[i]] += gal_n_i
        R_gal_n[v2d_n[j]] -= gal_n_i
        R_gal_p[v2d_p[i]] += gal_p_i
        R_gal_p[v2d_p[j]] -= gal_p_i
    return R_gal_n, R_gal_p


def _correction_residual_1d(
    psi_arr, n_arr, p_arr,
    cell_vertices, h_cells, v2d_psi, v2d_n, v2d_p,
    mu_n_hat, mu_p_hat, L0_sq,
):
    """`-sg_n - g_n` (and same for p) — the quantity the SG callback adds to b_M13."""
    sg_n, sg_p = assemble_sg_residual_1d(
        psi_arr, n_arr, p_arr,
        cell_vertices, h_cells, v2d_psi, v2d_n, v2d_p,
        mu_n_hat, mu_p_hat, L0_sq,
    )
    g_n, g_p = _galerkin_convdiff_residual_1d(
        psi_arr, n_arr, p_arr,
        cell_vertices, h_cells, v2d_psi, v2d_n, v2d_p,
        mu_n_hat, mu_p_hat, L0_sq,
    )
    return -sg_n - g_n, -sg_p - g_p


def test_assemble_sg_jacobian_correction_fd_verified():
    """
    Finite-difference verification of `assemble_sg_jacobian_correction_1d`
    against the corresponding correction residual `-sg_n - g_n` on a
    3-node 1D mesh at a non-equilibrium state.

    Why this test exists
    --------------------
    The 1D SG transient runner overrides the SNES residual so that

        b_corrected = b_M13_galerkin - g_n - sg_n   (n block; mirror for p)

    The Jacobian "correction" added to the auto-derived UFL Jacobian is
    `-d(g_n)/dx - d(sg_n)/dx`. A sign or factor-of-2 error in any of the
    16 entries per edge is invisible at zero field (where SG = Galerkin)
    but blows up Newton convergence once a depletion region with a
    non-trivial dpsi develops. This test is the regression guard that
    locks in the formulae in `assemble_sg_jacobian_correction_1d`.

    The script `scripts/verify_sg_jacobian_fd.py` is the same
    verification with prettier output for diagnostic use.
    """
    n_verts = 3
    cell_vertices = np.array([[0, 1], [1, 2]], dtype=np.int64)
    h_cells = np.array([5.0e-6, 5.0e-6], dtype=np.float64)
    v2d = np.arange(n_verts, dtype=np.int64)

    dof_offsets = (0, n_verts, 2 * n_verts)
    bc_dofs_n: set[int] = set()
    bc_dofs_p: set[int] = set()

    mu_n_hat = 1.0
    mu_p_hat = 0.32
    L0_sq = 1.0e-12  # ~ (1 um)^2

    # Non-equilibrium state with non-trivial dpsi (~ depletion-region scale).
    psi_arr = np.array([0.0, 1.5, 4.0], dtype=np.float64)
    n_arr = np.array([1.0e18, 5.0e16, 1.0e16], dtype=np.float64)
    p_arr = np.array([1.0e15, 5.0e16, 1.0e18], dtype=np.float64)

    # Analytic.
    rows, cols, vals = assemble_sg_jacobian_correction_1d(
        psi_arr, n_arr, p_arr,
        cell_vertices, h_cells, v2d, v2d, v2d,
        mu_n_hat, mu_p_hat, L0_sq,
        dof_offsets, bc_dofs_n, bc_dofs_p,
    )
    n_total = 3 * n_verts
    J_analytic = np.zeros((n_total, n_total), dtype=np.float64)
    for r, c, v in zip(rows, cols, vals):
        J_analytic[int(r), int(c)] += float(v)

    # Finite-difference reference.
    eps_rel = 1.0e-7
    J_fd = np.zeros((n_total, n_total), dtype=np.float64)
    for col in range(n_total):
        psi_p = psi_arr.copy(); psi_m = psi_arr.copy()
        n_p = n_arr.copy(); n_m = n_arr.copy()
        p_p = p_arr.copy(); p_m = p_arr.copy()

        if col < n_verts:
            v = col
            scale = max(abs(psi_arr[v]), 1.0)
            h = eps_rel * scale
            psi_p[v] += h; psi_m[v] -= h
        elif col < 2 * n_verts:
            v = col - n_verts
            scale = max(abs(n_arr[v]), 1.0)
            h = eps_rel * scale
            n_p[v] += h; n_m[v] -= h
        else:
            v = col - 2 * n_verts
            scale = max(abs(p_arr[v]), 1.0)
            h = eps_rel * scale
            p_p[v] += h; p_m[v] -= h

        Rp_n, Rp_p = _correction_residual_1d(
            psi_p, n_p, p_p,
            cell_vertices, h_cells, v2d, v2d, v2d,
            mu_n_hat, mu_p_hat, L0_sq,
        )
        Rm_n, Rm_p = _correction_residual_1d(
            psi_m, n_m, p_m,
            cell_vertices, h_cells, v2d, v2d, v2d,
            mu_n_hat, mu_p_hat, L0_sq,
        )
        Rp = np.zeros(n_total)
        Rp[n_verts:2 * n_verts] = Rp_n
        Rp[2 * n_verts:] = Rp_p
        Rm = np.zeros(n_total)
        Rm[n_verts:2 * n_verts] = Rm_n
        Rm[2 * n_verts:] = Rm_p
        J_fd[:, col] = (Rp - Rm) / (2.0 * h)

    # Compare entry-by-entry over the n,p row block (psi rows are zero).
    rel_thresh = 1.0e-5
    abs_floor = 1.0e-12  # treat anything smaller as zero
    mismatches = []
    for r in range(n_verts, 3 * n_verts):
        for c in range(3 * n_verts):
            j_a = J_analytic[r, c]
            j_f = J_fd[r, c]
            denom = max(abs(j_a), abs(j_f), abs_floor)
            rel = abs(j_a - j_f) / denom
            if rel > rel_thresh and max(abs(j_a), abs(j_f)) > abs_floor:
                mismatches.append((r, c, j_a, j_f, rel))

    assert not mismatches, (
        f"{len(mismatches)} Jacobian-correction entries failed FD verification "
        f"at rel tol {rel_thresh}:\n"
        + "\n".join(
            f"  J[{r},{c}] = {j_a:+.6e}  vs FD {j_f:+.6e}  (rel {rel:.3e})"
            for r, c, j_a, j_f, rel in mismatches
        )
    )
