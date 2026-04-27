"""
Finite-difference verification of `assemble_sg_jacobian_correction_1d`
against the corresponding correction residual `-sg_n - g_n` (the
quantity that the SG callback ADDS to the M13 Galerkin block residual).

Goal
----
The SG transient runner overrides the SNES residual so that

    b_corrected = b_M13_galerkin - g_n - sg_n   (n block; mirror for p)

where `b_M13_galerkin` is the standard UFL block residual, `g_n` is the
M13 Galerkin convection-diffusion contribution (auto-assembled by UFL,
also subtracted), and `sg_n` is the per-edge SG residual in the SG sign
convention. The "correction" we add to the UFL Jacobian is therefore
`-d(g_n)/dx - d(sg_n)/dx`, which is exactly what
`assemble_sg_jacobian_correction_1d` returns.

This script FD-verifies that returned correction by:

1. Building a 3-node 1D mesh on [0, 5e-6, 1e-5] m (two cells).
2. Setting a non-equilibrium state with |dpsi| ~ 2-5 (depletion-region scale).
3. Computing the analytic Jacobian correction.
4. Computing the same Jacobian by central-difference of `R_correction`.
5. Reporting any entry whose relative error exceeds the threshold.

If the Jacobian correction is correct, the FD and analytic matrices
agree to ~ 1e-6 relative on every nonzero entry. A sign or factor-of-2
error in any of the 16 entries per edge will surface as a column whose
analytic and FD values differ by > 1e-5 relative.
"""
from __future__ import annotations

import numpy as np

from semi.fem.sg_assembly import (
    assemble_sg_jacobian_correction_1d,
    assemble_sg_residual_1d,
)


def compute_galerkin_convdiff_residual_1d(
    psi_arr, n_arr, p_arr,
    cell_vertices, h_cells, v2d_psi, v2d_n, v2d_p,
    mu_n_hat, mu_p_hat, L0_sq,
):
    """
    Per-vertex Galerkin convection-diffusion residual on a 1D P1 mesh,
    matching the M13 form `inner(grad(*), grad(v)) -+ *grad(psi).grad(v)`
    (minus for n, plus for p).

    On a 1D cell K = (i, j) with hat-function basis v_i (1 at i, 0 at j,
    grad = -1/h on K), `grad(n) = (n_j - n_i)/h`, and lumped-mass
    n_avg = (n_i + n_j)/2. Element-by-element the contribution at vertex i
    (left of edge) integrated over K is:

        gal_n_i = -L0_sq * mu_n * ((n_j - n_i)/h - n_avg * dpsi / h)

    Mirror sign at vertex j. Same form for holes with the drift sign flipped.
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


def compute_correction_residual(
    psi_arr, n_arr, p_arr,
    cell_vertices, h_cells, v2d_psi, v2d_n, v2d_p,
    mu_n_hat, mu_p_hat, L0_sq,
):
    """The quantity that the SG callback ADDS to b_M13: `-sg_n - g_n` per block."""
    sg_n, sg_p = assemble_sg_residual_1d(
        psi_arr, n_arr, p_arr,
        cell_vertices, h_cells, v2d_psi, v2d_n, v2d_p,
        mu_n_hat, mu_p_hat, L0_sq,
    )
    g_n, g_p = compute_galerkin_convdiff_residual_1d(
        psi_arr, n_arr, p_arr,
        cell_vertices, h_cells, v2d_psi, v2d_n, v2d_p,
        mu_n_hat, mu_p_hat, L0_sq,
    )
    return -sg_n - g_n, -sg_p - g_p


def main() -> int:
    # 3-node mesh on [0, 5e-6, 1e-5] m -> 2 cells of length 5e-6.
    n_verts = 3
    cell_vertices = np.array([[0, 1], [1, 2]], dtype=np.int64)
    h_cells = np.array([5.0e-6, 5.0e-6], dtype=np.float64)
    v2d = np.arange(n_verts, dtype=np.int64)
    v2d_psi = v2d.copy()
    v2d_n = v2d.copy()
    v2d_p = v2d.copy()

    # Block-vector layout for the Jacobian: psi block, then n block, then p block.
    dof_offsets = (0, n_verts, 2 * n_verts)
    bc_dofs_n: set[int] = set()  # no BCs in the FD test
    bc_dofs_p: set[int] = set()

    # Physical-ish constants. Values are illustrative; the FD test is
    # algebraic — only the formulas in the assembly need to be self-consistent.
    mu_n_hat = 1.0
    mu_p_hat = 0.32
    L0_sq = 1.0e-12  # ~ (1 um)^2

    # Non-equilibrium state with non-trivial dpsi (depletion region scale).
    # Hand-tuned: psi varies over ~3 (scaled units, ~ 80 mV physical at V_t),
    # carrier densities ~1e16-1e18 across the junction.
    psi_arr = np.array([0.0, 1.5, 4.0], dtype=np.float64)
    n_arr = np.array([1.0e18, 5.0e16, 1.0e16], dtype=np.float64)
    p_arr = np.array([1.0e15, 5.0e16, 1.0e18], dtype=np.float64)

    print(f"State:")
    print(f"  psi = {psi_arr}")
    print(f"  n   = {n_arr}")
    print(f"  p   = {p_arr}")
    print(f"  dpsi (per cell) = "
          f"{psi_arr[1] - psi_arr[0]:.3f}, {psi_arr[2] - psi_arr[1]:.3f}")
    print()

    # 1. Analytic correction.
    rows, cols, vals = assemble_sg_jacobian_correction_1d(
        psi_arr, n_arr, p_arr,
        cell_vertices, h_cells, v2d_psi, v2d_n, v2d_p,
        mu_n_hat, mu_p_hat, L0_sq,
        dof_offsets, bc_dofs_n, bc_dofs_p,
    )
    n_total = 3 * n_verts  # 9 dofs total in block layout
    J_analytic = np.zeros((n_total, n_total), dtype=np.float64)
    for r, c, v in zip(rows, cols, vals):
        J_analytic[int(r), int(c)] += float(v)

    # 2. Finite-difference Jacobian of the correction residual `-sg - g`.
    # Use central differences. Step size scaled per variable: psi is O(1),
    # carrier densities O(1e18). Use relative perturbation eps_rel = 1e-7
    # for each variable (so the absolute step adapts to the variable scale).
    R0_n, R0_p = compute_correction_residual(
        psi_arr, n_arr, p_arr,
        cell_vertices, h_cells, v2d_psi, v2d_n, v2d_p,
        mu_n_hat, mu_p_hat, L0_sq,
    )
    R0 = np.zeros(n_total, dtype=np.float64)
    R0[dof_offsets[1]:dof_offsets[2]] = R0_n
    R0[dof_offsets[2]:] = R0_p
    # psi rows of the correction are identically zero (the correction only
    # touches the n and p continuity blocks).

    eps_rel = 1.0e-7
    J_fd = np.zeros((n_total, n_total), dtype=np.float64)
    for col in range(n_total):
        # Decide which array to perturb based on the block layout.
        block_psi = (0, n_verts)
        block_n = (n_verts, 2 * n_verts)
        block_p = (2 * n_verts, 3 * n_verts)

        psi_p = psi_arr.copy(); psi_m = psi_arr.copy()
        n_p = n_arr.copy(); n_m = n_arr.copy()
        p_p = p_arr.copy(); p_m = p_arr.copy()

        if block_psi[0] <= col < block_psi[1]:
            v = col - block_psi[0]
            scale = max(abs(psi_arr[v]), 1.0)
            h = eps_rel * scale
            psi_p[v] += h
            psi_m[v] -= h
        elif block_n[0] <= col < block_n[1]:
            v = col - block_n[0]
            scale = max(abs(n_arr[v]), 1.0)
            h = eps_rel * scale
            n_p[v] += h
            n_m[v] -= h
        else:
            v = col - block_p[0]
            scale = max(abs(p_arr[v]), 1.0)
            h = eps_rel * scale
            p_p[v] += h
            p_m[v] -= h

        Rp_n, Rp_p = compute_correction_residual(
            psi_p, n_p, p_p,
            cell_vertices, h_cells, v2d_psi, v2d_n, v2d_p,
            mu_n_hat, mu_p_hat, L0_sq,
        )
        Rm_n, Rm_p = compute_correction_residual(
            psi_m, n_m, p_m,
            cell_vertices, h_cells, v2d_psi, v2d_n, v2d_p,
            mu_n_hat, mu_p_hat, L0_sq,
        )
        Rp = np.zeros(n_total)
        Rp[dof_offsets[1]:dof_offsets[2]] = Rp_n
        Rp[dof_offsets[2]:] = Rp_p
        Rm = np.zeros(n_total)
        Rm[dof_offsets[1]:dof_offsets[2]] = Rm_n
        Rm[dof_offsets[2]:] = Rm_p

        J_fd[:, col] = (Rp - Rm) / (2.0 * h)

    # 3. Compare entry-by-entry. Only the n,p row blocks (3..8) have entries.
    block_labels = []
    for v in range(n_verts):
        block_labels.append(f"psi_{v}")
    for v in range(n_verts):
        block_labels.append(f"n_{v}")
    for v in range(n_verts):
        block_labels.append(f"p_{v}")

    rel_thresh = 1.0e-5
    abs_floor = 1.0e-12  # below this we consider the entry "zero"
    mismatches = []
    for r in range(n_verts, 3 * n_verts):  # skip psi rows
        for c in range(3 * n_verts):
            j_a = J_analytic[r, c]
            j_f = J_fd[r, c]
            denom = max(abs(j_a), abs(j_f), abs_floor)
            rel = abs(j_a - j_f) / denom
            if rel > rel_thresh and max(abs(j_a), abs(j_f)) > abs_floor:
                mismatches.append((r, c, j_a, j_f, rel))

    print(f"Analytic Jacobian correction (rows {n_verts}..{3*n_verts-1}, "
          f"all columns):")
    print(f"  shape (n,p rows) = ({2*n_verts}, {3*n_verts})")
    print(f"  Frobenius norm  = {np.linalg.norm(J_analytic[n_verts:]):.6e}")
    print(f"  FD Frobenius    = {np.linalg.norm(J_fd[n_verts:]):.6e}")
    print(f"  Frobenius diff  = "
          f"{np.linalg.norm(J_analytic[n_verts:] - J_fd[n_verts:]):.6e}")
    print()

    if not mismatches:
        print(f"FD verification CLEAN. All entries with magnitude > {abs_floor}"
              f" agree to relative tol {rel_thresh}.")
        # Print a few sample entries for documentation.
        print()
        print("Sample entries (analytic vs FD):")
        for r, c in [(3, 0), (3, 3), (3, 4), (4, 1), (4, 4), (5, 5),
                     (6, 0), (6, 6), (6, 7), (7, 1), (7, 7), (8, 8)]:
            print(f"  J[{block_labels[r]:>5} , {block_labels[c]:>5}] = "
                  f"{J_analytic[r, c]:+.6e}  vs FD {J_fd[r, c]:+.6e}")
        return 0

    print(f"FD MISMATCHES: {len(mismatches)} entries exceed rel tol {rel_thresh}")
    print()
    print(f"  {'row':>6}  {'col':>6}  {'analytic':>15}  {'FD':>15}  {'rel_err':>10}")
    for r, c, j_a, j_f, rel in mismatches:
        print(f"  {block_labels[r]:>6}  {block_labels[c]:>6}  "
              f"{j_a:+.6e}  {j_f:+.6e}  {rel:.3e}")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
