"""
Per-edge Scharfetter-Gummel residual assembly (M13.1 Path B1).

Why this lives outside UFL
--------------------------
UFL's symbolic differentiator and FFCx's kernel compiler do not handle
the SG `B(x) = x / (exp(x) - 1)` Bernoulli function reliably across the
full Peclet range; both the conditional and tanh-blended UFL Bernoulli
implementations produced NaN or line-search failures in the transient
runner once carriers entered the strong-injection regime around
t ~ tau_p. This module assembles the SG flux per edge in pure NumPy
on the DOF arrays, scatters into the residual vector, and lets PETSc's
SNES drive the Newton outer loop.

The Poisson, mass, recombination, and history blocks of the residual
remain in UFL (assembled via dolfinx). Only the convection-diffusion
blocks are replaced by per-edge SG.

SG flux (Sandia / Farrell convention, ADR 0012)
-----------------------------------------------
For each interior edge `e_ij` of the 1D mesh with vertices i, j (j on
the right, mesh oriented +x), edge length `h_ij`, and scaled potential
difference `dpsi = psi_j - psi_i`:

    F_n_ij = (mu_n_hat / h_ij) * ( n_j * B(+dpsi) - n_i * B(-dpsi) )
    F_p_ij = (mu_p_hat / h_ij) * ( p_i * B(+dpsi) - p_j * B(-dpsi) )

(All quantities scaled per the project's standard nondimensionalisation;
see docs/PHYSICS.md and `semi.scaling`.)

The residual contribution at vertex i from a single edge `e_ij`:
- For F_n: divergence of conventional current: residual at i gets
  `+L0_sq * F_n_ij`, residual at j gets `-L0_sq * F_n_ij`
  (factor ensures consistency with the UFL form's L0_sq prefix on the
  spatial operator).
- For F_p: same pattern with mu_p and the hole formula above.

This module exposes `assemble_sg_residual_1d` which takes the current
state vector and writes the per-vertex SG contributions into a PETSc
residual Vec.
"""
from __future__ import annotations

import numpy as np

from .scharfetter_gummel import bernoulli_array


def compute_edge_topology_1d(msh) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Return per-cell vertex pairs and lengths for a 1D dolfinx mesh.

    Returns
    -------
    cell_vertices : (n_cells, 2) int array
        cell_vertices[k] = (i, j) where i is the left vertex DOF and
        j is the right vertex DOF for cell k. Sorted by x so that
        psi_j - psi_i has the conventional sign.
    h_cells : (n_cells,) float array
        Cell lengths.
    cell_x : (n_cells, 2) float array
        Cell vertex x-coordinates.

    Notes
    -----
    For 1D dolfinx P1 meshes, vertex indices in cell connectivity are
    NOT guaranteed to be sorted by spatial position. We re-sort each
    cell's two vertices by their x-coordinate so the SG flux convention
    `dpsi = psi_j - psi_i` (j on the right) is consistent across the
    mesh.
    """
    tdim = msh.topology.dim
    if tdim != 1:
        raise NotImplementedError(
            "compute_edge_topology_1d only supports 1D meshes; "
            f"got tdim={tdim}."
        )
    msh.topology.create_connectivity(tdim, 0)
    c2v = msh.topology.connectivity(tdim, 0)
    n_cells = msh.topology.index_map(tdim).size_local

    coords = msh.geometry.x  # (n_vertices, 3)
    cell_vertices = np.zeros((n_cells, 2), dtype=np.int64)
    h_cells = np.zeros(n_cells, dtype=np.float64)
    cell_x = np.zeros((n_cells, 2), dtype=np.float64)
    for k in range(n_cells):
        verts = c2v.links(k)
        x0 = coords[verts[0], 0]
        x1 = coords[verts[1], 0]
        if x0 <= x1:
            cell_vertices[k] = (verts[0], verts[1])
            cell_x[k] = (x0, x1)
        else:
            cell_vertices[k] = (verts[1], verts[0])
            cell_x[k] = (x1, x0)
        h_cells[k] = abs(x1 - x0)
    return cell_vertices, h_cells, cell_x


def vertex_to_dof_map(V) -> np.ndarray:
    """
    Map from vertex index to (single-component) DOF index for a P1
    Lagrange function space on a 1D mesh.

    For a P1 space on a 1D dolfinx mesh, each vertex carries exactly
    one DOF. We need this mapping to scatter per-edge contributions
    into the right entries of the assembled residual vector.

    Returns
    -------
    v2d : (n_vertices,) int array
        v2d[v] = dof index of the DOF at vertex v.
    """
    msh = V.mesh
    tdim = msh.topology.dim
    msh.topology.create_connectivity(0, tdim)
    n_verts = msh.topology.index_map(0).size_local + msh.topology.index_map(0).num_ghosts
    # locate_dofs_topological returns the DOF that lies on each vertex.
    # For P1 in 1D with one component, this is a 1-1 mapping.
    from dolfinx.fem import locate_dofs_topological
    v2d = np.zeros(n_verts, dtype=np.int64)
    for v in range(n_verts):
        d = locate_dofs_topological(V, 0, np.array([v], dtype=np.int32))
        if len(d) == 0:
            v2d[v] = -1
        else:
            v2d[v] = int(d[0])
    return v2d


def solve_sg_block_1d(  # pragma: no cover - WIP integration, see ADR 0012
    F_list,
    u_list,
    bcs,
    prefix: str,
    *,
    msh,
    spaces,
    sc,
    mu_n_hat: float,
    mu_p_hat: float,
    petsc_options: dict | None = None,
):
    """
    Custom block SNES solve for the 1D SG transient (M13.1 Path B1).

    Drives PETSc SNES directly via petsc4py with custom F and J
    callbacks. The UFL forms in `F_list` provide the mass + recombination
    + history + Poisson contributions; the SG convection-diffusion
    contribution is added per-edge in NumPy via
    `assemble_sg_residual_1d` and scattered into the residual vector.

    Jacobian: assembled by PETSc finite-difference coloring
    (-snes_fd_color), since the SG block's analytic Jacobian is
    non-trivial to derive. The UFL parts of the Jacobian could be
    auto-derived but mixing analytic+FD in one Jacobian is fragile;
    the cost of full FD-color is acceptable for M13.1's 1D mesh sizes
    (under 1k DOFs the FD assembly is microseconds per Newton step).

    Parameters
    ----------
    F_list : sequence of ufl.Form
        [F_psi, F_n_no_convdiff, F_p_no_convdiff] residual forms.
    u_list : sequence of dolfinx.fem.Function
        [psi, n_hat, p_hat] unknowns. Updated in place.
    bcs : list of dolfinx.fem.DirichletBC
        Flat list of Dirichlet BCs.
    prefix : str
        PETSc options prefix.
    msh : dolfinx.mesh.Mesh
        The 1D mesh.
    spaces : DDBlockSpaces
    sc : Scaling
    mu_n_hat, mu_p_hat : float
        Scaled mobilities.
    petsc_options : dict, optional
        SNES tolerances (snes_rtol, snes_atol, snes_stol, snes_max_it).

    Returns
    -------
    dict
        {'iterations', 'reason', 'converged', 'problem': None}
    """
    import functools

    from dolfinx.fem import form as fem_form
    from dolfinx.fem.petsc import (
        assemble_jacobian,
        assemble_residual,
        create_matrix,
        create_vector,
    )
    from petsc4py import PETSc
    from ufl import derivative as ufl_derivative

    cell_v, h_cells, _ = compute_edge_topology_1d(msh)
    v2d_psi = vertex_to_dof_map(spaces.V_psi)
    v2d_n = vertex_to_dof_map(spaces.V_phi_n)
    v2d_p = vertex_to_dof_map(spaces.V_phi_p)
    L0_sq = sc.L0 ** 2

    psi, n_hat, p_hat = u_list
    F_psi, F_n, F_p = F_list

    # Compile residual forms (no convection-diffusion) and Jacobian
    # forms (UFL auto-derived). Block ordering: [psi, n, p].
    F_compiled = [fem_form(F_psi), fem_form(F_n), fem_form(F_p)]
    J_compiled = [
        [fem_form(ufl_derivative(F_psi, psi)),  fem_form(ufl_derivative(F_psi, n_hat)),  fem_form(ufl_derivative(F_psi, p_hat))],
        [fem_form(ufl_derivative(F_n, psi)),    fem_form(ufl_derivative(F_n, n_hat)),    fem_form(ufl_derivative(F_n, p_hat))],
        [fem_form(ufl_derivative(F_p, psi)),    fem_form(ufl_derivative(F_p, n_hat)),    fem_form(ufl_derivative(F_p, p_hat))],
    ]

    # Allocate block residual and Jacobian via dolfinx 0.10's
    # block-aware create_matrix / create_vector.
    A = create_matrix(J_compiled, kind=PETSc.Mat.Type.AIJ)
    b = create_vector(F_compiled, kind="mpi")
    x = b.duplicate()

    # Map (block index, local dof) -> global block-vector position.
    # dolfinx's block layout interleaves owned blocks: all psi DOFs,
    # then all n DOFs, then all p DOFs (within each rank).
    n_psi_dofs = spaces.V_psi.dofmap.index_map.size_local
    n_n_dofs = spaces.V_phi_n.dofmap.index_map.size_local
    n_p_dofs = spaces.V_phi_p.dofmap.index_map.size_local

    # Block-vector slice offsets (assumes single MPI rank for simplicity;
    # M13.1 1D test runs on a single rank).
    off_psi = 0
    off_n = n_psi_dofs
    off_p = n_psi_dofs + n_n_dofs

    def _x_to_functions(x_arr):
        """Update psi, n_hat, p_hat in place from a block-layout array."""
        psi.x.array[:n_psi_dofs] = x_arr[off_psi : off_psi + n_psi_dofs]
        n_hat.x.array[:n_n_dofs] = x_arr[off_n : off_n + n_n_dofs]
        p_hat.x.array[:n_p_dofs] = x_arr[off_p : off_p + n_p_dofs]
        psi.x.scatter_forward()
        n_hat.x.scatter_forward()
        p_hat.x.scatter_forward()

    def _functions_to_x(x_vec):
        x_arr = x_vec.getArray(readonly=False)
        x_arr[off_psi : off_psi + n_psi_dofs] = psi.x.array[:n_psi_dofs]
        x_arr[off_n : off_n + n_n_dofs] = n_hat.x.array[:n_n_dofs]
        x_arr[off_p : off_p + n_p_dofs] = p_hat.x.array[:n_p_dofs]
        x_vec.restoreArray(x_arr)

    # Build the standard dolfinx residual / Jacobian callbacks from the
    # UFL forms (mass + recombination + history + Poisson). We then
    # post-hook the SG contribution into the residual.
    _ufl_residual_callable = functools.partial(
        assemble_residual, u_list, F_compiled, J_compiled, bcs,
    )
    _ufl_jacobian_callable = functools.partial(
        assemble_jacobian, u_list, J_compiled, [], bcs,
    )

    def F_callback(snes, x_vec, b_vec):
        # Step 1: dolfinx assembles the UFL block residual into b_vec.
        # This call also updates the Function objects from x_vec, applies
        # BC lifting, and zeroes BC rows.
        _ufl_residual_callable(snes, x_vec, b_vec)
        # Step 2: add SG convection-diffusion contribution per edge.
        psi_full = psi.x.array
        n_full = n_hat.x.array
        p_full = p_hat.x.array
        R_n_sg, R_p_sg = assemble_sg_residual_1d(
            psi_full, n_full, p_full,
            cell_v, h_cells, v2d_psi, v2d_n, v2d_p,
            mu_n_hat, mu_p_hat, L0_sq,
        )
        b_arr = b_vec.getArray(readonly=False)
        # ZERO the SG contribution at BC rows (those rows are pinned
        # by the Dirichlet BC and must not be modified).
        bc_dofs_n = set()
        bc_dofs_p = set()
        for bc in bcs:
            if bc.function_space is spaces.V_phi_n:
                for d in bc.dof_indices()[0]:
                    bc_dofs_n.add(int(d))
            elif bc.function_space is spaces.V_phi_p:
                for d in bc.dof_indices()[0]:
                    bc_dofs_p.add(int(d))
        for d in range(n_n_dofs):
            if d not in bc_dofs_n:
                b_arr[off_n + d] += R_n_sg[d]
        for d in range(n_p_dofs):
            if d not in bc_dofs_p:
                b_arr[off_p + d] += R_p_sg[d]
        b_vec.restoreArray(b_arr)
        b_vec.assemble()

    def J_callback(snes, x_vec, J_mat, P_mat):
        # UFL Jacobian assembly. The SG block's Jacobian is added
        # implicitly via PETSc -snes_fd_color (set in petsc options
        # below): SNES will FD-perturb F_callback to fill the Jacobian
        # entries that the analytic UFL forms miss.
        _ufl_jacobian_callable(snes, x_vec, J_mat, P_mat)

    snes = PETSc.SNES().create(msh.comm)
    snes.setOptionsPrefix(prefix)

    # Push options into the prefix-namespaced PETSc options database.
    opts = PETSc.Options()
    opts_dict = {
        "snes_type": "newtonls",
        "snes_linesearch_type": "bt",
        "snes_monitor": None,
        "ksp_type": "preonly",
        "pc_type": "lu",
        "pc_factor_mat_solver_type": "mumps",
        "snes_fd_color": None,  # FD-color Jacobian for the SG block (B1)
        "snes_fd_color_use_mat": None,  # use Jacobian's sparsity for coloring
    }
    if petsc_options:
        opts_dict.update({k: v for k, v in petsc_options.items() if v is not None or k.startswith("snes_")})
    for k, v in opts_dict.items():
        if v is None:
            opts[f"{prefix}{k}"] = ""
        else:
            opts[f"{prefix}{k}"] = v

    snes.setFunction(F_callback, b)
    snes.setJacobian(J_callback, A, A)
    snes.setFromOptions()

    _functions_to_x(x)
    snes.solve(None, x)
    reason = snes.getConvergedReason()
    n_iter = snes.getIterationNumber()
    converged = reason > 0
    # Read solution back into Function objects.
    _x_to_functions(x.getArray(readonly=True).copy())

    # Cleanup options to avoid polluting other solves.
    for k in opts_dict.keys():
        full_key = f"{prefix}{k}"
        if full_key in opts:
            opts.delValue(full_key)

    snes.destroy()
    A.destroy()
    b.destroy()
    x.destroy()

    return {
        "iterations": int(n_iter),
        "reason": int(reason),
        "converged": bool(converged),
        "problem": None,
    }


def assemble_sg_residual_1d(
    psi_arr: np.ndarray,
    n_arr: np.ndarray,
    p_arr: np.ndarray,
    cell_vertices: np.ndarray,
    h_cells: np.ndarray,
    v2d_psi: np.ndarray,
    v2d_n: np.ndarray,
    v2d_p: np.ndarray,
    mu_n_hat: float,
    mu_p_hat: float,
    L0_sq: float,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute per-vertex SG contributions to the n and p continuity
    residuals on a 1D mesh.

    Parameters
    ----------
    psi_arr, n_arr, p_arr : (n_dofs,) float64 arrays
        Current values of the primary unknowns at all DOFs (interior +
        boundary).
    cell_vertices : (n_cells, 2) int64
    h_cells : (n_cells,) float64
    v2d_psi, v2d_n, v2d_p : (n_vertices,) int64
        Vertex-to-DOF maps for the three function spaces.
    mu_n_hat, mu_p_hat : float
        Scaled mobilities (mu / mu0).
    L0_sq : float
        Squared characteristic length, m^2 (matches the L0_sq prefix
        on the UFL spatial form).

    Returns
    -------
    R_n_sg : (n_dofs,) float64
        Per-DOF SG contribution to the electron continuity residual.
        ADD this to the UFL-assembled mass + recombination + history
        residual to get the total F_n.
    R_p_sg : (n_dofs,) float64
        Same for the hole continuity residual.
    """
    R_n = np.zeros_like(n_arr)
    R_p = np.zeros_like(p_arr)

    # Per-cell quantities, vectorised.
    v_i = cell_vertices[:, 0]
    v_j = cell_vertices[:, 1]
    psi_i = psi_arr[v2d_psi[v_i]]
    psi_j = psi_arr[v2d_psi[v_j]]
    n_i = n_arr[v2d_n[v_i]]
    n_j = n_arr[v2d_n[v_j]]
    p_i = p_arr[v2d_p[v_i]]
    p_j = p_arr[v2d_p[v_j]]
    h = h_cells

    dpsi = psi_j - psi_i
    B_pos = bernoulli_array(dpsi)
    B_neg = bernoulli_array(-dpsi)

    # SG flux scalars per cell (Sandia / Farrell convention).
    F_n = (mu_n_hat / h) * (n_j * B_pos - n_i * B_neg)
    F_p = (mu_p_hat / h) * (p_i * B_pos - p_j * B_neg)

    # Scatter contributions: residual at vertex i gets +L0_sq * F,
    # residual at vertex j gets -L0_sq * F. (Divergence of the flux
    # over the dual control volume, summing the per-edge contributions
    # to the per-vertex residual.)
    np.add.at(R_n, v2d_n[v_i],  L0_sq * F_n)
    np.add.at(R_n, v2d_n[v_j], -L0_sq * F_n)
    np.add.at(R_p, v2d_p[v_i],  L0_sq * F_p)
    np.add.at(R_p, v2d_p[v_j], -L0_sq * F_p)

    return R_n, R_p
