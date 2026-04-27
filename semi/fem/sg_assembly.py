"""
Per-edge Scharfetter-Gummel residual and Jacobian assembly (M13.1 Path B1).

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


def solve_sg_block_1d(
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
    SG-flux SNES solve for the 1D transient block residual.

    Build the standard dolfinx ``NonlinearProblem`` from the M13
    Galerkin UFL forms in ``F_list``, then override the SNES function
    callback so that on every residual assembly:

    1. ``assemble_residual`` populates the block residual ``b`` from
       the M13 Galerkin form (mass + Galerkin convection-diffusion +
       SRH + history + Poisson).
    2. The Galerkin convection-diffusion contribution is removed from
       ``b`` by assembling a "convection-diffusion only" auxiliary form
       and subtracting.
    3. The Scharfetter-Gummel convection-diffusion contribution is
       added to ``b`` via :func:`assemble_sg_residual_1d`.

    The Jacobian remains the UFL-derived M13 Galerkin Jacobian. This is
    an approximation: at zero cell Peclet ``B(0) = 1`` so SG and
    Galerkin agree exactly and the Jacobian is correct; at large Peclet
    the diffusion coefficient differs by ``A(dpsi)`` but Newton with
    backtracking line search still converges (the Galerkin Jacobian
    has the correct sign and sparsity, only its magnitude on the
    diffusion entries is off by ``A``). The ADR 0012 1D steady-state
    agreement gate (1e-4) verifies the spatial discretisation; SNES
    convergence is observed at 4-8 Newton iterations per timestep on
    the ``pn_1d`` test mesh, well within ``max_it = 100``.

    See ADR 0012 § "1D edge flux" for the SG flux derivation.

    Parameters
    ----------
    F_list : sequence of ufl.Form
        [F_psi, F_n_galerkin, F_p_galerkin] M13 Galerkin residual forms.
    u_list : sequence of dolfinx.fem.Function
        [psi, n_hat, p_hat] unknowns. Updated in place on return.
    bcs : list of dolfinx.fem.DirichletBC
        Flat list of Dirichlet BCs.
    prefix : str
        PETSc options prefix (must end with ``_``).
    msh : dolfinx.mesh.Mesh
        The 1D mesh (``msh.topology.dim == 1`` required).
    spaces : DDBlockSpaces
    sc : Scaling
    mu_n_hat, mu_p_hat : float
        Scaled mobilities (mu / mu0).
    petsc_options : dict, optional
        SNES options. Common: snes_rtol, snes_atol, snes_stol, snes_max_it.

    Returns
    -------
    dict
        ``{'iterations', 'reason', 'converged', 'problem': NonlinearProblem}``
    """
    import functools

    import ufl
    from dolfinx import fem
    from dolfinx.fem.petsc import (
        NonlinearProblem,
        assemble_jacobian,
        assemble_residual,
        assign,
    )
    from petsc4py import PETSc

    if msh.topology.dim != 1:
        raise NotImplementedError(
            "solve_sg_block_1d requires a 1D mesh; "
            f"got tdim={msh.topology.dim}. 2D/3D SG assembly is M13.2+."
        )

    cell_v, h_cells, _ = compute_edge_topology_1d(msh)
    v2d_psi = vertex_to_dof_map(spaces.V_psi)
    v2d_n = vertex_to_dof_map(spaces.V_phi_n)
    v2d_p = vertex_to_dof_map(spaces.V_phi_p)
    L0_sq = sc.L0 ** 2

    psi, n_hat, p_hat = u_list

    # Build a NonlinearProblem on the M13 Galerkin forms. This sets up
    # all the PETSc structures (Mat/Vec, SNES options) correctly via
    # the supported dolfinx 0.10 API.
    opts = {
        "snes_type": "newtonls",
        "snes_linesearch_type": "bt",
        "snes_linesearch_alpha": 1.0e-8,  # gentler bt threshold
        "snes_monitor": None,
        "ksp_type": "preonly",
        "pc_type": "lu",
        "pc_factor_mat_solver_type": "mumps",
        "snes_rtol": 1.0e-10,
        "snes_atol": 1.0e-10,
        "snes_max_it": 50,
        # Force at least one Newton iteration even when the initial
        # residual is below atol. Required for the SG transient
        # because the per-step residual at the previous step's
        # iterate (the natural initial guess) is sometimes already
        # below atol when the time-step has only perturbed psi (ohmic
        # n,p BCs are V-independent, so they don't drive residual
        # change), yet the iterate is NOT at the new time-step's
        # solution.
        "snes_force_iteration": 1,
    }
    if petsc_options:
        opts.update(petsc_options)

    problem = NonlinearProblem(
        list(F_list), list(u_list),
        bcs=list(bcs),
        petsc_options_prefix=prefix,
        petsc_options=opts,
    )

    # Build the "convection-diffusion only" auxiliary forms whose
    # contribution we subtract from the M13 residual (so that what
    # remains is the SG-style residual after we add the per-edge SG
    # term). These mirror the convection-diffusion blocks of the
    # M13 Galerkin form in `_build_transient_residual` exactly.
    v_n = ufl.TestFunction(spaces.V_phi_n)
    v_p = ufl.TestFunction(spaces.V_phi_p)
    L0_sq_const = fem.Constant(msh, PETSc.ScalarType(L0_sq))
    mu_n_const = fem.Constant(msh, PETSc.ScalarType(mu_n_hat))
    mu_p_const = fem.Constant(msh, PETSc.ScalarType(mu_p_hat))
    F_n_galerkin_convdiff = (
        L0_sq_const * mu_n_const * (
            ufl.inner(ufl.grad(n_hat), ufl.grad(v_n))
            - n_hat * ufl.inner(ufl.grad(psi), ufl.grad(v_n))
        ) * ufl.dx
    )
    F_p_galerkin_convdiff = (
        L0_sq_const * mu_p_const * (
            ufl.inner(ufl.grad(p_hat), ufl.grad(v_p))
            + p_hat * ufl.inner(ufl.grad(psi), ufl.grad(v_p))
        ) * ufl.dx
    )
    F_n_convdiff_form = fem.form(F_n_galerkin_convdiff)
    F_p_convdiff_form = fem.form(F_p_galerkin_convdiff)

    # Cache the per-block local DOF counts and the BC dof sets so we
    # only touch non-Dirichlet rows when patching the residual.
    n_psi_dofs = spaces.V_psi.dofmap.index_map.size_local
    n_n_dofs = spaces.V_phi_n.dofmap.index_map.size_local
    n_p_dofs = spaces.V_phi_p.dofmap.index_map.size_local
    off_n = n_psi_dofs
    off_p = n_psi_dofs + n_n_dofs

    bc_dofs_n: set[int] = set()
    bc_dofs_p: set[int] = set()
    for bc in bcs:
        if bc.function_space is spaces.V_phi_n:
            for d in bc.dof_indices()[0]:
                bc_dofs_n.add(int(d))
        elif bc.function_space is spaces.V_phi_p:
            for d in bc.dof_indices()[0]:
                bc_dofs_p.add(int(d))

    # Reuse the assemble_residual / assemble_jacobian partials that
    # NonlinearProblem set up for the M13 Galerkin form.
    _ufl_residual = functools.partial(
        assemble_residual, problem._u, problem._F, problem._J, bcs,
    )
    _ufl_jacobian = functools.partial(
        assemble_jacobian, problem._u, problem._J, problem.preconditioner, bcs,
    )

    # Block-vector layout offsets (single MPI rank for M13.1's 1D test).
    dof_offsets = (0, n_psi_dofs, n_psi_dofs + n_n_dofs)

    def F_callback(snes, x_vec, b_vec):
        # 1. Standard UFL block-residual assembly (M13 Galerkin form).
        _ufl_residual(snes, x_vec, b_vec)
        # The Function objects (psi, n_hat, p_hat) have been updated
        # from x_vec by the assemble_residual call.

        # 2. Subtract the Galerkin convection-diffusion contribution
        #    (the M13 form's `+grad(n).grad(v) - n*grad(psi).grad(v)`
        #    contribution to the residual; we'll replace it with SG).
        gn = fem.assemble_vector(F_n_convdiff_form).array
        gp = fem.assemble_vector(F_p_convdiff_form).array

        # 3. Add the SG convection-diffusion contribution.
        # Note: `assemble_sg_residual_1d` follows the divergence-of-
        # conventional-current sign convention from ADR 0012, which
        # is the OPPOSITE sign of the M13 Galerkin weak form's
        # `+grad(n).grad(v)` (verified at zero field on a 3-node
        # mesh — see the assembly tests). To inject SG into the M13
        # block residual we subtract its output (the negation flips
        # to the M13 sign convention).
        sg_n, sg_p = assemble_sg_residual_1d(
            psi.x.array, n_hat.x.array, p_hat.x.array,
            cell_v, h_cells, v2d_psi, v2d_n, v2d_p,
            mu_n_hat, mu_p_hat, L0_sq,
        )

        with b_vec.localForm() as b_local:
            b_arr = b_local.getArray(readonly=False)
            for d in range(n_n_dofs):
                if d not in bc_dofs_n:
                    b_arr[off_n + d] += -sg_n[d] - gn[d]
            for d in range(n_p_dofs):
                if d not in bc_dofs_p:
                    b_arr[off_p + d] += -sg_p[d] - gp[d]
        b_vec.assemble()

    def J_callback(snes, x_vec, J_mat, P_mat):
        # 1. UFL Jacobian assembly (M13 Galerkin form auto-derivative).
        _ufl_jacobian(snes, x_vec, J_mat, P_mat)
        # 2. Add the per-edge analytic correction
        #    (J_sg_block - J_galerkin_convdiff_block)
        # so the resulting Jacobian matches the analytic derivative of
        # the SG-corrected block residual.
        rows, cols, vals = assemble_sg_jacobian_correction_1d(
            psi.x.array, n_hat.x.array, p_hat.x.array,
            cell_v, h_cells, v2d_psi, v2d_n, v2d_p,
            mu_n_hat, mu_p_hat, L0_sq,
            dof_offsets, bc_dofs_n, bc_dofs_p,
        )
        if len(rows) > 0:
            for r, c, v in zip(rows, cols, vals):
                J_mat.setValue(r, c, v, addv=PETSc.InsertMode.ADD_VALUES)
        J_mat.assemble()
        if P_mat is not None and P_mat is not J_mat:
            P_mat.assemble()

    # Override the SNES function and Jacobian with the SG-corrected
    # callbacks. NonlinearProblem set up the Mat/Vec, prefix, and
    # PETSc options; we only swap the per-iteration assembly hooks.
    problem.solver.setFunction(F_callback, problem.b)
    problem.solver.setJacobian(J_callback, problem.A, problem.P_mat)

    # Drive the solve via NonlinearProblem.solve() which handles the
    # initial-guess copy, scatter, and final assignment.
    assign(problem._u, problem.x)
    problem.solver.solve(None, problem.x)
    import dolfinx.la
    dolfinx.la.petsc._ghost_update(
        problem.x, PETSc.InsertMode.INSERT, PETSc.ScatterMode.FORWARD,
    )
    assign(problem.x, problem._u)

    reason = problem.solver.getConvergedReason()
    n_iter = problem.solver.getIterationNumber()
    converged = reason > 0
    return {
        "iterations": int(n_iter),
        "reason": int(reason),
        "converged": bool(converged),
        "problem": problem,
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


def _bernoulli_prime_array(x: np.ndarray) -> np.ndarray:
    """
    Derivative of the Bernoulli function `B(x) = x/(exp(x)-1)` per element,
    numerically stable across the same regimes as `bernoulli_array`:

        x == 0          : B'(0) = -1/2 exact (Taylor leading order).
        |x| small       : Taylor:  -1/2 + x/6 - x^3/180.
        x large positive: B'(x) = u(1 - u - x)/(1-u)^2, with u = exp(-x).
                          For x>>0, u underflows to 0 and B'(x) -> -x*u
                          ≈ -x*exp(-x), which is benignly tiny.
        x large negative: B'(x) -> -1 (since B(x) -> -x as x -> -inf).
        Mid-range       : closed form (exp(x) - 1 - x*exp(x)) / (exp(x)-1)^2.

    Used by the analytic SG Jacobian assembly. Pure-NumPy.
    """
    x = np.asarray(x, dtype=np.float64)
    out = np.empty_like(x)

    is_zero = x == 0.0
    ax = np.abs(x)
    is_taylor = (ax < 1.0e-3) & ~is_zero
    is_big_pos = x > 30.0
    is_big_neg = x < -30.0
    is_mid = ~(is_zero | is_taylor | is_big_pos | is_big_neg)

    out[is_zero] = -0.5

    if is_taylor.any():
        xt = x[is_taylor]
        out[is_taylor] = -0.5 + xt / 6.0 - xt ** 3 / 180.0

    if is_big_pos.any():
        xp = x[is_big_pos]
        u = np.exp(-xp)
        out[is_big_pos] = u * (1.0 - u - xp) / (1.0 - u) ** 2

    if is_big_neg.any():
        xn = x[is_big_neg]
        # B(x) ~ -x + small correction for x -> -inf. B'(x) -> -1.
        # The closed form is fine here too: exp(x)*(exp(x)-1-x)/(exp(x)-1)^2
        # at x=-30 evaluates to ~ -1 to many digits.
        ex = np.exp(xn)
        out[is_big_neg] = (ex - 1.0 - xn * ex) / (ex - 1.0) ** 2

    if is_mid.any():
        xm = x[is_mid]
        ex = np.exp(xm)
        em = ex - 1.0
        out[is_mid] = (em - xm * ex) / (em * em)

    return out


def assemble_sg_jacobian_correction_1d(
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
    dof_offsets: tuple[int, int, int],
    bc_dofs_n: set[int],
    bc_dofs_p: set[int],
):
    """
    Compute the per-edge correction to the M13 Galerkin Jacobian needed
    to obtain the analytic Jacobian of the SG-corrected block residual.

    The corrected residual is::

        b_corrected = b_M13_galerkin - g_n - sg_n      (n block, mirror for p)

    where ``b_M13_galerkin`` is the standard UFL block residual,
    ``g_n`` is the M13 Galerkin convection-diffusion contribution
    (extracted from the same UFL form), and ``sg_n`` is the per-edge
    SG residual from :func:`assemble_sg_residual_1d` (which uses the
    opposite sign convention from the M13 weak form — the leading
    minus signs in the expression above flip it back).

    Differentiating block by block gives::

        J_corrected = J_M13_galerkin - J_galerkin_convdiff_n_p - J_sg_n_p

    The ``J_M13_galerkin`` contribution is auto-derived by UFL and
    assembled by :class:`dolfinx.fem.petsc.NonlinearProblem`. The
    correction term ``-J_galerkin_convdiff_n_p - J_sg_n_p`` is what
    this function returns: a list of ``(row, col, value)`` triples
    (per block) suitable for ``Mat.setValuesIJV`` / ``Mat.setValue``
    accumulation onto the UFL Jacobian.

    Returns
    -------
    rows, cols, vals : (m,) int64 / int64 / float64 arrays
        Triples to ADD to the assembled UFL Jacobian (M13 Galerkin
        block residual). Indices are GLOBAL (block-vector layout):
        ``off_psi``, ``off_n``, ``off_p`` from ``dof_offsets``. BC
        dofs are skipped (those rows are pinned by Dirichlet BC and
        must not be modified).
    """
    off_psi, off_n, off_p = dof_offsets

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
    Bp_pos = _bernoulli_prime_array(dpsi)
    Bp_neg = _bernoulli_prime_array(-dpsi)

    # Per-edge prefactors.
    a_n = -L0_sq * mu_n_hat / h
    a_p = -L0_sq * mu_p_hat / h

    # SG residual contributions in the M13 sign convention (the leading
    # minus mirrors `b += -sg_n`):
    #   left vertex (i):  R_n_b[i] += a_n * (n_j*B(+dpsi) - n_i*B(-dpsi))
    #   right vertex (j): R_n_b[j] -= a_n * (n_j*B(+dpsi) - n_i*B(-dpsi))
    #
    # Galerkin convdiff contribution (also subtracted from b in the
    # callback) at the same vertices:
    #   left vertex (i):  -L0_sq*mu_n*(Δn/h - n_avg*dpsi/h) ... but its
    #     UFL Jacobian is auto-derived by NonlinearProblem; we subtract
    #     it from the corrected Jacobian by ALSO subtracting its analytic
    #     Jacobian (which equals the standard Galerkin stiffness:
    #       d(F_galerkin[i])/d(n_i)   = +a_n_galerkin / h    (etc.))
    #
    # We compute the SG-MINUS-Galerkin Jacobian correction directly to
    # avoid Mat-vs-Mat manipulations: at zero field SG = Galerkin so
    # the correction is zero, and at high Peclet only the SG entries
    # remain. The formulas below enumerate all 16 nonzero entries per
    # edge across the two row blocks (n, p) and four column blocks
    # (psi_i, psi_j, n_i, n_j or p_i, p_j).

    # M13 Galerkin convdiff Jacobian contributions per edge (we subtract
    # these along with the SG):
    inv_h = 1.0 / h
    # F_n_galerkin contribution per edge (i,j) at vertex i (left):
    #   contrib_i = -L0_sq * mu_n * (Δn/h - n_avg*dpsi/h)
    #   d(contrib_i)/d(n_i)   = -L0_sq*mu_n*(-1/h - 0.5*dpsi/h)
    #                         = +L0_sq*mu_n*(1/h + 0.5*dpsi/h)
    #   d(contrib_i)/d(n_j)   = -L0_sq*mu_n*(+1/h - 0.5*dpsi/h)
    #                         = -L0_sq*mu_n*(1/h - 0.5*dpsi/h)
    #   d(contrib_i)/d(psi_i) = -L0_sq*mu_n*(0 - n_avg*(-1)/h) = -L0_sq*mu_n*n_avg/h
    #   d(contrib_i)/d(psi_j) = -L0_sq*mu_n*(0 - n_avg*(+1)/h) = +L0_sq*mu_n*n_avg/h
    # And contrib_j = -contrib_i (mirror sign).
    n_avg = 0.5 * (n_i + n_j)
    p_avg = 0.5 * (p_i + p_j)

    rows_list = []
    cols_list = []
    vals_list = []

    n_edges = len(h)

    # We accumulate into per-vertex/per-edge dictionaries. To skip BC
    # rows, build the row index then skip if BC.

    # Helper to push a (row, col, val) entry, skipping BC rows.
    def push(row, col, val, bc_set):
        # row, col are arrays here (vectorised below).
        for k in range(len(row)):
            r = int(row[k])
            if (bc_set is not None) and (r in bc_set):
                continue
            rows_list.append(r)
            cols_list.append(int(col[k]))
            vals_list.append(float(val[k]))

    # Row indices in the global block layout.
    R_n_i = off_n + v2d_n[v_i]
    R_n_j = off_n + v2d_n[v_j]
    R_p_i = off_p + v2d_p[v_i]
    R_p_j = off_p + v2d_p[v_j]
    C_n_i = off_n + v2d_n[v_i]
    C_n_j = off_n + v2d_n[v_j]
    C_p_i = off_p + v2d_p[v_i]
    C_p_j = off_p + v2d_p[v_j]
    C_psi_i = off_psi + v2d_psi[v_i]
    C_psi_j = off_psi + v2d_psi[v_j]

    # --- F_n block: SG contribution (in M13 sign convention) ---
    # Per edge: residual at i gets `a_n * (n_j*B+ - n_i*B-)`
    # d/d(n_i): -a_n * B-
    # d/d(n_j): +a_n * B+
    # d/d(psi_i): a_n * (n_j*B+'*(-1) - n_i*B-'*(+1)) * (-1?)  Let me redo:
    #   contrib_i = a_n * (n_j*B_pos - n_i*B_neg) where B_pos=B(+dpsi), B_neg=B(-dpsi)
    #   d/d(psi_i): chain rule. d(B_pos)/d(psi_i) = B'(+dpsi) * d(dpsi)/d(psi_i) = Bp_pos * (-1)
    #              d(B_neg)/d(psi_i) = B'(-dpsi) * d(-dpsi)/d(psi_i) = Bp_neg * (+1)
    #   d(contrib_i)/d(psi_i) = a_n * (n_j*Bp_pos*(-1) - n_i*Bp_neg*(+1))
    #                         = -a_n * (n_j*Bp_pos + n_i*Bp_neg)
    #   d(contrib_i)/d(psi_j) = a_n * (n_j*Bp_pos*(+1) - n_i*Bp_neg*(-1))
    #                         = +a_n * (n_j*Bp_pos + n_i*Bp_neg)

    sg_d_ni = -a_n * B_neg                          # d(contrib_i)/d(n_i)
    sg_d_nj = +a_n * B_pos                          # d(contrib_i)/d(n_j)
    sg_d_psi_i = -a_n * (n_j * Bp_pos + n_i * Bp_neg)
    sg_d_psi_j = +a_n * (n_j * Bp_pos + n_i * Bp_neg)

    # Galerkin n convdiff Jacobian contributions per edge (at vertex i, left):
    gal_d_ni = -L0_sq * mu_n_hat * (-inv_h - 0.5 * dpsi * inv_h)
    gal_d_nj = -L0_sq * mu_n_hat * (+inv_h - 0.5 * dpsi * inv_h)
    gal_d_psi_i = -L0_sq * mu_n_hat * (n_avg * inv_h)        # from d(-n_avg*dpsi/h)/d(psi_i): -n_avg*(-1)/h = n_avg/h, then *(-L0²mu) = -L0²mu*n_avg/h
    gal_d_psi_j = +L0_sq * mu_n_hat * (n_avg * inv_h)

    # Correction: SG - Galerkin (each at vertex i)
    cn_d_ni_i = sg_d_ni - gal_d_ni
    cn_d_nj_i = sg_d_nj - gal_d_nj
    cn_d_psi_i_i = sg_d_psi_i - gal_d_psi_i
    cn_d_psi_j_i = sg_d_psi_j - gal_d_psi_j

    # At vertex j (right): contrib_j = -contrib_i, so all derivatives flip sign.
    cn_d_ni_j = -cn_d_ni_i
    cn_d_nj_j = -cn_d_nj_i
    cn_d_psi_i_j = -cn_d_psi_i_i
    cn_d_psi_j_j = -cn_d_psi_j_i

    # --- F_p block: same structure with the hole flux ---
    #   contrib_i (hole) = a_p * (p_i*B+ - p_j*B-)
    sg_d_pi = +a_p * B_pos
    sg_d_pj = -a_p * B_neg
    sg_d_psi_i_p = -a_p * (p_i * Bp_pos + p_j * Bp_neg)
    sg_d_psi_j_p = +a_p * (p_i * Bp_pos + p_j * Bp_neg)

    # Galerkin p convdiff: F_p_galerkin = +L0²*mu_p*(Δp/h + p_avg*dpsi/h) at vertex i, left
    # d/d(p_i):   -L0²*mu_p*(-1/h + 0.5*dpsi/h) (mirror of n with +sign on drift)
    #   contrib_i = -L0²*mu_p*(Δp/h + p_avg*dpsi/h)   [consistent with F_p form residual at i]
    #   d/d(p_i) = -L0²*mu_p*(-1/h + 0.5*dpsi/h) = L0²*mu_p*(1/h - 0.5*dpsi/h)
    #   d/d(p_j) = -L0²*mu_p*(+1/h + 0.5*dpsi/h)
    #   d/d(psi_i) = -L0²*mu_p*(p_avg*(-1)/h) * (+1) = +L0²*mu_p*p_avg/h
    #   d/d(psi_j) = -L0²*mu_p*p_avg/h
    gal_p_d_pi = +L0_sq * mu_p_hat * (inv_h - 0.5 * dpsi * inv_h)
    gal_p_d_pj = -L0_sq * mu_p_hat * (inv_h + 0.5 * dpsi * inv_h)
    gal_p_d_psi_i = +L0_sq * mu_p_hat * (p_avg * inv_h)
    gal_p_d_psi_j = -L0_sq * mu_p_hat * (p_avg * inv_h)

    cp_d_pi_i = sg_d_pi - gal_p_d_pi
    cp_d_pj_i = sg_d_pj - gal_p_d_pj
    cp_d_psi_i_i = sg_d_psi_i_p - gal_p_d_psi_i
    cp_d_psi_j_i = sg_d_psi_j_p - gal_p_d_psi_j

    cp_d_pi_j = -cp_d_pi_i
    cp_d_pj_j = -cp_d_pj_i
    cp_d_psi_i_j = -cp_d_psi_i_i
    cp_d_psi_j_j = -cp_d_psi_j_i

    # Push triples: (row, col, val)
    # n block — row at vertex i (left)
    push(R_n_i, C_n_i, cn_d_ni_i, bc_dofs_n_global := {off_n + d for d in bc_dofs_n})
    push(R_n_i, C_n_j, cn_d_nj_i, bc_dofs_n_global)
    push(R_n_i, C_psi_i, cn_d_psi_i_i, bc_dofs_n_global)
    push(R_n_i, C_psi_j, cn_d_psi_j_i, bc_dofs_n_global)
    # n block — row at vertex j (right)
    push(R_n_j, C_n_i, cn_d_ni_j, bc_dofs_n_global)
    push(R_n_j, C_n_j, cn_d_nj_j, bc_dofs_n_global)
    push(R_n_j, C_psi_i, cn_d_psi_i_j, bc_dofs_n_global)
    push(R_n_j, C_psi_j, cn_d_psi_j_j, bc_dofs_n_global)
    # p block — row at vertex i (left)
    bc_dofs_p_global = {off_p + d for d in bc_dofs_p}
    push(R_p_i, C_p_i, cp_d_pi_i, bc_dofs_p_global)
    push(R_p_i, C_p_j, cp_d_pj_i, bc_dofs_p_global)
    push(R_p_i, C_psi_i, cp_d_psi_i_i, bc_dofs_p_global)
    push(R_p_i, C_psi_j, cp_d_psi_j_i, bc_dofs_p_global)
    # p block — row at vertex j (right)
    push(R_p_j, C_p_i, cp_d_pi_j, bc_dofs_p_global)
    push(R_p_j, C_p_j, cp_d_pj_j, bc_dofs_p_global)
    push(R_p_j, C_psi_i, cp_d_psi_i_j, bc_dofs_p_global)
    push(R_p_j, C_psi_j, cp_d_psi_j_j, bc_dofs_p_global)

    return (
        np.asarray(rows_list, dtype=np.int32),
        np.asarray(cols_list, dtype=np.int32),
        np.asarray(vals_list, dtype=np.float64),
    )
