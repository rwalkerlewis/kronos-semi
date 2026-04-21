"""
Solver driver.

Wraps dolfinx.fem.petsc.NonlinearProblem (PETSc SNES) with default
options tuned for semiconductor problems:
    - Direct LU (MUMPS) linear solve per Newton iteration — robust for
      the block-ill-conditioned Jacobians that show up here
    - Backtracking line search — helps when initial guess is far from the
      root, especially with exponentials in the residual
    - SNES monitor on by default so we can diagnose convergence issues

For large 3D problems, LU is overkill and we'd switch to GAMG or
fieldsplit preconditioning, but that's a Day 6+ concern.
"""
from __future__ import annotations

from typing import Any

DEFAULT_PETSC_OPTIONS = {
    "snes_type": "newtonls",
    "snes_linesearch_type": "bt",
    "snes_rtol": 1.0e-10,
    "snes_atol": 1.0e-12,
    "snes_max_it": 50,
    "snes_monitor": None,
    "ksp_type": "preonly",
    "pc_type": "lu",
    "pc_factor_mat_solver_type": "mumps",
}


def solve_nonlinear(F, u, bcs: list, prefix: str,
                    petsc_options: dict[str, Any] | None = None):
    """
    Solve a nonlinear variational problem F(u; v) = 0.

    Parameters
    ----------
    F : ufl.Form
        Residual form.
    u : dolfinx.fem.Function
        Unknown function (will be modified in place with the solution).
    bcs : list of dolfinx.fem.DirichletBC
        Dirichlet boundary conditions.
    prefix : str
        PETSc options prefix. Must be unique per problem (required in 0.10+).
    petsc_options : dict, optional
        Override DEFAULT_PETSC_OPTIONS.

    Returns
    -------
    dict
        {'iterations': int, 'reason': int, 'converged': bool}
    """
    from dolfinx.fem.petsc import NonlinearProblem

    opts = dict(DEFAULT_PETSC_OPTIONS)
    if petsc_options:
        opts.update(petsc_options)

    problem = NonlinearProblem(
        F, u,
        bcs=bcs,
        petsc_options_prefix=prefix,
        petsc_options=opts,
    )
    problem.solve()
    reason = problem.solver.getConvergedReason()
    n_iter = problem.solver.getIterationNumber()
    converged = reason > 0
    return {
        "iterations": int(n_iter),
        "reason": int(reason),
        "converged": bool(converged),
        "problem": problem,  # kept so caller can inspect further
    }


def solve_nonlinear_block(
    F_list,
    u_list,
    bcs: list,
    prefix: str,
    petsc_options: dict[str, Any] | None = None,
    kind: str | None = None,
):
    """
    Solve a coupled block nonlinear problem via SNES.

    Wraps dolfinx 0.10's `NonlinearProblem` with its blocked-mode
    signature (F and u passed as sequences, plus optional `kind`).
    Dirichlet BCs are passed as a single flat list; each BC is tied to
    the subspace it was built against.

    Parameters
    ----------
    F_list : sequence of ufl.Form
        One residual form per block (for example [F_psi, F_phi_n, F_phi_p]).
    u_list : sequence of dolfinx.fem.Function
        Unknown Functions for each block (updated in place on return).
    bcs : list of dolfinx.fem.DirichletBC
        Flat list of Dirichlet BCs.
    prefix : str
        PETSc options prefix (must end with `_`).
    petsc_options : dict, optional
        Override DEFAULT_PETSC_OPTIONS.
    kind : str, optional
        PETSc matrix/vector kind. `None` for monolithic blocked (default,
        works with the MUMPS LU solver). Use `"nest"` if you supply a
        fieldsplit preconditioner.

    Returns
    -------
    dict
        Same keys as :func:`solve_nonlinear`.
    """
    from dolfinx.fem.petsc import NonlinearProblem

    opts = dict(DEFAULT_PETSC_OPTIONS)
    if petsc_options:
        opts.update(petsc_options)

    problem = NonlinearProblem(
        list(F_list), list(u_list),
        bcs=list(bcs),
        petsc_options_prefix=prefix,
        petsc_options=opts,
        kind=kind,
    )
    problem.solve()
    reason = problem.solver.getConvergedReason()
    n_iter = problem.solver.getIterationNumber()
    converged = reason > 0
    return {
        "iterations": int(n_iter),
        "reason": int(reason),
        "converged": bool(converged),
        "problem": problem,
    }
