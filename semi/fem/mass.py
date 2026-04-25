"""
Lumped mass matrix assembly for the carrier continuity equations.

The Poisson equation does NOT get a mass matrix entry because its
unknown (psi) is not time-differentiated; only the carrier densities
n and p are time-dependent in the transient continuity equations.

Row-sum lumping technique: the consistent mass matrix M[i,j] =
integral(phi_i * phi_j) dx integrates to a diagonal D[i,i] = integral(phi_i)
dx when the row is summed, because sum_j phi_j = 1 for P1 elements
(partition of unity). We assemble this diagonal directly by assembling a
scalar functional with test function v and constant coefficient 1.

Only used at module scope for the dolfinx-dependent assembly inside
runner subprocesses; the dolfinx import is deferred inside the function.
"""
from __future__ import annotations

from typing import Any


def assemble_lumped_mass(V_n: Any, V_p: Any, dx: Any) -> tuple:
    """
    Assemble the lumped mass matrix diagonal for n and p spaces.

    The lumped mass diagonal M_diag[i] = integral(phi_i) dx is computed
    by assembling the linear functional

        M[v] = integral(1 * v) dx

    for each subspace. For P1 Lagrange elements this equals the row-sum
    of the consistent mass matrix (partition-of-unity property).

    The Poisson equation does NOT get a mass matrix entry because psi is
    not time-differentiated.

    Parameters
    ----------
    V_n : dolfinx.fem.FunctionSpace
        P1 function space for the electron density n_hat.
    V_p : dolfinx.fem.FunctionSpace
        P1 function space for the hole density p_hat.
    dx : ufl.Measure
        Volume integration measure on the mesh.

    Returns
    -------
    tuple of (petsc4py.PETSc.Vec, petsc4py.PETSc.Vec)
        (M_n_diag, M_p_diag) -- row-sum-lumped mass vectors. Each vector
        has the same layout as the corresponding function space's DOF
        vector; entry i holds the lumped mass for DOF i.
    """
    import ufl
    from dolfinx import fem
    from dolfinx.fem.petsc import assemble_vector
    from petsc4py import PETSc

    msh = V_n.mesh

    # Assemble M_n: integral(v_n) dx for each DOF in V_n
    v_n = ufl.TestFunction(V_n)
    L_n = fem.form(v_n * dx)
    M_n_diag = assemble_vector(L_n)
    M_n_diag.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)

    # Assemble M_p: integral(v_p) dx for each DOF in V_p
    v_p = ufl.TestFunction(V_p)
    L_p = fem.form(v_p * dx)
    M_p_diag = assemble_vector(L_p)
    M_p_diag.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)

    return M_n_diag, M_p_diag
