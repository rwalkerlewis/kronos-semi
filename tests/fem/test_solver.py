"""
Direct unit tests for the SNES wrappers in `semi.solver`.

The wrappers are exercised transitively by every benchmark and V&V
study, but a hands-on linear-Poisson sanity check makes any
infrastructure regression in the wrapper itself surface immediately
without having to re-run a full simulation.
"""
from __future__ import annotations

import numpy as np
import pytest

from semi.solver import solve_nonlinear, solve_nonlinear_block


def test_solve_nonlinear_recovers_linear_poisson_on_unit_interval():
    """
    Solve -u'' = 2 with u(0) = u(1) = 0; analytic solution u(x) = x (1 - x).
    Newton converges in one step (linear residual), proving the wrapper
    parses the form, applies BCs, and returns iteration metadata.
    """
    import ufl
    from dolfinx import fem, mesh
    from mpi4py import MPI
    from petsc4py import PETSc

    msh = mesh.create_interval(MPI.COMM_WORLD, 32, [0.0, 1.0])
    V = fem.functionspace(msh, ("Lagrange", 1))

    u = fem.Function(V, name="u")
    v = ufl.TestFunction(V)
    f = fem.Constant(msh, PETSc.ScalarType(2.0))
    F = ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx - f * v * ufl.dx

    # Homogeneous Dirichlet on both endpoints.
    fdim = msh.topology.dim - 1
    boundary = mesh.locate_entities_boundary(
        msh, fdim, lambda x: np.full(x.shape[1], True)
    )
    dofs = fem.locate_dofs_topological(V, fdim, boundary)
    bc = fem.dirichletbc(PETSc.ScalarType(0.0), dofs, V)

    info = solve_nonlinear(F, u, [bc], prefix="test_lin_poisson_")
    assert info["converged"]
    # A linear problem needs at most a couple of Newton iterations.
    assert info["iterations"] <= 2

    # Check the solution at midpoint matches the analytic value 0.25.
    x_dof = V.tabulate_dof_coordinates()[:, 0]
    mid_idx = int(np.argmin(np.abs(x_dof - 0.5)))
    assert u.x.array[mid_idx] == pytest.approx(0.25, abs=5.0e-4)


def test_solve_nonlinear_petsc_options_override_propagates():
    """
    Pass an explicit `snes_max_it` override and confirm the wrapper
    accepts it without raising. Convergence is checked indirectly by
    the iteration count staying below the (very loose) override.
    """
    import ufl
    from dolfinx import fem, mesh
    from mpi4py import MPI
    from petsc4py import PETSc

    msh = mesh.create_interval(MPI.COMM_WORLD, 16, [0.0, 1.0])
    V = fem.functionspace(msh, ("Lagrange", 1))
    u = fem.Function(V)
    v = ufl.TestFunction(V)
    F = ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx \
        - fem.Constant(msh, PETSc.ScalarType(1.0)) * v * ufl.dx

    fdim = msh.topology.dim - 1
    boundary = mesh.locate_entities_boundary(
        msh, fdim, lambda x: np.full(x.shape[1], True)
    )
    dofs = fem.locate_dofs_topological(V, fdim, boundary)
    bc = fem.dirichletbc(PETSc.ScalarType(0.0), dofs, V)

    info = solve_nonlinear(
        F, u, [bc],
        prefix="test_options_",
        petsc_options={"snes_max_it": 10},
    )
    assert info["converged"]
    assert info["iterations"] <= 10


def test_solve_nonlinear_block_decoupled_two_blocks():
    """
    Build a synthetic 2-block decoupled problem (`-u_1'' = 1` and
    `-u_2'' = 2`, each on its own scalar space) and solve as a block.
    Exercises the block wrapper machinery without DD physics.
    """
    import ufl
    from dolfinx import fem, mesh
    from mpi4py import MPI
    from petsc4py import PETSc

    msh = mesh.create_interval(MPI.COMM_WORLD, 16, [0.0, 1.0])
    V1 = fem.functionspace(msh, ("Lagrange", 1))
    V2 = fem.functionspace(msh, ("Lagrange", 1))

    u1 = fem.Function(V1, name="u1")
    u2 = fem.Function(V2, name="u2")
    v1 = ufl.TestFunction(V1)
    v2 = ufl.TestFunction(V2)

    f1 = fem.Constant(msh, PETSc.ScalarType(1.0))
    f2 = fem.Constant(msh, PETSc.ScalarType(2.0))
    F1 = ufl.dot(ufl.grad(u1), ufl.grad(v1)) * ufl.dx - f1 * v1 * ufl.dx
    F2 = ufl.dot(ufl.grad(u2), ufl.grad(v2)) * ufl.dx - f2 * v2 * ufl.dx

    fdim = msh.topology.dim - 1
    boundary = mesh.locate_entities_boundary(
        msh, fdim, lambda x: np.full(x.shape[1], True)
    )
    dofs1 = fem.locate_dofs_topological(V1, fdim, boundary)
    dofs2 = fem.locate_dofs_topological(V2, fdim, boundary)
    bc1 = fem.dirichletbc(PETSc.ScalarType(0.0), dofs1, V1)
    bc2 = fem.dirichletbc(PETSc.ScalarType(0.0), dofs2, V2)

    info = solve_nonlinear_block(
        [F1, F2], [u1, u2], [bc1, bc2],
        prefix="test_block_",
    )
    assert info["converged"]
    # Linear problems converge quickly even when posed as a block.
    assert info["iterations"] <= 3

    # u1 midpoint solution to -u'' = 1, u(0)=u(1)=0 is u(0.5) = 1/8.
    x_dof = V1.tabulate_dof_coordinates()[:, 0]
    mid_idx = int(np.argmin(np.abs(x_dof - 0.5)))
    assert u1.x.array[mid_idx] == pytest.approx(0.125, abs=5.0e-4)
    # u2 with f = 2 is u(0.5) = 0.25.
    assert u2.x.array[mid_idx] == pytest.approx(0.25, abs=5.0e-4)
