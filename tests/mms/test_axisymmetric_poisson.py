"""
Cylindrical-coordinate Poisson MMS convergence test (ADR 0015).

Method
------
Manufactured solution on the (r, z) "unit" half-cylinder
``r in [0, 1]``, ``z in [0, 1]``::

    phi(r, z) = (1 - r^2) * sin(pi * z)

* Vanishes on the outer boundary ``r = 1`` and on ``z = 0`` and
  ``z = 1`` (homogeneous Dirichlet on three sides; no BC
  interpolation error).
* Has zero radial derivative at ``r = 0`` (``d phi/d r = -2 r = 0``),
  so the natural Neumann zero-flux condition on the symmetry axis is
  satisfied exactly. The r-weighted measure enforces this BC for free
  (ADR 0015), so we leave that boundary unconstrained in the variational
  problem.

The cylindrical-coordinate Laplacian on an axially-symmetric ``phi``:

    Lap_cyl(phi) = (1/r) d/dr (r d phi/dr) + d2 phi/dz2
                 = -4 sin(pi z) - pi^2 (1 - r^2) sin(pi z)

so the manufactured source for ``-Lap_cyl(phi) = f`` is

    f(r, z) = [4 + pi^2 (1 - r^2)] sin(pi z).

Multiplying through by ``r`` (the axisymmetric weak form per ADR 0015)
gives the 2D mesh integral::

    int_Omega [grad(phi) . grad(v) - f * v] r dr dz = 0

which we solve on a sequence of CG1 triangle meshes
``h in {1/4, 1/8, 1/16, 1/32}``. The L2 rate (finest pair) must
exceed 1.9.

The test only exercises the dispatch helpers in
``semi.fem.coordinates``; it is independent of any drift-diffusion
runner so a regression here points at the cylindrical assembly
formulation directly.
"""
from __future__ import annotations

import math

import numpy as np


def _solve_axisymmetric_poisson(N: int) -> tuple[float, float]:
    """Solve the manufactured cylindrical Poisson at resolution N.

    Returns ``(h, e_L2)`` where ``h = 1/N``.
    """
    import ufl
    from dolfinx import fem
    from dolfinx.fem.petsc import LinearProblem
    from dolfinx.mesh import CellType, create_rectangle, locate_entities_boundary
    from mpi4py import MPI
    from petsc4py import PETSc

    from semi.fem.coordinates import get_volume_measure

    msh = create_rectangle(
        MPI.COMM_WORLD,
        [np.array([0.0, 0.0]), np.array([1.0, 1.0])],
        [N, N],
        cell_type=CellType.triangle,
    )

    V = fem.functionspace(msh, ("Lagrange", 1))
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)

    x = ufl.SpatialCoordinate(msh)
    r = x[0]
    z = x[1]
    phi_exact = (1.0 - r * r) * ufl.sin(ufl.pi * z)
    f = (4.0 + ufl.pi ** 2 * (1.0 - r * r)) * ufl.sin(ufl.pi * z)

    dx = get_volume_measure(msh, "axisymmetric")

    a = ufl.inner(ufl.grad(u), ufl.grad(v)) * dx
    L = f * v * dx

    # Homogeneous Dirichlet on r=1, z=0, z=1. The r=0 axis is
    # left as a natural Neumann zero-flux boundary, enforced
    # automatically by the r-weighted measure (ADR 0015).
    fdim = msh.topology.dim - 1
    msh.topology.create_connectivity(fdim, msh.topology.dim)

    def on_dirichlet(x_pts):
        return (
            np.isclose(x_pts[0], 1.0)
            | np.isclose(x_pts[1], 0.0)
            | np.isclose(x_pts[1], 1.0)
        )

    facets = locate_entities_boundary(msh, fdim, on_dirichlet)
    dofs = fem.locate_dofs_topological(V, fdim, facets)
    zero = fem.Constant(msh, PETSc.ScalarType(0.0))
    bcs = [fem.dirichletbc(zero, dofs, V)]

    uh = fem.Function(V, name="phi_h")
    problem = LinearProblem(
        a, L, u=uh, bcs=bcs,
        petsc_options_prefix="axisym_mms_",
        petsc_options={
            "ksp_type": "preonly",
            "pc_type": "lu",
            "pc_factor_mat_solver_type": "mumps",
        },
    )
    problem.solve()

    # L2 error using the cartesian (unweighted) measure: the
    # discretisation error rate is a property of the FE space on
    # the same 2D mesh and is invariant under the r-weighting.
    err_sq = ufl.inner(uh - phi_exact, uh - phi_exact) * ufl.dx
    err_form = fem.form(err_sq)
    e2_local = fem.assemble_scalar(err_form)
    e2 = msh.comm.allreduce(e2_local, op=MPI.SUM)
    e_L2 = math.sqrt(max(e2, 0.0))

    h = 1.0 / N
    return h, e_L2


def test_axisymmetric_poisson_l2_rate():
    """L2 convergence rate (finest pair) for cylindrical Poisson MMS.

    Theoretical rate for CG1 is 2; we gate at 1.9 to allow some
    degradation near the symmetry axis from the r-weighted measure.
    """
    Ns = [4, 8, 16, 32]
    hs: list[float] = []
    errs: list[float] = []
    for N in Ns:
        h, e = _solve_axisymmetric_poisson(N)
        hs.append(h)
        errs.append(e)

    for i in range(len(errs) - 1):
        assert errs[i + 1] < errs[i], (
            f"L2 error not decreasing: {errs}"
        )

    rate = math.log(errs[-2] / errs[-1]) / math.log(hs[-2] / hs[-1])
    assert rate >= 1.9, (
        f"axisymmetric Poisson MMS L2 rate {rate:.3f} < 1.9 "
        f"(hs={hs}, errs={errs})"
    )
