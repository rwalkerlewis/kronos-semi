"""
FEM-side smoke tests for the UFL `srh_rate` helper.

The drift-diffusion residual builder inlines its own SRH expression, so
this UFL helper has no in-package callers; the tests here keep it from
silently rotting and confirm it agrees with the NumPy reference at
representative carrier densities.
"""
from __future__ import annotations

import numpy as np
import pytest

from semi.physics.recombination import srh_rate, srh_rate_np


def _project_to_constant_value(expr, msh):
    """Assemble `expr * dx` and divide by the domain measure to get the (constant) value."""
    import ufl
    from dolfinx import fem

    val = fem.assemble_scalar(fem.form(expr * ufl.dx))
    area = fem.assemble_scalar(fem.form(1.0 * ufl.dx(domain=msh)))
    return float(val / area)


def test_srh_rate_ufl_matches_np_at_uniform_state():
    """
    Build constant n_hat, p_hat fields on a 1D interval and assemble
    `srh_rate * dx`; compare to the NumPy formula at the same densities.
    """
    from dolfinx import fem, mesh
    from mpi4py import MPI
    from petsc4py import PETSc

    msh = mesh.create_interval(MPI.COMM_WORLD, 8, [0.0, 1.0])

    n_hat_val = 5.0
    p_hat_val = 0.2
    n_i_hat_val = 1.0
    tau_n_hat_val = 1.0
    tau_p_hat_val = 1.0
    E_t = 0.0

    n_hat = fem.Constant(msh, PETSc.ScalarType(n_hat_val))
    p_hat = fem.Constant(msh, PETSc.ScalarType(p_hat_val))
    n_i_hat = fem.Constant(msh, PETSc.ScalarType(n_i_hat_val))
    tau_n_hat = fem.Constant(msh, PETSc.ScalarType(tau_n_hat_val))
    tau_p_hat = fem.Constant(msh, PETSc.ScalarType(tau_p_hat_val))

    expr = srh_rate(n_hat, p_hat, n_i_hat, tau_n_hat, tau_p_hat, E_t_over_Vt=E_t)
    R_ufl = _project_to_constant_value(expr, msh)

    R_np = srh_rate_np(
        np.array(n_hat_val), np.array(p_hat_val), n_i_hat_val,
        tau_n_hat_val, tau_p_hat_val, E_t,
    )
    assert R_ufl == pytest.approx(float(R_np), rel=1e-10)


def test_srh_rate_ufl_off_midgap_trap_matches_np():
    """Same comparison with E_t / V_t != 0 to exercise the n1 / p1 branch."""
    from dolfinx import fem, mesh
    from mpi4py import MPI
    from petsc4py import PETSc

    msh = mesh.create_interval(MPI.COMM_WORLD, 4, [0.0, 1.0])

    args = dict(n_hat_val=3.0, p_hat_val=0.5, n_i_hat_val=1.0,
                tau_n=2.0, tau_p=0.5, E_t=1.5)
    n_hat = fem.Constant(msh, PETSc.ScalarType(args["n_hat_val"]))
    p_hat = fem.Constant(msh, PETSc.ScalarType(args["p_hat_val"]))
    n_i_hat = fem.Constant(msh, PETSc.ScalarType(args["n_i_hat_val"]))
    tau_n = fem.Constant(msh, PETSc.ScalarType(args["tau_n"]))
    tau_p = fem.Constant(msh, PETSc.ScalarType(args["tau_p"]))

    expr = srh_rate(n_hat, p_hat, n_i_hat, tau_n, tau_p, E_t_over_Vt=args["E_t"])
    R_ufl = _project_to_constant_value(expr, msh)
    R_np = srh_rate_np(
        np.array(args["n_hat_val"]), np.array(args["p_hat_val"]),
        args["n_i_hat_val"], args["tau_n"], args["tau_p"], args["E_t"],
    )
    assert R_ufl == pytest.approx(float(R_np), rel=1e-10)
