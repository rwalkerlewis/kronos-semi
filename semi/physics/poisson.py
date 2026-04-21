"""
Equilibrium Poisson equation (Boltzmann statistics).

At thermal equilibrium with Boltzmann statistics and zero applied bias,
electrons and holes follow the potential directly:

    n(psi) = n_i exp(+psi / V_t)
    p(psi) = n_i exp(-psi / V_t)

Substituting into Poisson gives a single nonlinear PDE for psi. In
scaled units (psi_hat = psi / V_t, densities / C_0, length / L_0):

    -lambda^2 * eps_r * nabla^2 psi_hat = n_i_hat (exp(-psi_hat) - exp(psi_hat))
                                          + N_net_hat

where lambda^2 = eps_0 V_t / (q C_0 L_0^2).

This gives the built-in potential of pn junctions, MOS C-V curves at
zero bias, etc. Under bias, we'd need the full drift-diffusion system;
this module is Day 1 scope.
"""
from __future__ import annotations


def build_equilibrium_poisson_form(V, psi, N_hat_fn, sc, eps_r_value: float):
    """
    Build the UFL residual form for equilibrium Poisson.

    Parameters
    ----------
    V : dolfinx.fem.FunctionSpace
        P1 Lagrange space for psi_hat.
    psi : dolfinx.fem.Function
        The unknown scaled potential (will be solved for).
    N_hat_fn : dolfinx.fem.Function
        Scaled net doping profile interpolated into V (or a compatible space).
    sc : semi.scaling.Scaling
        Scaling object providing lambda2, n_i/C_0 ratio.
    eps_r_value : float
        Relative permittivity (assumed uniform here; pass a Function for
        multi-region once we add 2D/3D).

    Returns
    -------
    F : ufl.Form
        Residual form. Pass to dolfinx.fem.petsc.NonlinearProblem.
    """
    import ufl
    from dolfinx import fem
    from petsc4py import PETSc

    msh = V.mesh
    v = ufl.TestFunction(V)

    # The mesh coordinates are in physical meters (see semi/mesh.py), so the
    # stiffness coefficient is the dimensional squared Debye length
    #     L_D^2 = eps_0 V_0 / (q C_0) = lambda2 * L_0^2
    # rather than the bare scaled `lambda2`. The nondimensionalization only
    # affects psi (scaled by V_t) and densities (scaled by C_0); the spatial
    # coordinate is left in meters so that mesh extents and facet locations
    # from the JSON can be specified in SI units.
    L_D2   = fem.Constant(msh, PETSc.ScalarType(sc.lambda2 * sc.L0 ** 2))
    eps_r  = fem.Constant(msh, PETSc.ScalarType(eps_r_value))
    ni_hat = fem.Constant(msh, PETSc.ScalarType(sc.n_i / sc.C0))

    rho_hat = ni_hat * (ufl.exp(-psi) - ufl.exp(psi)) + N_hat_fn

    F = (
        L_D2 * eps_r * ufl.inner(ufl.grad(psi), ufl.grad(v)) * ufl.dx
        - rho_hat * v * ufl.dx
    )
    return F
