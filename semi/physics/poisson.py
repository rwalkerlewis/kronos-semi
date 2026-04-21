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


def build_equilibrium_poisson_form(V, psi, N_hat_fn, sc, eps_r):
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
    eps_r : float | dolfinx.fem.Function
        Relative permittivity. Scalar (single-region fast path,
        byte-identical with Day 2-5) or a cellwise DG0 Function
        (multi-region, Day 6). The scalar branch wraps the value in
        `fem.Constant`; the Function branch uses it directly in the
        bilinear form so coefficient-jump assembly picks up the
        piecewise eps_r at quadrature.

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
    L_D2 = fem.Constant(msh, PETSc.ScalarType(sc.lambda2 * sc.L0 ** 2))
    if isinstance(eps_r, (int, float)):
        eps_r_ufl = fem.Constant(msh, PETSc.ScalarType(float(eps_r)))
    else:
        eps_r_ufl = eps_r
    ni_hat = fem.Constant(msh, PETSc.ScalarType(sc.n_i / sc.C0))

    rho_hat = ni_hat * (ufl.exp(-psi) - ufl.exp(psi)) + N_hat_fn

    F = (
        L_D2 * eps_r_ufl * ufl.inner(ufl.grad(psi), ufl.grad(v)) * ufl.dx
        - rho_hat * v * ufl.dx
    )
    return F


def build_equilibrium_poisson_form_mr(
    V, psi, N_hat_fn, sc, eps_r_fn, cell_tags, semi_tag,
):
    """
    Build the UFL residual for multi-region equilibrium Poisson (MOS).

    Stiffness integrates over the full mesh with cellwise `eps_r_fn`;
    the Boltzmann space-charge term is restricted to semiconductor cells
    via a `dx(subdomain_id=semi_tag)` measure, so the oxide region
    carries only the Laplacian (no mobile carriers, no doping). The
    Si/SiO2 interface natural condition
    eps_r_Si grad psi . n = eps_r_ox grad psi . n
    is enforced automatically by the piecewise eps_r in the bilinear
    form (docs/mos_derivation.md section 3.1).
    """
    import ufl
    from dolfinx import fem
    from petsc4py import PETSc

    msh = V.mesh
    v = ufl.TestFunction(V)

    L_D2 = fem.Constant(msh, PETSc.ScalarType(sc.lambda2 * sc.L0 ** 2))
    ni_hat = fem.Constant(msh, PETSc.ScalarType(sc.n_i / sc.C0))

    dx_full = ufl.Measure("dx", domain=msh)
    dx_semi = ufl.Measure(
        "dx", domain=msh, subdomain_data=cell_tags, subdomain_id=int(semi_tag),
    )

    rho_hat = ni_hat * (ufl.exp(-psi) - ufl.exp(psi)) + N_hat_fn

    F = (
        L_D2 * eps_r_fn * ufl.inner(ufl.grad(psi), ufl.grad(v)) * dx_full
        - rho_hat * v * dx_semi
    )
    return F
