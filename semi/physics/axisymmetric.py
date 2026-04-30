"""
Axisymmetric (cylindrical) Poisson and drift-diffusion weak forms.

Geometry. We solve on the meridian half-plane (r, z) with r = x[0] >= 0
and z = x[1]. The physical 3D domain is recovered by revolving the
half-plane about the z-axis. The volume measure becomes

    dV_3D = 2*pi * r * dr * dz,

so every volume integrand of the cartesian weak forms picks up an
extra factor of r. The constant 2*pi cancels on both sides of the
residual.

This applies to both Poisson and the carrier continuity equations.
The Galerkin weak form for a scalar PDE  -div(A grad u) = f  becomes

    integral_Omega  A grad(u) . grad(v) * r dr dz
        =  integral_Omega  f * v * r dr dz                         (1)

with A possibly piecewise (e.g. eps_r jumping at the Si/SiO2 interface).

Boundary conditions.

- The symmetry axis r = 0 is a *natural* boundary: the r-weighted weak
  form has dV proportional to r, so test-function weights vanish there
  automatically and no Dirichlet row is needed (or allowed; see
  semi.schema._validate_coordinate_system).
- The outer radial wall r = R is treated as homogeneous Neumann (no
  flux). The user must choose R large enough that the solution near
  the axis is insensitive to the cutoff (>= ~5 * W_dmax for MOS).
- Top/bottom planar boundaries follow the usual Dirichlet/Neumann
  contact specification.

Sign and coefficient conventions match `semi.physics.poisson`:

    L_D^2 = lambda2 * L0^2     (scaled stiffness coefficient)

so the radial weighting is the only structural change from the
cartesian form.

References
----------
- docs/PHYSICS.md section 2 (weak forms in scaled units).
- Hu, Modern Semiconductor Devices for IC, chapter 5 (MOSCAP).
"""
from __future__ import annotations


def _radial_coord(msh):
    """Return UFL r = x[0] for use as the integrand weighting."""
    import ufl
    return ufl.SpatialCoordinate(msh)[0]


def build_equilibrium_poisson_form_axisym(V, psi, N_hat_fn, sc, eps_r):
    """
    Axisymmetric counterpart of
    :func:`semi.physics.poisson.build_equilibrium_poisson_form`.

    The single structural difference is the factor of r in both the
    stiffness and the source integrands; see equation (1) in the
    module docstring.

    Parameters
    ----------
    V : dolfinx.fem.FunctionSpace
        P1 Lagrange space for psi_hat on the meridian mesh.
    psi : dolfinx.fem.Function
        Scaled potential unknown.
    N_hat_fn : dolfinx.fem.Function
        Scaled net doping.
    sc : semi.scaling.Scaling
        Scaling object.
    eps_r : float | dolfinx.fem.Function
        Relative permittivity (scalar or DG0 cellwise).

    Returns
    -------
    F : ufl.Form
        Residual in scaled units, r-weighted.
    """
    import ufl
    from dolfinx import fem
    from petsc4py import PETSc

    msh = V.mesh
    v = ufl.TestFunction(V)
    r = _radial_coord(msh)

    L_D2 = fem.Constant(msh, PETSc.ScalarType(sc.lambda2 * sc.L0 ** 2))
    if isinstance(eps_r, (int, float)):
        eps_r_ufl = fem.Constant(msh, PETSc.ScalarType(float(eps_r)))
    else:
        eps_r_ufl = eps_r
    ni_hat = fem.Constant(msh, PETSc.ScalarType(sc.n_i / sc.C0))

    rho_hat = ni_hat * (ufl.exp(-psi) - ufl.exp(psi)) + N_hat_fn

    F = (
        L_D2 * eps_r_ufl * ufl.inner(ufl.grad(psi), ufl.grad(v)) * r * ufl.dx
        - rho_hat * v * r * ufl.dx
    )
    return F


def build_equilibrium_poisson_form_axisym_mr(
    V, psi, N_hat_fn, sc, eps_r_fn, cell_tags, semi_tag,
):
    """
    Multi-region axisymmetric Poisson, mirroring
    :func:`semi.physics.poisson.build_equilibrium_poisson_form_mr`.

    Stiffness integrates over the full meridian mesh with the cellwise
    DG0 ``eps_r_fn``; the Boltzmann space-charge term is restricted to
    semiconductor cells. Both integrals carry the radial weight r.
    """
    import ufl
    from dolfinx import fem
    from petsc4py import PETSc

    msh = V.mesh
    v = ufl.TestFunction(V)
    r = _radial_coord(msh)

    L_D2 = fem.Constant(msh, PETSc.ScalarType(sc.lambda2 * sc.L0 ** 2))
    ni_hat = fem.Constant(msh, PETSc.ScalarType(sc.n_i / sc.C0))

    dx_full = ufl.Measure("dx", domain=msh)
    dx_semi = ufl.Measure(
        "dx", domain=msh, subdomain_data=cell_tags, subdomain_id=int(semi_tag),
    )

    rho_hat = ni_hat * (ufl.exp(-psi) - ufl.exp(psi)) + N_hat_fn

    F = (
        L_D2 * eps_r_fn * ufl.inner(ufl.grad(psi), ufl.grad(v)) * r * dx_full
        - rho_hat * v * r * dx_semi
    )
    return F


def build_dd_block_residual_axisym(
    spaces,
    N_hat_fn,
    sc,
    eps_r,
    mu_n_over_mu0: float,
    mu_p_over_mu0: float,
    tau_n_hat: float,
    tau_p_hat: float,
    E_t_over_Vt: float = 0.0,
):
    """
    Axisymmetric counterpart of
    :func:`semi.physics.drift_diffusion.build_dd_block_residual`.

    Slotboom primary unknowns (psi, phi_n, phi_p). Every volume
    integrand carries the radial weight r:

        integral [ L_D^2 eps_r grad(psi).grad(v) - rho_hat v ] r dx
        integral [ L0^2 mu_n n_hat grad(phi_n).grad(v_n) - R v_n ] r dx
        integral [ L0^2 mu_p p_hat grad(phi_p).grad(v_p) + R v_p ] r dx

    Sign conventions and the SRH expression are identical to the
    cartesian path; only the measure is r-weighted.
    """
    import ufl
    from dolfinx import fem
    from petsc4py import PETSc

    from .slotboom import n_from_slotboom, p_from_slotboom

    psi = spaces.psi
    phi_n = spaces.phi_n
    phi_p = spaces.phi_p
    msh = spaces.V_psi.mesh
    r = _radial_coord(msh)

    v_psi = ufl.TestFunction(spaces.V_psi)
    v_n = ufl.TestFunction(spaces.V_phi_n)
    v_p = ufl.TestFunction(spaces.V_phi_p)

    L_D2 = fem.Constant(msh, PETSc.ScalarType(sc.lambda2 * sc.L0 ** 2))
    L0_sq = fem.Constant(msh, PETSc.ScalarType(sc.L0 ** 2))
    if isinstance(eps_r, (int, float)):
        eps_r_ufl = fem.Constant(msh, PETSc.ScalarType(float(eps_r)))
    else:
        eps_r_ufl = eps_r
    ni_hat = fem.Constant(msh, PETSc.ScalarType(sc.n_i / sc.C0))
    mu_n_hat = fem.Constant(msh, PETSc.ScalarType(mu_n_over_mu0))
    mu_p_hat = fem.Constant(msh, PETSc.ScalarType(mu_p_over_mu0))
    tau_n = fem.Constant(msh, PETSc.ScalarType(tau_n_hat))
    tau_p = fem.Constant(msh, PETSc.ScalarType(tau_p_hat))

    n_hat = n_from_slotboom(psi, phi_n, ni_hat)
    p_hat = p_from_slotboom(psi, phi_p, ni_hat)

    import math as _math
    n1 = ni_hat * _math.exp(E_t_over_Vt)
    p1 = ni_hat * _math.exp(-E_t_over_Vt)
    R = (n_hat * p_hat - ni_hat * ni_hat) / (
        tau_p * (n_hat + n1) + tau_n * (p_hat + p1)
    )

    rho_hat = p_hat - n_hat + N_hat_fn

    F_psi = (
        L_D2 * eps_r_ufl * ufl.inner(ufl.grad(psi), ufl.grad(v_psi)) * r * ufl.dx
        - rho_hat * v_psi * r * ufl.dx
    )
    F_phi_n = (
        L0_sq * mu_n_hat * n_hat * ufl.inner(ufl.grad(phi_n), ufl.grad(v_n)) * r * ufl.dx
        - R * v_n * r * ufl.dx
    )
    F_phi_p = (
        L0_sq * mu_p_hat * p_hat * ufl.inner(ufl.grad(phi_p), ufl.grad(v_p)) * r * ufl.dx
        + R * v_p * r * ufl.dx
    )
    return [F_psi, F_phi_n, F_phi_p]
