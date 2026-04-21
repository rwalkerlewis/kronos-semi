"""
Coupled drift-diffusion block residual (Slotboom form).

Primary unknowns: (psi_hat, phi_n_hat, phi_p_hat) in scaled units.
Densities are recovered pointwise via the Slotboom relations in
:mod:`semi.physics.slotboom`.

Scaled equations (see docs/PHYSICS.md section 2):

    Poisson:
        -div( L_D^2 eps_r grad psi_hat ) = p_hat - n_hat + N_hat

    Electron continuity:
        -div( L_0^2 mu_n_hat n_hat grad phi_n_hat ) = R_hat

    Hole continuity:
        -div( L_0^2 mu_p_hat p_hat grad phi_p_hat ) = -R_hat

where L_D^2 = lambda2 * L_0^2, mu_n_hat = mu_n / mu_0, and R_hat is the
scaled SRH rate (see :mod:`semi.physics.recombination`). The mesh stays
in physical meters (Invariant 3, see PLAN.md), so L_0^2 and L_D^2 appear
explicitly as UFL constants in the forms.

The residuals returned here are weak (Galerkin) forms with homogeneous
Neumann natural boundaries. Dirichlet contact data on psi, phi_n, phi_p
are applied via `dolfinx.fem.DirichletBC` by the caller.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass
class DDBlockSpaces:
    """Container for the three P1 subspaces and their unknowns."""
    V_psi: Any
    V_phi_n: Any
    V_phi_p: Any
    psi: Any
    phi_n: Any
    phi_p: Any


def make_dd_block_spaces(msh) -> DDBlockSpaces:
    """
    Create three independent P1 Lagrange spaces on `msh` and their
    unknown Functions.

    Using three separate scalar spaces (rather than a vector MixedElement
    on the same mesh) keeps the block residual assembly explicit and
    lets us re-use the Day 1 equilibrium Poisson solver for the initial
    guess.
    """
    from dolfinx import fem

    V_psi = fem.functionspace(msh, ("Lagrange", 1))
    V_phi_n = fem.functionspace(msh, ("Lagrange", 1))
    V_phi_p = fem.functionspace(msh, ("Lagrange", 1))
    psi = fem.Function(V_psi, name="psi_hat")
    phi_n = fem.Function(V_phi_n, name="phi_n_hat")
    phi_p = fem.Function(V_phi_p, name="phi_p_hat")
    return DDBlockSpaces(V_psi, V_phi_n, V_phi_p, psi, phi_n, phi_p)


def build_dd_block_residual(
    spaces: DDBlockSpaces,
    N_hat_fn,
    sc,
    eps_r_value: float,
    mu_n_over_mu0: float,
    mu_p_over_mu0: float,
    tau_n_hat: float,
    tau_p_hat: float,
    E_t_over_Vt: float = 0.0,
):
    """
    Build the three-block residual for the coupled drift-diffusion system.

    Parameters
    ----------
    spaces : DDBlockSpaces
        Function spaces and unknowns (psi, phi_n, phi_p in scaled units).
    N_hat_fn : dolfinx.fem.Function
        Scaled net doping interpolated in V_psi or a compatible space.
    sc : semi.scaling.Scaling
        Scaling object; provides lambda2, L0, n_i/C0 ratio.
    eps_r_value : float
        Uniform relative permittivity (1D assumption for Day 2).
    mu_n_over_mu0, mu_p_over_mu0 : float
        Ratios of carrier mobility to the scaling reference mobility.
    tau_n_hat, tau_p_hat : float
        Scaled lifetimes, tau / t0.
    E_t_over_Vt : float
        Trap level relative to the intrinsic level, divided by V_t.

    Returns
    -------
    list of ufl.Form
        [F_psi, F_phi_n, F_phi_p] -- pass to the blocked NonlinearProblem.
    """
    import ufl
    from dolfinx import fem
    from petsc4py import PETSc

    from .slotboom import n_from_slotboom, p_from_slotboom

    psi = spaces.psi
    phi_n = spaces.phi_n
    phi_p = spaces.phi_p
    msh = spaces.V_psi.mesh

    v_psi = ufl.TestFunction(spaces.V_psi)
    v_n = ufl.TestFunction(spaces.V_phi_n)
    v_p = ufl.TestFunction(spaces.V_phi_p)

    L_D2 = fem.Constant(msh, PETSc.ScalarType(sc.lambda2 * sc.L0 ** 2))
    L0_sq = fem.Constant(msh, PETSc.ScalarType(sc.L0 ** 2))
    eps_r = fem.Constant(msh, PETSc.ScalarType(eps_r_value))
    ni_hat = fem.Constant(msh, PETSc.ScalarType(sc.n_i / sc.C0))
    mu_n_hat = fem.Constant(msh, PETSc.ScalarType(mu_n_over_mu0))
    mu_p_hat = fem.Constant(msh, PETSc.ScalarType(mu_p_over_mu0))
    tau_n = fem.Constant(msh, PETSc.ScalarType(tau_n_hat))
    tau_p = fem.Constant(msh, PETSc.ScalarType(tau_p_hat))

    n_hat = n_from_slotboom(psi, phi_n, ni_hat)
    p_hat = p_from_slotboom(psi, phi_p, ni_hat)

    # SRH rate inlined here so we can share the same ni_hat Constant with
    # the Poisson block and avoid UFL-type surprises.
    import math as _math
    n1 = ni_hat * _math.exp(E_t_over_Vt)
    p1 = ni_hat * _math.exp(-E_t_over_Vt)
    R = (n_hat * p_hat - ni_hat * ni_hat) / (
        tau_p * (n_hat + n1) + tau_n * (p_hat + p1)
    )

    rho_hat = p_hat - n_hat + N_hat_fn

    F_psi = (
        L_D2 * eps_r * ufl.inner(ufl.grad(psi), ufl.grad(v_psi)) * ufl.dx
        - rho_hat * v_psi * ufl.dx
    )
    F_phi_n = (
        L0_sq * mu_n_hat * n_hat * ufl.inner(ufl.grad(phi_n), ufl.grad(v_n)) * ufl.dx
        - R * v_n * ufl.dx
    )
    F_phi_p = (
        L0_sq * mu_p_hat * p_hat * ufl.inner(ufl.grad(phi_p), ufl.grad(v_p)) * ufl.dx
        + R * v_p * ufl.dx
    )
    return [F_psi, F_phi_n, F_phi_p]
