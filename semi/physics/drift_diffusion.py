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
    lets us re-use the M1 equilibrium Poisson solver for the initial
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
    eps_r,
    mu_n_over_mu0: float,
    mu_p_over_mu0: float,
    tau_n_hat: float,
    tau_p_hat: float,
    E_t_over_Vt: float = 0.0,
    mobility_cfg: dict | None = None,
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
    eps_r : float | dolfinx.fem.Function
        Relative permittivity. Scalar (single-region fast path,
        byte-identical with M2-M5) or a DG0 cellwise Function on the
        parent mesh for the multi-region Poisson coefficient jump.
    mu_n_over_mu0, mu_p_over_mu0 : float
        Low-field mobility ratios (mu / sc.mu0). Under the
        caughey_thomas branch these are the `mu0` arguments to the
        closed form (the low-field reference).
    tau_n_hat, tau_p_hat : float
        Scaled lifetimes, tau / t0.
    E_t_over_Vt : float
        Trap level relative to the intrinsic level, divided by V_t.
    mobility_cfg : dict, optional
        The `cfg["physics"]["mobility"]` sub-dict. `None` (default) or
        `{"model": "constant"}` is bit-identical to pre-M16.1.
        `{"model": "caughey_thomas", ...}` substitutes a closed-form
        velocity-saturation expression (carrier-specific
        |grad(phi_n)| / |grad(phi_p)| under ADR 0004 Slotboom flux
        form; see docs/PHYSICS.md section 1.3 for the flux derivation
        and `semi.physics.mobility` for the closed form).

    Returns
    -------
    list of ufl.Form
        [F_psi, F_phi_n, F_phi_p] -- pass to the blocked NonlinearProblem.
    """
    import ufl
    from dolfinx import fem
    from petsc4py import PETSc

    from .mobility import build_mobility_expressions
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
    if isinstance(eps_r, (int, float)):
        eps_r_ufl = fem.Constant(msh, PETSc.ScalarType(float(eps_r)))
    else:
        eps_r_ufl = eps_r
    ni_hat = fem.Constant(msh, PETSc.ScalarType(sc.n_i / sc.C0))
    mu_n_hat, mu_p_hat, _mob_model = build_mobility_expressions(
        mobility_cfg, phi_n, phi_p, mu_n_over_mu0, mu_p_over_mu0, sc,
    )
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
        L_D2 * eps_r_ufl * ufl.inner(ufl.grad(psi), ufl.grad(v_psi)) * ufl.dx
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


@dataclass
class DDBlockSpacesMR:
    """
    Multi-region block spaces for the MOS capacitor and similar devices.

    V_psi lives on the parent mesh (spans silicon + oxide); V_phi_n and
    V_phi_p live on the semiconductor submesh only (Slotboom variables
    are ill-defined in an ideal insulator; see docs/mos_derivation.md
    section 2.2). `entity_map` is the parent<->submesh cell mapping
    returned by `dolfinx.mesh.create_submesh`; it threads through the
    mixed-domain form compilation (`fem.form(..., entity_maps=[em])`)
    and into the SNES solver (`NonlinearProblem(..., entity_maps=[em])`).
    """
    V_psi: Any
    V_phi_n: Any
    V_phi_p: Any
    psi: Any
    phi_n: Any
    phi_p: Any
    parent_mesh: Any
    submesh: Any
    entity_map: Any


def make_dd_block_spaces_mr(msh, submesh, entity_map) -> DDBlockSpacesMR:
    """Create the multi-region P1 spaces.

    V_psi lives on the parent mesh; V_phi_n and V_phi_p live on the
    semiconductor submesh. Unknown Functions are instantiated on their
    respective spaces.
    """
    from dolfinx import fem

    V_psi = fem.functionspace(msh, ("Lagrange", 1))
    V_phi_n = fem.functionspace(submesh, ("Lagrange", 1))
    V_phi_p = fem.functionspace(submesh, ("Lagrange", 1))
    psi = fem.Function(V_psi, name="psi_hat")
    phi_n = fem.Function(V_phi_n, name="phi_n_hat")
    phi_p = fem.Function(V_phi_p, name="phi_p_hat")
    return DDBlockSpacesMR(
        V_psi=V_psi, V_phi_n=V_phi_n, V_phi_p=V_phi_p,
        psi=psi, phi_n=phi_n, phi_p=phi_p,
        parent_mesh=msh, submesh=submesh, entity_map=entity_map,
    )


def build_dd_block_residual_mr(
    spaces: DDBlockSpacesMR,
    N_hat_fn,
    sc,
    eps_r,
    mu_n_over_mu0: float,
    mu_p_over_mu0: float,
    tau_n_hat: float,
    tau_p_hat: float,
    cell_tags,
    semi_tag: int,
    E_t_over_Vt: float = 0.0,
    mobility_cfg: dict | None = None,
):
    """Multi-region (submesh-based) block residual for the coupled DD system.

    Structural differences from the single-region path:

    - The Poisson stiffness term integrates over the full parent mesh
      (oxide included), using the cellwise DG0 eps_r(x) Function. The
      space-charge term integrates only over silicon via the parent
      mesh's `dx(subdomain_data=cell_tags, subdomain_id=semi_tag)`
      measure, which is restricted to semiconductor cells where
      (phi_n, phi_p) live on the submesh.
    - The continuity-equation residuals integrate on the submesh
      (`Measure('dx', domain=submesh)`), reading psi from the parent
      mesh.

    The mixed-domain mapping threads through `NonlinearProblem`'s
    `entity_maps` argument; callers must pass `entity_maps=[spaces.entity_map]`
    to `solve_nonlinear_block`. The derivation is in
    docs/mos_derivation.md section 4.3.

    Returns
    -------
    list of ufl.Form
        [F_psi, F_phi_n, F_phi_p] -- pass to `solve_nonlinear_block`
        with `entity_maps=[spaces.entity_map]`.
    """
    import ufl
    from dolfinx import fem
    from petsc4py import PETSc

    from .mobility import build_mobility_expressions
    from .slotboom import n_from_slotboom, p_from_slotboom

    psi = spaces.psi
    phi_n = spaces.phi_n
    phi_p = spaces.phi_p
    msh = spaces.parent_mesh
    submesh = spaces.submesh

    v_psi = ufl.TestFunction(spaces.V_psi)
    v_n = ufl.TestFunction(spaces.V_phi_n)
    v_p = ufl.TestFunction(spaces.V_phi_p)

    L_D2 = fem.Constant(msh, PETSc.ScalarType(sc.lambda2 * sc.L0 ** 2))
    L0_sq = fem.Constant(submesh, PETSc.ScalarType(sc.L0 ** 2))
    if isinstance(eps_r, (int, float)):
        eps_r_ufl = fem.Constant(msh, PETSc.ScalarType(float(eps_r)))
    else:
        eps_r_ufl = eps_r
    ni_hat_parent = fem.Constant(msh, PETSc.ScalarType(sc.n_i / sc.C0))
    ni_hat_sub = fem.Constant(submesh, PETSc.ScalarType(sc.n_i / sc.C0))
    mu_n_hat, mu_p_hat, _mob_model = build_mobility_expressions(
        mobility_cfg, phi_n, phi_p, mu_n_over_mu0, mu_p_over_mu0, sc,
    )
    tau_n = fem.Constant(submesh, PETSc.ScalarType(tau_n_hat))
    tau_p = fem.Constant(submesh, PETSc.ScalarType(tau_p_hat))

    # Measures
    dx_parent = ufl.Measure("dx", domain=msh, subdomain_data=cell_tags)
    dx_semi_parent = dx_parent(int(semi_tag))
    dx_sub = ufl.Measure("dx", domain=submesh)

    # Parent-mesh Slotboom (for the Poisson source term, which is integrated
    # over silicon cells of the parent mesh). phi_n, phi_p live on the
    # submesh; entity_maps threads them through at form compilation.
    n_hat_parent = n_from_slotboom(psi, phi_n, ni_hat_parent)
    p_hat_parent = p_from_slotboom(psi, phi_p, ni_hat_parent)
    rho_hat = p_hat_parent - n_hat_parent + N_hat_fn

    # Submesh-mesh Slotboom (for the continuity blocks). Here psi is the
    # parent-mesh Function and entity_maps pulls it onto submesh cells.
    n_hat_sub = n_from_slotboom(psi, phi_n, ni_hat_sub)
    p_hat_sub = p_from_slotboom(psi, phi_p, ni_hat_sub)

    import math as _math
    n1 = ni_hat_sub * _math.exp(E_t_over_Vt)
    p1 = ni_hat_sub * _math.exp(-E_t_over_Vt)
    R = (n_hat_sub * p_hat_sub - ni_hat_sub * ni_hat_sub) / (
        tau_p * (n_hat_sub + n1) + tau_n * (p_hat_sub + p1)
    )

    F_psi = (
        L_D2 * eps_r_ufl * ufl.inner(ufl.grad(psi), ufl.grad(v_psi)) * dx_parent
        - rho_hat * v_psi * dx_semi_parent
    )
    F_phi_n = (
        L0_sq * mu_n_hat * n_hat_sub
        * ufl.inner(ufl.grad(phi_n), ufl.grad(v_n)) * dx_sub
        - R * v_n * dx_sub
    )
    F_phi_p = (
        L0_sq * mu_p_hat * p_hat_sub
        * ufl.inner(ufl.grad(phi_p), ufl.grad(v_p)) * dx_sub
        + R * v_p * dx_sub
    )
    return [F_psi, F_phi_n, F_phi_p]
