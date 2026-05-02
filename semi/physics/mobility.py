"""
Field-dependent carrier mobility models (M16.1).

Layer 4 (FEM). Imports of dolfinx, ufl, and petsc4py are deferred to
function bodies so the pure-Python core stays dolfinx-free per
PLAN.md invariant 4 (covered by tests/test_lazy_imports.py).

Currently supported:
    constant       mu(F) = mu0 (no field dependence). The V&V reference;
                   bit-identical to pre-M16.1.
    caughey_thomas mu(F) = mu0 / (1 + (mu0 * F_par / vsat)^beta)^(1/beta)
                   Closed-form velocity saturation. F_par is the
                   magnitude of the gradient of the carrier-specific
                   quasi-Fermi potential; mu0 is the low-field
                   electron / hole mobility.

Asymptotics for caughey_thomas:
    F_par -> 0:    mu(F) -> mu0
    F_par -> inf:  mu(F) * F_par -> vsat / mu0 * mu0 = vsat (drift
                                                               saturation)

Scaling. The DD form integrates `mu_hat * carrier_density * grad(phi)`
against `grad(test)` in scaled units (Slotboom form, see
docs/PHYSICS.md section 2.5 and ADR 0004). The mesh stays in physical
meters (PLAN.md invariant 3), so `ufl.grad(phi_n_hat)` has units of
1/m; the field magnitude is `F = V_t * grad(phi_hat)` with units V/m,
and the SI velocity argument of the saturation formula is
`mu0 * F = mu0 * V_t * grad(phi_hat)`. Dividing by `vsat` gives a
dimensionless ratio when

    vsat_for_form = vsat_SI / (sc.mu0 * sc.V0)        [units 1/m]

so the closed form

    mu_hat(F) = mu0_hat / (1 + (mu0_hat * F_par / vsat_for_form)^beta)^(1/beta)

is dimensionless. `F_par` here is `sqrt(grad(phi).grad(phi) + eps)`
straight from UFL (units 1/m). See `caughey_thomas_vsat_for_form` for
the conversion helper.

The regularization epsilon is tiny (1e-30 in scaled-form units) so
that the `(F/vsat)^beta` term is below the floor of float precision
and shifts mu only at machine precision; near the equilibrium initial
guess where `grad(phi) = 0` it averts a `0^0` UFL nan.
"""
from __future__ import annotations

from typing import Any

# Regularization for sqrt(grad(phi).grad(phi)) at quadrature points where
# the gradient is identically zero (equilibrium initial guess). 1e-30 in
# (1/m)^2 units leaves the dimensionless ratio (mu0 * F / vsat)^beta
# below double-precision underflow.
_F_PAR_EPS_SQ: float = 1.0e-30


def caughey_thomas_mu(mu0, F_par, vsat, beta):
    """
    Caughey-Thomas field-dependent mobility, dimensionless or UFL form.

        mu(F) = mu0 / (1 + (mu0 * F_par / vsat)^beta)^(1/beta)

    Parameters
    ----------
    mu0 : float | ufl expression | fem.Constant
        Low-field mobility. Pass either a Python scalar (for the unit
        tests) or a UFL-compatible quantity in scaled units (for the
        DD form builder).
    F_par : float | ufl expression
        Magnitude of the parallel field. In SI use V/m and a SI vsat;
        in scaled units use the gradient of the scaled quasi-Fermi
        potential (1/m) and pass the corresponding scaled vsat (also
        1/m, see `caughey_thomas_vsat_for_form`).
    vsat : float | fem.Constant
        Saturation velocity. Units must match F_par * mu0.
    beta : float | fem.Constant
        Caughey-Thomas exponent (Si defaults: beta_n = 2, beta_p = 1).

    Returns
    -------
    Same kind as the inputs (Python scalar or UFL expression).
    """
    ratio = (mu0 * F_par) / vsat
    return mu0 / (1.0 + ratio ** beta) ** (1.0 / beta)


def constant_mu(mu0):
    """Identity wrapper. Keeps `build_mobility()` returning the same
    shape (UFL-compatible object) for either model branch. The constant
    branch reuses the input `fem.Constant` so byte-equivalence with
    pre-M16.1 is preserved at the form-compilation layer.
    """
    return mu0


def caughey_thomas_vsat_for_form(vsat_cm_per_s: float, sc) -> float:
    """
    Convert a saturation velocity in cm/s to the dimensionless form
    consumed by `caughey_thomas_mu` against `F_par = sqrt(grad(phi_hat)
    . grad(phi_hat))` in 1/m (the gradient of the scaled quasi-Fermi
    potential under the kronos-semi mesh-in-meters convention).

    Derivation (cite docs/PHYSICS.md section 2.5):
        physical field    F = V_t * grad(phi_hat)          [V/m]
        physical velocity v = mu0_SI * F                   [m/s]
        dimensionless arg (mu0_SI * F) / vsat_SI
                      = (mu0_hat * mu0_ref * V_t * grad(phi_hat)) / vsat_SI
                      = mu0_hat * grad(phi_hat) / [vsat_SI / (mu0_ref * V_t)]
        so the scaled vsat that closes the dimensionless ratio is
            vsat_for_form = vsat_SI / (sc.mu0 * sc.V0)     [1/m]

    `vsat_cm_per_s` is the JSON value in cm/s; convert to m/s via 1e-2
    before dividing.
    """
    vsat_SI = float(vsat_cm_per_s) * 1.0e-2
    return vsat_SI / (sc.mu0 * sc.V0)


def build_mobility_expressions(
    mobility_cfg: dict[str, Any] | None,
    phi_n,
    phi_p,
    mu_n_over_mu0: float,
    mu_p_over_mu0: float,
    sc,
) -> tuple[Any, Any, str]:
    """
    Dispatch on `mobility_cfg["model"]` and return UFL-compatible
    mobility expressions substitutable for the legacy
    `fem.Constant(mu_n_over_mu0)` / `fem.Constant(mu_p_over_mu0)` in
    the DD residual builders.

    Parameters
    ----------
    mobility_cfg : dict or None
        The validated `cfg["physics"]["mobility"]` sub-dict, or `None`
        which is equivalent to `{"model": "constant"}` (preserves
        pre-M16.1 byte-identity).
    phi_n, phi_p : dolfinx.fem.Function
        The carrier-specific scaled quasi-Fermi potentials. The
        Caughey-Thomas branch reads `ufl.grad(phi_n)` for the electron
        continuity row and `ufl.grad(phi_p)` for the hole continuity
        row (ADR 0004 Slotboom flux form).
    mu_n_over_mu0, mu_p_over_mu0 : float
        Low-field mobility ratios (mu / sc.mu0) already computed by
        the caller.
    sc : semi.scaling.Scaling
        Scaling object; used by `caughey_thomas_vsat_for_form`.

    Returns
    -------
    (mu_n_expr, mu_p_expr, model_name) : tuple
        UFL or fem.Constant objects substitutable for the legacy
        `fem.Constant(mu_n_over_mu0)` / `fem.Constant(mu_p_over_mu0)`
        in the DD residual. `model_name` echoes the dispatch for
        diagnostics.

    Branches
    --------
    constant (default):
        Returns `(fem.Constant(msh, mu_n_over_mu0),
                  fem.Constant(msh, mu_p_over_mu0))`. Bit-identical to
        pre-M16.1 since the existing DD form already wraps the inputs
        in fem.Constants of the same scalar values.
    caughey_thomas:
        Builds `caughey_thomas_mu(mu0_const, F_par, vsat_for_form,
        beta)` for each carrier where `F_par` is
        `sqrt(grad(phi) . grad(phi) + eps)` of the carrier-specific
        scaled quasi-Fermi potential.
    """
    import ufl
    from dolfinx import fem
    from petsc4py import PETSc

    mob = mobility_cfg or {}
    model = mob.get("model", "constant")

    # phi_n / phi_p may live on the parent mesh or on the semiconductor
    # submesh (multi-region variant); pull the mesh off each carrier's
    # function space directly so the constant fem.Constants land on
    # the same mesh as the legacy code did.
    msh_n = phi_n.function_space.mesh
    msh_p = phi_p.function_space.mesh

    if model == "constant":
        mu_n_expr = fem.Constant(msh_n, PETSc.ScalarType(mu_n_over_mu0))
        mu_p_expr = fem.Constant(msh_p, PETSc.ScalarType(mu_p_over_mu0))
        return mu_n_expr, mu_p_expr, "constant"

    if model != "caughey_thomas":
        raise ValueError(
            f"physics.mobility.model={model!r} not supported; "
            "expected 'constant' or 'caughey_thomas'"
        )

    vsat_n = caughey_thomas_vsat_for_form(mob.get("vsat_n", 1.0e7), sc)
    vsat_p = caughey_thomas_vsat_for_form(mob.get("vsat_p", 8.0e6), sc)
    beta_n = float(mob.get("beta_n", 2.0))
    beta_p = float(mob.get("beta_p", 1.0))

    mu0_n = fem.Constant(msh_n, PETSc.ScalarType(mu_n_over_mu0))
    mu0_p = fem.Constant(msh_p, PETSc.ScalarType(mu_p_over_mu0))
    vsat_n_c = fem.Constant(msh_n, PETSc.ScalarType(vsat_n))
    vsat_p_c = fem.Constant(msh_p, PETSc.ScalarType(vsat_p))
    beta_n_c = fem.Constant(msh_n, PETSc.ScalarType(beta_n))
    beta_p_c = fem.Constant(msh_p, PETSc.ScalarType(beta_p))
    eps_n = fem.Constant(msh_n, PETSc.ScalarType(_F_PAR_EPS_SQ))
    eps_p = fem.Constant(msh_p, PETSc.ScalarType(_F_PAR_EPS_SQ))

    grad_phi_n = ufl.grad(phi_n)
    grad_phi_p = ufl.grad(phi_p)
    F_par_n = ufl.sqrt(ufl.dot(grad_phi_n, grad_phi_n) + eps_n)
    F_par_p = ufl.sqrt(ufl.dot(grad_phi_p, grad_phi_p) + eps_p)

    mu_n_expr = caughey_thomas_mu(mu0_n, F_par_n, vsat_n_c, beta_n_c)
    mu_p_expr = caughey_thomas_mu(mu0_p, F_par_p, vsat_p_c, beta_p_c)
    return mu_n_expr, mu_p_expr, "caughey_thomas"
