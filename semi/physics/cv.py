"""
LF and HF small-signal capacitance for MOS C-V (M14.2).

Both curves are extracted from the same DC sweep:

* LF (quasi-static): both electrons and holes follow the AC signal.
  The full Boltzmann derivative ``-(q/V_t) (n + p)`` enters the
  linearised Poisson Jacobian, so the M14.1 ``mos_cap_ac`` sensitivity
  solve already produces the LF capacitance.

* HF: minority carriers are frozen at their DC value; only the
  majority-carrier Boltzmann response and the dopant background
  contribute to ``dQ/dV``. This is the "high-frequency" curve in
  Hu Fig. 5-18 that saturates at ``Cmin = Cox * Cdep,min / (Cox + Cdep,min)``
  in inversion. We assemble a custom Jacobian whose space-charge term
  drops the minority-carrier exponential and solve for the per-unit-V
  potential perturbation ``delta_psi`` against the same Dirichlet rows
  used by the DC solve.

Axisymmetric (r, z) models pass ``axisymmetric_axis`` through to all
forms so volume integrals carry the radial weight ``r``; the bare 2-pi
factor cancels in the residual but reappears when total charge is
integrated below.
"""
from __future__ import annotations


def _radial_weight(msh, axisymmetric_axis, L0):
    """Radial weight ``r / L0`` (dimensionless) or 1 for planar.

    Matches the normalisation used in
    :func:`semi.physics.poisson.build_equilibrium_poisson_form_mr` so
    the residual / Jacobian / charge forms all share the same scale.
    The factor of ``L0`` is reintroduced in the runner's ``area_factor``
    when the integral is converted into Q per gate area.
    """
    if axisymmetric_axis is None:
        return 1
    import ufl
    from dolfinx import fem
    from petsc4py import PETSc
    x = ufl.SpatialCoordinate(msh)
    return x[int(axisymmetric_axis)] / fem.Constant(
        msh, PETSc.ScalarType(float(L0))
    )


def build_hf_jacobian_form(
    V_psi, psi, N_hat_fn, sc, eps_r_fn, cell_tags, semi_tag,
    majority_carrier: str,
    *, axisymmetric_axis: int | None = None,
):
    """
    Build the bilinear Jacobian form a(d, v) = K_HF d for the HF
    small-signal solve.

    The HF Jacobian is the Frechet derivative of the equilibrium
    Poisson residual *with the minority-carrier Boltzmann term
    dropped*. For a P-body device (majority = holes), the space-charge
    derivative is ``+ ni_hat * exp(-psi)`` (no electron term). For an
    N-body device (majority = electrons), it is ``+ ni_hat * exp(psi)``.
    The dopant background ``N_hat_fn`` is psi-independent and drops
    out of the linearisation either way.
    """
    import ufl
    from dolfinx import fem
    from petsc4py import PETSc

    msh = V_psi.mesh
    v = ufl.TestFunction(V_psi)
    d = ufl.TrialFunction(V_psi)

    L_D2 = fem.Constant(msh, PETSc.ScalarType(sc.lambda2 * sc.L0 ** 2))
    ni_hat = fem.Constant(msh, PETSc.ScalarType(sc.n_i / sc.C0))

    dx_full = ufl.Measure("dx", domain=msh)
    dx_semi = ufl.Measure(
        "dx", domain=msh, subdomain_data=cell_tags, subdomain_id=int(semi_tag),
    )
    weight = _radial_weight(msh, axisymmetric_axis, sc.L0)

    if majority_carrier == "holes":
        # rho_hat = ni_hat (exp(-psi) - exp(psi)) + N_hat
        # d(rho_hat)/d(psi)|_HF = -ni_hat exp(-psi)   (holes only)
        # F = L*grad psi . grad v - rho_hat v ; differentiating in psi
        # gives K = L grad d . grad v - d(rho)/dpsi d v.
        drho_dpsi = -ni_hat * ufl.exp(-psi)
    elif majority_carrier == "electrons":
        # d(rho_hat)/d(psi)|_HF = -ni_hat exp(psi)
        drho_dpsi = -ni_hat * ufl.exp(psi)
    else:  # pragma: no cover - validated upstream
        raise ValueError(f"majority_carrier must be 'holes' or 'electrons', got {majority_carrier!r}")

    a = (
        L_D2 * eps_r_fn * ufl.inner(ufl.grad(d), ufl.grad(v)) * weight * dx_full
        - drho_dpsi * d * v * weight * dx_semi
    )
    return a


def build_lf_charge_sensitivity_form(
    V_psi, psi, sc, cell_tags, semi_tag, delta_psi,
    *, axisymmetric_axis: int | None = None,
):
    """
    Scalar form for the LF differential charge:
        integral_Si [ ni_hat (exp(-psi) + exp(psi)) ] * delta_psi * dx (* r)
    which is the linearisation of the semiconductor charge with both
    carriers responding. The factor of -q*C0 / A_gate is applied
    outside the form when the result is converted to F / m^2.
    """
    import ufl
    from dolfinx import fem
    from petsc4py import PETSc

    msh = V_psi.mesh
    ni_hat = fem.Constant(msh, PETSc.ScalarType(sc.n_i / sc.C0))
    dx_semi = ufl.Measure(
        "dx", domain=msh, subdomain_data=cell_tags, subdomain_id=int(semi_tag),
    )
    weight = _radial_weight(msh, axisymmetric_axis, sc.L0)
    return fem.form(ni_hat * (ufl.exp(-psi) + ufl.exp(psi)) * delta_psi * weight * dx_semi)


def build_hf_charge_sensitivity_form(
    V_psi, psi, sc, cell_tags, semi_tag, delta_psi, majority_carrier: str,
    *, axisymmetric_axis: int | None = None,
):
    """
    Scalar form for the HF differential charge:
        integral_Si [ ni_hat * exp(-psi) ] * delta_psi * dx (* r)   (P-body)
        integral_Si [ ni_hat * exp(psi)  ] * delta_psi * dx (* r)   (N-body)
    Only the majority-carrier Boltzmann derivative contributes; the
    minority carrier is frozen at the DC value.
    """
    import ufl
    from dolfinx import fem
    from petsc4py import PETSc

    msh = V_psi.mesh
    ni_hat = fem.Constant(msh, PETSc.ScalarType(sc.n_i / sc.C0))
    dx_semi = ufl.Measure(
        "dx", domain=msh, subdomain_data=cell_tags, subdomain_id=int(semi_tag),
    )
    weight = _radial_weight(msh, axisymmetric_axis, sc.L0)

    if majority_carrier == "holes":
        integrand = ni_hat * ufl.exp(-psi) * delta_psi
    elif majority_carrier == "electrons":
        integrand = ni_hat * ufl.exp(psi) * delta_psi
    else:  # pragma: no cover
        raise ValueError(f"majority_carrier must be 'holes' or 'electrons', got {majority_carrier!r}")

    return fem.form(integrand * weight * dx_semi)


def build_charge_form(
    V_psi, psi, N_hat_fn, sc, cell_tags, semi_tag,
    *, axisymmetric_axis: int | None = None,
):
    """Total semiconductor space-charge form: integral_Si rho_hat dx (* r)."""
    import ufl
    from dolfinx import fem
    from petsc4py import PETSc

    msh = V_psi.mesh
    ni_hat = fem.Constant(msh, PETSc.ScalarType(sc.n_i / sc.C0))
    dx_semi = ufl.Measure(
        "dx", domain=msh, subdomain_data=cell_tags, subdomain_id=int(semi_tag),
    )
    weight = _radial_weight(msh, axisymmetric_axis, sc.L0)
    rho_hat = ni_hat * (ufl.exp(-psi) - ufl.exp(psi)) + N_hat_fn
    return fem.form(rho_hat * weight * dx_semi)


def infer_majority_carrier(cfg) -> str:
    """
    Determine which carrier type is the majority in the semiconductor
    region from the doping profile. Used when cv_analysis.majority_carrier
    is 'auto' (the default).

    Net acceptor doping (N_A > N_D) -> 'holes' (P-body, n+ gate -> NMOS).
    Net donor doping (N_D > N_A)    -> 'electrons' (N-body, p+ gate -> PMOS).
    """
    n_d_total = 0.0
    n_a_total = 0.0
    for entry in cfg.get("doping", []):
        prof = entry["profile"]
        if prof["type"] == "uniform":
            n_d_total += float(prof.get("N_D", 0.0))
            n_a_total += float(prof.get("N_A", 0.0))
        elif prof["type"] == "step":
            n_d_total += float(prof.get("N_D_left", 0.0)) + float(prof.get("N_D_right", 0.0))
            n_a_total += float(prof.get("N_A_left", 0.0)) + float(prof.get("N_A_right", 0.0))
        elif prof["type"] == "gaussian":
            if prof.get("dopant") == "donor":
                n_d_total += float(prof.get("peak", 0.0))
            else:
                n_a_total += float(prof.get("peak", 0.0))
            n_d_total += float(prof.get("background_N_D", 0.0))
            n_a_total += float(prof.get("background_N_A", 0.0))
    if n_a_total >= n_d_total:
        return "holes"
    return "electrons"


# ----------------------------- analytical reference -----------------------------

def analytical_moscap_metrics(
    *,
    N_body: float,        # m^-3, magnitude (acceptor for P-body, donor for N-body)
    body_type: str,       # 'p' or 'n'
    t_ox: float,          # m
    eps_r_si: float,
    eps_r_ox: float,
    n_i: float,           # m^-3
    T: float,             # K
    phi_ms: float,        # V (work function difference, gate metal minus body intrinsic level)
):
    """
    Closed-form MOSCAP reference values. Returns a dict with:

        Cox          F/m^2  oxide capacitance per area
        W_dmax       m      maximum depletion width at threshold
        Cdep_min     F/m^2  minimum depletion capacitance per area
        Cmin         F/m^2  HF inversion-plateau capacitance per area
        phi_F        V      bulk Fermi potential magnitude
        V_FB         V      flat-band voltage (= phi_ms in this idealised model)
        V_T          V      threshold voltage

    Conventions follow Hu Ch. 5; signs are consistent with V_T > V_FB
    for an NMOS (P-body, electron inversion) and V_T < V_FB for a PMOS.
    """
    from ..constants import EPS0, Q, thermal_voltage

    V_t = thermal_voltage(T)
    eps_si = eps_r_si * EPS0
    eps_ox = eps_r_ox * EPS0
    Cox = eps_ox / t_ox

    phi_F = V_t * _safe_log(N_body / n_i)
    # Surface potential at threshold = 2 phi_F (magnitude).
    phi_st = 2.0 * phi_F
    W_dmax = (4.0 * eps_si * phi_F / (Q * N_body)) ** 0.5
    Cdep_min = eps_si / W_dmax
    Cmin = Cox * Cdep_min / (Cox + Cdep_min)
    Q_dep_max = Q * N_body * W_dmax  # |bulk depletion charge| per area

    if body_type == "p":
        V_T = phi_ms + phi_st + Q_dep_max / Cox
    elif body_type == "n":
        V_T = phi_ms - phi_st - Q_dep_max / Cox
    else:  # pragma: no cover
        raise ValueError(f"body_type must be 'p' or 'n', got {body_type!r}")

    return {
        "Cox": Cox,
        "W_dmax": W_dmax,
        "Cdep_min": Cdep_min,
        "Cmin": Cmin,
        "phi_F": phi_F,
        "V_FB": phi_ms,
        "V_T": V_T,
        "Q_dep_max": Q_dep_max,
        "V_t": V_t,
    }


def _safe_log(x):
    import math
    return math.log(max(x, 1.0e-300))
