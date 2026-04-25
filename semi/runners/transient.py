"""
Transient drift-diffusion runner (BDF1 / BDF2).

Solves the coupled (psi, phi_n, phi_p) system in **Slotboom form** with
backward differentiation formula (BDF) time integration. The carrier
densities are recovered pointwise from the Slotboom relations:

    n = n_i exp(psi - phi_n)        [scaled, V_t = 1]
    p = n_i exp(phi_p - psi)

Because n and p are exponentials of the primary unknowns they are
strictly positive at every Newton iterate, eliminating the negative-
density failure modes of the (psi, n, p) primary-density formulation.

Time integration
----------------
BDF order is configurable (1 = backward Euler, 2 = BDF2). For BDF2, a
single BDF1 step seeds the two-level history before switching. The
time derivative on the carrier density is expanded by the chain rule
into temporal derivatives of the Slotboom unknowns:

    dn/dt = n * (dpsi/dt - dphi_n/dt)
    dp/dt = p * (dphi_p/dt - dpsi/dt)

with each (d/dt) approximated by BDF on the corresponding primary
unknown. The history sources f_hist_psi, f_hist_phi_n, f_hist_phi_p
carry the past values: f_hist_u = sum_{k=1}^K alpha_k/dt * u^{n+1-k}.
See docs/adr/0010-bdf-time-integration.md.

Residual form (scaled units)
-----------------------------
Poisson:
    L_D^2 eps_r grad(psi) . grad(v_psi)  -  (p - n + N) * v_psi  =  0

Electron continuity (Slotboom + chain-rule):
    n * (BDF[psi] - BDF[phi_n]) * v_n
    +  L_0^2 mu_n n grad(phi_n) . grad(v_n)
    -  R v_n  =  0

Hole continuity (Slotboom + chain-rule):
    p * (BDF[phi_p] - BDF[psi]) * v_p
    +  L_0^2 mu_p p grad(phi_p) . grad(v_p)
    +  R v_p  =  0

with n = n_i exp(psi - phi_n), p = n_i exp(phi_p - psi), and
BDF[u] = alpha_0/dt * u + f_hist_u.

Boundary conditions (ohmic contacts)
--------------------------------------
At an ohmic contact under applied bias V (scaled V_hat = V / V_t):
    psi_bc   = arcsinh(N_net / (2 n_i)) + V_hat
    phi_n_bc = phi_p_bc = V_hat                 [Shockley boundary]

Because (psi - phi_n) and (phi_p - psi) are independent of V at an
ohmic contact, the carrier densities at the contact remain pinned to
their equilibrium majority/minority values throughout the transient.
"""
from __future__ import annotations

from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from ..physics.drift_diffusion import DDBlockSpaces
    from ..results import TransientResult


def run_transient(
    cfg: dict[str, Any],
    *,
    progress_callback=None,
) -> TransientResult:
    """
    Coupled transient drift-diffusion solver (Slotboom + chain rule).

    Runs BDF1 or BDF2 time integration on the (psi, phi_n, phi_p)
    Slotboom block system, starting from V=0 equilibrium at t=0 and
    applying any bias specified in ``cfg["contacts"]`` from t=0+.

    Parameters
    ----------
    cfg : dict
        Validated JSON config dict. ``solver.type`` must be
        ``"transient"``. Required solver sub-keys:
        - ``t_end`` (float): simulation end time in seconds.
        - ``dt`` (float): fixed time step in seconds.
        - ``order`` (int, default 2): BDF order (1 or 2).
        - ``max_steps`` (int, default 10000): safety cap on iterations.
    progress_callback : callable or None
        If provided, called after each successful timestep with a dict::

            {"type": "step_done", "step": n, "t": t_current,
             "iterations": n_snes_iters}

    Returns
    -------
    TransientResult
    """
    import copy

    import ufl
    from dolfinx import fem

    from ..bcs import build_dd_dirichlet_bcs, resolve_contacts
    from ..doping import build_profile
    from ..fem.mass import assemble_lumped_mass
    from ..mesh import build_mesh
    from ..physics.drift_diffusion import make_dd_block_spaces
    from ..postprocess import (
        evaluate_partial_currents,
        fmt_tag,
        resolve_contact_facets,
    )
    from ..results import TransientResult
    from ..runners.equilibrium import run_equilibrium
    from ..scaling import make_scaling_from_config
    from ..solver import solve_nonlinear_block
    from ..timestepping import BDFCoefficients
    from ._common import reference_material

    # ------------------------------------------------------------------
    # Parse solver parameters
    # ------------------------------------------------------------------
    solver_cfg = cfg.get("solver", {})
    t_end = float(solver_cfg["t_end"])
    dt_val = float(solver_cfg["dt"])
    order = int(solver_cfg.get("order", 2))
    max_steps = int(solver_cfg.get("max_steps", 10000))
    output_every = int(solver_cfg.get("output_every", 50))

    bdf = BDFCoefficients(order)
    alpha_0 = bdf.coeffs[0]

    # SNES tolerances (same as bias_sweep defaults from ADR 0008)
    snes_opts = solver_cfg.get("snes", {}) or {}
    snes_petsc_options = {
        "snes_rtol": float(snes_opts.get("rtol", 1.0e-10)),
        "snes_atol": float(snes_opts.get("atol", 1.0e-7)),
        "snes_stol": float(snes_opts.get("stol", 1.0e-14)),
        "snes_max_it": int(snes_opts.get("max_it", 100)),
    }

    # ------------------------------------------------------------------
    # Build mesh, scaling, doping
    # ------------------------------------------------------------------
    ref_mat = reference_material(cfg)
    sc = make_scaling_from_config(cfg, ref_mat)

    msh, _cell_tags, facet_tags = build_mesh(cfg)

    N_raw_fn = build_profile(cfg["doping"])

    phys = cfg.get("physics", {})
    mob = phys.get("mobility", {})
    mu_n_SI = float(mob.get("mu_n", 1400.0)) * 1.0e-4
    mu_p_SI = float(mob.get("mu_p", 450.0)) * 1.0e-4
    mu_n_hat = mu_n_SI / sc.mu0
    mu_p_hat = mu_p_SI / sc.mu0

    rec = phys.get("recombination", {})
    tau_n_s = float(rec.get("tau_n", 1.0e-7))
    tau_p_s = float(rec.get("tau_p", 1.0e-7))
    E_t_eV = float(rec.get("E_t", 0.0))
    tau_n_hat = tau_n_s / sc.t0
    tau_p_hat = tau_p_s / sc.t0
    E_t_over_Vt = E_t_eV / sc.V0

    # ------------------------------------------------------------------
    # Slotboom block spaces and unknowns
    # ------------------------------------------------------------------
    spaces = make_dd_block_spaces(msh)
    V_psi = spaces.V_psi
    V_phi_n = spaces.V_phi_n
    V_phi_p = spaces.V_phi_p
    psi = spaces.psi
    phi_n = spaces.phi_n
    phi_p = spaces.phi_p

    N_hat_fn = fem.Function(V_psi, name="N_net_hat")
    N_hat_fn.interpolate(lambda x: N_raw_fn(x) / sc.C0)

    # ------------------------------------------------------------------
    # Step 0: Solve equilibrium at t=0 to obtain initial psi
    # ------------------------------------------------------------------
    eq_cfg = copy.deepcopy(cfg)
    eq_cfg["solver"] = {"type": "equilibrium"}
    for c in eq_cfg["contacts"]:
        c["voltage"] = 0.0
        c.pop("voltage_sweep", None)
    eq_result = run_equilibrium(eq_cfg)

    psi.x.array[:] = eq_result.psi.x.array
    psi.x.scatter_forward()

    # At V=0 thermal equilibrium the quasi-Fermi levels coincide with
    # the (arbitrary) reference and we set both to zero. Carrier
    # densities are then n = n_i exp(psi), p = n_i exp(-psi) by the
    # Slotboom relations -- no separate carrier initialisation needed.
    phi_n.x.array[:] = 0.0
    phi_p.x.array[:] = 0.0
    phi_n.x.scatter_forward()
    phi_p.x.scatter_forward()

    # ------------------------------------------------------------------
    # Lumped mass diagonals (kept for diagnostic compatibility with the
    # density-form runner; not used by the Slotboom residual itself).
    # ------------------------------------------------------------------
    dx = ufl.Measure("dx", domain=msh)
    _M_n_diag, _M_p_diag = assemble_lumped_mass(V_phi_n, V_phi_p, dx)

    # ------------------------------------------------------------------
    # Build contact facet info for IV recording
    # ------------------------------------------------------------------
    contact_facet_infos: dict[str, Any] = {}
    for c in cfg["contacts"]:
        if c["type"] != "ohmic":
            continue
        info = resolve_contact_facets(cfg, msh, facet_tags, c["name"])
        contact_facet_infos[c["name"]] = info

    # ------------------------------------------------------------------
    # Static voltages (no sweep in transient)
    # ------------------------------------------------------------------
    static_voltages: dict[str, float] = {}
    for c in cfg["contacts"]:
        if c["type"] != "ohmic":
            continue
        static_voltages[c["name"]] = float(c.get("voltage", 0.0))

    # ------------------------------------------------------------------
    # Build the transient residual
    # ------------------------------------------------------------------
    f_hist_psi = fem.Function(V_psi, name="hist_psi")
    f_hist_phi_n = fem.Function(V_phi_n, name="hist_phi_n")
    f_hist_phi_p = fem.Function(V_phi_p, name="hist_phi_p")
    f_hist_psi.x.array[:] = 0.0
    f_hist_phi_n.x.array[:] = 0.0
    f_hist_phi_p.x.array[:] = 0.0

    # dt_hat = dt / t0 is the dimensionless time step consistent with the
    # scaled spatial coefficients (L0^2 * mu_hat, etc.) in the residual.
    dt_hat = dt_val / sc.t0

    from petsc4py import PETSc

    dt_const = fem.Constant(msh, PETSc.ScalarType(dt_hat))
    alpha0_const = fem.Constant(msh, PETSc.ScalarType(alpha_0))

    F_list = _build_transient_residual(
        psi, phi_n, phi_p, N_hat_fn,
        f_hist_psi, f_hist_phi_n, f_hist_phi_p,
        spaces, sc, ref_mat.epsilon_r,
        mu_n_hat, mu_p_hat, tau_n_hat, tau_p_hat, E_t_over_Vt,
        dt_const, alpha0_const,
    )

    # ------------------------------------------------------------------
    # Helper to build Slotboom (psi, phi_n, phi_p) BCs at t+dt
    # ------------------------------------------------------------------
    def _build_transient_bcs(voltages: dict[str, float]) -> list:
        contacts = resolve_contacts(cfg, facet_tags=facet_tags,
                                    voltages=voltages)
        return build_dd_dirichlet_bcs(
            spaces, msh, facet_tags, contacts, sc, ref_mat, N_raw_fn,
        )

    def _apply_bc_values(bcs: list):
        space_to_fn = {
            id(V_psi): psi,
            id(V_phi_n): phi_n,
            id(V_phi_p): phi_p,
        }
        for bc in bcs:
            fn = space_to_fn.get(id(bc.function_space))
            if fn is not None:
                bc.set(fn.x.array)
        for fn in (psi, phi_n, phi_p):
            fn.x.scatter_forward()

    # ------------------------------------------------------------------
    # History storage (numpy arrays, most-recent LAST). We track the
    # primary Slotboom unknowns; the chain-rule temporal term applies
    # BDF directly to (psi, phi_n, phi_p), so f_hist_u carries past
    # values of u itself.
    # ------------------------------------------------------------------
    psi_hist: list[np.ndarray] = [psi.x.array.copy()]
    phi_n_hist: list[np.ndarray] = [phi_n.x.array.copy()]
    phi_p_hist: list[np.ndarray] = [phi_p.x.array.copy()]

    def _trim_history():
        while len(psi_hist) > order:
            psi_hist.pop(0)
        while len(phi_n_hist) > order:
            phi_n_hist.pop(0)
        while len(phi_p_hist) > order:
            phi_p_hist.pop(0)

    # ------------------------------------------------------------------
    # IV recording helper
    # ------------------------------------------------------------------
    def _record_all_iv(t_val: float, iv_rows: list):
        for c in cfg["contacts"]:
            if c["type"] != "ohmic":
                continue
            cname = c["name"]
            V_applied = static_voltages.get(cname, 0.0)
            finfo = contact_facet_infos.get(cname)
            if finfo is None:
                continue

            J_n, J_p = evaluate_partial_currents(
                spaces, sc, ref_mat, finfo, mu_n_SI, mu_p_SI,
            )
            J_total = J_n + J_p

            iv_rows.append({
                "t": float(t_val),
                "contact": cname,
                "V": float(V_applied),
                "J": float(J_total),
                "J_n": float(J_n),
                "J_p": float(J_p),
            })

    # ------------------------------------------------------------------
    # Time loop
    # ------------------------------------------------------------------
    iv_rows: list[dict[str, Any]] = []
    t_vals: list[float] = []
    fields_out: dict[str, list[np.ndarray]] = {"psi": [], "n": [], "p": []}
    snap_t: list[float] = []
    n_steps_taken = 0
    step_count = 0

    t_current = 0.0
    _record_all_iv(t_current, iv_rows)
    t_vals.append(t_current)

    bdf1 = BDFCoefficients(1)

    while t_current < t_end - 0.5 * dt_val and step_count < max_steps:
        t_next = t_current + dt_val
        step_count += 1

        # Decide BDF order for this step (BDF1 seeds BDF2 history)
        if order == 2 and len(psi_hist) >= 2:
            effective_order = 2
            use_bdf = bdf
        else:
            effective_order = 1
            use_bdf = bdf1

        alpha_k_use = use_bdf.coeffs[0]
        alpha0_const.value = PETSc.ScalarType(alpha_k_use)

        # Update history source functions:
        #   f_hist_u[i] = sum_{k=1}^K alpha_k * u_hist[-k][i] / dt_hat
        hist_psi_arr = np.zeros(len(psi_hist[-1]))
        hist_phi_n_arr = np.zeros(len(phi_n_hist[-1]))
        hist_phi_p_arr = np.zeros(len(phi_p_hist[-1]))
        for k in range(1, effective_order + 1):
            idx = -k
            hist_psi_arr += use_bdf.coeffs[k] * psi_hist[idx]
            hist_phi_n_arr += use_bdf.coeffs[k] * phi_n_hist[idx]
            hist_phi_p_arr += use_bdf.coeffs[k] * phi_p_hist[idx]
        f_hist_psi.x.array[:] = hist_psi_arr / dt_hat
        f_hist_phi_n.x.array[:] = hist_phi_n_arr / dt_hat
        f_hist_phi_p.x.array[:] = hist_phi_p_arr / dt_hat
        f_hist_psi.x.scatter_forward()
        f_hist_phi_n.x.scatter_forward()
        f_hist_phi_p.x.scatter_forward()

        bcs = _build_transient_bcs(static_voltages)
        _apply_bc_values(bcs)

        tag = fmt_tag(t_next)
        info = solve_nonlinear_block(
            F_list, [psi, phi_n, phi_p], bcs,
            prefix=f"{cfg['name']}_tr_{tag}_",
            petsc_options=snes_petsc_options,
        )

        if not info["converged"]:
            raise RuntimeError(
                f"Transient SNES failed at t={t_next:.3e} s "
                f"(step {step_count}); reason={info['reason']}"
            )

        # Append the converged primary unknowns to history
        psi_hist.append(psi.x.array.copy())
        phi_n_hist.append(phi_n.x.array.copy())
        phi_p_hist.append(phi_p.x.array.copy())
        _trim_history()

        t_current = t_next
        n_steps_taken += 1

        _record_all_iv(t_current, iv_rows)
        t_vals.append(t_current)

        if step_count % output_every == 0 or abs(t_current - t_end) < 0.5 * dt_val:
            n_phys = ref_mat.n_i * np.exp(psi.x.array - phi_n.x.array)
            p_phys = ref_mat.n_i * np.exp(phi_p.x.array - psi.x.array)
            fields_out["psi"].append(psi.x.array.copy() * sc.V0)
            fields_out["n"].append(n_phys)
            fields_out["p"].append(p_phys)
            snap_t.append(t_current)

        if progress_callback is not None:
            progress_callback({
                "type": "step_done",
                "step": step_count,
                "t": float(t_current),
                "iterations": int(info.get("iterations", 0)),
            })

    x_dof = V_psi.tabulate_dof_coordinates()

    return TransientResult(
        t=t_vals,
        iv=iv_rows,
        fields=fields_out,
        meta={
            "order": order,
            "dt": dt_val,
            "n_steps_taken": n_steps_taken,
            "n_failed_steps": 0,
            "snap_t": snap_t,
        },
        x_dof=x_dof,
    )


def _build_transient_residual(
    psi, phi_n, phi_p, N_hat_fn,
    f_hist_psi, f_hist_phi_n, f_hist_phi_p,
    spaces: DDBlockSpaces,
    sc,
    eps_r,
    mu_n_over_mu0: float,
    mu_p_over_mu0: float,
    tau_n_hat: float,
    tau_p_hat: float,
    E_t_over_Vt: float,
    dt_const,
    alpha0_const,
):
    """
    Build the transient three-block residual in Slotboom form with
    chain-rule temporal discretisation.

    Parameters
    ----------
    psi, phi_n, phi_p : dolfinx.fem.Function
        Slotboom unknowns (will be solved for).
    N_hat_fn : dolfinx.fem.Function
        Scaled net doping.
    f_hist_psi, f_hist_phi_n, f_hist_phi_p : dolfinx.fem.Function
        Per-DOF history sources for the primary unknowns. Updated
        externally each timestep:
            f_hist_u = sum_{k=1}^K alpha_k / dt_hat * u^{n+1-k}.
    spaces : DDBlockSpaces
    sc : Scaling
    eps_r : float | fem.Function
    mu_n_over_mu0, mu_p_over_mu0 : float
    tau_n_hat, tau_p_hat : float
    E_t_over_Vt : float
    dt_const, alpha0_const : dolfinx.fem.Constant
        Dimensionless dt_hat = dt / t0 and BDF alpha_0 (updated when
        switching from BDF1 to BDF2 between steps).

    Returns
    -------
    list of ufl.Form
        [F_psi, F_phi_n, F_phi_p]
    """
    import math as _math

    import ufl
    from dolfinx import fem
    from petsc4py import PETSc

    from ..physics.slotboom import n_from_slotboom, p_from_slotboom

    V_psi = spaces.V_psi
    V_phi_n = spaces.V_phi_n
    V_phi_p = spaces.V_phi_p
    msh = V_psi.mesh

    v_psi = ufl.TestFunction(V_psi)
    v_n = ufl.TestFunction(V_phi_n)
    v_p = ufl.TestFunction(V_phi_p)

    L_D2 = fem.Constant(msh, PETSc.ScalarType(sc.lambda2 * sc.L0 ** 2))
    L0_sq = fem.Constant(msh, PETSc.ScalarType(sc.L0 ** 2))
    if isinstance(eps_r, (int, float)):
        eps_r_ufl = fem.Constant(msh, PETSc.ScalarType(float(eps_r)))
    else:
        eps_r_ufl = eps_r

    ni_hat_c = fem.Constant(msh, PETSc.ScalarType(sc.n_i / sc.C0))
    mu_n_c = fem.Constant(msh, PETSc.ScalarType(mu_n_over_mu0))
    mu_p_c = fem.Constant(msh, PETSc.ScalarType(mu_p_over_mu0))
    tau_n_c = fem.Constant(msh, PETSc.ScalarType(tau_n_hat))
    tau_p_c = fem.Constant(msh, PETSc.ScalarType(tau_p_hat))

    # Slotboom carrier densities (always positive)
    n_ufl = n_from_slotboom(psi, phi_n, ni_hat_c)
    p_ufl = p_from_slotboom(psi, phi_p, ni_hat_c)

    # SRH recombination
    n1 = ni_hat_c * _math.exp(E_t_over_Vt)
    p1 = ni_hat_c * _math.exp(-E_t_over_Vt)
    R = (n_ufl * p_ufl - ni_hat_c * ni_hat_c) / (
        tau_p_c * (n_ufl + n1) + tau_n_c * (p_ufl + p1)
    )

    # Poisson block: same as steady-state Slotboom
    rho_hat = p_ufl - n_ufl + N_hat_fn
    F_psi = (
        L_D2 * eps_r_ufl * ufl.inner(ufl.grad(psi), ufl.grad(v_psi)) * ufl.dx
        - rho_hat * v_psi * ufl.dx
    )

    # Electron continuity (Slotboom + chain-rule temporal):
    #   dn/dt = n * (dpsi/dt - dphi_n/dt)
    # so the discrete time term is n * (BDF[psi] - BDF[phi_n]) where
    # BDF[u] = alpha_0/dt * u + f_hist_u carries past values.
    bdf_psi = alpha0_const / dt_const * psi + f_hist_psi
    bdf_phi_n = alpha0_const / dt_const * phi_n + f_hist_phi_n
    bdf_phi_p = alpha0_const / dt_const * phi_p + f_hist_phi_p

    F_phi_n = (
        n_ufl * (bdf_psi - bdf_phi_n) * v_n * ufl.dx
        + L0_sq * mu_n_c * n_ufl
        * ufl.inner(ufl.grad(phi_n), ufl.grad(v_n)) * ufl.dx
        - R * v_n * ufl.dx
    )

    # Hole continuity (Slotboom + chain-rule temporal):
    #   dp/dt = p * (dphi_p/dt - dpsi/dt)
    F_phi_p = (
        p_ufl * (bdf_phi_p - bdf_psi) * v_p * ufl.dx
        + L0_sq * mu_p_c * p_ufl
        * ufl.inner(ufl.grad(phi_p), ufl.grad(v_p)) * ufl.dx
        + R * v_p * ufl.dx
    )

    return [F_psi, F_phi_n, F_phi_p]
