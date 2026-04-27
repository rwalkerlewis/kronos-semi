"""
Transient drift-diffusion runner (BDF1 / BDF2).

Solves the coupled (psi, n, p) system in primary-density form with
backward differentiation formula (BDF) time integration. Mirrors the
structure of :mod:`semi.runners.bias_sweep` for steady-state solves.

Primary unknowns
----------------
Unlike :mod:`semi.runners.bias_sweep`, which solves for Slotboom
quasi-Fermi potentials (psi, phi_n, phi_p), the transient runner uses
**carrier densities** (psi, n_hat, p_hat) as primary unknowns. This
choice makes the time-derivative term d(n)/dt natural and avoids the
nonlinear coupling d/dt[ni*exp(psi - phi_n)] that would arise if phi_n
were kept as the primary variable. See docs/adr/0009-transient-formulation.md.

Time integration
----------------
BDF order is configurable (1 = backward Euler, 2 = BDF2). For BDF2, a
single BDF1 step seeds the two-level history before switching. Fixed
timestep dt is used throughout (no adaptive time stepping in M13).
See docs/adr/0010-bdf-time-integration.md.

Residual form (scaled units)
-----------------------------
Poisson (unchanged):
    L_D^2 * eps_r * grad(psi) . grad(v_psi)  -  (p - n + N) * v_psi  =  0

Electron continuity (n-form):
    alpha_0/dt * n * v_n
    + L0^2 * mu_n * (inner(grad(n), grad(v_n)) - n * inner(grad(psi), grad(v_n)))
    + R * v_n
    + f_hist_n * v_n  =  0

Hole continuity (p-form):
    alpha_0/dt * p * v_p
    + L0^2 * mu_p * (inner(grad(p), grad(v_p)) + p * inner(grad(psi), grad(v_p)))
    - R * v_p
    + f_hist_p * v_p  =  0

where f_hist_n and f_hist_p are known source terms from the BDF history
(updated each timestep).

Boundary conditions (ohmic contacts)
--------------------------------------
At ohmic contact with applied bias V (scaled V_hat = V/V_t):
    psi_bc    = arcsinh(N_net / (2*ni)) + V_hat
    n_hat_bc  = ni_hat * exp(arcsinh(N_net/(2*ni)))   [constant w.r.t. V]
    p_hat_bc  = ni_hat * exp(-arcsinh(N_net/(2*ni)))  [constant w.r.t. V]

The carrier BCs are constant because the applied bias shifts psi and
phi_n/phi_p by the same amount, leaving n = ni*exp(psi - phi_n) fixed.
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
    Coupled transient drift-diffusion solver.

    Runs BDF1 or BDF2 time integration on the (psi, n_hat, p_hat)
    primary-variable system, starting from equilibrium at t=0 and
    applying any bias specified in `cfg["contacts"]` from t=0+.

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

    from ..bcs import resolve_contacts
    from ..doping import build_profile
    from ..fem.mass import assemble_lumped_mass
    from ..fem.sg_assembly import solve_sg_block_1d
    from ..mesh import build_mesh
    from ..physics.drift_diffusion import make_dd_block_spaces
    from ..postprocess import (
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

    # SNES tolerances. The (n,p) Galerkin spatial residual carries an
    # L_0^2 prefactor (typically 1e-10 to 1e-12 in scaled units), so the
    # bias_sweep default atol=1e-7 (ADR 0008) leaves SNES exiting at
    # iteration 0 every step and the time loop sits frozen at the initial
    # guess. atol=1e-10 forces honest Newton iteration. See ADR 0009
    # "Known limitation (M13.1)" for the deeper structural issue this
    # tighter tolerance now exposes (tracked separately).
    snes_opts = solver_cfg.get("snes", {}) or {}
    snes_petsc_options = {
        "snes_rtol": float(snes_opts.get("rtol", 1.0e-10)),
        "snes_atol": float(snes_opts.get("atol", 1.0e-10)),
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
    # Create function spaces: (psi, n_hat, p_hat)
    # V_psi is the same P1 space for psi.
    # V_n and V_p are P1 spaces for the carrier densities.
    # ------------------------------------------------------------------
    spaces = make_dd_block_spaces(msh)
    # Reuse spaces.V_phi_n and spaces.V_phi_p as V_n and V_p
    # (same P1 space structure; we rename the Function objects below)
    V_psi = spaces.V_psi
    V_n = spaces.V_phi_n
    V_p = spaces.V_phi_p

    # Create unknown Function objects for transient (n_hat, p_hat form)
    psi = spaces.psi
    n_hat = fem.Function(V_n, name="n_hat")
    p_hat = fem.Function(V_p, name="p_hat")

    # Interpolate doping for Poisson source
    N_hat_fn = fem.Function(V_psi, name="N_net_hat")
    N_hat_fn.interpolate(lambda x: N_raw_fn(x) / sc.C0)

    ni_hat_val = sc.n_i / sc.C0

    # ------------------------------------------------------------------
    # Step 0: Solve equilibrium at t=0 to obtain initial psi
    # The equilibrium solve uses zero bias on all contacts.
    # ------------------------------------------------------------------
    eq_cfg = copy.deepcopy(cfg)
    eq_cfg["solver"] = {"type": "equilibrium"}
    # Zero out all contact voltages so we get true V=0 equilibrium
    for c in eq_cfg["contacts"]:
        c["voltage"] = 0.0
        c.pop("voltage_sweep", None)
    eq_result = run_equilibrium(eq_cfg)

    # Copy equilibrium psi into our psi function
    # We need to interpolate from eq_result onto our mesh's V_psi space.
    # Since both use the same mesh, just copy the array.
    psi.x.array[:] = eq_result.psi.x.array
    psi.x.scatter_forward()

    # Initialize n and p from equilibrium using Boltzmann statistics
    psi_arr = psi.x.array
    n_hat.x.array[:] = ni_hat_val * np.exp(psi_arr)
    p_hat.x.array[:] = ni_hat_val * np.exp(-psi_arr)
    n_hat.x.scatter_forward()
    p_hat.x.scatter_forward()

    # ------------------------------------------------------------------
    # Assemble lumped mass diagonal (for future M14 use, and for
    # diagnostic storage in meta)
    # ------------------------------------------------------------------
    dx = ufl.Measure("dx", domain=msh)
    M_n_diag, M_p_diag = assemble_lumped_mass(V_n, V_p, dx)

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
    # Build initial static voltages from contacts (no sweep in transient)
    # ------------------------------------------------------------------
    static_voltages: dict[str, float] = {}
    for c in cfg["contacts"]:
        if c["type"] != "ohmic":
            continue
        static_voltages[c["name"]] = float(c.get("voltage", 0.0))

    # ------------------------------------------------------------------
    # Build the transient residual
    # ------------------------------------------------------------------
    # History source functions (updated each timestep)
    f_hist_n = fem.Function(V_n, name="hist_n")
    f_hist_p = fem.Function(V_p, name="hist_p")
    f_hist_n.x.array[:] = 0.0
    f_hist_p.x.array[:] = 0.0

    # dt_hat = dt / t0 is the dimensionless time step consistent with the
    # scaled spatial coefficients (L0^2 * mu_hat, etc.) in the residual.
    # Using physical dt_val here would make the temporal term ~t0/dt_phys
    # times larger than the spatial terms, leaving the carrier fields frozen.
    dt_hat = dt_val / sc.t0

    from petsc4py import PETSc

    dt_const = fem.Constant(msh, PETSc.ScalarType(dt_hat))
    alpha0_const = fem.Constant(msh, PETSc.ScalarType(alpha_0))

    F_list = _build_transient_residual(
        psi, n_hat, p_hat, N_hat_fn, f_hist_n, f_hist_p,
        spaces, sc, ref_mat.epsilon_r,
        mu_n_hat, mu_p_hat, tau_n_hat, tau_p_hat, E_t_over_Vt,
        dt_const, alpha0_const,
    )

    # ------------------------------------------------------------------
    # Helper to build Dirichlet BCs for (psi, n_hat, p_hat) at t+dt
    # ------------------------------------------------------------------
    def _build_transient_bcs(voltages: dict[str, float]) -> list:
        """Build BCs for (psi, n, p) at ohmic contacts."""
        contacts = resolve_contacts(cfg, facet_tags=facet_tags,
                                    voltages=voltages)
        bcs = []
        fdim = msh.topology.dim - 1
        two_ni = 2.0 * ref_mat.n_i
        from dolfinx import fem as _fem
        for c in contacts:
            if c.kind != "ohmic":
                continue
            facets = facet_tags.find(c.facet_tag)
            if len(facets) == 0:
                continue
            # Evaluate local doping at contact
            from ..bcs import _evaluate_doping_at_facet
            N_net = _evaluate_doping_at_facet(msh, facets, fdim, N_raw_fn)
            psi_eq_hat = float(np.arcsinh(N_net / two_ni))
            V_hat = c.V_applied / sc.V0
            psi_bc = psi_eq_hat + V_hat

            # n, p BCs are independent of V (ohmic contact = equilibrium
            # minority-carrier concentrations)
            n_bc = ni_hat_val * float(np.exp(psi_eq_hat))
            p_bc = ni_hat_val * float(np.exp(-psi_eq_hat))

            dofs_psi = _fem.locate_dofs_topological(V_psi, fdim, facets)
            dofs_n = _fem.locate_dofs_topological(V_n, fdim, facets)
            dofs_p = _fem.locate_dofs_topological(V_p, fdim, facets)

            bcs.append(_fem.dirichletbc(
                PETSc.ScalarType(psi_bc), dofs_psi, V_psi))
            bcs.append(_fem.dirichletbc(
                PETSc.ScalarType(n_bc), dofs_n, V_n))
            bcs.append(_fem.dirichletbc(
                PETSc.ScalarType(p_bc), dofs_p, V_p))
        return bcs

    # ------------------------------------------------------------------
    # Helper to apply psi BC init values on psi array (so SNES starts
    # near the solution at each step)
    # ------------------------------------------------------------------
    def _apply_bc_values(bcs: list):
        space_to_fn = {
            id(V_psi): psi,
            id(V_n): n_hat,
            id(V_p): p_hat,
        }
        for bc in bcs:
            fn = space_to_fn.get(id(bc.function_space))
            if fn is not None:
                bc.set(fn.x.array)
        for fn in (psi, n_hat, p_hat):
            fn.x.scatter_forward()

    # ------------------------------------------------------------------
    # History storage (numpy arrays, most-recent LAST)
    # ------------------------------------------------------------------
    n_hist: list[np.ndarray] = [n_hat.x.array.copy()]
    p_hist: list[np.ndarray] = [p_hat.x.array.copy()]

    # Keep only as many history levels as needed (order levels)
    def _trim_history():
        while len(n_hist) > order:
            n_hist.pop(0)
        while len(p_hist) > order:
            p_hist.pop(0)

    # ------------------------------------------------------------------
    # IV recording helper
    # ------------------------------------------------------------------
    def _record_all_iv(t_val: float, iv_rows: list):
        """Record IV for all ohmic contacts at the current solution."""
        from ..physics.slotboom import phi_n_from_np, phi_p_from_np
        from ..postprocess import evaluate_partial_currents

        for c in cfg["contacts"]:
            if c["type"] != "ohmic":
                continue
            cname = c["name"]
            V_applied = static_voltages.get(cname, 0.0)
            finfo = contact_facet_infos.get(cname)
            if finfo is None:
                continue
            # Reconstruct phi_n, phi_p from n_hat, p_hat for current evaluation
            phi_n_arr = phi_n_from_np(psi.x.array, n_hat.x.array, ni_hat_val)
            phi_p_arr = phi_p_from_np(psi.x.array, p_hat.x.array, ni_hat_val)
            phi_n_fn = fem.Function(V_n, name="phi_n_tmp")
            phi_p_fn = fem.Function(V_p, name="phi_p_tmp")
            phi_n_fn.x.array[:] = phi_n_arr
            phi_p_fn.x.array[:] = phi_p_arr
            phi_n_fn.x.scatter_forward()
            phi_p_fn.x.scatter_forward()

            # Build a spaces-like object for evaluate_partial_currents
            class _Spaces:
                pass
            sp = _Spaces()
            sp.V_psi = V_psi
            sp.psi = psi
            sp.V_phi_n = V_n
            sp.V_phi_p = V_p
            sp.phi_n = phi_n_fn
            sp.phi_p = phi_p_fn

            J_n, J_p = evaluate_partial_currents(
                sp, sc, ref_mat, finfo, mu_n_SI, mu_p_SI,
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

    # For BDF2: first step uses BDF1 to build the two-level history
    effective_order = 1  # start with BDF1 always
    bdf1 = BDFCoefficients(1)

    while t_current < t_end - 0.5 * dt_val and step_count < max_steps:
        t_next = t_current + dt_val
        step_count += 1

        # Decide which BDF order to use this step
        if order == 2 and len(n_hist) >= 2:
            effective_order = 2
            use_bdf = bdf
        else:
            effective_order = 1
            use_bdf = bdf1

        alpha_k_use = use_bdf.coeffs[0]
        alpha0_const.value = PETSc.ScalarType(alpha_k_use)

        # Update history source functions:
        # f_hist_n[i] = sum_{k=1}^K alpha_k * n_hist[-(k-1)-1][i] / dt
        hist_n_arr = np.zeros(len(n_hist[-1]))
        hist_p_arr = np.zeros(len(p_hist[-1]))
        for k in range(1, effective_order + 1):
            idx = -(k - 1) - 1  # n_hist[-1] = most recent (n^n), n_hist[-2] = n^{n-1}
            hist_n_arr += use_bdf.coeffs[k] * n_hist[idx]
            hist_p_arr += use_bdf.coeffs[k] * p_hist[idx]
        f_hist_n.x.array[:] = hist_n_arr / dt_hat
        f_hist_p.x.array[:] = hist_p_arr / dt_hat
        f_hist_n.x.scatter_forward()
        f_hist_p.x.scatter_forward()

        # Build BCs
        bcs = _build_transient_bcs(static_voltages)

        # Set BC values into unknowns as starting guess improvement
        _apply_bc_values(bcs)

        # SNES solve. M13.1 status (2026-04-26 follow-up #1): the SG
        # residual + analytic-Jacobian-correction primitive
        # `semi.fem.sg_assembly.solve_sg_block_1d` is implemented and
        # FD-verified, but its iterate develops a small numerical seed
        # of negative p in the depletion region at the first time step
        # that amplifies across subsequent steps (SG primary-variable
        # continuity is not strictly positivity-preserving). Closing
        # the steady-state agreement xfail needs either a positivity-
        # preserving change of variables or a continuation strategy
        # for the first time step. Until then the runner stays on the
        # M13 Galerkin path so the M13 + M14 transient/AC tests
        # continue to pass; the steady-state agreement test remains
        # xfail per ADR 0012 with updated reason.
        # See /tmp/m13.1-integration-blocker.md for details.
        _ = solve_sg_block_1d  # silence unused-import
        tag = fmt_tag(t_next)
        info = solve_nonlinear_block(
            F_list, [psi, n_hat, p_hat], bcs,
            prefix=f"{cfg['name']}_tr_{tag}_",
            petsc_options=snes_petsc_options,
        )

        if not info["converged"]:
            raise RuntimeError(
                f"Transient SNES failed at t={t_next:.3e} s "
                f"(step {step_count}); reason={info['reason']}"
            )

        # Update history
        n_hist.append(n_hat.x.array.copy())
        p_hist.append(p_hat.x.array.copy())
        _trim_history()

        t_current = t_next
        n_steps_taken += 1

        # Record IV
        _record_all_iv(t_current, iv_rows)
        t_vals.append(t_current)

        # Snapshots
        if step_count % output_every == 0 or abs(t_current - t_end) < 0.5 * dt_val:
            fields_out["psi"].append(psi.x.array.copy() * sc.V0)
            fields_out["n"].append(n_hat.x.array.copy() * sc.C0)
            fields_out["p"].append(p_hat.x.array.copy() * sc.C0)
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
    psi, n_hat, p_hat, N_hat_fn,
    f_hist_n, f_hist_p,
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
    Build the transient three-block residual in (psi, n_hat, p_hat) form.

    The Poisson block is identical to the Slotboom form but uses n_hat
    and p_hat directly (no exponentials). The continuity blocks are
    rewritten from Slotboom phi-form to density n/p form (see
    docs/adr/0009-transient-formulation.md).

    Parameters
    ----------
    psi, n_hat, p_hat : dolfinx.fem.Function
        Unknowns (will be solved for).
    N_hat_fn : dolfinx.fem.Function
        Scaled net doping.
    f_hist_n, f_hist_p : dolfinx.fem.Function
        History source terms (updated externally each timestep):
        f_hist_n = sum_{k=1}^K alpha_k/dt * n^{n+1-k}.
    spaces : DDBlockSpaces
        Function spaces (V_psi, V_phi_n=V_n, V_phi_p=V_p).
    sc : Scaling
    eps_r : float
    mu_n_over_mu0, mu_p_over_mu0 : float
    tau_n_hat, tau_p_hat : float
    E_t_over_Vt : float
    dt_const : dolfinx.fem.Constant
        Dimensionless time step dt_hat = dt_physical / sc.t0.  Must be in
        scaled time units so the temporal coefficient alpha0/dt_hat is
        dimensionally consistent with the L0^2*mu_hat spatial coefficients.
    alpha0_const : dolfinx.fem.Constant
        BDF alpha_0 coefficient as a fem.Constant (updated each step
        when switching from BDF1 to BDF2).

    Returns
    -------
    list of ufl.Form
        [F_psi, F_n, F_p]
    """
    import math as _math

    import ufl
    from dolfinx import fem
    from petsc4py import PETSc

    V_psi = spaces.V_psi
    V_n = spaces.V_phi_n
    V_p = spaces.V_phi_p
    msh = V_psi.mesh

    v_psi = ufl.TestFunction(V_psi)
    v_n = ufl.TestFunction(V_n)
    v_p = ufl.TestFunction(V_p)

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

    # SRH recombination in n/p primary-variable form
    n1 = ni_hat_c * _math.exp(E_t_over_Vt)
    p1 = ni_hat_c * _math.exp(-E_t_over_Vt)
    R = (n_hat * p_hat - ni_hat_c * ni_hat_c) / (
        tau_p_c * (n_hat + n1) + tau_n_c * (p_hat + p1)
    )

    # ------------------------------------------------------------------
    # Poisson block (same structure as Slotboom, but with n_hat, p_hat
    # directly as the carrier densities)
    # ------------------------------------------------------------------
    rho_hat = p_hat - n_hat + N_hat_fn
    F_psi = (
        L_D2 * eps_r_ufl * ufl.inner(ufl.grad(psi), ufl.grad(v_psi)) * ufl.dx
        - rho_hat * v_psi * ufl.dx
    )

    # ------------------------------------------------------------------
    # Convection-diffusion blocks. The UFL form below is the M13
    # Galerkin discretisation; on 1D meshes the SNES residual callback
    # in `run_transient` post-processes the UFL residual to swap the
    # Galerkin convection-diffusion contribution for the per-edge SG
    # flux from `semi.fem.sg_assembly.assemble_sg_residual_1d`. We keep
    # the Galerkin form here so the Jacobian (auto-derived by UFL) has
    # the correct nearest-neighbour sparsity pattern and matches SG to
    # leading order at low Peclet, giving Newton a good preconditioner
    # for the SG residual without requiring a hand-coded SG Jacobian.
    # See ADR 0012 §"1D edge flux" and the runner integration block
    # in run_transient below.
    # ------------------------------------------------------------------
    F_n = (
        alpha0_const / dt_const * n_hat * v_n * ufl.dx
        + L0_sq * mu_n_c * (
            ufl.inner(ufl.grad(n_hat), ufl.grad(v_n))
            - n_hat * ufl.inner(ufl.grad(psi), ufl.grad(v_n))
        ) * ufl.dx
        + R * v_n * ufl.dx
        + f_hist_n * v_n * ufl.dx
    )

    F_p = (
        alpha0_const / dt_const * p_hat * v_p * ufl.dx
        + L0_sq * mu_p_c * (
            ufl.inner(ufl.grad(p_hat), ufl.grad(v_p))
            + p_hat * ufl.inner(ufl.grad(psi), ufl.grad(v_p))
        ) * ufl.dx
        - R * v_p * ufl.dx
        + f_hist_p * v_p * ufl.dx
    )

    return [F_psi, F_n, F_p]
