"""
Coupled drift-diffusion runner with adaptive bias continuation.

Solves the three-block (psi, phi_n, phi_p) Slotboom residual under SRH
recombination, ramping the bias from V=0 to the final sweep endpoint
through `semi.continuation.AdaptiveStepController`. Per-step BC
construction is delegated to `semi.bcs`; per-step current evaluation
and IV recording are delegated to `semi.postprocess`.
"""
from __future__ import annotations

from typing import Any

import numpy as np


def run_bias_sweep(
    cfg: dict[str, Any],
    *,
    post_step_hook=None,
):
    """
    Coupled drift-diffusion solver with optional bias ramping.

    `post_step_hook`, if provided, is called after every successful
    bias step with the signature

        post_step_hook(V_applied, spaces, iv_row)

    where `iv_row` is the freshly-appended dict in `result.iv`. The
    hook is expected to mutate `iv_row` in place (e.g. to record
    Phase 3 V&V continuity samples on top of the V / J the sweep
    already records). No return value is inspected.
    """
    from dolfinx import fem

    from ..bcs import build_dd_dirichlet_bcs, resolve_contacts
    from ..continuation import AdaptiveStepController, StepTooSmall
    from ..doping import build_profile
    from ..mesh import build_mesh
    from ..physics.drift_diffusion import build_dd_block_residual, make_dd_block_spaces
    from ..postprocess import (
        fmt_tag,
        record_iv,
        resolve_contact_facets,
    )
    from ..run import SimulationResult
    from ..scaling import make_scaling_from_config
    from ..solver import solve_nonlinear_block
    from ._common import reference_material

    ref_mat = reference_material(cfg)
    sc = make_scaling_from_config(cfg, ref_mat)

    msh, _cell_tags, facet_tags = build_mesh(cfg)

    N_raw_fn = build_profile(cfg["doping"])

    spaces = make_dd_block_spaces(msh)
    V_psi = spaces.V_psi

    N_hat_fn = fem.Function(V_psi, name="N_net_hat")
    N_hat_fn.interpolate(lambda x: N_raw_fn(x) / sc.C0)

    two_ni = 2.0 * ref_mat.n_i
    spaces.psi.interpolate(lambda x: np.arcsinh(N_raw_fn(x) / two_ni))
    spaces.phi_n.x.array[:] = 0.0
    spaces.phi_p.x.array[:] = 0.0
    for fn in (spaces.psi, spaces.phi_n, spaces.phi_p):
        fn.x.scatter_forward()

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

    sweep_contact, sweep_values = _resolve_sweep(cfg)

    cont = cfg.get("solver", {}).get("continuation", {})
    max_halvings = int(cont.get("max_halvings", 6))
    min_step = float(cont.get("min_step", 1.0e-4))
    easy_iter_threshold = int(cont.get("easy_iter_threshold", 4))
    grow_factor = float(cont.get("grow_factor", 1.5))

    F_list = build_dd_block_residual(
        spaces, N_hat_fn, sc, ref_mat.epsilon_r,
        mu_n_hat, mu_p_hat, tau_n_hat, tau_p_hat, E_t_over_Vt,
    )

    static_voltages: dict[str, float] = {}
    for c in cfg["contacts"]:
        if c["type"] != "ohmic":
            continue
        if sweep_contact is not None and c["name"] == sweep_contact:
            continue
        static_voltages[c["name"]] = float(c.get("voltage", 0.0))

    if sweep_contact is None:
        v_sweep_list = [0.0]
        baked = float(next(
            (c.get("voltage", 0.0) for c in cfg["contacts"]
             if c["type"] == "ohmic" and c.get("voltage", 0.0) != 0.0),
            0.0,
        ))
        if baked != 0.0:
            v_sweep_list.append(baked)
    else:
        v_sweep_list = list(sweep_values)
        if not v_sweep_list or v_sweep_list[0] != 0.0:
            v_sweep_list = [0.0] + [v for v in v_sweep_list if v != 0.0]

    # Nominal step size from the sweep spacing (e.g., voltage_sweep.step).
    # continuation.max_step defaults to this, which preserves the Day 2
    # halving-only behaviour when growth cannot exceed the initial step.
    if len(v_sweep_list) >= 2:
        nominal_step_abs = abs(v_sweep_list[1] - v_sweep_list[0])
    else:
        nominal_step_abs = 0.0
    max_step_abs = float(cont.get("max_step", nominal_step_abs))

    def snapshot():
        return (
            spaces.psi.x.array.copy(),
            spaces.phi_n.x.array.copy(),
            spaces.phi_p.x.array.copy(),
        )

    def restore(snap):
        spaces.psi.x.array[:] = snap[0]
        spaces.phi_n.x.array[:] = snap[1]
        spaces.phi_p.x.array[:] = snap[2]
        for fn in (spaces.psi, spaces.phi_n, spaces.phi_p):
            fn.x.scatter_forward()

    def solve_at(V_by_contact: dict[str, float], tag: str):
        contacts = resolve_contacts(cfg, facet_tags=facet_tags, voltages=V_by_contact)
        bcs = build_dd_dirichlet_bcs(
            spaces, msh, facet_tags, contacts, sc, ref_mat, N_raw_fn,
        )
        space_to_fn = {
            id(spaces.V_psi): spaces.psi,
            id(spaces.V_phi_n): spaces.phi_n,
            id(spaces.V_phi_p): spaces.phi_p,
        }
        for bc in bcs:
            fn = space_to_fn.get(id(bc.function_space))
            if fn is not None:
                bc.set(fn.x.array)
        for fn in (spaces.psi, spaces.phi_n, spaces.phi_p):
            fn.x.scatter_forward()
        return solve_nonlinear_block(
            F_list, [spaces.psi, spaces.phi_n, spaces.phi_p],
            bcs, prefix=f"{cfg['name']}_dd_{tag}_",
            petsc_options={
                "snes_rtol": 1.0e-14,
                "snes_atol": 1.0e-14,
                "snes_stol": 1.0e-14,
                "snes_max_it": 60,
            },
        )

    sweep_facet_info = None
    if sweep_contact is not None:
        sweep_facet_info = resolve_contact_facets(
            cfg, msh, facet_tags, sweep_contact,
        )

    iv_rows: list[dict[str, float]] = []
    last_info: dict[str, Any] = {}

    # Seed solve at the first sweep entry (V=0 by construction above).
    V_seed = float(v_sweep_list[0])
    voltages = dict(static_voltages)
    if sweep_contact is not None:
        voltages[sweep_contact] = V_seed
    snap = snapshot()
    info = solve_at(voltages, fmt_tag(V_seed))
    if not info["converged"]:
        restore(snap)
        raise RuntimeError(f"SNES failed at seed bias V={V_seed:+.4f} V")
    last_info = info
    V_prev = V_seed
    record_iv(iv_rows, V_seed, spaces, sc, ref_mat,
              sweep_contact, sweep_facet_info, mu_n_SI, mu_p_SI)
    if post_step_hook is not None:
        post_step_hook(V_seed, spaces, iv_rows[-1])

    V_end = float(v_sweep_list[-1])
    if abs(V_end - V_seed) > 0.0 and max_step_abs > 0.0:
        sign = 1.0 if V_end > V_seed else -1.0
        # Start at the nominal sweep spacing (voltage_sweep.step); the
        # adaptive controller grows toward max_step_abs on consecutive
        # easy solves and halves on failure.
        initial_abs = nominal_step_abs if nominal_step_abs > 0.0 else max_step_abs
        initial_abs = min(initial_abs, max_step_abs)
        initial_step = sign * initial_abs

        controller = AdaptiveStepController(
            initial_step=initial_step,
            max_step_abs=max_step_abs,
            min_step_abs=min_step,
            easy_iter_threshold=easy_iter_threshold,
            grow_factor=grow_factor,
        )

        halvings = 0
        while abs(V_end - V_prev) > 1.0e-12:
            remaining = V_end - V_prev
            step = controller.clamp_to_endpoint(remaining)
            V_try = V_prev + step
            voltages = dict(static_voltages)
            if sweep_contact is not None:
                voltages[sweep_contact] = V_try
            snap = snapshot()
            try:
                info = solve_at(voltages, fmt_tag(V_try))
                converged = info["converged"]
            except Exception:
                converged = False
                info = {"converged": False, "iterations": -1, "reason": -99}

            if converged:
                last_info = info
                V_prev = V_try
                record_iv(iv_rows, V_try, spaces, sc, ref_mat,
                          sweep_contact, sweep_facet_info, mu_n_SI, mu_p_SI)
                if post_step_hook is not None:
                    post_step_hook(V_try, spaces, iv_rows[-1])
                halvings = 0
                controller.on_success(int(info.get("iterations", 0)))
                continue

            restore(snap)
            halvings += 1
            try:
                controller.on_failure()
            except StepTooSmall as exc:
                raise RuntimeError(
                    f"Bias ramp failed near V={V_try:+.4f} V: {exc}"
                ) from exc
            if halvings > max_halvings:
                raise RuntimeError(
                    f"Bias ramp failed near V={V_try:+.4f} V after "
                    f"{halvings} halvings (min_step={min_step})."
                )

    x_dof = V_psi.tabulate_dof_coordinates()
    psi_hat = spaces.psi.x.array
    phi_n_hat = spaces.phi_n.x.array
    phi_p_hat = spaces.phi_p.x.array
    n_hat = (ref_mat.n_i / sc.C0) * np.exp(psi_hat - phi_n_hat)
    p_hat = (ref_mat.n_i / sc.C0) * np.exp(phi_p_hat - psi_hat)

    return SimulationResult(
        cfg=cfg, mesh=msh, V=V_psi,
        psi=spaces.psi, phi_n=spaces.phi_n, phi_p=spaces.phi_p,
        psi_phys=psi_hat * sc.V0,
        phi_n_phys=phi_n_hat * sc.V0,
        phi_p_phys=phi_p_hat * sc.V0,
        n_phys=n_hat * sc.C0,
        p_phys=p_hat * sc.C0,
        x_dof=x_dof, N_hat=N_hat_fn, scaling=sc,
        solver_info=last_info, iv=iv_rows, bias_contact=sweep_contact,
    )


def _resolve_sweep(cfg):
    """Return (contact_name, list_of_voltage_values) or (None, []) if none."""
    for c in cfg["contacts"]:
        if c["type"] != "ohmic":
            continue
        sweep = c.get("voltage_sweep")
        if sweep is None:
            continue
        start = float(sweep["start"])
        stop = float(sweep["stop"])
        step = float(sweep["step"])
        if step <= 0.0:
            raise ValueError("voltage_sweep.step must be positive")
        direction = 1.0 if stop >= start else -1.0
        n = int(round((stop - start) / (direction * step))) + 1
        values = [start + i * direction * step for i in range(n)]
        return c["name"], values
    return None, []
