"""
Coupled drift-diffusion runner with adaptive bias continuation.

Solves the three-block (psi, phi_n, phi_p) Slotboom residual under SRH
recombination, ramping the bias from V=0 to the final sweep endpoint
through `semi.continuation.AdaptiveStepController`. Per-step BC
construction is delegated to `semi.bcs`; per-step current evaluation
and IV recording are delegated to `semi.postprocess`.

Dispatch. A config whose `regions` contain both a `semiconductor`-role
and an `insulator`-role region (MOS-style devices: MOSFET, SiGe HBT
with isolation oxide, etc.) routes through the multi-region machinery
from M6: a semiconductor submesh for phi_n/phi_p, a cellwise DG0
`eps_r(x)` on the parent mesh, and the `build_dd_block_residual_mr`
form with `entity_maps`. Single-region configs (pn junction, simple
resistor) continue on the M2-M5 single-region path unchanged.
"""
from __future__ import annotations

from typing import Any

import numpy as np


def run_bias_sweep(
    cfg: dict[str, Any],
    *,
    post_step_hook=None,
    progress_callback=None,
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

    `progress_callback`, if provided, is called with a plain dict per
    successful bias step so the M10 HTTP server can surface SNES
    progress on its WebSocket channel.
    """
    from dolfinx import fem

    from ..bcs import build_dd_dirichlet_bcs, build_psi_dirichlet_bcs, resolve_contacts
    from ..continuation import AdaptiveStepController, StepTooSmall
    from ..doping import build_profile
    from ..mesh import (
        build_eps_r_function,
        build_mesh,
        build_submesh_by_role,
        map_parent_facets_to_submesh,
    )
    from ..physics.drift_diffusion import (
        build_dd_block_residual,
        build_dd_block_residual_mr,
        make_dd_block_spaces,
        make_dd_block_spaces_mr,
    )
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

    msh, cell_tags, facet_tags = build_mesh(cfg)

    N_raw_fn = build_profile(cfg["doping"])

    regions = cfg.get("regions", {})
    roles = {r.get("role", "semiconductor") for r in regions.values()}
    is_multi_region = "insulator" in roles and "semiconductor" in roles

    if is_multi_region and cell_tags is None:
        raise RuntimeError(
            "bias_sweep: multi-region dispatch requires cell_tags, but "
            "the mesh carries none. Check regions_by_box or gmsh physical "
            "groups for cell tags."
        )

    # --- Function spaces and residual forms ---------------------------------
    if is_multi_region:
        semi_tag = _resolve_semi_tag(regions)
        submesh, em_cell, em_vert, _em_geom = build_submesh_by_role(
            msh, cell_tags, regions, role="semiconductor",
        )
        spaces = make_dd_block_spaces_mr(msh, submesh, em_cell)
        V_psi = spaces.V_psi

        eps_r_coeff = build_eps_r_function(msh, cell_tags, regions)

        N_hat_fn = fem.Function(V_psi, name="N_net_hat")
        N_hat_fn.interpolate(lambda x: N_raw_fn(x) / sc.C0)

        entity_maps: list[Any] | None = [em_cell]

        # Precompute the per-parent-tag -> submesh-facet lookup once so the
        # bias loop can build Dirichlet BCs on phi_n/phi_p cheaply.
        sub_facets_by_tag: dict[int, np.ndarray] = {}
        if facet_tags is not None:
            unique_tags = sorted({int(t) for t in np.unique(facet_tags.values)})
            for t in unique_tags:
                parent_facets = facet_tags.find(int(t))
                sub_f = map_parent_facets_to_submesh(
                    msh, submesh, em_vert, parent_facets,
                )
                sub_facets_by_tag[int(t)] = sub_f
    else:
        spaces = make_dd_block_spaces(msh)
        V_psi = spaces.V_psi
        eps_r_coeff = ref_mat.epsilon_r
        N_hat_fn = fem.Function(V_psi, name="N_net_hat")
        N_hat_fn.interpolate(lambda x: N_raw_fn(x) / sc.C0)
        entity_maps = None
        submesh = None
        sub_facets_by_tag = {}

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

    if is_multi_region:
        F_list = build_dd_block_residual_mr(
            spaces, N_hat_fn, sc, eps_r_coeff,
            mu_n_hat, mu_p_hat, tau_n_hat, tau_p_hat,
            cell_tags, semi_tag, E_t_over_Vt,
        )
    else:
        F_list = build_dd_block_residual(
            spaces, N_hat_fn, sc, eps_r_coeff,
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

    bipolar_legs = compute_bipolar_legs(v_sweep_list)

    # Nominal step size from the sweep spacing (e.g., voltage_sweep.step).
    # continuation.max_step defaults to this, which preserves the M2
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

    def build_bcs(contacts):
        if not is_multi_region:
            return build_dd_dirichlet_bcs(
                spaces, msh, facet_tags, contacts, sc, ref_mat, N_raw_fn,
            )
        # Multi-region: psi Dirichlets on the parent (both ohmic body and
        # gate); phi_n / phi_p Dirichlets only on the semiconductor submesh
        # facets corresponding to ohmic contacts. Gate contacts contribute
        # no (phi_n, phi_p) BCs -- their facets are not part of the
        # semiconductor submesh.
        from dolfinx import fem as _fem
        from petsc4py import PETSc as _PETSc

        psi_bcs = build_psi_dirichlet_bcs(
            spaces.V_psi, msh, facet_tags, contacts, sc, ref_mat, N_raw_fn,
        )
        fdim_sub = submesh.topology.dim - 1
        submesh.topology.create_connectivity(
            fdim_sub, submesh.topology.dim,
        )
        carrier_bcs = []
        for c in contacts:
            if c.kind != "ohmic":
                continue
            sub_f = sub_facets_by_tag.get(int(c.facet_tag))
            if sub_f is None or len(sub_f) == 0:
                raise RuntimeError(
                    f"Multi-region bias_sweep: ohmic contact {c.name!r} "
                    f"(parent facet tag {c.facet_tag}) does not lie on "
                    f"the semiconductor submesh."
                )
            V_hat = c.V_applied / sc.V0
            phi_bc = V_hat
            dofs_n = _fem.locate_dofs_topological(
                spaces.V_phi_n, fdim_sub, sub_f,
            )
            dofs_p = _fem.locate_dofs_topological(
                spaces.V_phi_p, fdim_sub, sub_f,
            )
            carrier_bcs.append(_fem.dirichletbc(
                _PETSc.ScalarType(phi_bc), dofs_n, spaces.V_phi_n))
            carrier_bcs.append(_fem.dirichletbc(
                _PETSc.ScalarType(phi_bc), dofs_p, spaces.V_phi_p))
        return psi_bcs + carrier_bcs

    def solve_at(V_by_contact: dict[str, float], tag: str):
        contacts = resolve_contacts(cfg, facet_tags=facet_tags, voltages=V_by_contact)
        bcs = build_bcs(contacts)
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
            entity_maps=entity_maps,
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
              sweep_contact, sweep_facet_info, mu_n_SI, mu_p_SI,
              entity_maps=entity_maps)
    if post_step_hook is not None:
        post_step_hook(V_seed, spaces, iv_rows[-1])
    if progress_callback is not None:
        progress_callback({
            "type": "step_done", "bias_step": 0, "V_applied": float(V_seed),
            "iterations": int(info.get("iterations", 0)),
        })

    def ramp_leg(V_start: float, V_target: float):
        """Adaptive single-direction ramp from V_start to V_target.

        Mutates `spaces` in place, appends one iv row per successful step,
        and returns the SNES info from the last accepted solve. The
        AdaptiveStepController is local to the leg, so a bipolar sweep
        re-seeds growth/halving state at each turning point.
        """
        nonlocal last_info, V_prev
        if abs(V_target - V_start) <= 1.0e-12 or max_step_abs <= 0.0:
            return
        sign = 1.0 if V_target > V_start else -1.0
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
        while abs(V_target - V_prev) > 1.0e-12:
            remaining = V_target - V_prev
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
                          sweep_contact, sweep_facet_info, mu_n_SI, mu_p_SI,
                          entity_maps=entity_maps)
                if post_step_hook is not None:
                    post_step_hook(V_try, spaces, iv_rows[-1])
                if progress_callback is not None:
                    progress_callback({
                        "type": "step_done",
                        "bias_step": len(iv_rows) - 1,
                        "V_applied": float(V_try),
                        "iterations": int(info.get("iterations", 0)),
                    })
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

    if bipolar_legs:
        # Walk: V_seed (=0) -> most-negative -> most-positive.
        for V_target in bipolar_legs:
            ramp_leg(V_prev, float(V_target))
    else:
        ramp_leg(V_prev, float(v_sweep_list[-1]))

    x_dof = V_psi.tabulate_dof_coordinates()
    psi_hat = spaces.psi.x.array
    phi_n_hat = spaces.phi_n.x.array
    phi_p_hat = spaces.phi_p.x.array
    if is_multi_region:
        # phi_n / phi_p live on the submesh, so their arrays have a different
        # length than psi. Recovering (n, p) at arbitrary parent-mesh DOFs
        # isn't meaningful for the MR runner -- consumers that need it
        # should interpolate through the submesh directly. Leave n/p as None
        # and let the result carry psi for plotting.
        n_phys = None
        p_phys = None
    else:
        n_hat = (ref_mat.n_i / sc.C0) * np.exp(psi_hat - phi_n_hat)
        p_hat = (ref_mat.n_i / sc.C0) * np.exp(phi_p_hat - psi_hat)
        n_phys = n_hat * sc.C0
        p_phys = p_hat * sc.C0

    return SimulationResult(
        cfg=cfg, mesh=msh, V=V_psi,
        psi=spaces.psi, phi_n=spaces.phi_n, phi_p=spaces.phi_p,
        psi_phys=psi_hat * sc.V0,
        phi_n_phys=phi_n_hat * sc.V0,
        phi_p_phys=phi_p_hat * sc.V0,
        n_phys=n_phys,
        p_phys=p_phys,
        x_dof=x_dof, N_hat=N_hat_fn, scaling=sc,
        solver_info=last_info, iv=iv_rows, bias_contact=sweep_contact,
    )


def compute_bipolar_legs(v_sweep_list: list[float]) -> list[float]:
    """Endpoints for a sign-spanning bias walk, empty for unipolar sweeps.

    A sweep list that contains both a negative entry and a positive entry
    (ignoring numerical-noise values within 1e-12 V of zero) is walked
    in two single-direction legs: V = 0 -> min(V) -> max(V). Unipolar
    sweeps return `[]`, preserving the original single-endpoint ramp
    path used by the pn-junction and MOS benchmarks.
    """
    has_neg = any(v < -1.0e-12 for v in v_sweep_list)
    has_pos = any(v > 1.0e-12 for v in v_sweep_list)
    if has_neg and has_pos:
        return [float(min(v_sweep_list)), float(max(v_sweep_list))]
    return []


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


def _resolve_semi_tag(regions_cfg: dict) -> int:
    """First region with role='semiconductor' wins. Raises if none."""
    for r in regions_cfg.values():
        if "tag" not in r:
            continue
        if r.get("role", "semiconductor") == "semiconductor":
            return int(r["tag"])
    raise ValueError(
        "bias_sweep multi-region dispatch: no region with role='semiconductor'."
    )
