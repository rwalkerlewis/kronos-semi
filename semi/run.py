"""
Top-level entry point: run a simulation from a validated config dict.

Supports:
    solver.type == "equilibrium"      equilibrium Poisson (Day 1).
    solver.type == "drift_diffusion"  coupled solve at baked biases.
    solver.type == "bias_sweep"       bias ramp: walk a contact's
                                      voltage_sweep and solve coupled.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np


@dataclass
class SimulationResult:
    """Container for simulation outputs."""
    cfg: dict[str, Any]
    mesh: Any = None
    V: Any = None
    psi: Any = None
    phi_n: Any = None
    phi_p: Any = None
    psi_phys: np.ndarray | None = None
    phi_n_phys: np.ndarray | None = None
    phi_p_phys: np.ndarray | None = None
    n_phys: np.ndarray | None = None
    p_phys: np.ndarray | None = None
    x_dof: np.ndarray | None = None
    N_hat: Any = None
    scaling: Any = None
    solver_info: dict[str, Any] = field(default_factory=dict)
    iv: list[dict[str, float]] = field(default_factory=list)
    bias_contact: str | None = None


def run(cfg: dict[str, Any]) -> SimulationResult:
    """Dispatch on solver.type."""
    stype = cfg.get("solver", {}).get("type", "equilibrium")
    if stype == "equilibrium":
        return run_equilibrium(cfg)
    if stype in ("drift_diffusion", "bias_sweep"):
        return run_bias_sweep(cfg)
    raise ValueError(f"Unknown solver.type {stype!r}")


def run_equilibrium(cfg: dict[str, Any]) -> SimulationResult:
    """Day 1 equilibrium Poisson solver (Boltzmann statistics)."""
    from dolfinx import fem

    from .doping import build_profile
    from .mesh import build_mesh
    from .physics.poisson import build_equilibrium_poisson_form
    from .scaling import make_scaling_from_config
    from .solver import solve_nonlinear

    ref_mat = _reference_material(cfg)
    sc = make_scaling_from_config(cfg, ref_mat)

    msh, _cell_tags, facet_tags = build_mesh(cfg)
    V = fem.functionspace(msh, ("Lagrange", 1))

    N_raw_fn = build_profile(cfg["doping"])

    def N_hat_expr(x):
        return N_raw_fn(x) / sc.C0

    N_hat_fn = fem.Function(V, name="N_net_hat")
    N_hat_fn.interpolate(N_hat_expr)

    psi = fem.Function(V, name="psi_hat")
    bcs = _build_ohmic_bcs_psi(cfg, V, msh, facet_tags, sc, ref_mat, N_raw_fn)

    two_ni = 2.0 * ref_mat.n_i

    def psi_init_expr(x):
        return np.arcsinh(N_raw_fn(x) / two_ni)

    psi.interpolate(psi_init_expr)
    for bc in bcs:
        bc.set(psi.x.array)
    psi.x.scatter_forward()

    F = build_equilibrium_poisson_form(V, psi, N_hat_fn, sc, ref_mat.epsilon_r)
    info = solve_nonlinear(F, psi, bcs, prefix=f"{cfg['name']}_")

    x_dof = V.tabulate_dof_coordinates()
    psi_phys = psi.x.array * sc.V0
    n_phys = ref_mat.n_i * np.exp(psi.x.array)
    p_phys = ref_mat.n_i * np.exp(-psi.x.array)

    return SimulationResult(
        cfg=cfg, mesh=msh, V=V, psi=psi,
        psi_phys=psi_phys, n_phys=n_phys, p_phys=p_phys,
        x_dof=x_dof, N_hat=N_hat_fn, scaling=sc, solver_info=info,
    )


def run_bias_sweep(cfg: dict[str, Any]) -> SimulationResult:
    """Coupled drift-diffusion solver with optional bias ramping."""
    from dolfinx import fem

    from .doping import build_profile
    from .mesh import build_mesh
    from .physics.drift_diffusion import build_dd_block_residual, make_dd_block_spaces
    from .scaling import make_scaling_from_config
    from .solver import solve_nonlinear_block

    ref_mat = _reference_material(cfg)
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
        bcs = _build_dd_ohmic_bcs(
            cfg, spaces, msh, facet_tags, sc, ref_mat, N_raw_fn, V_by_contact,
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
        sweep_facet_info = _resolve_contact_facets(
            cfg, msh, facet_tags, sweep_contact,
        )

    iv_rows: list[dict[str, float]] = []
    last_info: dict[str, Any] = {}

    from .continuation import AdaptiveStepController, StepTooSmall

    # Seed solve at the first sweep entry (V=0 by construction above).
    V_seed = float(v_sweep_list[0])
    voltages = dict(static_voltages)
    if sweep_contact is not None:
        voltages[sweep_contact] = V_seed
    snap = snapshot()
    info = solve_at(voltages, _fmt_tag(V_seed))
    if not info["converged"]:
        restore(snap)
        raise RuntimeError(f"SNES failed at seed bias V={V_seed:+.4f} V")
    last_info = info
    V_prev = V_seed
    _record_iv(iv_rows, V_seed, spaces, sc, ref_mat,
               sweep_contact, sweep_facet_info, mu_n_SI, mu_p_SI)

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
                info = solve_at(voltages, _fmt_tag(V_try))
                converged = info["converged"]
            except Exception:
                converged = False
                info = {"converged": False, "iterations": -1, "reason": -99}

            if converged:
                last_info = info
                V_prev = V_try
                _record_iv(iv_rows, V_try, spaces, sc, ref_mat,
                           sweep_contact, sweep_facet_info, mu_n_SI, mu_p_SI)
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


def _fmt_tag(v: float) -> str:
    return f"{v:+.4f}".replace("+", "p").replace("-", "m").replace(".", "d")


def _record_iv(iv_rows, V_applied, spaces, sc, ref_mat,
               sweep_contact, sweep_facet_info, mu_n_SI, mu_p_SI):
    if sweep_contact is None or sweep_facet_info is None:
        iv_rows.append({"V": float(V_applied), "J": 0.0})
        return
    J = _evaluate_current_at_contact(
        spaces, sc, ref_mat, sweep_facet_info, mu_n_SI, mu_p_SI,
    )
    iv_rows.append({"V": float(V_applied), "J": float(J)})


def _evaluate_current_at_contact(spaces, sc, ref_mat, facet_info,
                                 mu_n_SI, mu_p_SI) -> float:
    """
    Evaluate J = J_n + J_p at a contact via a UFL weak form:

        J_n = q mu_n n (grad phi_n . n_outward)
        J_p = q mu_p p (grad phi_p . n_outward)

    Integrated over the contact facet and divided by its measure.
    Positive J = current flowing outward (out of the device through
    the contact).
    """
    import ufl
    from dolfinx import fem
    from mpi4py import MPI
    from petsc4py import PETSc

    from .constants import Q
    from .physics.slotboom import n_from_slotboom, p_from_slotboom

    msh = spaces.V_psi.mesh
    tag = int(facet_info["tag"])
    facets = facet_info["facets"]
    fdim = msh.topology.dim - 1

    if len(facets) == 0:
        return 0.0

    import numpy as _np
    from dolfinx.mesh import meshtags
    values = _np.full(len(facets), tag, dtype=_np.int32)
    indices = _np.asarray(facets, dtype=_np.int32)
    sort = _np.argsort(indices)
    facet_mt = meshtags(msh, fdim, indices[sort], values[sort])
    ds = ufl.Measure("ds", domain=msh, subdomain_data=facet_mt, subdomain_id=tag)
    n_vec = ufl.FacetNormal(msh)

    ni_hat = fem.Constant(msh, PETSc.ScalarType(ref_mat.n_i / sc.C0))
    psi = spaces.psi
    phi_n = spaces.phi_n
    phi_p = spaces.phi_p

    n_ufl = n_from_slotboom(psi, phi_n, ni_hat) * sc.C0
    p_ufl = p_from_slotboom(psi, phi_p, ni_hat) * sc.C0
    grad_phi_n_phys = sc.V0 * ufl.grad(phi_n)
    grad_phi_p_phys = sc.V0 * ufl.grad(phi_p)

    Jn = Q * mu_n_SI * n_ufl * ufl.dot(grad_phi_n_phys, n_vec)
    Jp = Q * mu_p_SI * p_ufl * ufl.dot(grad_phi_p_phys, n_vec)
    current_form = fem.form((Jn + Jp) * ds)
    area_form = fem.form(1.0 * ds)

    I_local = fem.assemble_scalar(current_form)
    A_local = fem.assemble_scalar(area_form)
    I = msh.comm.allreduce(I_local, op=MPI.SUM)
    A = msh.comm.allreduce(A_local, op=MPI.SUM)

    if msh.topology.dim == 1:
        if A == 0.0:
            A = 1.0
        J_total = I / A
    else:
        if A == 0.0:
            return 0.0
        J_total = I / A

    return float(J_total)


def _resolve_contact_facets(cfg, msh, facet_tags, contact_name):
    from dolfinx import fem

    tag_by_name = {p["name"]: int(p["tag"])
                   for p in cfg["mesh"].get("facets_by_plane", [])}
    contact = next(c for c in cfg["contacts"] if c["name"] == contact_name)
    facet_ref = contact["facet"]
    tag = tag_by_name[facet_ref] if isinstance(facet_ref, str) else int(facet_ref)

    fdim = msh.topology.dim - 1
    facets = facet_tags.find(tag)

    outward_sign = 1
    extents = cfg["mesh"].get("extents", [])
    for p in cfg["mesh"].get("facets_by_plane", []):
        if int(p["tag"]) == tag:
            axis = int(p["axis"])
            if axis < len(extents):
                xmin, xmax = extents[axis]
                outward_sign = 1 if float(p["value"]) >= 0.5 * (xmin + xmax) else -1
            break

    V_ref = fem.functionspace(msh, ("Lagrange", 1))
    dofs = fem.locate_dofs_topological(V_ref, fdim, facets)
    return {
        "dofs": np.asarray(dofs, dtype=np.int64),
        "facets": facets,
        "outward_sign": int(outward_sign),
        "tag": int(tag),
    }


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


def _reference_material(cfg: dict[str, Any]):
    from .materials import get_material
    for _region_name, region in cfg["regions"].items():
        mat = get_material(region["material"])
        if mat.is_semiconductor():
            return mat
    raise ValueError("No semiconductor region found; nothing to solve.")


def _build_ohmic_bcs_psi(cfg, V, msh, facet_tags, sc, ref_mat, N_raw_fn):
    """Dirichlet BCs on psi only (equilibrium Poisson path)."""
    from dolfinx import fem
    from petsc4py import PETSc

    fdim = msh.topology.dim - 1
    tag_by_name = {p["name"]: int(p["tag"])
                   for p in cfg["mesh"].get("facets_by_plane", [])}

    bcs = []
    for contact in cfg["contacts"]:
        if contact["type"] != "ohmic":
            continue
        facet_ref = contact["facet"]
        tag = tag_by_name[facet_ref] if isinstance(facet_ref, str) else int(facet_ref)
        if facet_tags is None:
            raise RuntimeError("Mesh has no facet tags.")
        facets = facet_tags.find(tag)
        if len(facets) == 0:
            raise RuntimeError(
                f"No facets with tag {tag} for contact {contact['name']!r}"
            )
        N_net = _evaluate_doping_at_facet(msh, facets, fdim, N_raw_fn)
        psi_eq_hat = float(np.arcsinh(N_net / (2.0 * ref_mat.n_i)))
        V_applied = contact.get("voltage", 0.0)
        psi_bc = psi_eq_hat + V_applied / sc.V0
        dofs = fem.locate_dofs_topological(V, fdim, facets)
        bcs.append(fem.dirichletbc(PETSc.ScalarType(psi_bc), dofs, V))
    return bcs


def _build_dd_ohmic_bcs(cfg, spaces, msh, facet_tags, sc, ref_mat, N_raw_fn,
                        voltages: dict[str, float]):
    """
    Dirichlet BCs for the (psi, phi_n, phi_p) block system.

    psi_hat   = asinh(N_net / (2 n_i)) + V_app / V_t
    phi_n_hat = phi_p_hat = V_app / V_t
    """
    from dolfinx import fem
    from petsc4py import PETSc

    fdim = msh.topology.dim - 1
    tag_by_name = {p["name"]: int(p["tag"])
                   for p in cfg["mesh"].get("facets_by_plane", [])}

    bcs = []
    for contact in cfg["contacts"]:
        if contact["type"] != "ohmic":
            continue
        name = contact["name"]
        facet_ref = contact["facet"]
        tag = tag_by_name[facet_ref] if isinstance(facet_ref, str) else int(facet_ref)
        facets = facet_tags.find(tag)
        if len(facets) == 0:
            raise RuntimeError(f"No facets with tag {tag} for contact {name!r}")
        N_net = _evaluate_doping_at_facet(msh, facets, fdim, N_raw_fn)
        psi_eq_hat = float(np.arcsinh(N_net / (2.0 * ref_mat.n_i)))
        V_app = float(voltages.get(name, contact.get("voltage", 0.0)))
        V_hat = V_app / sc.V0
        psi_bc = psi_eq_hat + V_hat
        phi_bc = V_hat

        dofs_psi = fem.locate_dofs_topological(spaces.V_psi, fdim, facets)
        dofs_n = fem.locate_dofs_topological(spaces.V_phi_n, fdim, facets)
        dofs_p = fem.locate_dofs_topological(spaces.V_phi_p, fdim, facets)

        bcs.append(fem.dirichletbc(PETSc.ScalarType(psi_bc), dofs_psi, spaces.V_psi))
        bcs.append(fem.dirichletbc(PETSc.ScalarType(phi_bc), dofs_n, spaces.V_phi_n))
        bcs.append(fem.dirichletbc(PETSc.ScalarType(phi_bc), dofs_p, spaces.V_phi_p))
    return bcs


def _evaluate_doping_at_facet(msh, facets, fdim, N_raw_fn) -> float:
    tdim = msh.topology.dim
    msh.topology.create_connectivity(fdim, 0)
    f2v = msh.topology.connectivity(fdim, 0)
    coords = msh.geometry.x
    verts = f2v.links(int(facets[0]))
    centroid = coords[verts, :tdim].mean(axis=0)
    pt = centroid.reshape(tdim, 1)
    return float(N_raw_fn(pt)[0])
