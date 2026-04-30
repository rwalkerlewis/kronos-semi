"""
MOSCAP LF + HF C-V runner (M14.2).

Sweeps a gate contact over the requested bias range, solves
multi-region equilibrium Poisson at each step, then extracts both the
low-frequency (quasi-static) and high-frequency differential
capacitances via two linearised Poisson solves per bias point.

Supports planar 2D meshes and axisymmetric (r, z) meshes (set
``mesh.axisymmetric: true`` in the config). The axisymmetric path
reproduces a circular MOSCAP rotated about ``r = 0``; the radial
weight is folded into both the residual and the charge-extraction
forms via :mod:`semi.physics.cv`.

Each row of ``result.iv`` is

    {V, Q_gate, J: 0, C_LF, C_HF}

with capacitances in F/m^2 (axisymmetric: per unit gate area;
planar 2D: per unit transverse length, F/m).
"""
from __future__ import annotations

from typing import Any

import numpy as np


def run_moscap_lf_hf(cfg: dict[str, Any], *, progress_callback=None):
    """
    Run a MOSCAP gate-bias sweep and return both LF and HF C-V.

    Returns a SimulationResult whose ``iv`` rows carry the differential
    capacitances. ``solver_info`` is augmented with ``{"C_LF": [...],
    "C_HF": [...], "V_g": [...], "axisymmetric": bool}`` for downstream
    plotting and CI tolerance checks.
    """
    import math

    import ufl
    from dolfinx import fem
    from dolfinx.fem.petsc import LinearProblem
    from mpi4py import MPI
    from petsc4py import PETSc

    from ..bcs import build_psi_dirichlet_bcs, resolve_contacts
    from ..constants import Q
    from ..doping import build_profile
    from ..mesh import build_eps_r_function, build_mesh
    from ..physics.cv import (
        build_charge_form,
        build_hf_charge_sensitivity_form,
        build_hf_jacobian_form,
        build_lf_charge_sensitivity_form,
        infer_majority_carrier,
    )
    from ..physics.poisson import build_equilibrium_poisson_form_mr
    from ..run import SimulationResult
    from ..scaling import make_scaling_from_config
    from ..solver import solve_nonlinear
    from ._common import reference_material

    ref_mat = reference_material(cfg)
    sc = make_scaling_from_config(cfg, ref_mat)

    msh, cell_tags, facet_tags = build_mesh(cfg)
    if cell_tags is None:
        raise RuntimeError(
            "moscap_lf_hf requires regions_by_box: the multi-region "
            "assembly cannot run without cell tags."
        )

    mesh_cfg = cfg["mesh"]
    axisymmetric = bool(mesh_cfg.get("axisymmetric", False))
    axisymmetric_axis: int | None = (
        int(mesh_cfg.get("axisymmetric_axis", 0)) if axisymmetric else None
    )

    if axisymmetric:
        if int(cfg["dimension"]) != 2:
            raise ValueError(
                "mesh.axisymmetric is only supported for dimension == 2."
            )
        # Axis r=0 must be in the mesh: enforce that the radial extent starts at 0.
        rmin = float(mesh_cfg["extents"][axisymmetric_axis][0])
        if not math.isclose(rmin, 0.0, abs_tol=1.0e-12):
            raise ValueError(
                "Axisymmetric meshes must include the symmetry axis r=0 in their "
                f"extents along axis {axisymmetric_axis}; got rmin={rmin}."
            )

    semi_tag = _resolve_semi_tag(cfg["regions"])
    eps_r_fn = build_eps_r_function(msh, cell_tags, cfg["regions"])

    V_psi = fem.functionspace(msh, ("Lagrange", 1))
    N_raw_fn = build_profile(cfg["doping"])

    N_hat_fn = fem.Function(V_psi, name="N_net_hat")
    N_hat_fn.interpolate(lambda x: N_raw_fn(x) / sc.C0)

    psi = fem.Function(V_psi, name="psi_hat")
    two_ni = 2.0 * ref_mat.n_i
    psi.interpolate(lambda x: np.arcsinh(N_raw_fn(x) / two_ni))

    F = build_equilibrium_poisson_form_mr(
        V_psi, psi, N_hat_fn, sc, eps_r_fn, cell_tags, semi_tag,
        axisymmetric_axis=axisymmetric_axis,
    )

    sweep_contact, sweep_values = _resolve_gate_sweep(cfg)
    if sweep_contact is None:
        raise ValueError(
            "moscap_lf_hf requires a contact of type 'gate' with a voltage_sweep."
        )

    static_voltages: dict[str, float] = {}
    for c in cfg["contacts"]:
        if c["name"] == sweep_contact:
            continue
        if c["type"] not in ("ohmic", "gate"):
            continue
        static_voltages[c["name"]] = float(c.get("voltage", 0.0))

    cv_cfg = cfg.get("cv_analysis", {})
    modes = cv_cfg.get("modes", ["LF", "HF"])
    majority = cv_cfg.get("majority_carrier", "auto")
    if majority == "auto":
        majority = infer_majority_carrier(cfg)
    if majority not in ("holes", "electrons"):
        raise ValueError(
            f"cv_analysis.majority_carrier must be holes/electrons/auto, got {majority!r}"
        )

    # Charge extraction forms
    charge_form = build_charge_form(
        V_psi, psi, N_hat_fn, sc, cell_tags, semi_tag,
        axisymmetric_axis=axisymmetric_axis,
    )

    # Sensitivity solves: LF uses the existing DC Jacobian; HF uses a
    # custom Jacobian where the minority-carrier exponential is dropped.
    delta_psi_lf = fem.Function(V_psi, name="delta_psi_LF")
    delta_psi_hf = fem.Function(V_psi, name="delta_psi_HF")

    trial_psi = ufl.TrialFunction(V_psi)
    a_lf = ufl.derivative(F, psi, trial_psi)
    a_hf = build_hf_jacobian_form(
        V_psi, psi, N_hat_fn, sc, eps_r_fn, cell_tags, semi_tag,
        majority,
        axisymmetric_axis=axisymmetric_axis,
    )
    L_zero = fem.Constant(msh, PETSc.ScalarType(0.0)) * ufl.TestFunction(V_psi) * ufl.dx

    sens_form_lf = build_lf_charge_sensitivity_form(
        V_psi, psi, sc, cell_tags, semi_tag, delta_psi_lf,
        axisymmetric_axis=axisymmetric_axis,
    )
    sens_form_hf = build_hf_charge_sensitivity_form(
        V_psi, psi, sc, cell_tags, semi_tag, delta_psi_hf, majority,
        axisymmetric_axis=axisymmetric_axis,
    )

    extents = mesh_cfg["extents"]
    if axisymmetric:
        # The volume forms use weight r/L0 (dimensionless); the
        # actual physical r-integral is L0 * (computed integral).
        # Q_semi = 2 pi q C0 * L0 * computed_int   (full 3D)
        # Q_gate per gate area = -Q_semi / (pi R_gate^2)
        #                      = -(2 q C0 L0 / R_gate^2) * computed_int.
        R_gate = float(cv_cfg.get("gate_radius", extents[axisymmetric_axis][1]))
        area_factor = -(2.0 * Q * sc.C0 * sc.L0) / (R_gate ** 2)
        normalisation_label = "F/m^2"
    else:
        # Planar 2D: integral is dimensionless * m^2; W_lat is the
        # transverse width along the non-vertical axis (axis 0 here).
        W_lat = float(extents[0][1] - extents[0][0])
        area_factor = -(Q * sc.C0) / W_lat
        normalisation_label = "F/m" if int(cfg["dimension"]) == 2 else "F/m^2"

    # Sensitivity Dirichlet BCs: delta_psi = 1/V_t at swept gate, 0 elsewhere.
    sensitivity_bcs = _build_sensitivity_bcs(
        cfg, V_psi, msh, facet_tags, sc, sweep_contact,
    )

    lf_problem = LinearProblem(
        a_lf, L_zero, u=delta_psi_lf, bcs=sensitivity_bcs,
        petsc_options_prefix=f"{cfg['name']}_moscap_lf_",
        petsc_options={
            "ksp_type": "preonly",
            "pc_type": "lu",
            "pc_factor_mat_solver_type": "mumps",
        },
    )
    hf_problem = LinearProblem(
        a_hf, L_zero, u=delta_psi_hf, bcs=sensitivity_bcs,
        petsc_options_prefix=f"{cfg['name']}_moscap_hf_",
        petsc_options={
            "ksp_type": "preonly",
            "pc_type": "lu",
            "pc_factor_mat_solver_type": "mumps",
        },
    )

    iv_rows: list[dict[str, float]] = []
    last_info: dict[str, Any] = {}
    psi_backup = psi.x.array.copy()

    def _solve_at(V_gate: float) -> dict[str, Any]:
        voltages = dict(static_voltages)
        voltages[sweep_contact] = float(V_gate)
        contacts = resolve_contacts(cfg, facet_tags=facet_tags, voltages=voltages)
        bcs = build_psi_dirichlet_bcs(
            V_psi, msh, facet_tags, contacts, sc, ref_mat, N_raw_fn,
        )
        for bc in bcs:
            bc.set(psi.x.array)
        psi.x.scatter_forward()
        tag = _fmt_gate_tag(V_gate)
        return solve_nonlinear(
            F, psi, bcs, prefix=f"{cfg['name']}_moscap_lf_hf_{tag}_",
            petsc_options={
                "snes_rtol": 1.0e-12,
                "snes_atol": 1.0e-14,
                "snes_max_it": 50,
            },
        )

    for V_gate in sweep_values:
        info = _solve_at(V_gate)
        if not info["converged"]:  # pragma: no cover - sweep is well-conditioned
            psi.x.array[:] = psi_backup
            psi.x.scatter_forward()
            raise RuntimeError(
                f"moscap_lf_hf: SNES failed at V_gate={V_gate:+.4f} V "
                f"(reason={info.get('reason')})"
            )
        last_info = info
        psi_backup = psi.x.array.copy()

        rho_int_local = float(fem.assemble_scalar(charge_form))
        rho_int = msh.comm.allreduce(rho_int_local, op=MPI.SUM)
        Q_gate = area_factor * rho_int

        row: dict[str, float] = {
            "V": float(V_gate),
            "Q_gate": float(Q_gate),
            "J": 0.0,
        }

        if "LF" in modes:
            lf_problem.solve()
            delta_psi_lf.x.scatter_forward()
            sens_lf = msh.comm.allreduce(
                float(fem.assemble_scalar(sens_form_lf)), op=MPI.SUM,
            )
            # delta_psi was solved with delta_V_gate = 1 V (=> 1/V_t in
            # scaled units at the gate). The sensitivity integral is
            # already per-unit-V_gate. C = -dQ_gate/dV_gate; the area_factor
            # carries the leading minus sign convention.
            row["C_LF"] = float(area_factor * sens_lf * (-1.0))
        if "HF" in modes:
            hf_problem.solve()
            delta_psi_hf.x.scatter_forward()
            sens_hf = msh.comm.allreduce(
                float(fem.assemble_scalar(sens_form_hf)), op=MPI.SUM,
            )
            row["C_HF"] = float(area_factor * sens_hf * (-1.0))

        iv_rows.append(row)
        if progress_callback is not None:
            progress_callback({
                "type": "step_done",
                "bias_step": len(iv_rows) - 1,
                "V_applied": float(V_gate),
                "iterations": int(info.get("iterations", 0)),
            })

    # Aggregate per-bias capacitance arrays for downstream consumers.
    last_info = dict(last_info)
    last_info["V_g"] = [float(r["V"]) for r in iv_rows]
    if "LF" in modes:
        last_info["C_LF"] = [float(r["C_LF"]) for r in iv_rows]
    if "HF" in modes:
        last_info["C_HF"] = [float(r["C_HF"]) for r in iv_rows]
    last_info["axisymmetric"] = axisymmetric
    last_info["majority_carrier"] = majority
    last_info["capacitance_units"] = normalisation_label

    x_dof = V_psi.tabulate_dof_coordinates()
    psi_phys = psi.x.array * sc.V0
    n_phys = ref_mat.n_i * np.exp(psi.x.array)
    p_phys = ref_mat.n_i * np.exp(-psi.x.array)

    return SimulationResult(
        cfg=cfg, mesh=msh, V=V_psi, psi=psi,
        psi_phys=psi_phys, n_phys=n_phys, p_phys=p_phys,
        x_dof=x_dof, N_hat=N_hat_fn, scaling=sc,
        solver_info=last_info, iv=iv_rows, bias_contact=sweep_contact,
    )


def _build_sensitivity_bcs(cfg, V_psi, msh, facet_tags, sc, sweep_contact):
    """Dirichlet BCs for the sensitivity solve: delta_psi=1/V_t at sweep contact, 0 elsewhere."""
    from dolfinx import fem
    from petsc4py import PETSc

    from ..bcs import resolve_contacts

    contacts = resolve_contacts(cfg, facet_tags=facet_tags, voltages={})
    bcs: list = []
    fdim = msh.topology.dim - 1
    for c in contacts:
        if c.kind not in ("ohmic", "gate"):
            continue
        facets = facet_tags.find(c.facet_tag)
        if len(facets) == 0:  # pragma: no cover
            continue
        dofs = fem.locate_dofs_topological(V_psi, fdim, facets)
        value = 1.0 / sc.V0 if c.name == sweep_contact else 0.0
        bcs.append(fem.dirichletbc(PETSc.ScalarType(value), dofs, V_psi))
    return bcs


def _resolve_semi_tag(regions_cfg: dict) -> int:
    for r in regions_cfg.values():
        if "tag" not in r:  # pragma: no cover
            continue
        if r.get("role", "semiconductor") == "semiconductor":
            return int(r["tag"])
    raise ValueError("moscap_lf_hf: no region with role='semiconductor'.")  # pragma: no cover


def _resolve_gate_sweep(cfg) -> tuple[str | None, list[float]]:
    for c in cfg["contacts"]:
        if c["type"] != "gate":
            continue
        sweep = c.get("voltage_sweep")
        if sweep is None:  # pragma: no cover
            continue
        start = float(sweep["start"])
        stop = float(sweep["stop"])
        step = float(sweep["step"])
        if step <= 0.0:  # pragma: no cover
            raise ValueError("voltage_sweep.step must be positive")
        direction = 1.0 if stop >= start else -1.0
        n = int(round((stop - start) / (direction * step))) + 1
        values = [start + i * direction * step for i in range(n)]
        return c["name"], values
    return None, []


def _fmt_gate_tag(v: float) -> str:
    return f"{v:+.4f}".replace("+", "p").replace("-", "m").replace(".", "d")
