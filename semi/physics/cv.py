"""
C-V curve computation for axisymmetric MOS capacitor (LF and HF).

LF (quasi-static): equilibrium Poisson at each Vg — both majority and
minority carriers respond to the gate voltage.

HF (frozen minority): at each Vg, solve LF to get psi_lf, then freeze
the minority carrier density n_hat = ni*exp(psi_lf) and re-solve Poisson
so only majority holes respond to gate perturbations.

Both use central differences on Q_sub to get C:
    C(Vg) = -[Q_sub(Vg+dV) - Q_sub(Vg-dV)] / (2*dV)

Q_sub is the total semiconductor charge integrated over the 2D cross
section with r-weighting (the 2*pi azimuthal factor is included so that
Q_sub has units of Coulombs per gate area in C/m^2 after division by A_gate).
"""
from __future__ import annotations

import math
import ufl
from typing import Any

import numpy as np

# Guard against division-by-zero in central-difference denominator.
# Use a threshold safely above float64 machine epsilon (~2.2e-16).
_MIN_VOLTAGE_STEP = 1.0e-12


def compute_cv_curve(
    cfg: dict[str, Any],
    mode: str = "LF",
    Vg_sweep: list[float] | None = None,
    dV: float = 0.005,
):
    """
    Compute LF or HF C-V curves for the axisymmetric MOSCAP.

    Parameters
    ----------
    cfg : dict
        Validated config dict. Must have dimension='axisymmetric_2d' and
        mesh.source='builtin_axi'.
    mode : str
        'LF' for quasi-static (equilibrium Poisson) or 'HF' for frozen
        minority carrier (minority n is frozen at each Vg).
    Vg_sweep : list of float
        Gate voltages to evaluate. If None, uses np.arange(-1.5, 1.55, 0.1)
        giving voltages at 0.1 V intervals from -1.5 V to +1.5 V.
    dV : float
        Half-step for central-difference capacitance calculation.

    Returns
    -------
    dict with keys:
        'Vg'    : list of gate voltages
        'C'     : list of capacitances in F/m^2 (per unit gate area)
        'Q_sub' : list of semiconductor charge in C/m^2
    """
    import ufl
    from dolfinx import fem
    from mpi4py import MPI
    from petsc4py import PETSc

    from ..bcs import build_psi_dirichlet_bcs, resolve_contacts
    from ..constants import Q as Q_e
    from ..doping import build_profile
    from ..mesh import build_eps_r_function, build_mesh
    from ..physics.poisson import (
        build_equilibrium_poisson_form_mr_axi,
    )
    from ..run import SimulationResult
    from ..runners._common import reference_material
    from ..scaling import make_scaling_from_config
    from ..solver import solve_nonlinear

    ref_mat = reference_material(cfg)
    sc = make_scaling_from_config(cfg, ref_mat)

    msh, cell_tags, facet_tags = build_mesh(cfg)
    if cell_tags is None:
        raise RuntimeError("compute_cv_curve requires cell_tags (multi-region mesh).")

    semi_tag = _resolve_semi_tag(cfg["regions"])
    eps_r_fn = build_eps_r_function(msh, cell_tags, cfg["regions"])

    V = fem.functionspace(msh, ("Lagrange", 1))
    N_raw_fn = build_profile(cfg["doping"])

    N_hat_fn = fem.Function(V, name="N_net_hat")
    N_hat_fn.interpolate(lambda x: N_raw_fn(x) / sc.C0)

    psi = fem.Function(V, name="psi_hat")
    two_ni = 2.0 * ref_mat.n_i
    psi.interpolate(lambda x: np.arcsinh(N_raw_fn(x) / two_ni))

    # Build the LF (equilibrium) r-weighted Poisson form
    F_lf = build_equilibrium_poisson_form_mr_axi(
        V, psi, N_hat_fn, sc, eps_r_fn, cell_tags, semi_tag,
    )

    # Identify gate and body contacts
    gate_name = _find_contact_name_by_type(cfg, ("gate", "mos_gate"))
    R_gate = float(cfg["mesh"].get("R_gate", cfg["mesh"]["extents"]["r"][1]))
    A_gate = math.pi * R_gate ** 2

    # Default sweep
    if Vg_sweep is None:
        Vg_sweep = list(np.arange(-1.5, 1.55, 0.1))

    # ------------------------------------------------------------------ #
    # Charge integral helper
    # ------------------------------------------------------------------ #
    x_coord = ufl.SpatialCoordinate(msh)
    r_coord = x_coord[0]
    ni_hat_const = fem.Constant(msh, PETSc.ScalarType(sc.n_i / sc.C0))
    dx_semi = ufl.Measure(
        "dx", domain=msh, subdomain_data=cell_tags, subdomain_id=int(semi_tag)
    )

    # LF: rho = ni*(exp(-psi) - exp(psi)) + N_hat  (both carriers respond)
    rho_hat_lf = ni_hat_const * (ufl.exp(-psi) - ufl.exp(psi)) + N_hat_fn
    charge_form_lf = fem.form(rho_hat_lf * r_coord * dx_semi)

    def integrate_charge_lf():
        """Integrate rho * r * dx over Si, multiply 2*pi*q*C0 -> Coulombs."""
        local_val = float(fem.assemble_scalar(charge_form_lf))
        total = msh.comm.allreduce(local_val, op=MPI.SUM)
        return 2.0 * math.pi * Q_e * sc.C0 * total

    def solve_lf_at(Vg: float, prefix_tag: str):
        """Solve equilibrium Poisson at gate voltage Vg."""
        voltages = {gate_name: Vg}
        contacts = resolve_contacts(cfg, facet_tags=facet_tags, voltages=voltages)
        bcs = build_psi_dirichlet_bcs(
            V, msh, facet_tags, contacts, sc, ref_mat, N_raw_fn
        )
        for bc in bcs:
            bc.set(psi.x.array)
        psi.x.scatter_forward()
        return solve_nonlinear(
            F_lf, psi, bcs,
            prefix=f"cv_{prefix_tag}_",
            petsc_options={
                "snes_rtol": 1.0e-12,
                "snes_atol": 1.0e-14,
                "snes_max_it": 50,
            },
        )

    # ------------------------------------------------------------------ #
    # LF sweep
    # ------------------------------------------------------------------ #
    if mode == "LF":
        Q_list = []
        Vg_list = []
        psi_backup = psi.x.array.copy()

        for Vg in Vg_sweep:
            tag = _vg_tag(Vg)
            info = solve_lf_at(Vg, f"lf_{tag}")
            if not info["converged"]:
                import warnings
                warnings.warn(
                    f"LF solve failed to converge at Vg={Vg:.4f} V", stacklevel=2
                )
                psi.x.array[:] = psi_backup
                psi.x.scatter_forward()
            else:
                psi_backup = psi.x.array.copy()

            Q_total = integrate_charge_lf()
            Q_list.append(Q_total)
            Vg_list.append(Vg)

        C_list = _central_diff_C(Vg_list, Q_list, A_gate)
        Q_per_area = [q / A_gate for q in Q_list]
        return {"Vg": Vg_list, "C": C_list, "Q_sub": Q_per_area}

    # ------------------------------------------------------------------ #
    # HF sweep: at each Vg, solve LF, freeze n_hat, re-solve Poisson
    # ------------------------------------------------------------------ #
    elif mode == "HF":
        C_list = []
        Q_list = []
        Vg_list = []
        psi_backup = psi.x.array.copy()

        for Vg in Vg_sweep:
            tag = _vg_tag(Vg)

            # 1. Solve LF at Vg
            info = solve_lf_at(Vg, f"hf_lf_{tag}")
            if not info["converged"]:
                import warnings
                warnings.warn(f"HF base solve failed at Vg={Vg:.4f} V", stacklevel=2)
                psi.x.array[:] = psi_backup
                psi.x.scatter_forward()
            else:
                psi_backup = psi.x.array.copy()

            # Record Q at this Vg (LF charge)
            Q_total = integrate_charge_lf()
            Q_list.append(Q_total)
            Vg_list.append(Vg)

            # 2. Freeze minority carrier density using dolfinx Expression
            #    n_hat_frozen = ni_hat * exp(psi_lf)  evaluated via DG0 interpolation
            V_DG0 = fem.functionspace(msh, ("DG", 0))
            n_hat_frozen_fn = fem.Function(V_DG0, name="n_hat_frozen")
            ni_hat_val = sc.n_i / sc.C0
            n_expr = fem.Expression(
                fem.Constant(msh, PETSc.ScalarType(ni_hat_val)) * ufl.exp(psi),
                V_DG0.element.interpolation_points(),
            )
            n_hat_frozen_fn.interpolate(n_expr)
            n_hat_frozen_fn.x.scatter_forward()

            # 3. Build HF frozen Poisson form
            psi_hf = fem.Function(V, name=f"psi_hf_{tag}")
            psi_hf.x.array[:] = psi.x.array.copy()
            psi_hf.x.scatter_forward()

            def _solve_hf_frozen(psi_fn, Vg_val, ptag):
                voltages = {gate_name: Vg_val}
                contacts = resolve_contacts(
                    cfg, facet_tags=facet_tags, voltages=voltages
                )
                bcs_h = build_psi_dirichlet_bcs(
                    V, msh, facet_tags, contacts, sc, ref_mat, N_raw_fn
                )
                for bc in bcs_h:
                    bc.set(psi_fn.x.array)
                psi_fn.x.scatter_forward()
                F_h = _build_frozen_n_poisson_axi(
                    V, psi_fn, n_hat_frozen_fn, N_hat_fn, sc, eps_r_fn,
                    cell_tags, semi_tag
                )
                info_h = solve_nonlinear(
                    F_h, psi_fn, bcs_h,
                    prefix=f"cv_hf_frz_{ptag}_",
                    petsc_options={
                        "snes_rtol": 1.0e-12,
                        "snes_atol": 1.0e-14,
                        "snes_max_it": 50,
                    },
                )
                if not info_h["converged"]:
                    import warnings
                    warnings.warn(
                        f"HF frozen Poisson failed at Vg={Vg_val:.4f}", stacklevel=2
                    )
                return psi_fn

            def _q_hf(psi_fn):
                """Integrate HF charge: rho_HF = p - n_frozen + N."""
                p_hat_expr = ni_hat_const * ufl.exp(-psi_fn)
                rho_hf = p_hat_expr - n_hat_frozen_fn + N_hat_fn
                cf = fem.form(rho_hf * r_coord * dx_semi)
                local = float(fem.assemble_scalar(cf))
                total = msh.comm.allreduce(local, op=MPI.SUM)
                return 2.0 * math.pi * Q_e * sc.C0 * total

            # Solve at Vg+dV
            psi_p = fem.Function(V, name=f"psi_hf_p_{tag}")
            psi_p.x.array[:] = psi.x.array.copy()
            psi_p.x.scatter_forward()
            _solve_hf_frozen(psi_p, Vg + dV, f"{tag}_p")

            # Solve at Vg-dV
            psi_m = fem.Function(V, name=f"psi_hf_m_{tag}")
            psi_m.x.array[:] = psi.x.array.copy()
            psi_m.x.scatter_forward()
            _solve_hf_frozen(psi_m, Vg - dV, f"{tag}_m")

            Q_p = _q_hf(psi_p)
            Q_m = _q_hf(psi_m)
            C_hf = -(Q_p - Q_m) / (2.0 * dV) / A_gate
            C_list.append(C_hf)

        Q_per_area = [q / A_gate for q in Q_list]
        return {"Vg": Vg_list, "C": C_list, "Q_sub": Q_per_area}

    else:
        raise ValueError(f"Unknown C-V mode {mode!r}; use 'LF' or 'HF'.")


def _build_frozen_n_poisson_axi(
    V, psi, n_hat_frozen_fn, N_hat_fn, sc, eps_r_fn, cell_tags, semi_tag
):
    """
    Axisymmetric Poisson with frozen minority electron density.

    rho = p_hat - n_hat_frozen + N_hat
        = ni*exp(-psi) - n_hat_frozen + N_hat

    n_hat_frozen is a precomputed DG0 Function (does not depend on psi).
    """
    import ufl
    from dolfinx import fem
    from petsc4py import PETSc

    msh = V.mesh
    x = ufl.SpatialCoordinate(msh)
    r = x[0]
    v = ufl.TestFunction(V)

    L_D2 = fem.Constant(msh, PETSc.ScalarType(sc.lambda2 * sc.L0 ** 2))
    ni_hat = fem.Constant(msh, PETSc.ScalarType(sc.n_i / sc.C0))

    dx_full = ufl.Measure("dx", domain=msh)
    dx_semi = ufl.Measure(
        "dx", domain=msh, subdomain_data=cell_tags, subdomain_id=int(semi_tag)
    )

    # p responds to psi; n is frozen
    p_hat = ni_hat * ufl.exp(-psi)
    rho_hat = p_hat - n_hat_frozen_fn + N_hat_fn

    F = (
        L_D2 * eps_r_fn * ufl.inner(ufl.grad(psi), ufl.grad(v)) * r * dx_full
        - rho_hat * v * r * dx_semi
    )
    return F


def _central_diff_C(Vg_list, Q_list, A_gate):
    """
    C = -dQ/dVg per unit gate area, using central differences.
    For endpoints, use one-sided differences.
    """
    n = len(Vg_list)
    C_list = []
    for i in range(n):
        if i == 0:
            dQ = Q_list[1] - Q_list[0]
            dv = Vg_list[1] - Vg_list[0]
        elif i == n - 1:
            dQ = Q_list[-1] - Q_list[-2]
            dv = Vg_list[-1] - Vg_list[-2]
        else:
            dQ = Q_list[i + 1] - Q_list[i - 1]
            dv = Vg_list[i + 1] - Vg_list[i - 1]
        C_list.append(-dQ / dv / A_gate if abs(dv) > _MIN_VOLTAGE_STEP else 0.0)
    return C_list


def _resolve_semi_tag(regions_cfg: dict) -> int:
    """First region with role='semiconductor'."""
    for r in regions_cfg.values():
        if "tag" not in r:
            continue
        if r.get("role", "semiconductor") == "semiconductor":
            return int(r["tag"])
    raise ValueError("No region with role='semiconductor'.")


def _find_contact_name_by_type(cfg, kinds) -> str:
    """Return the name of the first contact whose type is in kinds."""
    for c in cfg["contacts"]:
        if c["type"] in kinds:
            return c["name"]
    raise ValueError(f"No contact of type {', '.join(kinds)} found in config.")


def _vg_tag(v: float) -> str:
    """PETSc-prefix-safe encoding of a gate voltage."""
    return f"{v:+.4f}".replace("+", "p").replace("-", "m").replace(".", "d")
