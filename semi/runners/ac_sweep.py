"""
Small-signal AC sweep runner (M14).

Linearises the steady-state drift-diffusion system around a converged
DC operating point and solves the complex linear system

    (J + j*omega*M) * delta_u(omega) = -dF/dV * delta_V

at each frequency in a user-specified sweep. The Jacobian J is the
steady-state DD operator at u_0; the mass matrix M is the lumped P1
diagonal on the carrier-density rows (psi has no time derivative —
Poisson's equation is instantaneous in the quasi-electrostatic limit
used throughout this engine; see ADR 0011).

Primary unknowns (matching the M13 transient runner): (psi, n_hat, p_hat).
The DC operating point is obtained by running :func:`run_bias_sweep`
in Slotboom variables (psi, phi_n, phi_p) and converting via Boltzmann
statistics. See ADR 0011 for the design choices and ADR 0009 for why
the primary-density form is used in time-dependent contexts.

PETSc complex-vs-real-block: the dolfinx-real PETSc build (the only
build available in CI today) uses real scalars. M14 therefore
implements the standard real 2x2 block reformulation:

    [  J   -omega*M  ] [ Re delta_u ]   [ Re b ]
    [ omega*M   J    ] [ Im delta_u ] = [ Im b ]

The block size is 2 * (3N) = 6N where N is the per-block DOF count.
A direct LU (MUMPS) factorisation is used at each frequency. For the
benchmark mesh sizes shipped with M14 (~80-200 cells in 1D), this is
sub-second per frequency. M16+ may switch to a complex-PETSc build,
which would shrink the matrix to 3N complex entries; the runner
contract (the AcSweepResult shape) does not depend on this choice.
"""
from __future__ import annotations

from typing import Any

import numpy as np


def _resolve_frequencies(freq_spec: dict) -> np.ndarray:
    """Expand the schema's frequency-sweep object to a numpy array of Hz."""
    kind = freq_spec["type"]
    if kind == "list":
        return np.asarray(freq_spec["values"], dtype=float)
    if kind == "logspace":
        start = float(freq_spec["start"])
        stop = float(freq_spec["stop"])
        n = int(freq_spec["n_points"])
        return np.logspace(np.log10(start), np.log10(stop), n)
    if kind == "linspace":
        start = float(freq_spec["start"])
        stop = float(freq_spec["stop"])
        n = int(freq_spec["n_points"])
        return np.linspace(start, stop, n)
    raise ValueError(f"Unknown ac.frequencies.type {kind!r}")


def _build_dc_subcfg(cfg: dict, contact_name: str, voltage: float) -> dict:
    """Build a bias_sweep cfg pinned at a single DC bias point.

    Returns a deep copy of `cfg` with solver.type set to "bias_sweep" and
    a single-point voltage_sweep on `contact_name` at `voltage`. Any
    pre-existing voltage_sweep on other contacts is removed (the sweep
    runner picks the first ohmic contact with a voltage_sweep, so we
    have to be the only one).
    """
    import copy

    sub = copy.deepcopy(cfg)
    # Inherit the user's solver options (snes, continuation) but switch
    # the type to bias_sweep and pin the sweep at the AC operating point.
    # We must override max_step / min_step to a sane pair: the AC-side
    # eps_V perturbation can drive the nominal_step far below the user's
    # min_step default (1e-4 V), which would trip the AdaptiveStepController
    # min_step >= max_step guard.
    snes_overrides = (cfg.get("solver", {}).get("snes") or {})
    user_cont = (cfg.get("solver", {}).get("continuation") or {})
    target_abs = abs(float(voltage))
    # Cap nominal_step so the DC ramp from V=0 to V_DC uses several
    # continuation steps rather than one giant jump (essential for high
    # reverse bias on diodes; see pn_1d_bias_reverse benchmark, which
    # uses voltage_sweep.step = 0.1 V for a -2 V endpoint). 0.1 V is a
    # safe upper bound for SRH-limited DD; the user can override
    # cfg.solver.continuation.max_step to grow it.
    nominal_step = max(min(target_abs, 0.1), 1.0e-3)
    cont = dict(user_cont)
    cont["max_step"] = max(float(cont.get("max_step", nominal_step)), nominal_step)
    # min_step must not exceed max_step.  Cap at nominal_step * 0.1 so a
    # tiny eps_V perturbation can subdivide aggressively if needed.
    cont["min_step"] = min(
        float(cont.get("min_step", 1.0e-4)),
        cont["max_step"] * 0.1,
        max(nominal_step * 0.01, 1.0e-7),
    )
    cont.setdefault("max_halvings", 10)
    cont.setdefault("easy_iter_threshold", 4)
    cont.setdefault("grow_factor", 1.5)
    sub["solver"] = {
        "type": "bias_sweep",
        "continuation": cont,
    }
    if snes_overrides:
        sub["solver"]["snes"] = dict(snes_overrides)
    found = False
    for c in sub["contacts"]:
        if c["name"] == contact_name:
            c["voltage"] = float(voltage)
            c["voltage_sweep"] = {
                "start": float(voltage),
                "stop": float(voltage),
                "step": nominal_step,
            }
            found = True
        else:
            c.pop("voltage_sweep", None)
    if not found:
        raise ValueError(
            f"ac_sweep dc_bias.contact {contact_name!r} not found in cfg.contacts"
        )
    return sub


def run_ac_sweep(cfg: dict[str, Any], *, progress_callback=None):
    """
    Small-signal AC sweep around a DC operating point.

    Parameters
    ----------
    cfg : dict
        Validated config dict. ``solver.type`` must be ``"ac_sweep"`` and
        ``solver.dc_bias`` and ``solver.ac`` must be provided per the
        schema (see ``schemas/input.v1.json``).
    progress_callback : callable or None
        If provided, called with ``{"type": "ac_step", "freq_index": k,
        "frequency": f}`` after each frequency point.

    Returns
    -------
    AcSweepResult
        See :class:`semi.results.AcSweepResult`.
    """
    import math

    import ufl
    from dolfinx import fem
    from dolfinx.fem.petsc import assemble_matrix

    from ..doping import build_profile
    from ..fem.mass import assemble_lumped_mass
    from ..mesh import build_mesh
    from ..physics.drift_diffusion import make_dd_block_spaces
    from ..postprocess import resolve_contact_facets
    from ..results import AcSweepResult
    from ..runners.bias_sweep import run_bias_sweep
    from ..scaling import make_scaling_from_config
    from ._common import reference_material

    solver_cfg = cfg.get("solver", {})
    if solver_cfg.get("type") != "ac_sweep":
        raise ValueError("run_ac_sweep requires solver.type == 'ac_sweep'")
    dc_cfg = solver_cfg["dc_bias"]
    ac_cfg = solver_cfg["ac"]
    contact_name = dc_cfg["contact"]
    V_DC = float(dc_cfg["voltage"])
    ac_amplitude = float(ac_cfg.get("amplitude", 1.0))
    frequencies = _resolve_frequencies(ac_cfg["frequencies"])

    # ------------------------------------------------------------------
    # 1. Solve the DC operating point at V_DC and at V_DC + eps_V.
    #    The pair lets us form the DC sensitivity delta_u_DC by finite
    #    difference, which gives J * delta_u_DC = -dF/dV at u_0 (since
    #    F is zero at both points and the linearisation is exact to
    #    O(eps_V) in the steady-state response). The bias_sweep runner
    #    is reused so we do not duplicate the multi-region BC and SNES
    #    plumbing.
    # ------------------------------------------------------------------
    # Fixed-magnitude perturbation. Smaller than V_t = 0.026 V (so the
    # response is genuinely small-signal) but large enough to stay well
    # above the SNES tolerance floor (atol ~1e-7 in scaled units).
    eps_V = 1.0e-3

    dc_sub = _build_dc_subcfg(cfg, contact_name, V_DC)
    res_dc = run_bias_sweep(dc_sub)
    dc_iter = int(res_dc.solver_info.get("iterations", 0))

    dc_sub_eps = _build_dc_subcfg(cfg, contact_name, V_DC + eps_V)
    res_dc_eps = run_bias_sweep(dc_sub_eps)
    dc_iter = max(dc_iter, int(res_dc_eps.solver_info.get("iterations", 0)))

    # ------------------------------------------------------------------
    # 2. Build a fresh (psi, n, p) function-space stack on the same
    #    mesh and project the Slotboom DC solution into primary-density
    #    form. We keep our own mesh build so the runner is independent
    #    of the dc-subcfg's plotter / artifact side effects.
    # ------------------------------------------------------------------
    ref_mat = reference_material(cfg)
    sc = make_scaling_from_config(cfg, ref_mat)
    msh, _cell_tags, facet_tags = build_mesh(cfg)

    N_raw_fn = build_profile(cfg["doping"])
    N_hat_fn = fem.Function(make_dd_block_spaces(msh).V_psi, name="N_net_hat")
    N_hat_fn.interpolate(lambda x: N_raw_fn(x) / sc.C0)

    spaces = make_dd_block_spaces(msh)
    V_psi = spaces.V_psi
    V_n = spaces.V_phi_n
    V_p = spaces.V_phi_p

    psi = spaces.psi
    n_hat = fem.Function(V_n, name="n_hat")
    p_hat = fem.Function(V_p, name="p_hat")

    ni_over_C0 = sc.n_i / sc.C0

    def _slotboom_to_density(res):
        psi_arr = res.psi.x.array
        phi_n_arr = res.phi_n.x.array
        phi_p_arr = res.phi_p.x.array
        n_arr = ni_over_C0 * np.exp(psi_arr - phi_n_arr)
        p_arr = ni_over_C0 * np.exp(phi_p_arr - psi_arr)
        return psi_arr.copy(), n_arr, p_arr

    psi_0, n_0, p_0 = _slotboom_to_density(res_dc)
    psi_eps, n_eps, p_eps = _slotboom_to_density(res_dc_eps)

    psi.x.array[:] = psi_0
    n_hat.x.array[:] = n_0
    p_hat.x.array[:] = p_0
    for fn in (psi, n_hat, p_hat):
        fn.x.scatter_forward()

    # DC sensitivity delta_u_DC per unit dV (concatenated [psi, n, p]).
    delta_psi_DC = (psi_eps - psi_0) / eps_V
    delta_n_DC = (n_eps - n_0) / eps_V
    delta_p_DC = (p_eps - p_0) / eps_V

    N_psi = psi.x.array.size
    N_n = n_hat.x.array.size

    # ------------------------------------------------------------------
    # 3. Material constants and physics-model coefficients (same as
    #    bias_sweep / transient).
    # ------------------------------------------------------------------
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
    tau_n_hat_v = tau_n_s / sc.t0
    tau_p_hat_v = tau_p_s / sc.t0
    E_t_over_Vt = E_t_eV / sc.V0

    # ------------------------------------------------------------------
    # 4. Build the steady-state (psi, n, p) residual and its Jacobian
    #    bilinear forms a[i][j] = derivative(F[i], u[j]) for i,j in
    #    {psi, n, p}. The 3x3 block layout is what we'll assemble below.
    # ------------------------------------------------------------------
    F_list = _build_steady_density_residual(
        psi, n_hat, p_hat, N_hat_fn,
        spaces, sc, ref_mat.epsilon_r,
        mu_n_hat, mu_p_hat, tau_n_hat_v, tau_p_hat_v, E_t_over_Vt,
    )
    a_blocks = [
        [ufl.derivative(F_list[i], u_j) for u_j in (psi, n_hat, p_hat)]
        for i in range(3)
    ]

    # ------------------------------------------------------------------
    # 5. Build BCs at the DC operating point. The block Jacobian is
    #    assembled with these BCs so BC rows are zeroed with diag=1
    #    (standard PETSc lifting). Carrier BCs are V-independent for
    #    ohmic contacts (n_bc, p_bc are set by the local doping); only
    #    psi shifts with the applied bias.
    # ------------------------------------------------------------------
    static_voltages = {
        c["name"]: float(c.get("voltage", 0.0))
        for c in cfg["contacts"]
        if c["type"] == "ohmic"
    }
    static_voltages[contact_name] = V_DC

    bcs_dc = _build_npp_bcs(
        cfg, msh, facet_tags, static_voltages, spaces, sc, ref_mat, N_raw_fn,
    )

    # Compiled forms for assembly.
    a_forms = [[fem.form(a_blocks[i][j]) for j in range(3)] for i in range(3)]

    # Block Jacobian as a single PETSc Mat (size 3N x 3N, monolithic).
    A_J = assemble_matrix(a_forms, bcs=bcs_dc, diag=1.0)
    A_J.assemble()

    # ------------------------------------------------------------------
    # 6. Lumped mass diagonal on the (n, p) blocks. The Poisson row has
    #    no time derivative -> M_psi = 0. We embed M as a 3N-long
    #    diagonal vector aligned with the same DOF ordering as A_J.
    # ------------------------------------------------------------------
    dx = ufl.Measure("dx", domain=msh)
    M_n_diag, M_p_diag = assemble_lumped_mass(V_n, V_p, dx)

    # Convert to numpy local arrays for the block solver. The transient
    # residual (semi/runners/transient.py) uses the time-derivative term
    # (alpha_0 / dt_hat) * n_hat * v dx with dt_hat = dt_phys/t0. In
    # frequency domain d/dt -> j*omega, so (1/dt_phys) -> j*omega and the
    # bilinear coefficient becomes j*omega*t0*M_lumped. Pre-multiplying
    # the diagonal by t0 here lets the solver use omega in physical
    # rad/s when forming the 2x2 real block. See ADR 0011.
    M_n_arr = M_n_diag.array.copy() * sc.t0
    M_p_arr = M_p_diag.array.copy() * sc.t0
    M_psi_arr = np.zeros(N_psi, dtype=float)

    # Zero-out BC rows in M (BCs are imposed as identity at all freqs;
    # the BC perturbation is the source, not subject to dt dynamics).
    bc_n_dofs, bc_p_dofs, bc_psi_dofs = _collect_bc_dof_indices(bcs_dc, V_psi, V_n, V_p)
    M_n_arr[bc_n_dofs] = 0.0
    M_p_arr[bc_p_dofs] = 0.0
    M_psi_arr[bc_psi_dofs] = 0.0

    # Concatenated 3N mass diagonal in the same block order as A_J:
    # [psi-block | n-block | p-block].
    M_diag_3N = np.concatenate([M_psi_arr, M_n_arr, M_p_arr])

    # ------------------------------------------------------------------
    # 7. RHS direction b = J @ delta_u_DC. delta_u_DC is the per-unit-V
    #    DC sensitivity. By construction, J @ delta_u_DC = -dF/dV at u_0.
    #    We use a PETSc mat-vec to stay consistent with the block
    #    layout (same ordering as A_J).
    # ------------------------------------------------------------------
    delta_u_DC_3N = np.concatenate([delta_psi_DC, delta_n_DC, delta_p_DC])
    b_RHS = np.zeros_like(delta_u_DC_3N)
    _matvec_into(A_J, delta_u_DC_3N, b_RHS)

    # At BC rows, A_J has identity, so (A_J @ delta_u_DC)[bc] = delta_u_DC[bc].
    # That is exactly the BC perturbation per-unit-V (e.g. 1.0 at the
    # swept-contact psi BC), which is the right RHS for the lifted
    # AC linear system. No further patching required.

    # ------------------------------------------------------------------
    # 8. Frequency loop. At each omega build the 6N x 6N real 2x2 block
    #    system, solve via direct LU, and compute terminal current
    #    perturbation Y(omega).
    # ------------------------------------------------------------------
    Y_list: list[complex] = []
    Z_list: list[complex] = []
    C_list: list[float] = []
    G_list: list[float] = []

    sweep_facet_info = resolve_contact_facets(cfg, msh, facet_tags, contact_name)

    psi_fn_pert = fem.Function(V_psi, name="delta_psi")
    n_fn_pert = fem.Function(V_n, name="delta_n")
    p_fn_pert = fem.Function(V_p, name="delta_p")

    for k, f in enumerate(frequencies):
        omega = 2.0 * math.pi * float(f)

        x_real, x_imag = _solve_2x2_real_block(
            A_J, M_diag_3N, omega, b_RHS,
        )

        delta_psi_re = x_real[:N_psi]
        delta_n_re = x_real[N_psi:N_psi + N_n]
        delta_p_re = x_real[N_psi + N_n:]
        delta_psi_im = x_imag[:N_psi]
        delta_n_im = x_imag[N_psi:N_psi + N_n]
        delta_p_im = x_imag[N_psi + N_n:]

        # Conduction current perturbation at the contact (linearised):
        #   delta_J_n = -q*mu_n*(delta_n*grad(psi_0) + n_0*grad(delta_psi))
        #               + q*mu_n*V_t * grad(delta_n)
        # and similarly for J_p (with sign flip on the drift coupling).
        # The contact-facet integral converts this to a current density
        # in A/m^2 (per the existing per-unit-area convention).
        I_cond_re = _eval_linearised_conduction_current(
            psi, n_hat, p_hat,
            delta_psi_re, delta_n_re, delta_p_re,
            psi_fn_pert, n_fn_pert, p_fn_pert,
            sc, mu_n_SI, mu_p_SI, sweep_facet_info, msh,
        )
        I_cond_im = _eval_linearised_conduction_current(
            psi, n_hat, p_hat,
            delta_psi_im, delta_n_im, delta_p_im,
            psi_fn_pert, n_fn_pert, p_fn_pert,
            sc, mu_n_SI, mu_p_SI, sweep_facet_info, msh,
        )

        # Displacement current at the contact:
        #   I_disp = j*omega * eps * integral(grad(delta_psi).n_outward) dS
        # In phasor form: I_total = I_cond + j*omega*eps*Q_psi(delta_psi).
        Q_psi_re = _eval_displacement_charge(
            delta_psi_re, psi_fn_pert, sc, ref_mat.epsilon_r,
            sweep_facet_info, msh,
        )
        Q_psi_im = _eval_displacement_charge(
            delta_psi_im, psi_fn_pert, sc, ref_mat.epsilon_r,
            sweep_facet_info, msh,
        )

        # delta_D = -eps * grad(delta_psi); the displacement-evaluator
        # returns Q_psi = +eps * grad(delta_psi).n_outward, so the
        # contact-outward displacement field is -Q_psi. Therefore
        #   I_disp_out = j*omega * delta_D_outward = -j*omega * Q_psi
        # I_total_out = I_cond_out - j*omega*(Q_re + j*Q_im)
        #             = (I_cond_re + omega*Q_im) + j*(I_cond_im - omega*Q_re)
        # The conduction-current evaluator and the displacement evaluator
        # both contract with the FE FacetNormal (outward), so I_total_out
        # is the small-signal terminal current flowing OUT of the device.
        # bias_sweep / postprocess.evaluate_current_at_contact reports
        # dI/dV with the opposite sign (positive at forward bias);
        # algebraically, in (n,p)-primary form the linearised conduction
        # current d(J_bs)/dV = -[delta_J_n_out + delta_J_p_out] (the two
        # forms differ by a global sign because Slotboom phi_n carries the
        # opposite sign to the (n,p)-primary outward flux; see ADR 0011
        # "Errata"). To produce a Y in the same convention as
        # bias_sweep -- so that Re(Y(omega->0)) == dI/dV at the same V_DC
        # -- we negate the assembled OUT-convention current. With this
        # convention an ideal capacitor has Y = +j*omega*C and the
        # corresponding capacitance read-out is C = +Im(Y)/(2*pi*f).
        I_total_re = -(I_cond_re + omega * Q_psi_im)
        I_total_im = -(I_cond_im - omega * Q_psi_re)

        Y = complex(I_total_re, I_total_im)
        Y_list.append(Y)

        if abs(Y) < 1.0e-300:
            Z_list.append(complex(0.0, 0.0))
        else:
            Z_list.append(1.0 / Y)

        if f > 0.0:
            # Convention: Y is reported with the "positive = current INTO
            # device" sign used by bias_sweep / postprocess
            # .evaluate_current_at_contact. An ideal capacitor terminal
            # has Y = +j*omega*C, so a positive depletion capacitance
            # maps to a positive Im(Y) and C = +Im(Y)/(2*pi*f).
            C_list.append(Y.imag / (2.0 * math.pi * float(f)))
        else:
            C_list.append(0.0)
        G_list.append(Y.real)

        if progress_callback is not None:
            progress_callback({
                "type": "ac_step",
                "freq_index": k,
                "frequency": float(f),
            })

    return AcSweepResult(
        frequencies=[float(f) for f in frequencies],
        Y=Y_list,
        Z=Z_list,
        C=C_list,
        G=G_list,
        dc_bias={"contact": contact_name, "voltage": V_DC},
        meta={
            "n_freqs": len(frequencies),
            "perturbation_voltage": float(eps_V),
            "ac_amplitude": ac_amplitude,
            "dc_solver_iterations": int(dc_iter),
            "backend": "real-2x2-block",
            "displacement_current_included": True,
        },
    )


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------


def _build_steady_density_residual(
    psi, n_hat, p_hat, N_hat_fn,
    spaces, sc, eps_r,
    mu_n_over_mu0: float,
    mu_p_over_mu0: float,
    tau_n_hat: float,
    tau_p_hat: float,
    E_t_over_Vt: float,
):
    """Steady-state DD residual in (psi, n_hat, p_hat) primary-density form.

    Identical to the spatial part of `_build_transient_residual` from
    `semi.runners.transient`, with the time-derivative term dropped.
    The Jacobian of this form (taken below via ufl.derivative) is the
    operator J that drives the AC solve.
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

    n1 = ni_hat_c * _math.exp(E_t_over_Vt)
    p1 = ni_hat_c * _math.exp(-E_t_over_Vt)
    R = (n_hat * p_hat - ni_hat_c * ni_hat_c) / (
        tau_p_c * (n_hat + n1) + tau_n_c * (p_hat + p1)
    )

    rho_hat = p_hat - n_hat + N_hat_fn

    F_psi = (
        L_D2 * eps_r_ufl * ufl.inner(ufl.grad(psi), ufl.grad(v_psi)) * ufl.dx
        - rho_hat * v_psi * ufl.dx
    )
    F_n = (
        L0_sq * mu_n_c * (
            ufl.inner(ufl.grad(n_hat), ufl.grad(v_n))
            - n_hat * ufl.inner(ufl.grad(psi), ufl.grad(v_n))
        ) * ufl.dx
        + R * v_n * ufl.dx
    )
    F_p = (
        L0_sq * mu_p_c * (
            ufl.inner(ufl.grad(p_hat), ufl.grad(v_p))
            + p_hat * ufl.inner(ufl.grad(psi), ufl.grad(v_p))
        ) * ufl.dx
        - R * v_p * ufl.dx
    )
    return [F_psi, F_n, F_p]


def _build_npp_bcs(cfg, msh, facet_tags, voltages: dict[str, float],
                   spaces, sc, ref_mat, N_raw_fn) -> list:
    """Build Dirichlet BCs for (psi, n_hat, p_hat) at every ohmic contact.

    Mirrors the BC builder inside the transient runner. The carrier BC
    values are V-independent (n_bc, p_bc come from the local doping at
    the contact); psi shifts by V_hat = V/V_t.
    """
    from dolfinx import fem as _fem
    from petsc4py import PETSc

    from ..bcs import _evaluate_doping_at_facet, resolve_contacts

    V_psi = spaces.V_psi
    V_n = spaces.V_phi_n
    V_p = spaces.V_phi_p

    contacts = resolve_contacts(cfg, facet_tags=facet_tags, voltages=voltages)
    bcs: list = []
    fdim = msh.topology.dim - 1
    two_ni = 2.0 * ref_mat.n_i
    ni_hat_val = sc.n_i / sc.C0

    for c in contacts:
        if c.kind != "ohmic":
            continue
        facets = facet_tags.find(c.facet_tag)
        if len(facets) == 0:
            continue
        N_net = _evaluate_doping_at_facet(msh, facets, fdim, N_raw_fn)
        psi_eq_hat = float(np.arcsinh(N_net / two_ni))
        V_hat = c.V_applied / sc.V0
        psi_bc = psi_eq_hat + V_hat
        n_bc = ni_hat_val * float(np.exp(psi_eq_hat))
        p_bc = ni_hat_val * float(np.exp(-psi_eq_hat))

        dofs_psi = _fem.locate_dofs_topological(V_psi, fdim, facets)
        dofs_n = _fem.locate_dofs_topological(V_n, fdim, facets)
        dofs_p = _fem.locate_dofs_topological(V_p, fdim, facets)

        bcs.append(_fem.dirichletbc(PETSc.ScalarType(psi_bc), dofs_psi, V_psi))
        bcs.append(_fem.dirichletbc(PETSc.ScalarType(n_bc), dofs_n, V_n))
        bcs.append(_fem.dirichletbc(PETSc.ScalarType(p_bc), dofs_p, V_p))
    return bcs


def _collect_bc_dof_indices(bcs, V_psi, V_n, V_p):
    """Pull the local DOF indices each Dirichlet BC constrains, by space."""
    psi_dofs: list[int] = []
    n_dofs: list[int] = []
    p_dofs: list[int] = []
    for bc in bcs:
        space_id = id(bc.function_space)
        if space_id == id(V_psi):
            psi_dofs.extend(bc._cpp_object.dof_indices()[0].tolist())
        elif space_id == id(V_n):
            n_dofs.extend(bc._cpp_object.dof_indices()[0].tolist())
        elif space_id == id(V_p):
            p_dofs.extend(bc._cpp_object.dof_indices()[0].tolist())
    return (
        np.asarray(n_dofs, dtype=np.int64),
        np.asarray(p_dofs, dtype=np.int64),
        np.asarray(psi_dofs, dtype=np.int64),
    )


def _matvec_into(A, x_arr: np.ndarray, out_arr: np.ndarray) -> None:
    """y = A @ x for monolithic block PETSc Mats, via local arrays."""
    x_vec = A.createVecRight()
    y_vec = A.createVecLeft()
    x_vec.array[:] = x_arr
    A.mult(x_vec, y_vec)
    out_arr[:] = y_vec.array
    x_vec.destroy()
    y_vec.destroy()


def _solve_2x2_real_block(A_J, M_diag_3N: np.ndarray,
                          omega: float, b_RHS: np.ndarray):
    """Solve the real 2x2 block AC system [J -wM; wM J][x;y]=[b;0].

    A_J is a PETSc Mat (3N x 3N) representing the steady-state Jacobian
    with BC rows replaced by identity. M_diag_3N is a numpy diagonal of
    length 3N, zero on BC rows and on the psi block. omega is the
    angular frequency in rad/s scaled to be consistent with J (M was
    pre-divided by t0 in run_ac_sweep, so omega is in physical rad/s).

    Returns a tuple (x_arr, y_arr) of length-3N numpy arrays, the real
    and imaginary parts of the solution. Uses a direct LU solve via
    PETSc's nested KSP. For the M14 benchmark mesh sizes (sub-1k DOF
    in 1D) this is sub-second per frequency.
    """
    from petsc4py import PETSc

    N3 = b_RHS.size
    n_local = N3
    big_size = 2 * N3

    # Build a 6N x 6N CSR matrix in PETSc by extracting J as CSR and
    # composing. For benchmark-sized meshes a Python-level CSR build is
    # cheap; the dominant cost is the per-frequency LU below.
    J_csr = A_J.getValuesCSR()
    J_indptr, J_indices, J_data = J_csr

    # Top-left and bottom-right blocks: J. Top-right: -omega*M (diag).
    # Bottom-left: +omega*M (diag).
    # Build COO triples then CSR.
    from scipy.sparse import bmat, csr_matrix, diags

    J_sp = csr_matrix((J_data, J_indices, J_indptr), shape=(N3, N3))
    wM = diags(omega * M_diag_3N, 0, shape=(N3, N3), format="csr")

    # bmat builds the big block matrix.
    big = bmat([[J_sp, -wM], [wM, J_sp]], format="csr")

    # Build PETSc Mat from CSR.
    A_big = PETSc.Mat().createAIJ(
        size=(big_size, big_size),
        csr=(big.indptr.astype(PETSc.IntType),
             big.indices.astype(PETSc.IntType),
             big.data.astype(PETSc.RealType)),
        comm=A_J.getComm(),
    )
    A_big.assemble()

    b_big = A_big.createVecLeft()
    x_big = A_big.createVecRight()
    b_big.array[:n_local] = b_RHS
    b_big.array[n_local:] = 0.0

    ksp = PETSc.KSP().create(A_big.getComm())
    ksp.setOperators(A_big)
    ksp.setType("preonly")
    pc = ksp.getPC()
    pc.setType("lu")
    pc.setFactorSolverType("mumps")
    ksp.solve(b_big, x_big)

    x_arr_full = x_big.array.copy()
    x_real = x_arr_full[:N3].copy()
    x_imag = x_arr_full[N3:].copy()

    ksp.destroy()
    A_big.destroy()
    b_big.destroy()
    x_big.destroy()

    return x_real, x_imag


def _eval_linearised_conduction_current(
    psi, n_hat, p_hat,
    delta_psi_arr, delta_n_arr, delta_p_arr,
    psi_pert_fn, n_pert_fn, p_pert_fn,
    sc, mu_n_SI, mu_p_SI, facet_info, msh,
) -> float:
    """Evaluate delta_I_conduction at the contact for one part (real or imag).

    Linearised current density (per the n,p primary-variable form, with
    Einstein relation D = mu*V_t):
        delta_J_n = q*mu_n*V_t*grad(delta_n) - q*mu_n*(delta_n*grad(psi_0)
                                                        + n_0*grad(delta_psi))
        delta_J_p = -q*mu_p*V_t*grad(delta_p) - q*mu_p*(delta_p*grad(psi_0)
                                                          + p_0*grad(delta_psi))

    All density-like quantities are scaled internally; the final result
    is multiplied back by C_0 to recover physical A/m^2.
    """
    import ufl
    from dolfinx import fem
    from dolfinx.mesh import meshtags
    from mpi4py import MPI

    from ..constants import Q

    psi_pert_fn.x.array[:] = delta_psi_arr
    n_pert_fn.x.array[:] = delta_n_arr
    p_pert_fn.x.array[:] = delta_p_arr
    for fn in (psi_pert_fn, n_pert_fn, p_pert_fn):
        fn.x.scatter_forward()

    tag = int(facet_info["tag"])
    facets = facet_info["facets"]
    fdim = msh.topology.dim - 1
    if len(facets) == 0:
        return 0.0
    values = np.full(len(facets), tag, dtype=np.int32)
    indices = np.asarray(facets, dtype=np.int32)
    sort = np.argsort(indices)
    facet_mt = meshtags(msh, fdim, indices[sort], values[sort])
    ds = ufl.Measure("ds", domain=msh, subdomain_data=facet_mt, subdomain_id=tag)
    n_vec = ufl.FacetNormal(msh)

    # Convert to physical units inside the form. n_hat -> n_phys = n_hat*C_0,
    # grad(n_hat) -> grad(n_hat) since the mesh is in meters (Invariant 3).
    grad_psi_phys = sc.V0 * ufl.grad(psi)
    grad_dpsi_phys = sc.V0 * ufl.grad(psi_pert_fn)
    n_phys = n_hat * sc.C0
    p_phys = p_hat * sc.C0
    dn_phys = n_pert_fn * sc.C0
    dp_phys = p_pert_fn * sc.C0
    grad_dn_phys = sc.C0 * ufl.grad(n_pert_fn)
    grad_dp_phys = sc.C0 * ufl.grad(p_pert_fn)

    V_t = sc.V0  # k_B T / q in V

    # Conduction current density (linearised) dotted with the FE
    # FacetNormal, which already points OUT of the computational domain
    # at boundary facets. This produces the OUT-of-device current; the
    # runner negates Y at assembly to convert to the INTO-device
    # convention used by bias_sweep / postprocess (see ADR 0011 Errata
    # and the assembly site in run_ac_sweep).
    Jn_lin = (
        Q * mu_n_SI * V_t * ufl.dot(grad_dn_phys, n_vec)
        - Q * mu_n_SI * (
            dn_phys * ufl.dot(grad_psi_phys, n_vec)
            + n_phys * ufl.dot(grad_dpsi_phys, n_vec)
        )
    )
    Jp_lin = (
        - Q * mu_p_SI * V_t * ufl.dot(grad_dp_phys, n_vec)
        - Q * mu_p_SI * (
            dp_phys * ufl.dot(grad_psi_phys, n_vec)
            + p_phys * ufl.dot(grad_dpsi_phys, n_vec)
        )
    )

    form_J = fem.form((Jn_lin + Jp_lin) * ds)
    form_A = fem.form(1.0 * ds)
    I_local = fem.assemble_scalar(form_J)
    A_local = fem.assemble_scalar(form_A)
    I_val = msh.comm.allreduce(I_local, op=MPI.SUM)
    A_val = msh.comm.allreduce(A_local, op=MPI.SUM)
    if A_val == 0.0:
        A_val = 1.0
    return float(I_val / A_val)


def _eval_displacement_charge(
    delta_psi_arr, psi_pert_fn, sc, eps_r,
    facet_info, msh,
) -> float:
    """Evaluate Q_psi = eps * integral(grad(delta_psi) . n) dS / area.

    The displacement-current contribution at frequency omega is
    j*omega*Q_psi (in A/m^2 when divided by the contact area, as is
    done here for consistency with the conduction-current evaluator).
    """
    import ufl
    from dolfinx import fem
    from dolfinx.mesh import meshtags
    from mpi4py import MPI

    from ..constants import EPS0

    psi_pert_fn.x.array[:] = delta_psi_arr
    psi_pert_fn.x.scatter_forward()

    tag = int(facet_info["tag"])
    facets = facet_info["facets"]
    fdim = msh.topology.dim - 1
    if len(facets) == 0:
        return 0.0
    values = np.full(len(facets), tag, dtype=np.int32)
    indices = np.asarray(facets, dtype=np.int32)
    sort = np.argsort(indices)
    facet_mt = meshtags(msh, fdim, indices[sort], values[sort])
    ds = ufl.Measure("ds", domain=msh, subdomain_data=facet_mt, subdomain_id=tag)
    n_vec = ufl.FacetNormal(msh)

    # eps_r may be a Function or a float; either is fine in UFL.
    if isinstance(eps_r, (int, float)):
        eps_r_v = float(eps_r)
    else:
        eps_r_v = eps_r

    grad_dpsi_phys = sc.V0 * ufl.grad(psi_pert_fn)
    # Q_psi = eps0 * eps_r * grad(delta_psi).n_FE / area. n_FE points
    # outward from the computational domain at boundary facets. The
    # physical displacement field is D = -eps*grad(psi), so the
    # outward-D-field is -Q_psi and the OUT-of-device displacement
    # current at the contact is I_disp_out = -j*omega*Q_psi. As with
    # the conduction-current evaluator, the runner converts the
    # OUT-convention total current to the bias_sweep INTO convention by
    # negating Y at assembly.
    D_n_form = EPS0 * eps_r_v * ufl.dot(grad_dpsi_phys, n_vec)
    form_Q = fem.form(D_n_form * ds)
    form_A = fem.form(1.0 * ds)
    Q_local = fem.assemble_scalar(form_Q)
    A_local = fem.assemble_scalar(form_A)
    Q_val = msh.comm.allreduce(Q_local, op=MPI.SUM)
    A_val = msh.comm.allreduce(A_local, op=MPI.SUM)
    if A_val == 0.0:
        A_val = 1.0
    return float(Q_val / A_val)
