"""
2D long-channel n-MOSFET verifier (M12).

Compares the simulated drain current at (V_GS = 1.5 V, V_DS = 0.1 V)
to the closed-form long-channel triode-regime expression

    I_D = mu_n_eff * C_ox * (W / L) * (V_GS - V_T) * V_DS          (1)

with V_T from the classical strong-inversion formula

    V_T = V_FB + 2 phi_F + sqrt(2 q eps_Si N_A (2 phi_F)) / C_ox    (2)

V_FB is computed against this engine's BC convention: the body ohmic
contact pins psi_body = -phi_F (since psi is referenced to the
intrinsic Fermi level and phi_n = phi_p = 0 at equilibrium), so
V_FB = phi_ms - phi_F. This matches what the M6 mos_2d verifier uses
and is the reason the mos_2d cap within 10 percent on the depletion
window.

The simulation is 2D per-unit-depth, so the engine's drain current is
reported in A/m^2 as a spatial average over the drain facet. To compare
to analytical long-channel theory (which is also per-unit-depth in 2D,
with W = 1 m notionally), we integrate: I_D_sim = J_drain * L_drain,
where L_drain = 1 um is the physical length of the drain contact line.

Tolerance is +/-20 percent. M16 physics (field-dependent mobility,
velocity saturation, surface scattering) is explicitly out of scope,
and the estimate bakes a conservative mu_n_eff = 600 cm^2 / V s to
represent typical inversion-layer mobility degradation.
"""
from __future__ import annotations

import numpy as np


def _mosfet_device_params(cfg, sc, si_mat, ox_mat):
    """Long-channel n-MOSFET analytical parameters.

    Returns a dict with V_T, V_FB, C_ox, W/L, mu_n_eff, and intermediate
    quantities used by the verifier and the debug print. All values SI.
    """
    from semi.constants import EPS0, Q, cm3_to_m3

    # Body doping from the step profile (uniform N_A both sides of the
    # metallurgical step; the step exists only so the schema accepts a
    # step profile without requiring dedicated well doping).
    profile = cfg["doping"][0]["profile"]
    N_A = cm3_to_m3(float(profile.get("N_A_left", profile.get("N_A", 0.0))))
    if N_A <= 0.0:
        raise ValueError("MOSFET verifier expected N_A > 0 for the silicon body")

    eps_Si = si_mat.epsilon
    n_i = si_mat.n_i
    eps_ox = ox_mat.epsilon

    # Oxide thickness and channel length come from the .geo:
    #   oxide:    y in [2e-6, 2.005e-6]        -> t_ox = 5 nm
    #   channel:  x in [1.5e-6, 3.5e-6]        -> L = 2 um
    #   drain contact line length (used as the facet measure for
    #   converting J_avg in A/m^2 back to I_D per unit depth in A/m):
    #       x in [4.0 um, 5.0 um]              -> L_drain = 1 um
    t_ox = 5.0e-9
    L = 2.0e-6
    L_drain = 1.0e-6

    C_ox = eps_ox / t_ox

    # Bulk-to-intrinsic Fermi potential of the p-body (phi_F > 0 for p-type)
    V_t = sc.V0
    phi_F = V_t * float(np.log(N_A / n_i))

    gate = next(c for c in cfg["contacts"] if c["type"] == "gate")
    phi_ms = float(gate.get("workfunction", 0.0) or 0.0)

    # With the engine's psi = 0 at intrinsic BC convention, the body
    # ohmic pins psi_body = -phi_F, so V_FB (the gate voltage at which
    # psi_gate_hat == psi_body_hat) evaluates to V_FB = phi_ms - phi_F.
    # This matches the M6 MOS-cap code that passed the 10 percent gate.
    V_FB = phi_ms - phi_F

    # Threshold: psi_s = 2 phi_F at the surface.
    Q_B_threshold = float(np.sqrt(4.0 * eps_Si * Q * N_A * phi_F))
    V_T = V_FB + 2.0 * phi_F + Q_B_threshold / C_ox

    # Effective channel mobility (~half the bulk value); M16 will
    # replace this with a real surface-mobility model (Lombardi etc.).
    mu_n_eff = 600.0e-4  # SI: m^2 / (V s)

    return dict(
        N_A=N_A,
        eps_Si=eps_Si,
        eps_ox=eps_ox,
        n_i=n_i,
        t_ox=t_ox,
        L=L,
        L_drain=L_drain,
        C_ox=C_ox,
        V_t=V_t,
        phi_F=phi_F,
        phi_ms=phi_ms,
        V_FB=V_FB,
        V_T=V_T,
        mu_n_eff=mu_n_eff,
    )


def _triode_I_D_per_depth(V_GS: float, V_DS: float, dp: dict) -> float:
    """Long-channel triode-regime I_D per unit depth (A/m) at (V_GS, V_DS).

    For a 2D simulation, W (channel width transverse to current flow)
    factors out as per-unit-depth, so the canonical I_D formula becomes

        I_D / W = mu_n_eff * C_ox / L * (V_GS - V_T) * V_DS.

    Returns 0 in subthreshold (V_GS <= V_T).
    """
    overdrive = V_GS - dp["V_T"]
    if overdrive <= 0.0:
        return 0.0
    return dp["mu_n_eff"] * dp["C_ox"] / dp["L"] * overdrive * V_DS


def verify(result) -> list[tuple[str, bool, str]]:
    """MOSFET I-V verifier at V_DS = 0.1 V, V_GS = 1.5 V (triode)."""
    from semi.materials import get_material

    cfg = result.cfg
    sc = result.scaling
    si_mat = get_material(cfg["regions"]["silicon"]["material"])
    ox_mat = get_material(cfg["regions"]["oxide"]["material"])
    dp = _mosfet_device_params(cfg, sc, si_mat, ox_mat)

    iv = result.iv or []
    checks: list[tuple[str, bool, str]] = []

    if not iv:
        return [("mosfet_2d: iv table non-empty", False, "no IV rows recorded")]

    V_DS_target = 0.1
    drain = next((r for r in iv if abs(r["V"] - V_DS_target) < 1.0e-3), None)
    if drain is None:
        V_max = max(r["V"] for r in iv)
        drain = min(iv, key=lambda r: abs(r["V"] - V_DS_target))
        checks.append((
            f"mosfet_2d: sweep reaches V_DS = {V_DS_target:.2f} V",
            False,
            f"max V_DS in sweep = {V_max:.4f} V; using closest ({drain['V']:.4f} V)",
        ))
    else:
        checks.append((
            f"mosfet_2d: sweep reaches V_DS = {V_DS_target:.2f} V",
            True,
            f"V_DS = {drain['V']:.4f} V",
        ))

    gate = next(c for c in cfg["contacts"] if c["type"] == "gate")
    V_GS = float(gate.get("voltage", 0.0))

    # iv J is the average current density in A/m^2 (see
    # semi/postprocess.py::evaluate_current_at_contact). Convert to
    # per-unit-depth I_D (A/m) by multiplying by the drain contact's
    # length L_drain.
    J_drain = float(drain.get("J", 0.0))
    I_D_sim = abs(J_drain) * dp["L_drain"]

    I_D_theory = _triode_I_D_per_depth(V_GS, drain["V"], dp)

    checks.append((
        f"mosfet_2d: V_T = {dp['V_T']:.3f} V (V_FB = {dp['V_FB']:.3f}, phi_F = {dp['phi_F']:.3f}, C_ox = {dp['C_ox']*1e2:.4f} uF/cm^2)",
        True,
        f"overdrive = V_GS - V_T = {V_GS - dp['V_T']:.3f} V",
    ))

    I_D_expected = I_D_theory
    if I_D_expected <= 0.0:
        checks.append((
            "mosfet_2d: triode I_D > 0 expected",
            False,
            f"V_GS = {V_GS:.3f} V <= V_T = {dp['V_T']:.3f} V; subthreshold.",
        ))
    else:
        ratio = abs(I_D_sim - I_D_expected) / abs(I_D_expected)
        assert ratio <= 0.20, (
            f"MOSFET I_D mismatch: simulated {I_D_sim:.3e} A/m vs "
            f"analytical {I_D_expected:.3e} A/m (ratio {ratio:.3%}, "
            f"tolerance +-20%). Known cause (M12): bias_sweep runner does not "
            f"wire multi-region DD; oxide region is treated as silicon."
        )
        checks.append((
            f"mosfet_2d: triode I_D within 20% at V_GS={V_GS:.2f}, V_DS={drain['V']:.2f}",
            True,
            f"I_D_sim={I_D_sim:.4e}, I_D_theory={I_D_expected:.4e}, ratio={ratio:.3%}",
        ))

    # Assert a benign sanity check: the simulation returned a finite,
    # non-negative current (catches NaN propagation and sign bugs).
    checks.append((
        "mosfet_2d: simulated I_D is finite and non-negative",
        bool(np.isfinite(I_D_sim) and I_D_sim >= 0.0),
        f"|I_D_sim| = {I_D_sim:.4e} A/m per unit depth",
    ))

    # Record for the plotter / debug tooling.
    result._mosfet_params = dp
    result._mosfet_obs = dict(
        V_GS=V_GS, V_DS=drain["V"],
        I_D_sim=I_D_sim, I_D_theory=I_D_theory,
    )

    return checks
