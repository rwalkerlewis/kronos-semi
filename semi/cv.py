"""
MOSCAP C-V analytics and FEM postprocessing helpers.

Two layers live here:

1. Pure-analytical textbook formulas (Hu chapter 5): flat-band voltage,
   surface-potential equation, threshold voltage, maximum depletion
   width, oxide capacitance, and the high-frequency C-V curve under the
   depletion approximation. These functions take only scalar inputs and
   are used both by the verifier and the regression tests; they do not
   import dolfinx.

2. FEM postprocessing helpers (`compute_lf_cv`, `compute_hf_cv_from_fem`)
   that consume a sweep of FEM solutions. The LF helper differentiates
   the total semiconductor charge with a centered finite difference
   (Hu Eq. 5.6.1). The HF helper uses the depletion-approximation
   clamp: once the FEM-extracted surface potential reaches 2*phi_B,
   the depletion width is frozen at W_dmax and C_HF is computed from
   the series oxide-plus-depletion stack. This is an explicit choice
   stated in the docstring, justified for a benchmark that targets
   Hu Fig. 5-18 visual reproduction; a fully principled small-signal
   solve with frozen minority carriers is a future extension.

References
----------
- Hu, Modern Semiconductor Devices for Integrated Circuits, ch. 5.
- Sze and Ng, Physics of Semiconductor Devices, ch. 4.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .constants import EPS0, Q, cm3_to_m3, thermal_voltage

EPS_R_SI = 11.7
EPS_R_SIO2 = 3.9


@dataclass(frozen=True)
class MoscapAnalytic:
    """Bundle of analytical MOSCAP parameters at a fixed bias = 0 reference.

    All quantities are SI: F/m^2 for capacitance per unit area, V for
    voltages, m^-3 for densities, m for lengths.
    """
    Vt: float          # thermal voltage kT/q
    phi_B: float       # bulk Fermi potential, |phi_B| = Vt * ln(N_a/n_i)
    V_fb: float        # flat-band voltage, V_fb = phi_ms - Q_f/C_ox
    V_t: float         # threshold voltage, n-MOS on p-body
    W_dmax: float      # maximum depletion width at strong inversion
    C_ox_per_area: float   # eps_ox / T_ox, F/m^2
    C_min_per_area: float  # series of C_ox and C_dep,max, F/m^2
    body_dopant: str   # "p" or "n"

    def C_min_normalized(self) -> float:
        return self.C_min_per_area / self.C_ox_per_area


def analytical_moscap_params(
    *,
    body_dopant: str = "p",
    N_body_cm3: float = 5.0e16,
    T_ox_m: float = 10.0e-9,
    phi_ms: float = -0.95,
    Q_f_per_area: float = 0.0,
    T: float = 300.0,
    n_i_cm3: float = 1.0e10,
    eps_r_semi: float = EPS_R_SI,
    eps_r_ox: float = EPS_R_SIO2,
) -> MoscapAnalytic:
    """
    Closed-form MOSCAP parameters from textbook formulas (Hu ch. 5).

    The defaults match the prompt's reference parameters: p-Si body,
    Na = 5e16 cm^-3, Tox = 10 nm SiO2, N+ poly gate
    (phi_ms ~ -0.95 V for Na = 5e16). Qf is taken to be zero unless
    overridden.

    Parameters
    ----------
    body_dopant : "p" | "n"
        Sign of the body's majority carrier; controls the sign of phi_B
        and which threshold (n- or p-channel) is computed.
    N_body_cm3 : float
        Body net doping magnitude, cm^-3.
    T_ox_m : float
        Gate oxide thickness, m.
    phi_ms : float
        Gate-to-body workfunction difference, V.
    Q_f_per_area : float
        Fixed oxide-charge sheet density (per area, C/m^2). Defaults 0.
    T : float
        Temperature, K.
    n_i_cm3 : float
        Intrinsic density at T, cm^-3.

    Returns
    -------
    MoscapAnalytic
        Frozen dataclass with all derived quantities in SI units.

    Notes
    -----
    - Hu Eq. 5.4.3:   V_fb = phi_ms - Q_f / C_ox
    - Hu Eq. 5.5.1:   V_t  = V_fb + 2 phi_B + sqrt(2 eps_s q N_a (2 phi_B)) / C_ox
    - Hu Eq. 5.6.4:   W_dmax = sqrt(2 eps_s (2 phi_B) / (q N_a))
    """
    Vt = thermal_voltage(T)
    N_body = cm3_to_m3(N_body_cm3)
    n_i = cm3_to_m3(n_i_cm3)
    eps_s = eps_r_semi * EPS0
    eps_o = eps_r_ox * EPS0

    if body_dopant not in ("p", "n"):
        raise ValueError(f"body_dopant must be 'p' or 'n', got {body_dopant!r}")

    # |phi_B| = Vt * ln(N_body / n_i); positive magnitude. Sign convention
    # for the bulk-Fermi level relative to mid-gap: p-body -> +|phi_B|,
    # n-body -> -|phi_B|. We carry the magnitude internally.
    phi_B_mag = Vt * math.log(N_body / n_i)
    phi_B_signed = +phi_B_mag if body_dopant == "p" else -phi_B_mag

    C_ox = eps_o / T_ox_m
    V_fb = phi_ms - Q_f_per_area / C_ox

    # Strong-inversion surface potential phi_s = +/- 2|phi_B|
    psi_s_inv = 2.0 * phi_B_mag
    W_dmax = math.sqrt(2.0 * eps_s * psi_s_inv / (Q * N_body))

    # Threshold voltage: for p-body / nMOS the bulk charge term is positive,
    # for n-body / pMOS it is negative. Sign mirrors phi_B_signed.
    Q_dep_max = math.sqrt(2.0 * eps_s * Q * N_body * psi_s_inv)
    if body_dopant == "p":
        V_t = V_fb + 2.0 * phi_B_mag + Q_dep_max / C_ox
    else:
        V_t = V_fb - 2.0 * phi_B_mag - Q_dep_max / C_ox

    C_dep_max_per_area = eps_s / W_dmax
    C_min_per_area = 1.0 / (1.0 / C_ox + 1.0 / C_dep_max_per_area)

    # Phi_B is reported with its (signed) physical convention, so
    # downstream callers can use it directly in formulas like phi_s = 2 phi_B.
    _ = phi_B_signed  # documented sign; magnitude returned below
    return MoscapAnalytic(
        Vt=Vt,
        phi_B=phi_B_mag,
        V_fb=V_fb,
        V_t=V_t,
        W_dmax=W_dmax,
        C_ox_per_area=C_ox,
        C_min_per_area=C_min_per_area,
        body_dopant=body_dopant,
    )


def hf_cv_depletion_approximation(
    Vg: np.ndarray,
    params: MoscapAnalytic,
    N_body_cm3: float = 5.0e16,
    eps_r_semi: float = EPS_R_SI,
) -> np.ndarray:
    """
    Closed-form high-frequency C-V under the depletion approximation.

    For accumulation (Vg < V_fb on a p-body), C tends to C_ox.
    For depletion (V_fb < Vg < V_t), the surface potential phi_s solves

        Vg - V_fb = phi_s + sqrt(2 q eps_s N_body phi_s) / C_ox       (Hu 5.6.2)

    and W_dep = sqrt(2 eps_s phi_s / (q N_body)), giving
    C = (1/C_ox + W_dep/eps_s)^-1.

    For inversion (Vg >= V_t on a p-body), W_dep clamps at W_dmax and
    C stays at C_min.

    Returns C / C_ox normalized.
    """
    Vg = np.asarray(Vg, dtype=float)
    N_body = cm3_to_m3(N_body_cm3)
    eps_s = eps_r_semi * EPS0
    C_ox = params.C_ox_per_area

    out = np.empty_like(Vg)
    sign = +1.0 if params.body_dopant == "p" else -1.0
    Vox_eff = sign * (Vg - params.V_fb)

    for i, dV in enumerate(Vox_eff):
        if dV <= 0.0:
            out[i] = C_ox  # accumulation
            continue
        # Solve quadratic for sqrt(phi_s):
        #   C_ox * phi_s + sqrt(2 q eps_s N) sqrt(phi_s) - C_ox * dV = 0
        a = C_ox
        b = math.sqrt(2.0 * Q * eps_s * N_body)
        c = -C_ox * dV
        disc = b * b - 4.0 * a * c
        sqrt_phi = (-b + math.sqrt(disc)) / (2.0 * a)
        phi_s = sqrt_phi * sqrt_phi
        if phi_s >= 2.0 * params.phi_B:
            phi_s = 2.0 * params.phi_B  # inversion clamp
        W_dep = math.sqrt(2.0 * eps_s * phi_s / (Q * N_body))
        out[i] = 1.0 / (1.0 / C_ox + W_dep / eps_s)

    return out / C_ox


def lf_cv_quasistatic(
    Vg: np.ndarray,
    params: MoscapAnalytic,
    N_body_cm3: float = 5.0e16,
    eps_r_semi: float = EPS_R_SI,
    n_i_cm3: float = 1.0e10,
) -> np.ndarray:
    """
    Closed-form low-frequency / quasi-static C-V (Hu Fig. 5-18 upper curve).

    In strong inversion the inversion charge follows the AC signal, so

        Q_s(phi_s) = - sqrt(2 eps_s kT N_body) * F(phi_s, n_i^2/N_body^2)

    with F the standard MOS charge function

        F(u, np_ratio) = sqrt( exp(-u) + u - 1
                              + np_ratio * ( exp(u) - u - 1 ) )       (Hu 5.6.3)

    and u = phi_s / Vt. The body-side normal-field continuity gives

        Vg - V_fb = phi_s + sqrt(2 eps_s kT N_body) F(phi_s) / C_ox

    so we invert for phi_s by Brent root-finding then differentiate
    Q_s with respect to Vg numerically. Returns C/C_ox.

    Sign convention: for p-body, phi_s is positive when surface bends
    upward (depletion of holes, inversion of electrons).
    """
    from scipy.optimize import brentq

    Vg = np.asarray(Vg, dtype=float)
    N_body = cm3_to_m3(N_body_cm3)
    n_i = cm3_to_m3(n_i_cm3)
    eps_s = eps_r_semi * EPS0
    Vt = params.Vt
    C_ox = params.C_ox_per_area
    sign = +1.0 if params.body_dopant == "p" else -1.0
    np_ratio = (n_i / N_body) ** 2  # very small for typical doping

    Ldebye_kT = math.sqrt(2.0 * eps_s * Vt * Q * N_body)
    # equals sqrt(2 eps_s kT N) since kT = q*Vt; carries units C/m^2.
    pref = Ldebye_kT

    def F_func(u: float) -> float:
        # u = phi_s / Vt (signed). The MOS body-charge function is
        #   F^2(u) = (e^{-u} + u - 1) + (n_i^2/N^2) (e^{u} - u - 1)
        # which is asymmetric in u: e^{-u} dominates in accumulation
        # (u<0 for p-body) and the n_i^2/N^2 e^{u} term dominates in
        # strong inversion (u>2*phi_B/Vt). Cap |u| at 50 to avoid
        # overflow at biases far outside the physically interesting
        # range. The sqrt argument is clamped at 0 from below to absorb
        # roundoff for tiny |u|.
        u_c = max(-50.0, min(50.0, u))
        arg = math.expm1(-u_c) + u_c + np_ratio * (math.expm1(u_c) - u_c)
        if arg < 0.0:
            arg = 0.0
        return math.copysign(math.sqrt(arg), u_c)

    def residual(phi_s: float, dV: float) -> float:
        u = phi_s / Vt
        Q_s = -pref * F_func(u)  # C/m^2
        return phi_s - sign * Q_s / C_ox - dV

    Q_s_arr = np.empty_like(Vg)
    for i, Vgi in enumerate(Vg):
        dV = sign * (Vgi - params.V_fb)
        # Wide static bracket; F_func is monotone in phi_s so brentq is robust.
        phi_s = brentq(residual, -1.5, 1.5, args=(dV,), maxiter=400, xtol=1e-9)
        u = phi_s / Vt
        Q_s_arr[i] = -pref * F_func(u)

    # Centered FD on Q_s vs Vg; Q_s is the body-side charge so
    # C_LF = -dQ_s/dVg per unit area (Hu Eq. 5.6.1 with sign convention).
    C_LF = np.gradient(-Q_s_arr, Vg)
    return C_LF / C_ox


def compute_lf_cv_fem(Vg: np.ndarray, Q_s: np.ndarray) -> np.ndarray:
    """
    Centered finite-difference low-frequency capacitance from a FEM sweep.

    Parameters
    ----------
    Vg : array shape (N,)
        Gate voltages in the sweep, monotonic.
    Q_s : array shape (N,)
        Total semiconductor charge per unit area at each Vg (C/m^2).
        Caller is responsible for the r-weighted volume integral over
        the meridian times 2*pi divided by the gate area.

    Returns
    -------
    C_LF : array shape (N,)
        -dQ_s/dVg in F/m^2.
    """
    Vg = np.asarray(Vg, dtype=float)
    Q_s = np.asarray(Q_s, dtype=float)
    if Vg.shape != Q_s.shape:
        raise ValueError("Vg and Q_s must have the same shape")
    return -np.gradient(Q_s, Vg)


def compute_hf_cv_depletion_clamp(
    Vg: np.ndarray,
    phi_s: np.ndarray,
    params: MoscapAnalytic,
    N_body_cm3: float = 5.0e16,
    eps_r_semi: float = EPS_R_SI,
) -> np.ndarray:
    """
    High-frequency capacitance from FEM-extracted surface potentials,
    using the depletion-approximation clamp.

    HF method (explicit): for each Vg, take phi_s extracted from the
    FEM solution at the Si/SiO2 interface (centerline r=0). Compute the
    depletion width W_dep = sqrt(2 eps_s |phi_s| / (q N_body)) capped at
    W_dmax; then C_HF = (1/C_ox + W_dep/eps_s)^-1.

    Rationale. Once phi_s reaches 2 phi_B the inversion layer dominates
    the DC charge but cannot follow MHz AC; freezing W at W_dmax matches
    the textbook HF curve. Hu chapter 5 derives this explicitly. A
    fully principled small-signal solve with delta_n=0 in p-body
    inversion is a future extension; this clamp suffices for the visual
    Fig. 5-18 reproduction.

    Returns C_HF in F/m^2.
    """
    Vg = np.asarray(Vg, dtype=float)
    phi_s = np.asarray(phi_s, dtype=float)
    if Vg.shape != phi_s.shape:
        raise ValueError("Vg and phi_s must have the same shape")

    N_body = cm3_to_m3(N_body_cm3)
    eps_s = eps_r_semi * EPS0
    C_ox = params.C_ox_per_area
    W_dmax = params.W_dmax

    out = np.empty_like(Vg)
    for i, ps in enumerate(phi_s):
        # In accumulation phi_s on a p-body has the opposite sign;
        # depletion charge is zero so C -> C_ox.
        ps_eff = ps if params.body_dopant == "p" else -ps
        if ps_eff <= 0.0:
            out[i] = C_ox
            continue
        ps_clamped = min(ps_eff, 2.0 * params.phi_B)
        W = math.sqrt(2.0 * eps_s * ps_clamped / (Q * N_body))
        W = min(W, W_dmax)
        out[i] = 1.0 / (1.0 / C_ox + W / eps_s)
    return out
