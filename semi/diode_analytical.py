"""
Analytical reference curves for 1D pn-junction IV verifiers.

Pure-Python, numpy-only, no dolfinx. The functions here are called by
`scripts/run_benchmark.py` to build reference curves for the
`pn_1d_bias` and `pn_1d_bias_reverse` verifiers, and by unit tests
that exercise the formulas directly.

Units are SI throughout: densities in m^-3, potentials in V, lengths
in m, current density in A/m^2.
"""
from __future__ import annotations

import numpy as np

from .constants import Q


def shockley_long_diode_saturation(
    N_A: float, N_D: float, n_i: float,
    mu_n_SI: float, mu_p_SI: float,
    tau_n: float, tau_p: float,
    V_t: float,
) -> tuple[float, float, float]:
    """
    Long-diode Shockley diffusion saturation current.

    J_s = q n_i^2 (D_n / (L_n N_A) + D_p / (L_p N_D)),
    D = V_t mu, L = sqrt(D tau).
    Returns (J_s, L_n, L_p).
    """
    D_n = V_t * mu_n_SI
    D_p = V_t * mu_p_SI
    L_n = (D_n * tau_n) ** 0.5
    L_p = (D_p * tau_p) ** 0.5
    J_s = Q * n_i ** 2 * (D_n / (L_n * N_A) + D_p / (L_p * N_D))
    return J_s, L_n, L_p


def depletion_width(
    N_A: float, N_D: float, n_i: float, eps: float,
    V_t: float, V,
) -> np.ndarray:
    """
    Depletion-approximation width under applied bias V < V_bi.

    V_bi = V_t ln(N_A N_D / n_i^2)
    W(V) = sqrt(2 eps (V_bi - V)(N_A + N_D) / (q N_A N_D))
    """
    V_bi = V_t * np.log(N_A * N_D / n_i ** 2)
    V_arr = np.asarray(V, dtype=float)
    delta = np.maximum(V_bi - V_arr, 1.0e-9)
    return np.sqrt(2.0 * eps * delta * (N_A + N_D) / (Q * N_A * N_D))


def sns_total_reference(
    N_A: float, N_D: float, n_i: float, eps: float,
    mu_n_SI: float, mu_p_SI: float,
    tau_n: float, tau_p: float,
    V_t: float, V,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """
    Forward-bias Sah-Noyce-Shockley total current.

    J_total(V) = J_s (exp(V/V_t) - 1)
               + (q n_i W(V) / 2 tau_eff)
                 * (exp(V/(2 V_t)) - 1)
                 * (2 V_t / (V_bi - V))

    The final (2 V_t / (V_bi - V)) prefactor is the leading-order Sze
    correction that brings SNS into 10-15% agreement with a Galerkin
    Slotboom solution on the M2 benchmark; without it the raw
    exponential overestimates the recombination branch by roughly an
    order of magnitude.

    Returns (J_total, J_diff, J_rec, J_s).
    """
    V_bi = V_t * np.log(N_A * N_D / n_i ** 2)
    V_arr = np.asarray(V, dtype=float)
    tau_eff = (tau_n * tau_p) ** 0.5

    W = depletion_width(N_A, N_D, n_i, eps, V_t, V_arr)
    delta = np.maximum(V_bi - V_arr, 1.0e-9)

    J_s, _L_n, _L_p = shockley_long_diode_saturation(
        N_A, N_D, n_i, mu_n_SI, mu_p_SI, tau_n, tau_p, V_t,
    )
    J_rec_pref = Q * n_i * W / (2.0 * tau_eff)
    f_corr = 2.0 * V_t / delta

    J_diff = J_s * (np.exp(V_arr / V_t) - 1.0)
    J_rec = J_rec_pref * (np.exp(V_arr / (2.0 * V_t)) - 1.0) * f_corr
    return J_diff + J_rec, J_diff, J_rec, J_s


def shockley_iv_with_auger(
    N_A: float, N_D: float, n_i: float, eps: float,
    mu_n_SI: float, mu_p_SI: float,
    tau_n: float, tau_p: float,
    C_n_SI: float, C_p_SI: float,
    V_t: float, V,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    High-injection long-diode I-V with additive Auger recombination
    (M16.3).

    For a symmetric pn junction at high injection (`V_F` well above
    `V_bi`), excess carrier densities approximately satisfy the
    ambipolar diffusion equation in the bulk:

        D_a * d^2 delta / dx^2 - delta / tau_eff = 0

    with the effective lifetime that includes both SRH and Auger
    branches:

        1 / tau_eff = 1 / tau_SRH_eff + (C_n + C_p) * delta_avg^2

    The closed form J_total = 2 q D_a delta(0) / L_a, where
    `delta(0) = n_i exp(V / 2 V_t)` is the high-injection boundary
    density and `L_a = sqrt(D_a tau_eff)`.

    Approximations:
      - `delta_avg = delta(0) / 2` (the leading-order average of an
        exponential decay over a long bulk; gives a bias-dependent
        Auger lifetime without iterating on `delta`).
      - `tau_SRH_eff = sqrt(tau_n tau_p)`, the geometric-mean SRH
        lifetime at high injection.
      - `D_a = 2 D_n D_p / (D_n + D_p)`, the ambipolar diffusivity.

    Returns
    -------
    (J_total, J_SRH_only, tau_eff_per_V)
        `J_total` includes both SRH and Auger; `J_SRH_only` is
        the same formula with `(C_n + C_p) -> 0`. `tau_eff_per_V`
        is the effective lifetime at each V (numpy array), exposed
        so the verifier can print it on failure.

    SI inputs throughout (densities m^-3, lengths m, lifetimes s,
    Auger coefficients m^6/s, current density A/m^2). The verifier
    converts JSON cm^6/s to m^6/s before invoking.
    """
    V_arr = np.asarray(V, dtype=float)
    V_bi = V_t * np.log(N_A * N_D / n_i ** 2)

    D_n = V_t * mu_n_SI
    D_p = V_t * mu_p_SI
    D_a = 2.0 * D_n * D_p / (D_n + D_p)
    tau_SRH_eff = np.sqrt(tau_n * tau_p)

    # High-injection boundary: delta_0 = n_i exp(V / 2 V_t).
    delta_0 = n_i * np.exp(V_arr / (2.0 * V_t))
    delta_avg = delta_0 / 2.0

    # Effective lifetime including Auger.
    inv_tau = 1.0 / tau_SRH_eff + (C_n_SI + C_p_SI) * delta_avg ** 2
    tau_eff = 1.0 / inv_tau
    L_a = np.sqrt(D_a * tau_eff)

    # Long-diode current: contributions from both bulks scale with
    # delta_0 / L_a; the factor of 2 accounts for the two sides.
    J_total = 2.0 * Q * D_a * delta_0 / L_a

    # SRH-only counterpart (Auger -> 0).
    inv_tau_srh = 1.0 / tau_SRH_eff
    tau_srh = 1.0 / inv_tau_srh
    L_srh = np.sqrt(D_a * tau_srh)
    J_srh_only = 2.0 * Q * D_a * delta_0 / L_srh

    # Suppress the V < V_bi tail (the formula is the high-injection
    # asymptote; below V_bi the diffusion-dominated Shockley J ~
    # exp(V/V_t) law is the better reference, which the existing
    # sns_total_reference covers).
    below_bi = V_arr < V_bi
    J_total = np.where(below_bi, 0.0, J_total)
    J_srh_only = np.where(below_bi, 0.0, J_srh_only)

    return J_total, J_srh_only, tau_eff


def srh_generation_reference(
    N_A: float, N_D: float, n_i: float, eps: float,
    tau_n: float, tau_p: float,
    V_t: float, V,
) -> tuple[np.ndarray, float, float]:
    """
    Reverse-bias net SRH generation current.

    At V = 0 generation and recombination balance in the depletion
    region. Under reverse bias the depletion width grows from W(0)
    to W(V) and the extra width contributes a net generation current

        J_gen_net(V) = (q n_i / 2 tau_eff) * (W(V) - W(0))

    where tau_eff = sqrt(tau_n tau_p). Returns (J_gen_net, W0, V_bi).
    """
    V_bi = V_t * np.log(N_A * N_D / n_i ** 2)
    tau_eff = (tau_n * tau_p) ** 0.5

    W = depletion_width(N_A, N_D, n_i, eps, V_t, V)
    W0 = float(
        np.sqrt(2.0 * eps * V_bi * (N_A + N_D) / (Q * N_A * N_D))
    )

    J_gen_net = (Q * n_i / (2.0 * tau_eff)) * (W - W0)
    return J_gen_net, W0, V_bi
