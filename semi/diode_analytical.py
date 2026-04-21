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
    Slotboom solution on the Day 2 benchmark; without it the raw
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
