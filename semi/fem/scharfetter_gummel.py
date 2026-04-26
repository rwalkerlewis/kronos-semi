"""
Scharfetter-Gummel edge-flux primitives (M13.1).

This module is the load-bearing primitive for the SG transient
discretisation specified in ADR 0012. It provides:

    * `bernoulli(x)`            : numerically stable B(x) = x / (exp(x) - 1).
    * `bernoulli_array(x)`      : vectorised B over a NumPy array.
    * `sg_edge_flux_n(...)`     : 1D electron edge flux per ADR 0012.
    * `sg_edge_flux_p(...)`     : 1D hole edge flux per ADR 0012.

The 2D simplicial assembly is layered on top of these primitives in a
follow-up commit; per the M13.1 plan we land 1D first, gate on
steady-state agreement, then build 2D. This module is pure-Python /
NumPy with no dolfinx import at module scope so it can be imported
and unit-tested without a FEM toolchain (Invariant 4).

Sign convention (ADR 0012 § "Sign convention")
----------------------------------------------
Sandia / Farrell-et-al convention: `F_n_ij` is the electron contribution
to the conventional current density on the edge from i to j. At zero
field with n_j > n_i, F > 0 (current in +∇n direction); at uniform n
with psi_j > psi_i, F < 0 (electrons drift toward j, conventional
current opposite). This matches the textbook
`J_n = q mu_n n E + q D_n grad(n)` once the kronos-semi phi_n sign
convention is unwound (see notes in `semi/postprocess.py`,
`evaluate_partial_currents`). The user-prompt formula
`F = (mu V_t/h)[B(-dpsi) n_j - B(+dpsi) n_i]` had the B-arguments
swapped relative to the standard convention; the swap was caught by
the midpoint-Galerkin cross-check guard before any 1D residual code
was wired. See ADR 0012 "Sources actually used" for the trail.

Bernoulli function regimes
--------------------------
- `x == 0`        : exact 1.0.
- `|x| < taylor_eps` (default 1e-3): 4th-order Taylor 1 - x/2 + x^2/12 - x^4/720.
- `x > big_x` (default 30): rewrite as `x exp(-x) / (1 - exp(-x))` to keep
                            the denominator bounded away from zero and let
                            `x exp(-x)` underflow benignly to 0.
- `x < -big_x` (default -30): closed form `x / (exp(x) - 1)` is fine; `exp(x)`
                              underflows to 0 and the denominator goes to -1,
                              so `B(x) -> -x` cleanly.
- mid-range      : closed form `x / (exp(x) - 1)`.

The identity `B(x) - B(-x) = -x` is verified by unit test across the full
sampling range. (The earlier prompt draft stated `B(x) + B(-x) = -x`, which
is wrong; the sum is `x coth(x/2)`. The corrected identity is the one
used in tests.)
"""
from __future__ import annotations

import math

import numpy as np

_TAYLOR_EPS = 1.0e-3
_BIG_X = 30.0


def bernoulli(x: float) -> float:
    """
    Scalar Bernoulli `B(x) = x / (exp(x) - 1)` with numerically stable
    branches in all four regimes.

    Parameters
    ----------
    x : float
        Argument (dimensionless; for SG flux this is the scaled potential
        difference `psi_j - psi_i`).

    Returns
    -------
    float
        `B(x)`. `B(0) = 1.0` exactly.
    """
    if x == 0.0:
        return 1.0
    ax = abs(x)
    if ax < _TAYLOR_EPS:
        x2 = x * x
        return 1.0 - 0.5 * x + x2 / 12.0 - (x2 * x2) / 720.0
    if x > _BIG_X:
        ex = math.exp(-x)
        return (x * ex) / (1.0 - ex)
    return x / (math.expm1(x))


def bernoulli_array(x: np.ndarray) -> np.ndarray:
    """
    Vectorised Bernoulli over a NumPy array of arguments.

    Same regimes as :func:`bernoulli`; uses `np.where` rather than Python
    branching so it composes with vectorised edge assembly without a
    Python-level loop.
    """
    x = np.asarray(x, dtype=np.float64)
    out = np.empty_like(x)

    is_zero = x == 0.0
    ax = np.abs(x)
    is_taylor = (ax < _TAYLOR_EPS) & ~is_zero
    is_big_pos = x > _BIG_X
    is_big_neg = x < -_BIG_X
    is_mid = ~(is_zero | is_taylor | is_big_pos | is_big_neg)

    # Zero branch.
    out[is_zero] = 1.0

    # Taylor branch.
    if is_taylor.any():
        xt = x[is_taylor]
        xt2 = xt * xt
        out[is_taylor] = 1.0 - 0.5 * xt + xt2 / 12.0 - (xt2 * xt2) / 720.0

    # Large positive branch.
    if is_big_pos.any():
        xp = x[is_big_pos]
        ex = np.exp(-xp)
        out[is_big_pos] = (xp * ex) / (1.0 - ex)

    # Large negative branch and mid-range branch share the closed form
    # x / expm1(x); expm1 is accurate near 0 and saturates safely for
    # large |x|.
    if is_big_neg.any():
        xn = x[is_big_neg]
        out[is_big_neg] = xn / np.expm1(xn)
    if is_mid.any():
        xm = x[is_mid]
        out[is_mid] = xm / np.expm1(xm)

    return out


def sg_edge_flux_n(
    n_i: float, n_j: float, dpsi: float,
    h_ij: float, mu_n: float, V_t: float,
) -> float:
    """
    1D Scharfetter-Gummel electron flux on edge (i, j) per ADR 0012.

        F_n_ij = (mu_n V_t / h_ij) * ( n_j B(+dpsi) - n_i B(-dpsi) )

    with `dpsi = psi_j - psi_i` in scaled units (i.e., psi already divided
    by V_t). All quantities in physical SI, except `dpsi` which is
    dimensionless because the SG argument is the scaled potential drop.

    Sign convention: Sandia / Farrell-et-al. F is the electron
    contribution to the conventional current density on edge i->j.
    F > 0 at dpsi=0 with n_j > n_i (diffusion in +grad(n) direction).
    F < 0 at uniform n with dpsi > 0 (electrons drift to j; conventional
    current is opposite). Source: arXiv 1911.00377 Eq (30); Sandia OSTI
    2011-3865 Eq (18) (after a_ij = -dpsi/(2 beta) substitution).
    """
    return (mu_n * V_t / h_ij) * (n_j * bernoulli(dpsi) - n_i * bernoulli(-dpsi))


def sg_edge_flux_p(
    p_i: float, p_j: float, dpsi: float,
    h_ij: float, mu_p: float, V_t: float,
) -> float:
    """
    1D Scharfetter-Gummel hole flux on edge (i, j) per ADR 0012.

        F_p_ij = (mu_p V_t / h_ij) * ( p_i B(+dpsi) - p_j B(-dpsi) )

    Hole formula in the same Sandia / Farrell convention as electrons.
    The i, j positions are swapped relative to the electron formula
    because the device-equation symmetry is `(n <-> p, psi -> -psi)`:
    starting from electrons, applying the symmetry, and re-grouping by
    the (i, j) positions of the swapped flux yields the spelling above.
    Cross-checked against Spevak Section 3.2.6 hole expression after
    sign-mapping for the convention difference (B_S(x) = -B(-x)):
    Spevak's hole `phi_p` equals `-F_p_ij` here (Spevak uses the
    electron-particle-flux convention, opposite of Sandia).

    F_p > 0 at dpsi=0 with p_i > p_j (diffusion in +grad(p) direction
    for the hole conventional current). F_p > 0 at uniform p with
    dpsi > 0 (holes drift from i to j; conventional current matches).

    The test guard 2b verifies the device-equation symmetry under
    `(n <-> p, psi -> -psi)` numerically to 1e-10 relative.
    """
    return (mu_p * V_t / h_ij) * (p_i * bernoulli(dpsi) - p_j * bernoulli(-dpsi))


def sg_edge_flux_n_array(
    n_i: np.ndarray, n_j: np.ndarray, dpsi: np.ndarray,
    h_ij: np.ndarray, mu_n: float, V_t: float,
) -> np.ndarray:
    """Vectorised electron edge flux over an array of edges."""
    return (mu_n * V_t / h_ij) * (
        n_j * bernoulli_array(dpsi) - n_i * bernoulli_array(-dpsi)
    )


def sg_edge_flux_p_array(
    p_i: np.ndarray, p_j: np.ndarray, dpsi: np.ndarray,
    h_ij: np.ndarray, mu_p: float, V_t: float,
) -> np.ndarray:
    """Vectorised hole edge flux over an array of edges."""
    return (mu_p * V_t / h_ij) * (
        p_i * bernoulli_array(dpsi) - p_j * bernoulli_array(-dpsi)
    )


def midpoint_galerkin_flux_n(
    n_i: float, n_j: float, dpsi: float,
    h_ij: float, mu_n: float, V_t: float,
) -> float:
    """
    Independent reference flux for the SG sign+scaling guard test
    (ADR 0012 test plan #2).

        F_ref = mu_n V_t (Delta n / h - n_mid * Delta psi_phys / (V_t h))
              = (mu_n V_t / h) (n_j - n_i)
              + (mu_n / h) * (-n_mid) * V_t * dpsi      [drift, scaled]

    Equivalent to a centred-difference / standard P1 Galerkin
    discretisation of `mu_n n grad psi + D_n grad n` with `n_mid = (n_i + n_j)/2`
    and `D_n = mu_n V_t` (Einstein relation). Sign convention matches
    `sg_edge_flux_n`.

    Used **only** in unit tests as an independent (different-derivation)
    cross-check of the SG flux at small cell Peclet number. The Slotboom
    2-node closed form is **not** used as the reference because it derives
    from the same `B(x)` and only catches typos.
    """
    n_mid = 0.5 * (n_i + n_j)
    return (mu_n * V_t / h_ij) * ((n_j - n_i) - n_mid * dpsi)
