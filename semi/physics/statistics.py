"""
Carrier statistics dispatch (M16.4).

Boltzmann is the M16.x default and the V&V reference. Fermi-Dirac is
gated by `physics.statistics: "fermi_dirac"` in the JSON contract and
substitutes the Blakemore approximation for the order-1/2 Fermi-Dirac
integral inside the generalized-Slotboom helpers (ADR 0004) and the
Poisson source.

Production approximation (Blakemore 1982 basic form):

    F_{1/2}(eta) ~ 1 / (exp(-eta) + 0.27)

Accuracy: this is the "good in the non-degenerate-to-onset" form, not
the "good everywhere" form. Numerical envelope vs the full integral:
< 1 % for `eta < -1`, ~3 % at `eta = 0`, ~1 % at `eta = 1`, ~3 % at
`eta = 1.25`, ~13 % at `eta = 2`, growing further into the deep-
degenerate regime. The M16.4 acceptance benchmark sits at `eta ~ 1.25`
(the n+ region with `N_D = 1e20 cm^-3` against Si `N_C = 2.86e19
cm^-3`), inside the sub-5 % regime. The improved Blakemore form (with
the `eta`-dependent prefactor `zeta(eta)` from the same 1982 paper) is
good to < 1 % across `[-inf, +inf]`; the basic form is the production
choice because the closed-form Slotboom substitution simplifies
cleanly and the M16.4 benchmark range stays inside the accurate
window. The acceptance gate that matters is the V_bi divergence
(>15 % between Boltzmann and FD) and the analytical-match comparison
against the scipy-driven full integral, both of which use the
reference helper (:func:`fermi_dirac_half_reference`), not the
production residual.

Verification reference uses the full integral

    F_{1/2}(eta) = -Li_{3/2}(-exp(eta))

via `mpmath.polylog(1.5, -exp(eta))` (the standard polylogarithm
identity for the Fermi-Dirac integral). Older scipy releases shipped a
`scipy.special.fdk` helper that took the same argument; it is not
available on the docker-fem image, so the production fallback is
mpmath. The reference is only called from the unit tests and from the
benchmark verifier, never from the production residual.

Generalized-Slotboom substitution. ADR 0004 keeps Slotboom primary
unknowns under FD: rather than rewriting the continuity rows, the
substitution rule for `n` and `p` in terms of `(psi, phi_n, phi_p)`
gains a smooth bounded prefactor:

    n = n_i * gamma_n(eta_n) * exp(psi - phi_n)         (scaled)
    p = n_i * gamma_p(eta_p) * exp(phi_p - psi)         (scaled)

with reduced Fermi levels

    eta_n = (psi - phi_n) / V_t + ln(n_i / N_C)
    eta_p = (phi_p - psi) / V_t + ln(n_i / N_V)

and prefactors

    gamma_n(eta) = F_{1/2}(eta) / exp(eta)
    gamma_p(eta) = F_{1/2}(eta) / exp(eta)

so that `n = N_C F_{1/2}(eta_n)` and `p = N_V F_{1/2}(eta_p)` (the
direct FD definitions) reduce to the textbook Slotboom result in the
non-degenerate limit (gamma -> 1). Inserting the Blakemore basic form
into the prefactor gives the closed expression

    gamma_blakemore(eta) = 1 / (1 + 0.27 * exp(eta))

which is what the production residual evaluates. The Einstein factor
`g(eta) = F_{1/2}(eta) / F_{-1/2}(eta)` does not appear in the
Slotboom-form current `J = -q mu n grad(phi)` because the Einstein
correction cancels against the FD prefactor in the substitution; this
is exactly why generalized Slotboom is the standard production-FD
path. See Schenk 1998 and the Sentaurus device manual.

This module is pure-Python (NumPy only) and ships the Blakemore
approximation, the full-integral reference, and the gamma_n / gamma_p
prefactors. The UFL counterparts of gamma_n and gamma_p live in
`semi.physics.slotboom` (Layer 4); imports are deferred to the function
bodies there. Keeping the closed-form math here means the Layer-3
invariant (no dolfinx in pure-Python core) holds.
"""
from __future__ import annotations

import math

import numpy as np

_BLAKEMORE_BASIC_OFFSET = 0.27


def fermi_dirac_half_blakemore(eta):
    """
    Basic Blakemore approximation for the order-1/2 Fermi-Dirac integral.

    F_{1/2}(eta) ~ 1 / (exp(-eta) + 0.27)

    Pure NumPy. Accepts a Python scalar, NumPy scalar, or NumPy array;
    returns the same shape. The basic form is what the production
    residual evaluates (the simpler-and-good-enough form per the M16.4
    Blakemore citation; the M16.4 acceptance benchmark sits at
    eta ~ 1.25 where the basic form is sub-1 % accurate. See the
    module docstring for the full accuracy envelope vs eta).
    """
    eta_arr = np.asarray(eta, dtype=float)
    return 1.0 / (np.exp(-eta_arr) + _BLAKEMORE_BASIC_OFFSET)


def fermi_dirac_half_reference(eta):
    """
    Full-integral reference for F_{1/2}(eta).

    Uses the polylogarithm identity F_{1/2}(eta) = -Li_{3/2}(-exp(eta))
    via :func:`mpmath.polylog`. Accepts a Python scalar, NumPy scalar,
    or NumPy array; returns the same shape (always real-valued because
    the imaginary part of the polylog vanishes for negative real
    argument).

    Used by :func:`tests.test_statistics` and the
    `diode_fermi_dirac_1d` benchmark verifier; never called from the
    production residual. Older scipy shipped `scipy.special.fdk` which
    computes the same quantity, but that function is not on the
    docker-fem image; mpmath is the supported fallback.
    """
    import mpmath

    eta_arr = np.atleast_1d(np.asarray(eta, dtype=float))
    out = np.empty_like(eta_arr)
    for i, e in enumerate(eta_arr.flat):
        z = -mpmath.exp(float(e))
        # Li_{3/2}(z) is real for real z <= 0; coerce to float.
        out.flat[i] = float(mpmath.re(-mpmath.polylog(1.5, z)))
    if np.isscalar(eta) or np.asarray(eta).ndim == 0:
        return float(out[0])
    return out.reshape(np.asarray(eta).shape)


def gamma_n_blakemore(psi_minus_phi_n_over_Vt, eta_offset):
    """
    FD-correction prefactor for the generalized-Slotboom electron form

        n = n_i * gamma_n * exp(psi - phi_n)         (scaled)
        gamma_n(eta) = F_{1/2}(eta) / exp(eta)
                     = 1 / (1 + 0.27 * exp(eta))    (Blakemore basic)

    where `eta = (psi - phi_n)/V_t + ln(n_i / N_C) = eta_n` in the
    standard convention. The first argument is the Slotboom drive
    (positive in n-type / forward-active electron-rich regions). The
    second argument absorbs the band-edge offset `(E_i - E_C)/kT` in
    the intrinsic-Fermi convention; it is a per-material constant the
    caller computes once at scaling-build time via
    :func:`_eta_offset_for_material`.

    Pure NumPy. As `(psi - phi_n)/V_t -> -inf` (deep p-region or
    negligible electron density), `eta -> -inf`, the exponential dies,
    and gamma_n -> 1 (Boltzmann limit). As `eta -> +inf` (degenerate
    n+ region), gamma_n decays to zero so the FD density remains below
    the Boltzmann prediction (Pauli blocking).
    """
    eta = np.asarray(psi_minus_phi_n_over_Vt, dtype=float) + float(eta_offset)
    return 1.0 / (1.0 + _BLAKEMORE_BASIC_OFFSET * np.exp(eta))


def gamma_p_blakemore(phi_p_minus_psi_over_Vt, eta_offset):
    """
    Hole counterpart of :func:`gamma_n_blakemore`.

    For holes the Slotboom drive is `(phi_p - psi)/V_t` and the
    corresponding reduced Fermi level is

        eta_p = (phi_p - psi)/V_t + ln(n_i / N_V).

    The closed Blakemore form is identical to the electron case once
    `eta_p` has been formed, so the implementation reuses the same
    scalar reduction.
    """
    eta = np.asarray(phi_p_minus_psi_over_Vt, dtype=float) + float(eta_offset)
    return 1.0 / (1.0 + _BLAKEMORE_BASIC_OFFSET * np.exp(eta))


def einstein_factor_blakemore(eta):
    """
    g(eta) = F_{1/2}(eta) / F_{-1/2}(eta), the FD Einstein correction.

    Used only by the unit tests to numerically witness the Einstein-
    factor cancellation in the generalized-Slotboom current expression
    (see ADR 0004 derivation note in `docs/PHYSICS.md` 1.3). The
    cancellation is exact in closed form when both F_{1/2} and F_{-1/2}
    are expressed through the basic Blakemore identity

        F_{-1/2}(eta) = d/d eta F_{1/2}(eta).

    Differentiating the basic Blakemore form yields

        F_{-1/2}(eta) ~ exp(-eta) / (exp(-eta) + 0.27)^2

    so g(eta) reduces to (exp(-eta) + 0.27) / exp(-eta)
    = 1 + 0.27 * exp(eta), which is exactly 1 / gamma_n_blakemore(eta).
    The unit tests verify this identity numerically.
    """
    eta_arr = np.asarray(eta, dtype=float)
    return 1.0 + _BLAKEMORE_BASIC_OFFSET * np.exp(eta_arr)


def _eta_offset_for_material(N_C, n_i):
    """
    Per-material reduced-Fermi-level offset.

    Under Boltzmann the intrinsic level satisfies
    `n_i = N_C exp((E_i - E_C)/kT)`, so the offset that turns the
    Slotboom drive `(psi - phi_n)/V_t` into the absolute reduced Fermi
    level `eta_n = (E_F_n - E_C)/kT` is `ln(n_i / N_C)`. Same shape
    for holes with N_C replaced by N_V.

    Pure-Python; the caller (typically `semi.scaling.Scaling`) stores
    the result so the closed form does not run on every form-builder
    invocation.
    """
    return math.log(float(n_i) / float(N_C))
