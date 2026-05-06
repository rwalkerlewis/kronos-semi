"""
Net recombination kernels.

Shockley-Read-Hall recombination (always available).

Dimensional form (PHYSICS.md 1.4):

    R_SRH = (n p - n_i^2) / (tau_p (n + n1) + tau_n (p + p1))

with n1 = n_i exp(E_t / V_t), p1 = n_i exp(-E_t / V_t), and E_t the
trap energy measured from the intrinsic level (zero for a mid-gap
trap).

Scaled form. Substituting n = C0 n_hat, p = C0 p_hat, n_i = C0 n_i_hat,
and dividing numerator and denominator by C0:

    R / (C0 / t0)  =  (n_hat p_hat - n_i_hat^2)
                      / ( tau_p_hat (n_hat + n1_hat)
                        + tau_n_hat (p_hat + p1_hat) )

where the scaled lifetimes are tau_hat = tau / t0. The return value is
scaled by the density rate C0/t0, which is what the scaled continuity
equation expects.

The pure-Python helper :func:`srh_rate_np` evaluates the same formula
on numpy arrays and is used by the unit tests. The UFL helper
:func:`srh_rate` consumes UFL expressions (densities built via
Slotboom) and returns a UFL expression suitable for a form builder.

Auger recombination (M16.3).

Dimensional form (PHYSICS.md 1.4, Auger paragraph):

    R_Auger = (C_n n + C_p p) (n p - n_i^2)

Two-carrier band-to-band Auger; same drive (n p - n_i^2) as SRH but a
density-weighted prefactor (C_n n + C_p p) instead of the SRH lifetime
denominator. Units: SI `C_n`, `C_p` in m^6/s; the JSON contract uses
cm^6/s (Si Dziewior-Schmid: C_n = 2.8e-31, C_p = 9.9e-32). At high
injection (n ~ p >> n_i, N) the kernel limits to

    R_Auger -> (C_n + C_p) n^3,

cubic in carrier density and dominant over the bilinear-over-linear
SRH form. Auger is additive to SRH: the form builder evaluates
R = R_SRH + R_Auger when the auger flag is on (see ADR 0001 for the
JSON contract; see ADR 0002 for the scaling convention).

Scaled form. Substituting n = C0 n_hat into the Auger expression and
dividing by C0/t0 (the density rate the scaled continuity equation
expects):

    R / (C0/t0) = C0^2 t0 (C_n n_hat + C_p p_hat) (n_hat p_hat
                                                    - n_i_hat^2)

so the dimensionless Auger coefficient is

    C_hat = C_SI * C0^2 * t0

with C_SI in m^6/s. SI conversion from cm^6/s to m^6/s is 1e-12
(cm^6 = 1e-12 m^6). The pure-Python helper :func:`scaled_auger_C`
performs the dimensional reduction and returns the dimensionless
ratio; runners convert cm^6/s -> m^6/s before calling.
"""
from __future__ import annotations

import math

import numpy as np


def srh_rate(n_hat, p_hat, n_i_hat, tau_n_hat, tau_p_hat, E_t_over_Vt=0.0):
    """
    UFL expression for scaled SRH recombination rate.

    Parameters
    ----------
    n_hat, p_hat : UFL expressions
        Scaled carrier densities (for example built with
        :mod:`semi.physics.slotboom`).
    n_i_hat : UFL expression or float
        Scaled intrinsic density.
    tau_n_hat, tau_p_hat : UFL expression or float
        Scaled electron and hole lifetimes (tau / t0).
    E_t_over_Vt : float
        Trap level referenced to the intrinsic Fermi level, divided by
        V_t. Zero for a mid-gap trap.

    Returns
    -------
    UFL expression
        Scaled recombination rate (per C0/t0).
    """
    import ufl
    n1_hat = n_i_hat * math.exp(E_t_over_Vt)
    p1_hat = n_i_hat * math.exp(-E_t_over_Vt)
    num = n_hat * p_hat - n_i_hat * n_i_hat
    den = tau_p_hat * (n_hat + n1_hat) + tau_n_hat * (p_hat + p1_hat)
    # Explicit UFL call clarifies evaluation when n_i_hat is a fem.Constant
    return ufl.as_ufl(1.0) * num / den


def srh_rate_np(n, p, n_i, tau_n, tau_p, E_t_over_Vt=0.0):
    """
    NumPy counterpart of :func:`srh_rate`.

    Accepts any units (dimensional or scaled) so long as all six
    arguments are consistent. Returns the rate in the matching unit
    system (for example m^-3 s^-1 if the inputs are SI).
    """
    n1 = n_i * np.exp(E_t_over_Vt)
    p1 = n_i * np.exp(-E_t_over_Vt)
    num = n * p - n_i ** 2
    den = tau_p * (n + n1) + tau_n * (p + p1)
    return num / den


def scaled_tau(tau_seconds, t0_seconds):
    """Dimensionless lifetime, tau_hat = tau / t0."""
    return float(tau_seconds) / float(t0_seconds)


# ---------------------------------------------------------------------------
# M16.3 Auger recombination.
# ---------------------------------------------------------------------------


def auger_rate(n_hat, p_hat, n_i_hat, C_n_hat, C_p_hat):
    """
    UFL expression for the scaled Auger recombination rate

        R_hat_Auger = (C_n_hat n_hat + C_p_hat p_hat)
                      * (n_hat p_hat - n_i_hat^2).

    The factor `(n_hat p_hat - n_i_hat^2)` is shared with the SRH
    kernel; an inline assembly in the form builder reuses it via UFL's
    automatic common subexpression elimination (see
    `semi/physics/drift_diffusion.py`).

    Parameters
    ----------
    n_hat, p_hat : UFL expressions
        Scaled carrier densities (Slotboom).
    n_i_hat : UFL expression or float
        Scaled intrinsic density.
    C_n_hat, C_p_hat : UFL expression or float
        Dimensionless Auger coefficients (C_SI * C0^2 * t0); see
        :func:`scaled_auger_C`.

    Returns
    -------
    UFL expression
        Scaled Auger recombination rate (per C0/t0).
    """
    import ufl
    np_minus_nieq = n_hat * p_hat - n_i_hat * n_i_hat
    prefactor = C_n_hat * n_hat + C_p_hat * p_hat
    return ufl.as_ufl(1.0) * prefactor * np_minus_nieq


def auger_rate_np(n, p, n_i, C_n, C_p):
    """
    NumPy counterpart of :func:`auger_rate`.

    Accepts any consistent unit system (dimensional or scaled). With
    SI inputs (n, p, n_i in m^-3 and C_n, C_p in m^6/s), the return
    value is in m^-3 s^-1.
    """
    return (C_n * n + C_p * p) * (n * p - n_i ** 2)


def scaled_auger_C(C_si, C0, t0):
    """
    Convert a dimensional Auger coefficient (m^6/s) to the
    dimensionless ratio consumed by the scaled form:

        C_hat = C_SI * C0^2 * t0

    where C0 has units m^-3 and t0 has units s. Callers reading from
    the JSON contract must first convert cm^6/s to m^6/s (factor 1e-12)
    before invoking this helper; see ADR 0002.
    """
    return float(C_si) * float(C0) ** 2 * float(t0)
