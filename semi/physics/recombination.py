"""
Shockley-Read-Hall recombination kernel.

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
