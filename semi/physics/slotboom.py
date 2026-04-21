"""
Slotboom variable helpers.

Under Boltzmann statistics the carrier densities are

    n = n_i exp((psi - phi_n) / V_t)
    p = n_i exp((phi_p - psi) / V_t)

We solve for the triple (psi, phi_n, phi_p) instead of (psi, n, p) so the
continuity equations become coercive on the quasi-Fermi potentials (see
ADR 0004). Inside `semi/physics/` everything is in scaled units where
psi, phi_n, phi_p are divided by V_t and densities are divided by C_0;
`n_i_hat = n_i / C_0` is the only material-specific constant that shows
up here.

UFL expressions in this module work with either dolfinx Functions or
fem.Constants for `psi`, `phi_n`, `phi_p`. The NumPy helpers are pure
Python and are imported for tests and post-processing. This module must
remain independent of the FE space construction; it builds expressions,
not Forms.
"""
from __future__ import annotations

import numpy as np


def n_from_slotboom(psi_hat, phi_n_hat, n_i_hat):
    """
    Scaled electron density, n_hat = n_i_hat * exp(psi_hat - phi_n_hat).

    All arguments are UFL-compatible (Functions, Constants, expressions)
    or plain floats / numpy arrays. Returns an expression of the same
    flavor.
    """
    import ufl
    return n_i_hat * ufl.exp(psi_hat - phi_n_hat)


def p_from_slotboom(psi_hat, phi_p_hat, n_i_hat):
    """Scaled hole density, p_hat = n_i_hat * exp(phi_p_hat - psi_hat)."""
    import ufl
    return n_i_hat * ufl.exp(phi_p_hat - psi_hat)


def n_from_slotboom_np(psi_hat, phi_n_hat, n_i_hat):
    """NumPy counterpart of :func:`n_from_slotboom` for post-processing."""
    return n_i_hat * np.exp(psi_hat - phi_n_hat)


def p_from_slotboom_np(psi_hat, phi_p_hat, n_i_hat):
    """NumPy counterpart of :func:`p_from_slotboom` for post-processing."""
    return n_i_hat * np.exp(phi_p_hat - psi_hat)


def phi_n_from_np(psi_hat, n_hat, n_i_hat):
    """
    Recover the electron quasi-Fermi potential from (psi, n).

        phi_n = psi - log(n / n_i)  (all scaled by V_t and C_0).

    Guards against non-positive densities by flooring at a tiny value.
    """
    n_safe = np.maximum(n_hat, 1.0e-300)
    return psi_hat - np.log(n_safe / n_i_hat)


def phi_p_from_np(psi_hat, p_hat, n_i_hat):
    """Recover the hole quasi-Fermi potential from (psi, p)."""
    p_safe = np.maximum(p_hat, 1.0e-300)
    return psi_hat + np.log(p_safe / n_i_hat)


def equilibrium_psi_hat(N_net, n_i):
    """
    Scaled equilibrium potential for local charge neutrality, psi/V_t.

        psi_eq = V_t * asinh(N_net / (2 n_i))

    Works on scalars or numpy arrays. `N_net` and `n_i` must be in the
    same (physical) units; the returned quantity is dimensionless
    (already scaled by V_t).
    """
    return np.arcsinh(N_net / (2.0 * n_i))


def contact_phi_hat(V_applied_volts, V_t_volts):
    """
    Scaled quasi-Fermi potential at an ohmic contact under applied bias.

    Both phi_n and phi_p equal the applied contact bias (PHYSICS.md 3.1).
    In scaled units this is simply V_applied / V_t.
    """
    return float(V_applied_volts) / float(V_t_volts)
