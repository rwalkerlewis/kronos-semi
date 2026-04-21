"""Tests for the Slotboom variable helpers (pure-Python path, no dolfinx)."""
from __future__ import annotations

import numpy as np
import pytest

from semi.physics import slotboom

# The UFL helpers are covered indirectly by the drift-diffusion Form
# builder tests that run in Docker; here we exercise the pure-NumPy
# pointwise helpers that do not require dolfinx to be importable.


@pytest.mark.parametrize(
    "psi_hat, phi_n_hat, phi_p_hat",
    [
        (0.0, 0.0, 0.0),           # intrinsic, equilibrium
        (10.0, 0.0, 0.0),          # heavy n-type, equilibrium
        (-10.0, 0.0, 0.0),         # heavy p-type, equilibrium
        (2.0, 5.0, 5.0),           # split quasi-Fermi potentials under bias
        (-3.0, -1.0, -1.0),
    ],
)
def test_slotboom_round_trip(psi_hat, phi_n_hat, phi_p_hat):
    """Recover (phi_n, phi_p) from (psi, n, p) built by the Slotboom law."""
    n_i_hat = 1.0e-6
    n_hat = slotboom.n_from_slotboom_np(psi_hat, phi_n_hat, n_i_hat)
    p_hat = slotboom.p_from_slotboom_np(psi_hat, phi_p_hat, n_i_hat)
    phi_n_rec = slotboom.phi_n_from_np(psi_hat, n_hat, n_i_hat)
    phi_p_rec = slotboom.phi_p_from_np(psi_hat, p_hat, n_i_hat)
    assert phi_n_rec == pytest.approx(phi_n_hat, rel=0, abs=1e-12)
    assert phi_p_rec == pytest.approx(phi_p_hat, rel=0, abs=1e-12)


def test_slotboom_equilibrium_mass_action():
    """At equilibrium (phi_n = phi_p = 0) the product n p = n_i^2 exactly."""
    n_i_hat = 1.0e-6
    psi_hat = np.linspace(-15.0, 15.0, 25)
    n_hat = slotboom.n_from_slotboom_np(psi_hat, 0.0, n_i_hat)
    p_hat = slotboom.p_from_slotboom_np(psi_hat, 0.0, n_i_hat)
    np.testing.assert_allclose(n_hat * p_hat, n_i_hat ** 2, rtol=1e-12)


def test_slotboom_vectorized_round_trip():
    """Vectorized round trip over an array of potentials."""
    n_i_hat = 1.0e-4
    rng = np.random.default_rng(0)
    psi_hat = rng.uniform(-10.0, 10.0, size=50)
    phi_n_hat = rng.uniform(-2.0, 2.0, size=50)
    phi_p_hat = rng.uniform(-2.0, 2.0, size=50)
    n_hat = slotboom.n_from_slotboom_np(psi_hat, phi_n_hat, n_i_hat)
    p_hat = slotboom.p_from_slotboom_np(psi_hat, phi_p_hat, n_i_hat)
    np.testing.assert_allclose(
        slotboom.phi_n_from_np(psi_hat, n_hat, n_i_hat), phi_n_hat, atol=1e-10,
    )
    np.testing.assert_allclose(
        slotboom.phi_p_from_np(psi_hat, p_hat, n_i_hat), phi_p_hat, atol=1e-10,
    )


def test_equilibrium_psi_hat_limits():
    """Extrinsic p-type and n-type limits: psi -> +- log(|N|/n_i)."""
    n_i = 1.0e16      # m^-3
    N_n = 1.0e22      # heavy n-type (N_D effective)
    N_p = -1.0e22     # heavy p-type
    psi_n = slotboom.equilibrium_psi_hat(N_n, n_i)
    psi_p = slotboom.equilibrium_psi_hat(N_p, n_i)
    # For |N| >> 2 n_i, asinh(x) ~ sign(x) * log(2|x|)
    assert psi_n == pytest.approx(np.log(2.0 * N_n / (2.0 * n_i)), rel=1e-12)
    assert psi_p == pytest.approx(-np.log(2.0 * abs(N_p) / (2.0 * n_i)), rel=1e-12)


def test_equilibrium_psi_hat_intrinsic_zero():
    """Intrinsic sample (N_net = 0) has psi_hat = 0."""
    assert slotboom.equilibrium_psi_hat(0.0, 1.0e16) == 0.0


def test_contact_phi_hat_scales_by_thermal_voltage():
    """phi_hat at an ohmic contact is just V_applied / V_t."""
    V_t = 0.02585
    assert slotboom.contact_phi_hat(0.0, V_t) == 0.0
    assert slotboom.contact_phi_hat(0.6, V_t) == pytest.approx(0.6 / V_t)
    assert slotboom.contact_phi_hat(-0.3, V_t) == pytest.approx(-0.3 / V_t)


def test_phi_n_clips_nonpositive_density():
    """phi_n recovery must not explode on zero density (guarded by floor)."""
    val = slotboom.phi_n_from_np(0.0, 0.0, 1.0e-6)
    assert np.isfinite(val)
