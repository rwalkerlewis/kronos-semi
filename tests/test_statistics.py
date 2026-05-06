"""
Pure-Python tests for the M16.4 Fermi-Dirac statistics module.

Covers:
  - Blakemore-vs-full-integral accuracy across the eta range the
    production benchmark exercises.
  - gamma_n / gamma_p limits (non-degenerate -> 1, degenerate -> 0).
  - Closed-form Einstein-factor cancellation with the basic Blakemore
    form (the algebraic identity behind ADR 0004's preservation under
    FD).
  - Default-fill bit-identity: the Slotboom NumPy helpers with
    statistics_cfg=None produce results bit-identical to the bare
    Boltzmann form.
  - Scaling.eta_offset_n / eta_offset_p surface a clear error when
    N_C / N_V are unset (Boltzmann-only solves stay green).
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from semi.physics.slotboom import (
    n_from_slotboom_np,
    p_from_slotboom_np,
)
from semi.physics.statistics import (
    _BLAKEMORE_BASIC_OFFSET,
    _eta_offset_for_material,
    einstein_factor_blakemore,
    fermi_dirac_half_blakemore,
    fermi_dirac_half_reference,
    gamma_n_blakemore,
    gamma_p_blakemore,
)
from semi.scaling import Scaling

# ---------------------------------------------------------------------------
# Blakemore vs full-integral reference accuracy.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("eta", [-5.0, -3.0, -1.0, 0.0, 1.0, 1.25])
def test_blakemore_within_5_percent_of_reference(eta):
    """The M16.4 acceptance benchmark operates at eta ~ 1.25 (n+
    region of the diode_fermi_dirac_1d benchmark with N_D = 1e20
    cm^-3, Si N_C = 2.86e19). Across the non-degenerate-to-onset
    range eta in [-5, +1.25] the basic Blakemore form is < 5 %
    accurate vs the full integral; the M16.4 acceptance gate uses
    the analytical scipy reference, not the production residual,
    so this 5 % envelope is documentation of the production-form
    accuracy rather than a gate on the benchmark."""
    approx = float(fermi_dirac_half_blakemore(eta))
    ref = float(fermi_dirac_half_reference(eta))
    rel_err = abs(approx - ref) / ref
    assert rel_err < 0.05, (
        f"Blakemore deviates {rel_err:.3f} from reference at eta={eta}: "
        f"approx={approx:.4e}, ref={ref:.4e}"
    )


def test_blakemore_at_benchmark_operating_point():
    """Documents the basic Blakemore form's accuracy at the actual
    M16.4 acceptance benchmark operating point eta ~ 1.25 (~3 % off
    the full integral). Protects against regressions in the closed
    form (e.g. a stray sign or a missing factor would push this
    well outside 5 %)."""
    eta_bench = 1.25
    approx = float(fermi_dirac_half_blakemore(eta_bench))
    ref = float(fermi_dirac_half_reference(eta_bench))
    rel_err = abs(approx - ref) / ref
    assert rel_err < 0.05


def test_blakemore_envelope_grows_in_deep_degenerate():
    """Documents the basic-form failure mode for future readers: in
    the deep-degenerate regime (eta = 5) the basic Blakemore form is
    materially off (~50 % low), which is why the M16.4 benchmark stays
    at eta ~ 1.25 and the acceptance gate is met by the analytical
    scipy reference, not by the production residual extrapolated."""
    eta_deep = 5.0
    approx = float(fermi_dirac_half_blakemore(eta_deep))
    ref = float(fermi_dirac_half_reference(eta_deep))
    rel_err = abs(approx - ref) / ref
    assert rel_err > 0.30, (
        "If this test fails the basic Blakemore form has been "
        "replaced; update the production-vs-reference accuracy "
        "documentation in semi/physics/statistics.py."
    )


def test_blakemore_array_input():
    eta = np.linspace(-5.0, 5.0, 11)
    approx = fermi_dirac_half_blakemore(eta)
    assert approx.shape == eta.shape
    assert np.all(approx > 0.0)
    # Monotonic in eta.
    assert np.all(np.diff(approx) > 0.0)


def test_reference_array_input():
    eta = np.linspace(-3.0, 3.0, 7)
    ref = fermi_dirac_half_reference(eta)
    assert ref.shape == eta.shape
    assert np.all(ref > 0.0)
    assert np.all(np.diff(ref) > 0.0)


# ---------------------------------------------------------------------------
# gamma_n / gamma_p limits.
# ---------------------------------------------------------------------------


def test_gamma_n_non_degenerate_limit_is_one():
    """As (psi - phi_n)/V_t -> -inf the Boltzmann limit is recovered
    and gamma_n -> 1 within 1e-4."""
    eta_offset_n = -22.0  # Si-like
    drive = -10.0  # well below the FD onset for Si parameters
    g = float(gamma_n_blakemore(drive, eta_offset_n))
    assert abs(g - 1.0) < 1.0e-4


def test_gamma_n_degenerate_value_below_one():
    """At eta = 5 (deeply degenerate), gamma_n is materially < 1.
    The closed Blakemore form gives 1 / (1 + 0.27 * exp(5)) ~ 0.024."""
    eta_offset_n = 0.0
    drive = 5.0
    g = float(gamma_n_blakemore(drive, eta_offset_n))
    expected = 1.0 / (1.0 + _BLAKEMORE_BASIC_OFFSET * math.exp(5.0))
    assert g == pytest.approx(expected, rel=1.0e-12)
    assert g < 0.05


def test_gamma_p_mirrors_gamma_n():
    """Hole-side prefactor uses identical closed form once eta_p has
    been formed; verify with the same reduced argument."""
    eta_offset = -5.0
    drive_n = 2.0
    drive_p = 2.0
    gn = float(gamma_n_blakemore(drive_n, eta_offset))
    gp = float(gamma_p_blakemore(drive_p, eta_offset))
    assert gn == pytest.approx(gp, rel=1.0e-12)


# ---------------------------------------------------------------------------
# Einstein-factor cancellation (ADR 0004 invariance under FD).
# ---------------------------------------------------------------------------


def test_einstein_factor_is_inverse_of_gamma():
    """The closed Blakemore form gives the algebraic identity
    g(eta) * gamma_blakemore(eta) = 1, which is exactly the
    cancellation that keeps the generalized-Slotboom current
    `J = -q mu n grad(phi)` independent of the FD correction."""
    for eta in (-3.0, -1.0, 0.0, 1.0, 3.0, 5.0):
        g = float(einstein_factor_blakemore(eta))
        gamma = 1.0 / (1.0 + _BLAKEMORE_BASIC_OFFSET * math.exp(eta))
        assert g * gamma == pytest.approx(1.0, rel=1.0e-12)


def test_fd_density_consistent_two_ways():
    """At five (psi, phi_n) Slotboom samples, the FD electron density
    built via `n_from_slotboom_np` must agree with the direct FD form
    `N_C * F_{1/2}(eta)` evaluated through the same Blakemore basic
    closed form (this is an algebraic-identity test, not an accuracy
    test against the full integral)."""
    n_i = 1.0e10
    N_C = 2.86e19
    n_i_hat = 1.0  # work in Slotboom-scaled units, n_i = C0
    eta_offset_n = math.log(n_i / N_C)
    cfg = {"statistics": "fermi_dirac"}
    pairs = [(-3.0, 0.0), (-1.0, 0.0), (0.0, 0.0), (2.0, 0.0), (4.0, 0.0)]
    for psi, phi_n in pairs:
        # Path 1: production residual (Blakemore-based generalized Slotboom).
        n_residual = float(n_from_slotboom_np(
            psi, phi_n, n_i_hat,
            statistics_cfg=cfg, eta_offset_n=eta_offset_n,
        ))
        # Path 2: direct definition n / n_i = (N_C / n_i) F_{1/2}(eta).
        eta = (psi - phi_n) + eta_offset_n
        n_direct = (N_C / n_i) * float(fermi_dirac_half_blakemore(eta))
        assert n_residual == pytest.approx(n_direct, rel=1.0e-12), (
            f"FD substitution disagrees with direct Blakemore at "
            f"(psi={psi}, phi_n={phi_n}): residual={n_residual:.6e}, "
            f"direct={n_direct:.6e}"
        )


# ---------------------------------------------------------------------------
# Default-fill bit-identity.
# ---------------------------------------------------------------------------


def test_n_from_slotboom_np_default_fill_is_boltzmann():
    """statistics_cfg=None must reproduce the bare Boltzmann result
    bit-identically (np.array_equal, not np.allclose)."""
    n_i_hat = 1.0
    samples = [(-1.0, 0.0), (0.0, 0.0), (0.5, 0.2), (1.0, 0.0), (3.0, 1.0)]
    for psi, phi_n in samples:
        boltzmann = n_i_hat * np.exp(psi - phi_n)
        result = n_from_slotboom_np(psi, phi_n, n_i_hat)
        assert np.array_equal(np.asarray(boltzmann), np.asarray(result))
        result_explicit = n_from_slotboom_np(
            psi, phi_n, n_i_hat,
            statistics_cfg={"statistics": "boltzmann"},
        )
        assert np.array_equal(np.asarray(boltzmann), np.asarray(result_explicit))


def test_p_from_slotboom_np_default_fill_is_boltzmann():
    n_i_hat = 1.0
    samples = [(-1.0, 0.0), (0.0, 0.0), (-0.5, 0.2), (-1.0, 0.0)]
    for psi, phi_p in samples:
        boltzmann = n_i_hat * np.exp(phi_p - psi)
        result = p_from_slotboom_np(psi, phi_p, n_i_hat)
        assert np.array_equal(np.asarray(boltzmann), np.asarray(result))


def test_n_from_slotboom_np_fd_requires_eta_offset():
    """The FD path raises a clear ValueError when eta_offset_n is
    missing rather than silently falling back to the Boltzmann form."""
    with pytest.raises(ValueError, match="eta_offset_n"):
        n_from_slotboom_np(0.0, 0.0, 1.0,
                           statistics_cfg={"statistics": "fermi_dirac"})


def test_p_from_slotboom_np_fd_requires_eta_offset():
    with pytest.raises(ValueError, match="eta_offset_p"):
        p_from_slotboom_np(0.0, 0.0, 1.0,
                           statistics_cfg={"statistics": "fermi_dirac"})


# ---------------------------------------------------------------------------
# Scaling.eta_offset_n / eta_offset_p.
# ---------------------------------------------------------------------------


def test_scaling_eta_offset_si_values():
    """Si reference: eta_offset_n = ln(1e16 / 2.86e25) ~ -21.78,
    eta_offset_p = ln(1e16 / 3.10e25) ~ -21.86."""
    sc = Scaling(
        L0=1.0e-6, C0=1.0e22, T=300.0, mu0=0.135, n_i=1.0e16,
        N_C=2.86e25, N_V=3.10e25,
    )
    assert sc.eta_offset_n == pytest.approx(math.log(1.0e16 / 2.86e25), rel=1.0e-12)
    assert sc.eta_offset_p == pytest.approx(math.log(1.0e16 / 3.10e25), rel=1.0e-12)


def test_scaling_eta_offset_raises_when_unset():
    """Pre-M16.4 callers (and Boltzmann-only tests) instantiate
    Scaling without N_C / N_V; the FD-only properties must surface a
    clear error rather than a silent log(None)."""
    sc = Scaling(L0=1.0e-6, C0=1.0e22, T=300.0, mu0=0.135, n_i=1.0e16)
    with pytest.raises(ValueError, match="N_C"):
        _ = sc.eta_offset_n
    with pytest.raises(ValueError, match="N_V"):
        _ = sc.eta_offset_p


def test_eta_offset_for_material_helper():
    """The pure-function reduction matches log(n_i / N_C)."""
    val = _eta_offset_for_material(2.86e25, 1.0e16)
    assert val == pytest.approx(math.log(1.0e16 / 2.86e25), rel=1.0e-12)
