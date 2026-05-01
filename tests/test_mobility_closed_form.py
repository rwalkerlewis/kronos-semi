"""
Pure-Python unit tests for the Caughey-Thomas closed-form mobility
math (M16.1).

Covers `semi/physics/mobility.py::caughey_thomas_mu` against Python
scalars; the FEM wiring is exercised by the MMS variant
(tests/fem/test_mms_caughey_thomas.py) and the diode_velsat_1d
benchmark.

Tests
-----
1. Low-field limit: mu(0) = mu0 to within 1e-12.
2. Saturation limit: mu(F_large) * F_large -> vsat to within 1%.
3. beta = 2 at the inflection F = vsat / mu0: mu = mu0 / sqrt(2).
4. beta = 1 at the inflection F = vsat / mu0: mu = mu0 / 2.
5. Three sample points for the Si electron defaults
   (mu0 = 1400 cm^2/Vs, vsat = 1e7 cm/s, beta = 2) at
   F in {1e3, 1e4, 1e5} V/cm match a hand calculation.

All tests run without dolfinx (the module's UFL-bearing functions are
imported via `_lazy_import`-style local imports inside the FEM
builder, not at module scope).
"""
from __future__ import annotations

import math

import pytest

from semi.physics.mobility import caughey_thomas_mu, constant_mu


def _hand_caughey_thomas(mu0: float, F: float, vsat: float, beta: float) -> float:
    """Reference closed-form, computed in Python doubles."""
    if F == 0.0:
        return mu0
    return mu0 / (1.0 + (mu0 * F / vsat) ** beta) ** (1.0 / beta)


def test_constant_mu_is_identity():
    assert constant_mu(1.5) == 1.5
    obj = object()
    assert constant_mu(obj) is obj


def test_low_field_limit_returns_mu0():
    """F_par -> 0: mu(F) = mu0 to within 1e-12 (machine-floor cancel)."""
    mu0 = 1400.0
    vsat = 1.0e7
    for beta in (1.0, 1.5, 2.0):
        mu = caughey_thomas_mu(mu0, 0.0, vsat, beta)
        assert mu == pytest.approx(mu0, rel=1.0e-12, abs=0.0), (
            f"low-field limit failed at beta={beta}: mu={mu}, mu0={mu0}"
        )


def test_high_field_drift_velocity_saturates():
    """F_par -> infinity: mu(F) * F -> vsat to within 1%."""
    mu0 = 1400.0          # cm^2/(V s)
    vsat = 1.0e7          # cm/s
    beta = 2.0
    F_huge = 1.0e10       # V/cm; ensures mu0*F/vsat = 1.4e6 >> 1

    mu = caughey_thomas_mu(mu0, F_huge, vsat, beta)
    v = mu * F_huge

    # Asymptote: (mu0 * F / vsat)^beta dominates 1, so
    # mu ~ mu0 / (mu0 * F / vsat) = vsat / F.
    rel_err = abs(v - vsat) / vsat
    assert rel_err < 0.01, (
        f"saturation limit failed: v={v} cm/s, vsat={vsat} cm/s, "
        f"rel_err={rel_err}"
    )


def test_beta_two_at_inflection_F_is_mu0_over_sqrt2():
    """At F = vsat / mu0 with beta = 2: closed-form mu = mu0 / sqrt(2)."""
    mu0 = 1400.0
    vsat = 1.0e7
    F = vsat / mu0
    mu = caughey_thomas_mu(mu0, F, vsat, 2.0)
    expected = mu0 / math.sqrt(2.0)
    assert mu == pytest.approx(expected, rel=1.0e-12, abs=0.0)


def test_beta_one_at_inflection_F_is_mu0_over_two():
    """At F = vsat / mu0 with beta = 1: closed-form mu = mu0 / 2."""
    mu0 = 450.0
    vsat = 8.0e6
    F = vsat / mu0
    mu = caughey_thomas_mu(mu0, F, vsat, 1.0)
    expected = mu0 / 2.0
    assert mu == pytest.approx(expected, rel=1.0e-12, abs=0.0)


@pytest.mark.parametrize("F_V_per_cm", [1.0e3, 1.0e4, 1.0e5])
def test_si_electron_sample_points_match_hand_calc(F_V_per_cm: float):
    """For Si electron defaults (mu0=1400, vsat=1e7, beta=2), three
    intermediate field samples must reproduce a hand calculation.

    The hand reference is computed in this same module; the test
    pins the closed-form algebra against a literal restatement so
    drift in the implementation surfaces immediately.
    """
    mu0 = 1400.0
    vsat = 1.0e7
    beta = 2.0
    mu = caughey_thomas_mu(mu0, F_V_per_cm, vsat, beta)
    ref = _hand_caughey_thomas(mu0, F_V_per_cm, vsat, beta)
    assert mu == pytest.approx(ref, rel=1.0e-14, abs=0.0)
    assert mu <= mu0 + 1.0e-12
    assert mu > 0.0
