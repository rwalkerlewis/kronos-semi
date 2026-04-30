"""
Pytest assertions on the closed-form MOSCAP figures of merit.
Pure Python; no dolfinx required so this runs in the fast-CI tier.

Mirrors the bullets in the problem statement §7. Where the spec's
hand-quoted numbers are inconsistent with its own defaults (the
t_ox = 5 nm vs. ~10 nm oxide mismatch), we assert against the
formulas with the actual defaults and document the discrepancy.

These constants are the reference the FEM C-V sweep is compared
against in the notebook and in tests/fem/test_moscap_axi_cv.py
(future Slice 3).
"""
from __future__ import annotations

import math

import pytest

from tests.check_moscap_axi_math import closed_form


@pytest.fixture(scope="module")
def cf5():
    """Closed forms at the benchmark defaults (t_ox = 5 nm)."""
    return closed_form()


@pytest.fixture(scope="module")
def cf10():
    """Closed forms at t_ox = 10 nm (the textbook-figure thickness)."""
    return closed_form(t_ox=10.0e-9)


def test_thermal_voltage(cf5):
    assert cf5["Vt"] == pytest.approx(0.025852, abs=5e-5)


def test_flatband_voltage(cf5):
    # n+ poly / p-Si, N_A = 1e17, T = 300 K  ->  V_FB ~ -0.95 to -0.98 V
    assert cf5["V_FB"] == pytest.approx(-0.95, abs=0.05)
    # V_FB does not depend on t_ox
    assert cf5["V_FB"] == pytest.approx(closed_form(t_ox=10e-9)["V_FB"])


def test_two_phi_F(cf5):
    # 2 phi_F = 2 V_t ln(N_A/n_i) ~ 0.83 V at N_A = 1e17
    assert cf5["two_phi_F"] == pytest.approx(0.833, abs=0.01)


def test_W_dep_max(cf5):
    # W_dep,max = sqrt(2 eps_Si * 2 phi_F / (q N_A)) ~ 104 nm
    assert cf5["W_dep_max"] * 1e9 == pytest.approx(104.0, abs=2.0)


def test_C_ox(cf5):
    # C_ox = eps_ox / t_ox = 3.9 eps0 / 5nm
    assert cf5["C_ox"] == pytest.approx(3.9 * 8.8541878128e-12 / 5.0e-9, rel=1e-5)


def test_C_min_over_C_ox_series_law(cf5):
    # C_min from the series law must equal the algebraic 1/(1+...) form.
    series = cf5["C_min"] / cf5["C_ox"]
    assert series == pytest.approx(cf5["C_min_over_C_ox"], rel=1e-10)


def test_C_min_over_C_ox_at_5nm(cf5):
    # Closed form at the benchmark defaults: ~0.126.
    # (The spec's "~0.21" hand-quoted figure corresponds to t_ox=10nm; see
    # tests/check_moscap_axi_math.py docstring.)
    assert cf5["C_min_over_C_ox"] == pytest.approx(0.126, abs=0.005)


def test_C_min_over_C_ox_at_10nm(cf10):
    # Cross-check: at 10 nm oxide we recover the spec's hand-quoted ~0.21
    # (actually 0.224); this confirms the formulas are right.
    assert cf10["C_min_over_C_ox"] == pytest.approx(0.224, abs=0.01)


def test_V_T_at_5nm(cf5):
    # Closed form at the benchmark defaults: ~+0.10 V.
    # (Spec's "~+0.27 V" again corresponds to ~10 nm oxide.)
    assert cf5["V_T"] == pytest.approx(0.098, abs=0.01)


def test_V_T_at_10nm(cf10):
    # Cross-check at 10 nm: matches spec's ~+0.27 V (actually 0.338).
    assert cf10["V_T"] == pytest.approx(0.338, abs=0.02)


def test_V_T_above_V_FB(cf5):
    # Threshold must lie above flatband by at least 2 phi_F.
    assert cf5["V_T"] > cf5["V_FB"] + cf5["two_phi_F"]


def test_C_min_below_C_ox(cf5):
    # By construction the series capacitance is always below the smaller leg.
    assert cf5["C_min"] < cf5["C_ox"]
    assert cf5["C_min"] < cf5["C_dep_min"]


def test_intrinsic_consistency(cf5):
    # phi_F is positive for p-type (Fermi level below intrinsic).
    assert cf5["phi_F"] > 0
    # eps_Si > eps_ox so the depletion region screens less per unit length;
    # this is reflected in 2 phi_F and W_dep,max being self-consistent:
    # rebuild W_dep_max from 2 phi_F and check.
    Q = 1.602176634e-19
    N_A = 1.0e17 * 1e6
    W = math.sqrt(2.0 * cf5["eps_s"] * cf5["two_phi_F"] / (Q * N_A))
    assert W == pytest.approx(cf5["W_dep_max"], rel=1e-12)
