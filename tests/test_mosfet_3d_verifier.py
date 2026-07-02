"""M19: pure-Python assertions on the 3D MOSFET analytical helpers.

Exercises `semi.diode_analytical.mosfet_3d_paosah_iv` (linear regime)
and `mosfet_3d_saturation_iv` (velocity-saturation regime). No dolfinx;
runs in the pure-Python CI suite alongside the FEM smoke test in
`tests/fem/test_mosfet_3d.py`.
"""
from __future__ import annotations

import numpy as np

from semi.diode_analytical import mosfet_3d_paosah_iv, mosfet_3d_saturation_iv

# Representative device operating point (matches the mosfet_3d benchmark:
# N_A = 1e16 cm^-3 body, 5 nm SiO2, L = 250 nm, W = 1 um).
_V_T = 0.4277
_C_OX = 0.006906266   # F/m^2 (eps_SiO2 / 5 nm)
_L = 250.0e-9
_W = 1.0e-6
_MU = 0.14            # m^2/(V s)  (1400 cm^2/Vs)
_VSAT = 1.0e5         # m/s        (1e7 cm/s)


def test_threshold_zero_below_and_at_vt():
    """Both currents are exactly zero at and below threshold."""
    V_GS = np.array([-0.5, 0.0, _V_T - 0.05, _V_T])
    I_lin = mosfet_3d_paosah_iv(V_GS, 0.05, _MU, _C_OX, _L, _W, _V_T, _VSAT)
    I_sat = mosfet_3d_saturation_iv(V_GS, _MU, _C_OX, _L, _W, _V_T, _VSAT)
    assert np.all(I_lin == 0.0)
    assert np.all(I_sat == 0.0)


def test_linear_regime_matches_closed_form():
    """Linear-regime current equals the Pao-Sah formula with the
    velocity-saturation mobility correction, to machine precision."""
    V_GS = 1.0
    V_DS = 0.05
    I = float(mosfet_3d_paosah_iv(V_GS, V_DS, _MU, _C_OX, _L, _W, _V_T, _VSAT))
    mu_corr = _MU / (1.0 + _MU * V_DS / (_VSAT * _L))
    expected = (_W / _L) * mu_corr * _C_OX * (V_GS - _V_T) * V_DS
    assert I == expected
    assert I > 0.0
    # The velocity-saturation correction must lower the current relative
    # to the uncorrected square law.
    uncorrected = (_W / _L) * _MU * _C_OX * (V_GS - _V_T) * V_DS
    assert I < uncorrected


def test_saturation_regime_matches_closed_form():
    """Saturation current equals the velocity-saturation I_DSAT formula."""
    V_GS = 1.5
    I = float(mosfet_3d_saturation_iv(V_GS, _MU, _C_OX, _L, _W, _V_T, _VSAT))
    V_ov = V_GS - _V_T
    denom = 1.0 + V_ov / (2.0 * _VSAT * _L / _MU)
    expected = (_W / (2.0 * _L)) * _MU * _C_OX * V_ov ** 2 / denom
    assert I == expected
    assert I > 0.0
    # Saturation current at fixed overdrive exceeds the small-V_DS linear
    # current (V_DS = 0.05 V), a basic sanity ordering.
    I_lin = float(mosfet_3d_paosah_iv(V_GS, 0.05, _MU, _C_OX, _L, _W, _V_T, _VSAT))
    assert I > I_lin


def test_monotonic_above_threshold():
    """Both currents increase strictly monotonically with V_GS above V_T."""
    V_GS = np.linspace(_V_T + 0.05, 2.0, 25)
    I_lin = mosfet_3d_paosah_iv(V_GS, 0.05, _MU, _C_OX, _L, _W, _V_T, _VSAT)
    I_sat = mosfet_3d_saturation_iv(V_GS, _MU, _C_OX, _L, _W, _V_T, _VSAT)
    assert np.all(np.diff(I_lin) > 0.0)
    assert np.all(np.diff(I_sat) > 0.0)
