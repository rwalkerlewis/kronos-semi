"""
Pure-Python unit tests for the Caughey-Thomas (M16.1) and Lombardi
(M16.2) closed-form mobility math.

Covers `semi/physics/mobility.py::caughey_thomas_mu`, `lombardi_mu_AC`,
`lombardi_mu_sr`, `lombardi_compose`, and `lombardi_unit_conversions`
against Python scalars. The FEM wiring is exercised by the MMS variants
(tests/fem/test_mms_caughey_thomas.py for M16.1; Phase D adds
tests/fem/test_mms_lombardi.py for M16.2) and the
diode_velsat_1d / mosfet_2d benchmarks.

Tests (M16.1, Caughey-Thomas)
-----------------------------
1. Low-field limit: mu(0) = mu0 to within 1e-12.
2. Saturation limit: mu(F_large) * F_large -> vsat to within 1%.
3. beta = 2 at the inflection F = vsat / mu0: mu = mu0 / sqrt(2).
4. beta = 1 at the inflection F = vsat / mu0: mu = mu0 / 2.
5. Three sample points for the Si electron defaults
   (mu0 = 1400 cm^2/Vs, vsat = 1e7 cm/s, beta = 2) at
   F in {1e3, 1e4, 1e5} V/cm match a hand calculation.

Tests (M16.2, Lombardi)
-----------------------
6. Low-doping low-field limit of `lombardi_mu_AC`: B / E_perp
   dominates; the C-term contribution is small.
7. `lombardi_mu_sr` returns a numerically large value as
   E_perp -> 0 (regularized via the same eps as caughey_thomas).
8. `lombardi_compose` reduces to mu_bulk when mu_AC, mu_sr -> +inf.
9. `lombardi_compose` reduces to the parallel combination of
   mu_AC and mu_sr when mu_bulk -> +inf.
10. The Si-electron sample point at N_total = 1e17 cm^-3,
    E_perp = 1e6 V/m, T = 300 K matches a hand-computed restatement
    of the Lombardi closed form.
11. `lombardi_unit_conversions` round-trips a known SI mu_AC at the
    sample point to a dimensionless mu_AC_hat consistent with the
    SI calculation.

All tests run without dolfinx (the module's UFL-bearing functions are
imported via local imports inside the FEM builder, not at module
scope).
"""
from __future__ import annotations

import math

import pytest

from semi.physics.mobility import (
    caughey_thomas_mu,
    constant_mu,
    lombardi_compose,
    lombardi_mu_AC,
    lombardi_mu_sr,
    lombardi_unit_conversions,
)


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


# ---------------------------------------------------------------------------
# M16.2 Lombardi surface mobility
# ---------------------------------------------------------------------------


def _hand_lombardi_mu_AC(B, C, N_total, E_perp, lam, T):
    """Reference closed-form for mu_AC, computed in Python doubles."""
    return B / E_perp + C * N_total ** lam / (E_perp ** (1.0 / 3.0) * T)


def _hand_lombardi_mu_sr(delta, E_perp):
    return delta / (E_perp ** 2)


def _hand_lombardi_compose(mu_bulk, mu_AC, mu_sr):
    return 1.0 / (1.0 / mu_bulk + 1.0 / mu_AC + 1.0 / mu_sr)


@pytest.mark.parametrize(
    "B, C, N_total, E_perp, lam, T",
    [
        # Three low-doping low-field tuples engineered so B / E_perp
        # dominates the C * N^lam / (E_perp^(1/3) * T) Coulomb-phonon
        # term by at least an order of magnitude. Units: SI m, V, s, K.
        (1.0e6,  1.0e2, 1.0e15, 1.0e4, 0.125, 300.0),
        (4.75e5, 1.0e2, 1.0e15, 1.0e5, 0.125, 300.0),
        (1.0e7,  5.0e1, 1.0e15, 1.0e5, 0.10,  77.0),
    ],
)
def test_lombardi_mu_AC_low_doping_B_term_dominates(B, C, N_total, E_perp, lam, T):
    """In the low-doping low-field limit, the B / E_perp acoustic
    phonon term must dominate the C * N^lam / (E_perp^(1/3) * T)
    Coulomb-phonon term, and the closed-form result must match a
    literal hand restatement to machine precision."""
    mu_AC = lombardi_mu_AC(B, C, N_total, E_perp, lam, T)
    ref = _hand_lombardi_mu_AC(B, C, N_total, E_perp, lam, T)
    assert mu_AC == pytest.approx(ref, rel=1.0e-14, abs=0.0)
    B_term = B / E_perp
    C_term = C * N_total ** lam / (E_perp ** (1.0 / 3.0) * T)
    assert B_term > 5.0 * C_term, (
        f"B_term={B_term} should dominate C_term={C_term} by >5x in "
        "the low-doping low-field limit"
    )


def test_lombardi_mu_sr_diverges_as_E_perp_zero():
    """mu_sr = delta / E_perp^2 must grow without bound as E_perp -> 0
    (so 1/mu_sr -> 0 and the resistor sum reduces to mu_bulk in the
    bulk far from the interface)."""
    delta = 5.82e10
    near_zero = 1.0e-3
    huge = 1.0e3
    mu_sr_near_zero = lombardi_mu_sr(delta, near_zero)
    mu_sr_huge = lombardi_mu_sr(delta, huge)
    assert mu_sr_near_zero > 1.0e9 * mu_sr_huge


def test_lombardi_compose_reduces_to_bulk_when_surface_terms_diverge():
    """1/mu = 1/mu_bulk + 1/mu_AC + 1/mu_sr; with mu_AC, mu_sr -> +inf
    the composite collapses to mu_bulk."""
    mu_bulk = 1400.0
    mu_AC = 1.0e30
    mu_sr = 1.0e30
    mu = lombardi_compose(mu_bulk, mu_AC, mu_sr)
    assert mu == pytest.approx(mu_bulk, rel=1.0e-20, abs=0.0)


def test_lombardi_compose_reduces_to_parallel_surface_terms_when_bulk_diverges():
    """When mu_bulk -> +inf the composite collapses to the parallel
    combination of the two surface terms,
    mu_AC * mu_sr / (mu_AC + mu_sr)."""
    mu_bulk = 1.0e30
    mu_AC = 200.0
    mu_sr = 600.0
    mu = lombardi_compose(mu_bulk, mu_AC, mu_sr)
    expected = mu_AC * mu_sr / (mu_AC + mu_sr)
    assert mu == pytest.approx(expected, rel=1.0e-12, abs=0.0)


def test_lombardi_si_electron_sample_point_matches_hand_calc():
    """At N_total = 1e17 cm^-3 (= 1e23 m^-3), E_perp = 1e6 V/m,
    T = 300 K with the Si electron Lombardi defaults
    (B_n = 4.75e7 cm/s = 4.75e5 m/s; C_n = 1.74e5 cm^(5/3)/V^(2/3)/s
    converted; lambda_n = 0.125; mu0_bulk = 1400 cm^2/V/s = 0.14 m^2/V/s;
    delta_n = 5.82e14 cm^2 V/s = 5.82e10 m^2 V/s), the Lombardi
    composite mobility matches a literal restatement of the formula
    in this module to machine precision.

    Reference values are SI; the test exercises the algebra, not the
    Lombardi physics validity (the latter belongs to the mosfet_2d
    benchmark Pao-Sah comparison)."""
    cm_to_m = 1.0e-2
    cm_to_m_5_3 = cm_to_m ** (5.0 / 3.0)
    cm_to_m_2 = cm_to_m ** 2

    # SI inputs derived from the Lombardi 1988 / Sentaurus Si electron
    # defaults shipped in schemas/input.v2.json.
    B_n = 4.75e7 * cm_to_m
    C_n = 1.74e5 * cm_to_m_5_3
    lam_n = 0.125
    delta_n = 5.82e14 * cm_to_m_2
    mu0_bulk = 1400.0 * cm_to_m_2  # m^2/V/s
    N_total = 1.0e17 * 1.0e6        # cm^-3 -> m^-3
    E_perp = 1.0e6                  # V/m
    T = 300.0

    mu_AC = lombardi_mu_AC(B_n, C_n, N_total, E_perp, lam_n, T)
    mu_sr = lombardi_mu_sr(delta_n, E_perp)
    mu = lombardi_compose(mu0_bulk, mu_AC, mu_sr)

    ref_AC = _hand_lombardi_mu_AC(B_n, C_n, N_total, E_perp, lam_n, T)
    ref_sr = _hand_lombardi_mu_sr(delta_n, E_perp)
    ref = _hand_lombardi_compose(mu0_bulk, ref_AC, ref_sr)

    assert mu_AC == pytest.approx(ref_AC, rel=1.0e-14, abs=0.0)
    assert mu_sr == pytest.approx(ref_sr, rel=1.0e-14, abs=0.0)
    assert mu == pytest.approx(ref, rel=1.0e-14, abs=0.0)
    # Lombardi composite must be no greater than mu_bulk and strictly
    # positive (resistor-sum monotonicity).
    assert 0.0 < mu <= mu0_bulk + 1.0e-12


class _StubScaling:
    """Minimal stand-in for `semi.scaling.Scaling` exposing only the
    fields read by `lombardi_unit_conversions`."""

    def __init__(self, V0: float, mu0: float):
        self.V0 = V0
        self.mu0 = mu0


def test_lombardi_unit_conversions_returns_dimensionless_constants():
    """`lombardi_unit_conversions` converts JSON cm-based parameters
    into the dimensionless ratios consumed by the UFL builder. The
    test checks that the returned `B_n_for_form / E_perp_for_form`
    reproduces the SI mu_AC_hat = mu_AC_SI / sc.mu0 to machine
    precision when E_perp_for_form = grad(psi_hat) . n_hat is the
    scaled gradient (units 1/m) and E_perp_SI = V_t * E_perp_for_form.
    """
    cfg = {
        "B_n": 4.75e7,
        "B_p": 9.93e6,
        "C_n": 1.74e5,
        "C_p": 8.84e5,
        "lambda_n": 0.125,
        "lambda_p": 0.0317,
        "delta_n": 5.82e14,
        "delta_p": 2.0546e14,
    }
    sc = _StubScaling(V0=0.0259, mu0=0.14)  # V_t at 300 K, mu0 in m^2/V/s

    out = lombardi_unit_conversions(cfg, sc)

    # B_n_for_form = B_n_SI / (V_t * sc.mu0).
    B_n_SI = 4.75e7 * 1.0e-2
    expected_B_n = B_n_SI / (sc.V0 * sc.mu0)
    assert out["B_n_for_form"] == pytest.approx(expected_B_n, rel=1.0e-14, abs=0.0)
    assert out["lambda_n"] == pytest.approx(0.125, rel=1.0e-14)
    assert out["lambda_p"] == pytest.approx(0.0317, rel=1.0e-14)

    # Symmetry: B_n / B_p in SI is preserved in the for-form values.
    B_p_SI = 9.93e6 * 1.0e-2
    expected_B_p = B_p_SI / (sc.V0 * sc.mu0)
    assert out["B_p_for_form"] == pytest.approx(expected_B_p, rel=1.0e-14, abs=0.0)

    # Strict positivity of every for-form constant.
    for key, value in out.items():
        if key in ("lambda_n", "lambda_p"):
            continue
        assert value > 0.0, f"{key} must be positive in for-form units"


def test_lombardi_unit_conversions_uses_defaults_for_missing_keys():
    """An empty config produces the Lombardi 1988 / Sentaurus
    defaults; equivalent to the schema default-fill path."""
    sc = _StubScaling(V0=0.0259, mu0=0.14)
    out_empty = lombardi_unit_conversions({}, sc)
    out_full = lombardi_unit_conversions(
        {
            "B_n": 4.75e7,
            "B_p": 9.93e6,
            "C_n": 1.74e5,
            "C_p": 8.84e5,
            "lambda_n": 0.125,
            "lambda_p": 0.0317,
            "delta_n": 5.82e14,
            "delta_p": 2.0546e14,
        },
        sc,
    )
    for key in out_full:
        assert out_empty[key] == pytest.approx(
            out_full[key], rel=1.0e-14, abs=0.0
        ), key


def test_build_mobility_expressions_lombardi_raises_pending_phase_c():
    """The Phase B dispatch raises a clear NotImplementedError on the
    lombardi branch so callers see the missing facet_tags wiring (the
    pure-Python helpers above remain usable directly). Phase C
    replaces this stub with the actual UFL builder."""
    from semi.physics.mobility import build_mobility_expressions

    class _StubFunctionSpace:
        def __init__(self):
            self.mesh = None

    class _StubFunction:
        def __init__(self):
            self.function_space = _StubFunctionSpace()

    cfg = {"model": "lombardi", "interface_facet_tag": 1}
    with pytest.raises(NotImplementedError, match="Phase C"):
        build_mobility_expressions(
            cfg,
            _StubFunction(),
            _StubFunction(),
            1.0,
            1.0,
            sc=_StubScaling(V0=0.0259, mu0=0.14),
        )
