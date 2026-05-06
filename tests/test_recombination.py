"""Tests for the SRH recombination kernel (pure-Python path)."""
from __future__ import annotations

import numpy as np
import pytest

from semi.physics import recombination

N_I = 1.0e16              # m^-3, Si at 300 K (order of magnitude)
TAU_N = 1.0e-7            # s
TAU_P = 1.0e-7            # s


def _np_rate(n, p, E_t_over_Vt=0.0):
    return recombination.srh_rate_np(n, p, N_I, TAU_N, TAU_P, E_t_over_Vt)


def test_srh_zero_at_equilibrium_intrinsic():
    """R = 0 when n p = n_i^2 at the intrinsic point."""
    n = N_I
    p = N_I
    assert _np_rate(n, p) == pytest.approx(0.0, abs=1e-30)


def test_srh_zero_at_equilibrium_extrinsic():
    """R = 0 for any (n, p) obeying mass action n p = n_i^2."""
    n = 1.0e22
    p = N_I ** 2 / n
    assert _np_rate(n, p) == pytest.approx(0.0, rel=0, abs=1.0e-10)


def test_srh_positive_under_forward_bias():
    """Above equilibrium (n p > n_i^2) recombination is strictly positive."""
    n = 1.0e18
    p = 1.0e16
    assert n * p > N_I ** 2
    assert _np_rate(n, p) > 0.0


def test_srh_negative_under_reverse_bias():
    """Below equilibrium (n p < n_i^2) the kernel models net generation."""
    n = 1.0e6
    p = 1.0e6
    assert n * p < N_I ** 2
    assert _np_rate(n, p) < 0.0


def test_srh_low_level_limit_n_type():
    """
    Low-level injection in n-type (n0 >> p, delta_p << n0): R ~ delta_p / tau_p.
    """
    n0 = 1.0e22
    p0 = N_I ** 2 / n0
    delta_p = 1.0e12              # << n0
    n = n0                        # delta_n ~ delta_p << n0 so n stays
    p = p0 + delta_p
    R = _np_rate(n, p)
    assert R == pytest.approx(delta_p / TAU_P, rel=1e-3)


def test_srh_low_level_limit_p_type():
    """Low-level injection in p-type: R ~ delta_n / tau_n."""
    p0 = 1.0e22
    n0 = N_I ** 2 / p0
    delta_n = 1.0e12
    n = n0 + delta_n
    p = p0
    R = _np_rate(n, p)
    assert R == pytest.approx(delta_n / TAU_N, rel=1e-3)


def test_srh_high_level_saturation():
    """
    High-level injection (delta_n = delta_p >> doping): R ~ delta_n / (tau_n + tau_p).
    """
    delta = 1.0e22
    n = delta
    p = delta
    R = _np_rate(n, p)
    # At delta >> n_i, num ~ delta^2 and den ~ (tau_p + tau_n) delta.
    assert R == pytest.approx(delta / (TAU_N + TAU_P), rel=1e-6)


def test_srh_deep_trap_suppresses_rate():
    """Traps far from mid-gap make recombination less efficient."""
    n = 1.0e18
    p = 1.0e16
    R_mid = _np_rate(n, p, E_t_over_Vt=0.0)
    R_deep = _np_rate(n, p, E_t_over_Vt=10.0)
    assert 0.0 < R_deep < R_mid


def test_srh_asymmetric_lifetimes_limit_on_slower_tau():
    """
    With tau_n >> tau_p in strongly n-type material (n >> p), R is set by
    the minority lifetime tau_p, not tau_n.
    """
    n0 = 1.0e22
    delta_p = 1.0e14
    p = N_I ** 2 / n0 + delta_p
    R = recombination.srh_rate_np(n0, p, N_I, tau_n=1.0e-3, tau_p=1.0e-7)
    assert R == pytest.approx(delta_p / 1.0e-7, rel=2e-2)


def test_srh_vector_input_matches_scalar():
    """Vectorized evaluation equals element-wise scalar evaluation."""
    rng = np.random.default_rng(1)
    n = rng.uniform(1.0e16, 1.0e20, size=20)
    p = rng.uniform(1.0e10, 1.0e16, size=20)
    R_vec = _np_rate(n, p)
    R_loop = np.array([_np_rate(float(nn), float(pp)) for nn, pp in zip(n, p, strict=True)])
    np.testing.assert_allclose(R_vec, R_loop, rtol=1e-12)


def test_scaled_tau_conversion():
    t0 = 1.0e-9
    assert recombination.scaled_tau(1.0e-7, t0) == pytest.approx(100.0)


# ---------------------------------------------------------------------------
# M16.3 Auger schema-surface tests (Phase A).
# ---------------------------------------------------------------------------

import glob  # noqa: E402
import json  # noqa: E402
import warnings  # noqa: E402
from pathlib import Path  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[1]
V2_SCHEMA_PATH = REPO_ROOT / "schemas" / "input.v2.json"
BENCHMARKS_DIR = REPO_ROOT / "benchmarks"


def _v2_validator():
    import jsonschema

    schema = json.loads(V2_SCHEMA_PATH.read_text())
    return jsonschema.Draft7Validator(schema)


def _base_v2_cfg(schema_version: str = "2.3.0") -> dict:
    return {
        "schema_version": schema_version,
        "name": "auger_test",
        "dimension": 1,
        "mesh": {
            "source": "builtin",
            "extents": [[0.0, 1.0e-6]],
            "resolution": [100],
            "facets_by_plane": [
                {"name": "L", "tag": 1, "axis": 0, "value": 0.0},
                {"name": "R", "tag": 2, "axis": 0, "value": 1.0e-6},
            ],
        },
        "regions": {"si": {"material": "Si", "tag": 1, "role": "semiconductor"}},
        "doping": [
            {"region": "si", "profile": {"type": "uniform", "N_D": 1e17, "N_A": 0.0}},
        ],
        "contacts": [
            {"name": "L", "facet": "L", "type": "ohmic", "voltage": 0.0},
            {"name": "R", "facet": "R", "type": "ohmic", "voltage": 0.0},
        ],
    }


def test_every_benchmark_validates_under_v2_3_0():
    """Existing v2.0.0 / v2.1.0 / v2.2.0 benchmark JSONs validate
    unchanged against the v2.3.0 schema (additive minor)."""
    v = _v2_validator()
    paths = sorted(glob.glob(str(BENCHMARKS_DIR / "*" / "*.json")))
    assert paths, "no benchmark JSONs found"
    for p in paths:
        cfg = json.loads(Path(p).read_text())
        errors = sorted(v.iter_errors(cfg), key=lambda e: list(e.path))
        assert not errors, (
            f"{p} fails v2.3.0 validation:\n"
            + "\n".join(
                f"  at {'.'.join(str(x) for x in e.absolute_path) or '<root>'}: {e.message}"
                for e in errors[:5]
            )
        )


def test_auger_true_with_default_C_n_C_p_validates():
    """auger=true with default Si C_n / C_p validates under v2.3.0."""
    v = _v2_validator()
    cfg = _base_v2_cfg()
    cfg["physics"] = {
        "recombination": {
            "auger": True,
            "C_n": 2.8e-31,
            "C_p": 9.9e-32,
        }
    }
    errors = sorted(v.iter_errors(cfg), key=lambda e: list(e.path))
    assert not errors, [e.message for e in errors]


def test_negative_C_n_is_rejected():
    """Schema declares minimum: 0 on C_n / C_p; a negative coefficient
    must be rejected."""
    v = _v2_validator()
    cfg = _base_v2_cfg()
    cfg["physics"] = {
        "recombination": {
            "auger": True,
            "C_n": -1.0e-30,
        }
    }
    errors = sorted(v.iter_errors(cfg), key=lambda e: list(e.path))
    assert errors, "negative C_n should fail validation"
    messages = " | ".join(e.message for e in errors)
    assert "C_n" in messages or "minimum" in messages


def test_negative_C_p_is_rejected():
    v = _v2_validator()
    cfg = _base_v2_cfg()
    cfg["physics"] = {
        "recombination": {
            "auger": True,
            "C_p": -1.0e-30,
        }
    }
    errors = sorted(v.iter_errors(cfg), key=lambda e: list(e.path))
    assert errors


def test_recombination_extra_property_is_rejected():
    """Strict-v2 additionalProperties:false on physics.recombination."""
    v = _v2_validator()
    cfg = _base_v2_cfg()
    cfg["physics"] = {"recombination": {"auger": True, "C_unknown": 1.0e-30}}
    errors = sorted(v.iter_errors(cfg), key=lambda e: list(e.path))
    assert errors
    messages = " | ".join(e.message for e in errors)
    assert "C_unknown" in messages or "additional" in messages.lower()


def test_default_fill_auger_off_with_C_defaults():
    """validate() default-fills auger=False and Si C_n / C_p."""
    from semi import schema as schema_mod

    cfg = _base_v2_cfg()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        cfg = schema_mod.validate(cfg)
    rec = cfg["physics"]["recombination"]
    assert rec["auger"] is False
    assert rec["C_n"] == pytest.approx(2.8e-31)
    assert rec["C_p"] == pytest.approx(9.9e-32)


def test_default_fill_does_not_override_user_C_n():
    """User-supplied C_n / C_p survive default-fill."""
    from semi import schema as schema_mod

    cfg = _base_v2_cfg()
    cfg["physics"] = {"recombination": {"auger": True, "C_n": 1.0e-30}}
    cfg = schema_mod.validate(cfg)
    rec = cfg["physics"]["recombination"]
    assert rec["auger"] is True
    assert rec["C_n"] == pytest.approx(1.0e-30)
    # C_p was unset; default fill applies.
    assert rec["C_p"] == pytest.approx(9.9e-32)


def test_v2_2_0_input_still_validates_after_minor_bump():
    """v2.2.0 inputs (no auger / C_n / C_p) continue to validate
    against the bumped schema; this is the additive-minor invariant."""
    v = _v2_validator()
    cfg = _base_v2_cfg(schema_version="2.2.0")
    errors = sorted(v.iter_errors(cfg), key=lambda e: list(e.path))
    assert not errors, [e.message for e in errors]


def test_schema_supported_minor_is_3():
    """Sanity: SCHEMA_SUPPORTED_MINOR was bumped to 3 in M16.3."""
    from semi import schema as schema_mod

    assert schema_mod.SCHEMA_SUPPORTED_MINOR == 3


# ---------------------------------------------------------------------------
# M16.3 Auger pure-Python tests (Phase B).
# ---------------------------------------------------------------------------

# Si Dziewior-Schmid coefficients in SI (m^6/s); cm^6/s = 1e-12 * m^6/s.
C_N_SI = 2.8e-31 * 1.0e-12
C_P_SI = 9.9e-32 * 1.0e-12


def test_auger_zero_at_equilibrium_intrinsic():
    """R_Auger = 0 when n p = n_i^2 at the intrinsic point."""
    n = N_I
    p = N_I
    R = recombination.auger_rate_np(n, p, N_I, C_N_SI, C_P_SI)
    assert R == pytest.approx(0.0, abs=1e-30)


def test_auger_zero_at_equilibrium_extrinsic():
    """R_Auger = 0 for any (n, p) obeying mass action n p = n_i^2."""
    n = 1.0e22
    p = N_I ** 2 / n
    R = recombination.auger_rate_np(n, p, N_I, C_N_SI, C_P_SI)
    # Numerical zero on a moderately wide dynamic range.
    assert abs(R) < 1.0e-15 * (C_N_SI * n + C_P_SI * p) * n * p


def test_auger_positive_under_forward_bias():
    """Above equilibrium (n p > n_i^2) Auger rate is strictly positive."""
    n = 1.0e22
    p = 1.0e21
    R = recombination.auger_rate_np(n, p, N_I, C_N_SI, C_P_SI)
    assert R > 0.0


def test_auger_high_injection_cubic_limit():
    """At n ~ p >> n_i, R_Auger -> (C_n + C_p) X^3 within 1 %."""
    X = 1.0e24                      # well above N_I = 1e16
    n = X
    p = X
    R = recombination.auger_rate_np(n, p, N_I, C_N_SI, C_P_SI)
    R_ref = (C_N_SI + C_P_SI) * X ** 3
    assert R == pytest.approx(R_ref, rel=0.01)


def test_auger_low_injection_n_type_dominated_by_C_n_n():
    """In strongly n-type material at modest injection, R_Auger is
    dominated by the C_n n (n p - n_i^2) branch (electron lifetime
    sets the scale)."""
    n0 = 1.0e25                     # heavy n-type
    p0 = N_I ** 2 / n0
    delta_p = 1.0e15
    n = n0
    p = p0 + delta_p
    R = recombination.auger_rate_np(n, p, N_I, C_N_SI, C_P_SI)
    # Expected leading-order: C_n n^2 delta_p (since C_n n >> C_p p
    # when n >> p, and n p - n_i^2 ~ n delta_p).
    R_ref = C_N_SI * n ** 2 * delta_p
    assert R == pytest.approx(R_ref, rel=1.0e-3)


def test_auger_negative_under_reverse_bias():
    """Below equilibrium (n p < n_i^2) Auger models net generation."""
    n = 1.0e6
    p = 1.0e6
    assert n * p < N_I ** 2
    R = recombination.auger_rate_np(n, p, N_I, C_N_SI, C_P_SI)
    assert R < 0.0


def test_auger_vector_input_matches_scalar():
    rng = np.random.default_rng(7)
    n = rng.uniform(1.0e16, 1.0e22, size=20)
    p = rng.uniform(1.0e10, 1.0e16, size=20)
    R_vec = recombination.auger_rate_np(n, p, N_I, C_N_SI, C_P_SI)
    R_loop = np.array([
        recombination.auger_rate_np(float(nn), float(pp), N_I, C_N_SI, C_P_SI)
        for nn, pp in zip(n, p, strict=True)
    ])
    np.testing.assert_allclose(R_vec, R_loop, rtol=1e-12)


def test_scaled_auger_C_round_trip():
    """C_hat = C_SI * C0^2 * t0; recover C_SI by dividing back out."""
    C_si = C_N_SI
    C0 = 1.0e23                     # m^-3
    t0 = 1.0e-12                    # s
    C_hat = recombination.scaled_auger_C(C_si, C0, t0)
    assert C_hat == pytest.approx(C_si * C0 ** 2 * t0)
    # And the inverse round-trip:
    assert C_hat / (C0 ** 2 * t0) == pytest.approx(C_si)


def test_scaled_auger_C_dimensionless_check():
    """For Si-like (C0 ~ 1e23 m^-3, t0 ~ 1e-12 s) the scaled C is O(1e3)
    or so; sanity check that we are not off by orders of magnitude."""
    C_hat = recombination.scaled_auger_C(C_N_SI, 1.0e23, 1.0e-12)
    # C_N_SI = 2.8e-43 m^6/s, C0^2 = 1e46 m^-6, t0 = 1e-12 s.
    # Product = 2.8e-43 * 1e46 * 1e-12 = 2.8e-9. Should be small but
    # finite and positive.
    assert 0 < C_hat < 1.0


def test_auger_rate_scaled_form_consistency():
    """Round-trip: R_SI computed via auger_rate_np, then divided by
    C0/t0, should equal auger_rate_np applied to scaled inputs with
    C_hat = C_SI * C0^2 * t0."""
    C0 = 1.0e23
    t0 = 1.0e-12

    n = 1.0e22                      # m^-3
    p = 5.0e21
    n_i = 1.0e16
    R_SI = recombination.auger_rate_np(n, p, n_i, C_N_SI, C_P_SI)
    R_hat_via_scaling = R_SI / (C0 / t0)

    n_hat = n / C0
    p_hat = p / C0
    n_i_hat = n_i / C0
    C_n_hat = recombination.scaled_auger_C(C_N_SI, C0, t0)
    C_p_hat = recombination.scaled_auger_C(C_P_SI, C0, t0)
    R_hat_direct = recombination.auger_rate_np(
        n_hat, p_hat, n_i_hat, C_n_hat, C_p_hat
    )
    assert R_hat_direct == pytest.approx(R_hat_via_scaling, rel=1.0e-12)
