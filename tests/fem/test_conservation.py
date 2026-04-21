"""
FEM-level conservation tests (Phase 3 V&V).

These tests run the full Poisson / drift-diffusion solver on the
standard pn_1d benchmarks and assert:

    1) charge neutrality (pn_1d equilibrium):
       |Q_net| < 1e-10 * q * max|N_net| * L

    2) current continuity on the final forward bias (V=0.6 V):
       max|J_total(x) - mean(J_total)| / |mean(J_total)| < 5%

    3) current continuity on the final reverse bias (V=-2.0 V):
       max|J_total(x) - mean(J_total)| / |mean(J_total)| < 15%

Skipped at collection time when dolfinx is not importable (see
tests/fem/conftest.py). These tests are slow (~1 min each for the
bias sweeps) so they are kept in tests/fem/ rather than tests/.
"""
from __future__ import annotations

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]


def test_pn_1d_charge_neutrality_at_equilibrium():
    """
    Integrate q * (p - n + N_net) over the pn_1d device and require
    the result to be at least ten orders of magnitude smaller than
    q * max|N_net| * L. That matches the snes_rtol ~ 1e-14 used by
    run_equilibrium and is the tightest honest threshold.
    """
    from semi import run as semi_run
    from semi import schema
    from semi.verification.conservation import charge_conservation_from_result

    cfg = schema.load(REPO_ROOT / "benchmarks" / "pn_1d" / "pn_junction.json")
    result = semi_run.run(cfg)
    m = charge_conservation_from_result(result)
    assert m.Q_ref > 0.0, f"empty Q_ref: {m}"
    assert m.rel_error < 1.0e-10, (
        f"|Q_net|={m.Q_net:.3e} C/m^2 fails 1e-10 * Q_ref "
        f"({m.Q_ref:.3e}); rel={m.rel_error:.3e}"
    )


def test_pn_1d_bias_forward_current_continuity():
    """
    At the end of the pn_1d_bias forward sweep (V=0.6 V) the interior
    current continuity must be within 5%.
    """
    from semi import run as semi_run
    from semi import schema
    from semi.verification.conservation import current_continuity_from_result

    cfg_path = REPO_ROOT / "benchmarks" / "pn_1d_bias" / "pn_junction_bias.json"
    cfg = schema.load(cfg_path)
    result = semi_run.run(cfg)
    cc = current_continuity_from_result(result, n_samples=10)
    assert cc.xs.size >= 8, f"expected ~10 samples, got {cc.xs.size}"
    assert cc.max_rel_dev < 0.05, (
        f"forward J continuity: max rel dev {cc.max_rel_dev*100:.2f}% "
        f"(mean={cc.mean_J:.3e}, max_dev={cc.max_abs_dev:.3e})"
    )


def test_pn_1d_bias_reverse_current_continuity():
    """
    At the end of the pn_1d_bias_reverse sweep (V=-2.0 V) the interior
    current continuity must be within 15%.
    """
    from semi import run as semi_run
    from semi import schema
    from semi.verification.conservation import current_continuity_from_result

    cfg_path = (
        REPO_ROOT / "benchmarks" / "pn_1d_bias_reverse" / "pn_junction_bias_reverse.json"
    )
    cfg = schema.load(cfg_path)
    result = semi_run.run(cfg)
    cc = current_continuity_from_result(result, n_samples=10)
    assert cc.xs.size >= 8, f"expected ~10 samples, got {cc.xs.size}"
    assert cc.max_rel_dev < 0.15, (
        f"reverse J continuity: max rel dev {cc.max_rel_dev*100:.2f}% "
        f"(mean={cc.mean_J:.3e}, max_dev={cc.max_abs_dev:.3e})"
    )
