"""
Pure-Python unit tests for the mosfet_2d Pao-Sah verifier helpers.

The full FEM benchmark exercises this code in CI via
`python scripts/run_benchmark.py mosfet_2d`. These tests cover the
pure-Python helpers (device-parameter extraction, Pao-Sah formula
evaluation, verifier pass/fail logic) by mocking a `SimulationResult`
with a fabricated iv table.
"""
from __future__ import annotations

import json
import sys
import types
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = REPO_ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))


def _ensure_matplotlib_stub() -> None:
    """Stub matplotlib only when it is not already importable.

    The verifier helpers do not use matplotlib at runtime; only the
    plotter does. Loading the script unconditionally imports
    matplotlib.pyplot, so a minimal stub lets the tests run on systems
    that do not have matplotlib. When matplotlib *is* installed, leave
    sys.modules untouched so other tests in the same pytest session
    continue to see the real package.
    """
    try:
        import matplotlib  # noqa: F401
        import matplotlib.pyplot  # noqa: F401
        import matplotlib.tri  # noqa: F401
        return
    except ImportError:
        pass

    mpl = types.ModuleType("matplotlib")
    mpl.use = lambda *a, **kw: None
    sys.modules["matplotlib"] = mpl

    pyplot = types.ModuleType("matplotlib.pyplot")
    def _noop(*a, **kw):
        return None, None
    pyplot.subplots = _noop
    pyplot.close = lambda *a, **kw: None
    sys.modules["matplotlib.pyplot"] = pyplot

    tri = types.ModuleType("matplotlib.tri")
    tri.Triangulation = lambda *a, **kw: None
    sys.modules["matplotlib.tri"] = tri


_ensure_matplotlib_stub()


@pytest.fixture(scope="module")
def benchmark_cfg() -> dict:
    path = REPO_ROOT / "benchmarks" / "mosfet_2d" / "mosfet_2d.json"
    with path.open() as fh:
        return json.load(fh)


def _import_run_benchmark():
    import importlib

    return importlib.import_module("run_benchmark")


def test_mosfet_device_params_extracts_V_T(benchmark_cfg):
    rb = _import_run_benchmark()
    from semi.materials import get_material

    mat = get_material(benchmark_cfg["regions"]["silicon"]["material"])
    dp = rb._mosfet_device_params(benchmark_cfg, mat)

    # Long-channel V_T for N_A = 1e16 cm^-3, t_ox = 5 nm, ideal gate
    # (phi_ms = 0) in the kronos-semi BC convention (psi_intrinsic = 0
    # at the ohmic body in equilibrium, V_FB = phi_ms - phi_F). For
    # N_A = 1e16, phi_F ~ 0.358 V at 300 K and the textbook V_T ~ 0.787
    # V drops to V_T_kronos ~ 0.429 V after the convention shift.
    assert 0.35 <= dp["V_T"] <= 0.50, f"V_T = {dp['V_T']:.3f} V outside expected band"
    assert dp["V_DS"] == pytest.approx(0.05, abs=1e-9)
    assert dp["L_ch"] == pytest.approx(3.0e-6, abs=1e-9)
    assert dp["L_dc"] == pytest.approx(2.005e-6, abs=1e-9)
    assert dp["t_ox"] == pytest.approx(5.0e-9, abs=1e-12)
    assert dp["mu_n_SI"] == pytest.approx(0.14, rel=1e-9)
    # C_ox = 3.9 * eps_0 / 5 nm = 6.906 mF/m^2 (~ 0.69 uF/cm^2)
    assert dp["C_ox"] == pytest.approx(6.91e-3, rel=1e-2)


def test_pao_sah_linear_formula_above_threshold(benchmark_cfg):
    rb = _import_run_benchmark()
    from semi.materials import get_material

    mat = get_material(benchmark_cfg["regions"]["silicon"]["material"])
    dp = rb._mosfet_device_params(benchmark_cfg, mat)

    V_GS = np.array([dp["V_T"] - 0.1, dp["V_T"] + 0.4])
    I_D_per_W = rb._pao_sah_linear_I_D_per_W(V_GS, dp)

    assert np.isnan(I_D_per_W[0]), "sub-threshold V_GS must yield NaN"
    expected = (dp["mu_n_SI"] / dp["L_ch"]) * dp["C_ox"] * 0.4 * dp["V_DS"]
    assert I_D_per_W[1] == pytest.approx(expected, rel=1e-9)


def _fabricated_result(cfg: dict, *, scale: float = 1.0):
    """Build a fake SimulationResult whose recorded J_drain matches Pao-Sah
    times `scale`. scale = 1.0 -> verifier passes; scale far from 1.0 ->
    verifier fails on the relative-error check.
    """
    rb = _import_run_benchmark()
    from semi.materials import get_material

    mat = get_material(cfg["regions"]["silicon"]["material"])
    dp = rb._mosfet_device_params(cfg, mat)

    V_GS_arr = np.linspace(0.0, 1.5, 16)
    iv = []
    for V in V_GS_arr:
        if V <= dp["V_T"]:
            J_drain = 0.0
        else:
            I_per_W = (dp["mu_n_SI"] / dp["L_ch"]) * dp["C_ox"] * (V - dp["V_T"]) * dp["V_DS"]
            J_drain = scale * I_per_W / dp["L_dc"]
        iv.append({"V": float(V), "J": 0.0, "J_drain": float(J_drain)})
    return SimpleNamespace(cfg=cfg, iv=iv)


def test_verifier_passes_when_simulation_matches_analytical(benchmark_cfg):
    rb = _import_run_benchmark()
    result = _fabricated_result(benchmark_cfg, scale=1.0)
    checks = rb.verify_mosfet_2d(result)

    assert checks, "verifier returned no checks"
    failed = [c for c in checks if not c[1]]
    assert not failed, f"verifier failed on a perfect match: {failed}"


def test_verifier_fails_when_simulation_diverges(benchmark_cfg):
    rb = _import_run_benchmark()
    result = _fabricated_result(benchmark_cfg, scale=1.5)
    checks = rb.verify_mosfet_2d(result)

    rel_err_check = next(
        c for c in checks if "I_D_PaoSah" in c[0]
    )
    assert rel_err_check[1] is False, (
        "verifier should reject 50 % over-estimate as outside 20 % tolerance"
    )


def test_verifier_reports_missing_J_drain(benchmark_cfg):
    rb = _import_run_benchmark()
    iv = [{"V": 0.5, "J": 0.0}]  # no J_drain key
    result = SimpleNamespace(cfg=benchmark_cfg, iv=iv)
    checks = rb.verify_mosfet_2d(result)

    assert any("J_drain" in c[0] for c in checks)
    assert all(not c[1] for c in checks if "J_drain" in c[0])


def test_verifier_handles_empty_iv(benchmark_cfg):
    rb = _import_run_benchmark()
    result = SimpleNamespace(cfg=benchmark_cfg, iv=[])
    checks = rb.verify_mosfet_2d(result)

    assert len(checks) == 1
    assert checks[0][1] is False
