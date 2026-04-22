"""
End-to-end tests for the Day 7 3D resistor benchmark.

Two scopes:

1. Smoke / wiring: a coarsened version of the production resistor JSON
   runs through the bias-sweep solver, populates `result.iv` with one row
   per sweep point, and the verifier returns at least one PASS line.

2. Theory match on the production JSON: builtin and gmsh fixtures both
   solve to within 1% of the analytical resistance R = L / (q N_D mu_n A).
   This is the same gate the runner enforces, exercised here so a CI
   regression caught at the unit-test level is easier to localize than a
   benchmark-script failure.
"""
from __future__ import annotations

import copy
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
RESISTOR_JSON = REPO_ROOT / "benchmarks" / "resistor_3d" / "resistor.json"
RESISTOR_GMSH_JSON = REPO_ROOT / "benchmarks" / "resistor_3d" / "resistor_gmsh.json"
GMSH_FIXTURE = (
    REPO_ROOT / "benchmarks" / "resistor_3d" / "fixtures" / "box.msh"
)


def _load_cfg(path: Path) -> dict:
    """Load a config exactly the way scripts/run_benchmark.py does."""
    from semi import schema

    return schema.load(str(path))


def _coarse_builtin_cfg() -> dict:
    """Coarsen the production builtin JSON for fast smoke testing."""
    cfg = _load_cfg(RESISTOR_JSON)
    cfg = copy.deepcopy(cfg)
    cfg["mesh"]["resolution"] = [16, 4, 4]
    # Keep the sweep small so the smoke runs in seconds.
    sweep = cfg["contacts"][1]["voltage_sweep"]
    sweep["start"] = -0.005
    sweep["stop"] = 0.005
    sweep["step"] = 0.005
    return cfg


@pytest.mark.skipif(not RESISTOR_JSON.exists(), reason="resistor.json missing")
def test_resistor_3d_smoke_runs_and_records_iv():
    """The bias_sweep solver runs end-to-end on a coarsened 3D bar."""
    from semi.run import run

    cfg = _coarse_builtin_cfg()
    result = run(cfg)

    assert result.solver_info.get("converged"), result.solver_info
    iv = result.iv
    assert iv is not None and len(iv) >= 3, f"got iv={iv!r}"
    assert all("V" in row and "J" in row for row in iv)
    # The last sweep point must equal the configured endpoint.
    Vs = sorted(float(r["V"]) for r in iv)
    assert Vs[0] == pytest.approx(-0.005, abs=1.0e-9)
    assert Vs[-1] == pytest.approx(0.005, abs=1.0e-9)


@pytest.mark.skipif(not RESISTOR_JSON.exists(), reason="resistor.json missing")
def test_resistor_3d_verifier_passes_on_builtin_production_json():
    """Production builtin config: verifier returns all PASS lines."""
    from scripts.run_benchmark import verify_resistor_3d
    from semi.run import run

    cfg = _load_cfg(RESISTOR_JSON)
    result = run(cfg)
    assert result.solver_info.get("converged"), result.solver_info

    checks = verify_resistor_3d(result)
    failures = [(name, msg) for name, ok, msg in checks if not ok]
    assert not failures, f"verifier failures: {failures}"

    # Spot-check the headline tolerance: max |R_sim - R_theory|/R_theory < 1%.
    curve = getattr(result, "_resistor_curve", None)
    assert curve is not None
    rel = np.abs(curve["R_sim_at_nonzero"] - curve["R_theory"]) / curve["R_theory"]
    assert rel.max() < 0.01, f"max rel err {rel.max():.4%}"


@pytest.mark.skipif(
    not (RESISTOR_GMSH_JSON.exists() and GMSH_FIXTURE.exists()),
    reason="resistor_gmsh.json or box.msh missing",
)
def test_resistor_3d_verifier_passes_on_gmsh_production_json():
    """Production gmsh config: verifier returns all PASS lines.

    Also asserts R_sim from the gmsh path agrees with R_theory within 1%,
    which by transitivity matches the builtin run's R_sim within ~1%.
    """
    from scripts.run_benchmark import verify_resistor_3d
    from semi.run import run

    cfg = _load_cfg(RESISTOR_GMSH_JSON)
    result = run(cfg)
    assert result.solver_info.get("converged"), result.solver_info

    checks = verify_resistor_3d(result)
    failures = [(name, msg) for name, ok, msg in checks if not ok]
    assert not failures, f"verifier failures: {failures}"

    curve = getattr(result, "_resistor_curve", None)
    assert curve is not None
    rel = np.abs(curve["R_sim_at_nonzero"] - curve["R_theory"]) / curve["R_theory"]
    assert rel.max() < 0.01, f"max rel err {rel.max():.4%}"
