"""
Mesh convergence for the M12 MOSFET benchmark.

Walks `geometry.characteristic_length` through `[1e-7, 5e-8, 2.5e-8]`
and checks that (a) the M12 geometry pipeline (gmsh subprocess -> cached
.msh -> dolfinx load) produces a monotonically finer mesh each time and
(b) the bias_sweep solver converges end-to-end at every resolution.

What this test does NOT assert: a tight Cauchy-sequence rate on
|I_D(V_DS = 0.1)|. The current at the drain contact is subthreshold
leakage dominated by numerical noise from the Poisson solver, because
`semi.runners.bias_sweep` today treats the Si/SiO2 multi-region mesh
as a single-material Si mesh (`build_dd_block_residual` uses a scalar
`ref_mat.epsilon_r`). Wiring multi-region coupled drift-diffusion
through `bias_sweep` (mirroring what `run_mos_cv` does for equilibrium
Poisson) is a separate milestone; this test intentionally stops at the
geometry-pipeline gate.

Skipped when gmsh is not on PATH -- the M12 geometry pipeline needs
gmsh as a subprocess, so the test has nothing to exercise without it.
"""
from __future__ import annotations

import copy
import json
import shutil
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
BENCHMARK_JSON = REPO_ROOT / "benchmarks" / "mosfet_2d" / "mosfet_2d.json"

CHAR_LENGTHS = [1.0e-7, 5.0e-8, 2.5e-8]

pytestmark = pytest.mark.skipif(
    shutil.which("gmsh") is None, reason="gmsh not on PATH"
)


def _solve_at_clmax(clmax: float) -> tuple[float, int, int]:
    """Run the MOSFET benchmark at a given characteristic_length.

    Returns (|I_D(V_DS=0.1 V)|, n_cells, n_iter_last).
    """
    from semi import run as semi_run
    from semi import schema

    cfg = json.loads(BENCHMARK_JSON.read_text())
    cfg["geometry"] = copy.deepcopy(cfg["geometry"])
    cfg["geometry"]["characteristic_length"] = float(clmax)
    cfg["_source_dir"] = str(BENCHMARK_JSON.parent.resolve())
    cfg = schema.validate(cfg)
    cfg["_source_dir"] = str(BENCHMARK_JSON.parent.resolve())

    result = semi_run.run(cfg)
    iv = result.iv or []
    assert iv, f"no IV rows at clmax={clmax}"
    target = min(iv, key=lambda r: abs(r["V"] - 0.1))
    I_D = abs(float(target.get("J", 0.0))) * 1.0e-6  # A/m per unit depth
    n_cells = int(result.mesh.topology.index_map(result.mesh.topology.dim).size_global)
    n_iter = int((result.solver_info or {}).get("iterations", -1))
    return I_D, n_cells, n_iter


def test_mosfet_geometry_pipeline_runs_at_three_resolutions(tmp_path, monkeypatch):
    """The M12 pipeline produces a monotonically finer mesh at each clmax.

    Each solve must converge end-to-end. The Cauchy rate on the
    resulting drain-current sequence is recorded for inspection only;
    see the module docstring for why this is not a gated assertion.
    """
    monkeypatch.setenv("KRONOS_MESH_CACHE", str(tmp_path))

    currents: list[float] = []
    cell_counts: list[int] = []

    for cl in CHAR_LENGTHS:
        I_D, n_cells, n_iter = _solve_at_clmax(cl)
        currents.append(I_D)
        cell_counts.append(n_cells)
        print(
            f"[mosfet_convergence] clmax={cl:.2e}  n_cells={n_cells}  "
            f"|I_D|={I_D:.4e}  last_snes_iters={n_iter}",
            file=sys.stderr,
        )
        assert n_iter > 0, f"SNES did not report iterations at clmax={cl}"

    # Mesh must refine with smaller clmax.
    assert cell_counts[0] < cell_counts[1] < cell_counts[2], (
        f"cell counts not strictly increasing: {cell_counts}"
    )

    # All currents are finite (non-NaN).
    for cl, I in zip(CHAR_LENGTHS, currents, strict=True):
        assert I == I, f"NaN current at clmax={cl}"

    # Report the Cauchy-sequence ratio for the record; not gated.
    if len(currents) >= 3:
        d_coarse = abs(currents[1] - currents[0])
        d_fine = abs(currents[2] - currents[1])
        ratio = d_coarse / d_fine if d_fine > 0.0 else float("inf")
        print(
            f"[mosfet_convergence] Cauchy ratio |I(h1)-I(h0)|/|I(h2)-I(h1)| "
            f"= {ratio:.3f} (not asserted; see module docstring)",
            file=sys.stderr,
        )
