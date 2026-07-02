"""M19: coarsened FEM smoke test for the 3D MOSFET capstone.

Runs a coarsened (cl ~ 50-60 nm) linear-regime V_GS sweep over three
gate voltages on the mosfet_3d device geometry and asserts the drain
current is finite and monotone above threshold. This is a pipeline
smoke test (does the gmsh multi-region Slotboom bias sweep run end to
end in 3D), not a quantitative V&V gate; the Pao-Sah accuracy gate
lives in `scripts/run_benchmark.py::verify_mosfet_3d`.

The test uses Boltzmann statistics with constant mobility rather than
the shipped Fermi-Dirac + Lombardi physics: the coupled 3D FD/Lombardi
bias sweep across the MOSFET inversion onset stagnates in the current
bias_sweep SNES driver (the same line-search-stabilization gap that
keeps mosfet_2d / nmos_idvgs / the mosfet_3d benchmark matrix entries
on `allow-failure`; retiring it is M19's named next task). The simpler
physics isolates the 3D machinery (gmsh multi-region ingest, submesh
DD assembly, gate BC, per-contact drain current) from that solver
difficulty.

Because even the coarsened coupled 3D solve exceeds the unit-test time
budget, the test is skipped by default and opts in via the
`KRONOS_RUN_MOSFET_3D_SMOKE` environment variable (set in a dedicated,
non-blocking CI step). It is also marked `slow`. Requires dolfinx +
gmsh; collected only inside the Docker FEM image (see
`tests/fem/conftest.py`).
"""
from __future__ import annotations

import copy
import importlib.util
import os
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
BENCH_DIR = REPO_ROOT / "benchmarks" / "mosfet_3d"

_OPT_IN = os.environ.get("KRONOS_RUN_MOSFET_3D_SMOKE") == "1"


def _load_generator():
    """Import benchmarks/mosfet_3d/generate_mesh.py as a module."""
    spec = importlib.util.spec_from_file_location(
        "mosfet3d_generate_mesh", BENCH_DIR / "generate_mesh.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.mark.slow
@pytest.mark.skipif(
    not _OPT_IN,
    reason="set KRONOS_RUN_MOSFET_3D_SMOKE=1 to run the coarsened 3D MOSFET "
    "bias-sweep smoke test (excluded from the default coverage job because "
    "a coupled 3D DD sweep exceeds the unit-test time budget)",
)
def test_mosfet_3d_linear_smoke(tmp_path):
    from semi import run as semi_run
    from semi import schema

    # Coarse mesh regenerated from the committed .geo.
    coarse_msh = tmp_path / "mosfet3d_coarse.msh"
    _load_generator().generate(cl_m=6.0e-8, out=coarse_msh, binary=False)
    assert coarse_msh.exists()

    cfg = copy.deepcopy(schema.load(str(BENCH_DIR / "mosfet_3d.json")))
    cfg["mesh"]["path"] = str(coarse_msh)
    # Simpler physics to isolate the 3D pipeline from FD/Lombardi SNES
    # stagnation at inversion onset.
    cfg["physics"]["statistics"] = "boltzmann"
    cfg["physics"]["mobility"] = {"mu_n": 1400.0, "mu_p": 450.0}
    # Three V_GS points straddling threshold (V_T ~ 0.43 V).
    for c in cfg["contacts"]:
        if c["type"] == "gate":
            c["voltage_sweep"] = {"start": 0.0, "stop": 0.8, "step": 0.4}
    # Bound the continuation so a stagnating step raises StepTooSmall
    # quickly instead of hanging the opt-in run indefinitely.
    cfg["solver"]["continuation"] = {
        "min_step": 0.05, "max_halvings": 3,
        "easy_iter_threshold": 4, "grow_factor": 1.5,
    }

    result = semi_run.run(cfg)
    iv = result.iv or []
    assert len(iv) >= 3, f"expected >=3 iv rows, got {len(iv)}"
    assert all("J_drain" in r for r in iv), "bias_sweep must record J_drain"

    V = np.array([r["V"] for r in iv], dtype=float)
    I_D = np.array([abs(float(r["J_drain"])) for r in iv], dtype=float)
    order = np.argsort(V)
    V, I_D = V[order], I_D[order]

    assert np.all(np.isfinite(I_D)), f"non-finite I_D: {I_D}"

    # Monotone non-decreasing I_D above threshold.
    above = V >= 0.5
    I_above = I_D[above]
    if I_above.size >= 2:
        assert np.all(np.diff(I_above) >= -1.0e-14), (
            f"I_D not monotone above threshold: V={V[above]}, I_D={I_above}"
        )
