"""
Shared utilities for the audit suite.

Every audit test consumes one or more existing benchmark JSON configs,
runs two or more runners on a shared problem, computes a comparison
metric, and writes the results to `/tmp/audit/<case>.csv` and a
markdown fragment to `/tmp/audit/<case>.md`. The conftest aggregator
collects the markdown fragments into `docs/PHYSICS_AUDIT.md`.
"""
from __future__ import annotations

import copy
import csv
import importlib
import json
import math
from pathlib import Path
from typing import Any

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
BENCH_DIR = REPO_ROOT / "benchmarks"
AUDIT_DIR = Path("/tmp/audit")


def require_dolfinx():
    """Skip the calling test if dolfinx is not importable.

    The audit tests exercise the FEM runners; without dolfinx they
    cannot run. The host environment routinely lacks dolfinx; CI runs
    this suite inside the dolfinx-env Docker image.
    """
    try:
        importlib.import_module("dolfinx")
    except Exception as exc:  # pragma: no cover - environment-dependent
        pytest.skip(f"dolfinx not available: {exc}")


def load_benchmark(rel_path: str) -> dict[str, Any]:
    """Load a benchmark JSON config relative to `benchmarks/`."""
    with open(BENCH_DIR / rel_path) as fh:
        return json.load(fh)


def write_csv(case: str, header: list[str], rows: list[list[Any]]) -> Path:
    """Write a CSV under /tmp/audit/<case>.csv and return the path."""
    AUDIT_DIR.mkdir(parents=True, exist_ok=True)
    path = AUDIT_DIR / f"{case}.csv"
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)
    return path


def write_markdown(case: str, title: str, body: str) -> Path:
    """Write a markdown fragment that the conftest aggregator will pick up."""
    AUDIT_DIR.mkdir(parents=True, exist_ok=True)
    path = AUDIT_DIR / f"{case}.md"
    path.write_text(f"## {title}\n\n{body.strip()}\n")
    return path


def relative_l2(a, b) -> float:
    """Relative L2 disagreement, nan-safe."""
    import numpy as np
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    denom = np.linalg.norm(b)
    if denom < 1e-300:
        denom = 1.0
    return float(np.linalg.norm(a - b) / denom)


def make_bias_sweep_cfg(cfg_base: dict, contact_name: str, V_target: float,
                        relaxed_snes: bool = True) -> dict:
    """Configure `cfg_base` for a `run_bias_sweep` ramp from 0 to V_target.

    `_resolve_sweep` reads the swept contact's `voltage_sweep` block
    rather than `solver.bias_ramp`; this helper rewrites the contact's
    `voltage_sweep` to land on V_target and overrides the contact's
    `voltage` field to match. Optionally injects the M12 (ADR 0008)
    relaxed SNES tolerances for configs whose solver block was authored
    for a tighter linear path (e.g., rc_ac_sweep).
    """
    cfg = copy.deepcopy(cfg_base)
    found = False
    for c in cfg.get("contacts", []):
        if c.get("name") != contact_name:
            continue
        c["voltage"] = float(V_target)
        if abs(float(V_target)) < 1e-12:
            c.pop("voltage_sweep", None)
        else:
            # Choose step so that V_target / step is an integer; this
            # avoids the rounding artifact in `_resolve_sweep` where a
            # `voltage_sweep.stop` that is not an exact multiple of
            # `.step` overshoots to the next step boundary.
            target_step = 0.05
            n = max(1, int(math.ceil(abs(float(V_target)) / target_step)))
            step = abs(float(V_target)) / n
            c["voltage_sweep"] = {
                "start": 0.0,
                "stop": float(V_target),
                "step": float(step),
            }
        found = True
        break
    if not found:
        raise KeyError(f"contact {contact_name!r} not found in cfg")

    solver = cfg.setdefault("solver", {})
    solver["type"] = "bias_sweep"
    if relaxed_snes:
        solver["snes"] = {
            "rtol": 1e-10,
            "atol": 1e-7,
            "stol": 1e-14,
            "max_it": 100,
        }
    # Drop AC-only sub-objects that are irrelevant for bias_sweep.
    for k in ("ac", "dc_bias"):
        solver.pop(k, None)
    return cfg


def final_field(result, name: str):
    """Return the final-state nodal array for `name` from any runner result.

    The runners disagree on output shape: `bias_sweep` and `equilibrium`
    return a `SimulationResult` with `psi_phys`/`n_phys`/`p_phys` arrays
    (one final state); `transient` exposes a `.fields` dict mapping
    field name to a list of snapshot arrays. This helper hides the
    difference.
    """
    import numpy as np
    # Transient-style: `.fields` dict of snapshot lists.
    fields = getattr(result, "fields", None)
    if fields:
        keys = [name]
        if name == "psi":
            keys.append("potential")
        for k in keys:
            if k in fields and len(fields[k]) > 0:
                return np.asarray(fields[k][-1])
    # SimulationResult-style: per-field arrays as attributes.
    attr_map = {"psi": "psi_phys", "n": "n_phys", "p": "p_phys"}
    attr = attr_map.get(name)
    if attr is not None and hasattr(result, attr):
        val = getattr(result, attr)
        if val is not None:
            return np.asarray(val)
    raise KeyError(
        f"field {name!r} not found on {type(result).__name__}; "
        f"fields keys = {list(fields) if fields else []}, "
        f"attrs = {[a for a in ('psi_phys','n_phys','p_phys') if hasattr(result,a)]}"
    )
