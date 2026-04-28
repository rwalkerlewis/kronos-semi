"""
Shared utilities for the audit suite.

Every audit test consumes one or more existing benchmark JSON configs,
runs two or more runners on a shared problem, computes a comparison
metric, and writes the results to `/tmp/audit/<case>.csv` and a
markdown fragment to `/tmp/audit/<case>.md`. The conftest aggregator
collects the markdown fragments into `docs/PHYSICS_AUDIT.md`.
"""
from __future__ import annotations

import csv
import importlib
import json
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
