"""
Tests for the M9 result artifact writer.

FEM-heavy tests (requiring dolfinx) are skipped when dolfinx is
unavailable, matching the pattern in tests/fem/conftest.py.
Pure-Python tests (schema validation round-trip) always run.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import pytest

try:
    import dolfinx  # noqa: F401
    HAS_DOLFINX = True
except ImportError:
    HAS_DOLFINX = False

REPO_ROOT = Path(__file__).parent.parent
SCHEMAS_DIR = REPO_ROOT / "schemas"
BENCHMARKS_DIR = REPO_ROOT / "benchmarks"

BENCHMARK_JSONS = [
    BENCHMARKS_DIR / "pn_1d" / "pn_junction.json",
    BENCHMARKS_DIR / "pn_1d_bias" / "pn_junction_bias.json",
    BENCHMARKS_DIR / "pn_1d_bias_reverse" / "pn_junction_bias_reverse.json",
    BENCHMARKS_DIR / "mos_2d" / "mos_cap.json",
    BENCHMARKS_DIR / "resistor_3d" / "resistor.json",
]


# ---- Pure-Python tests (always run) ----

def test_manifest_schema_valid_draft7():
    import jsonschema
    schema = json.loads((SCHEMAS_DIR / "manifest.v1.json").read_text())
    jsonschema.Draft7Validator.check_schema(schema)


def test_manifest_schema_round_trip():
    """A minimal valid manifest validates without error."""
    import jsonschema
    schema = json.loads((SCHEMAS_DIR / "manifest.v1.json").read_text())
    manifest = {
        "schema_version": "1.0.0",
        "engine": {"name": "kronos-semi", "version": "0.9.0", "commit": "abc1234abc1234"},
        "run_id": "2026-04-23T12-00-00Z_test_abc1234",
        "status": "completed",
        "wall_time_s": 1.5,
        "input_sha256": "a" * 64,
        "solver": {
            "type": "equilibrium",
            "backend": "petsc-cpu-mumps",
            "n_dofs": 401,
            "n_steps": 1,
            "converged": True,
        },
        "fields": [
            {"name": "psi", "units": "V", "path": "fields/psi.bp", "rank": 0},
        ],
        "mesh": {
            "path": "mesh/mesh.xdmf",
            "dimension": 1,
            "n_cells": 400,
            "n_vertices": 401,
            "regions": [{"name": "silicon", "material": "Si", "tag": 1}],
        },
        "warnings": [],
    }
    jsonschema.Draft7Validator(schema).validate(manifest)


def test_manifest_schema_optional_sweeps():
    """Manifest with sweeps validates without error."""
    import jsonschema
    schema = json.loads((SCHEMAS_DIR / "manifest.v1.json").read_text())
    manifest = {
        "schema_version": "1.0.0",
        "engine": {"name": "kronos-semi", "version": "0.9.0", "commit": "abc1234"},
        "run_id": "2026-04-23T12-00-00Z_test_abc1234",
        "status": "completed",
        "wall_time_s": 5.0,
        "input_sha256": "b" * 64,
        "solver": {
            "type": "bias_sweep",
            "backend": "petsc-cpu-mumps",
            "n_dofs": 801,
            "n_steps": 13,
            "converged": True,
        },
        "fields": [
            {"name": "psi", "units": "V", "path": "fields/psi.bp", "rank": 0},
        ],
        "mesh": {
            "path": "mesh/mesh.xdmf",
            "dimension": 1,
            "n_cells": 800,
            "n_vertices": 801,
            "regions": [{"name": "silicon", "material": "Si", "tag": 1}],
        },
        "sweeps": [
            {
                "kind": "voltage",
                "contact": "anode",
                "path": "iv/anode.csv",
                "n_steps": 13,
                "bipolar": False,
            }
        ],
        "warnings": ["test warning"],
    }
    jsonschema.Draft7Validator(schema).validate(manifest)


# ---- FEM-heavy tests (require dolfinx) ----

@pytest.mark.skipif(not HAS_DOLFINX, reason="dolfinx not available")
@pytest.mark.parametrize(
    "bench_json",
    BENCHMARK_JSONS,
    ids=[p.parent.name for p in BENCHMARK_JSONS],
)
def test_write_read_artifact(bench_json, tmp_path):
    """Write artifact for a benchmark and assert manifest validity."""
    import jsonschema

    from semi.io.artifact import write_artifact
    from semi.io.reader import read_manifest
    from semi.run import run
    from semi.schema import load as schema_load

    schema = json.loads((SCHEMAS_DIR / "manifest.v1.json").read_text())

    cfg = schema_load(str(bench_json))
    result = run(cfg)
    run_dir = write_artifact(result, tmp_path / "runs", input_json_path=bench_json)

    # 1. Manifest validates against the schema
    manifest = read_manifest(run_dir)
    jsonschema.Draft7Validator(schema).validate(manifest)

    # 2. All field files listed in manifest exist on disk
    for field_entry in manifest["fields"]:
        field_path = run_dir / field_entry["path"]
        assert field_path.exists(), f"Field file not found: {field_entry['path']}"

    # 3. Sweep CSVs have n_steps + 1 rows (header + data)
    for sweep in manifest.get("sweeps", []):
        sweep_path = run_dir / sweep["path"]
        assert sweep_path.exists(), f"Sweep CSV not found: {sweep['path']}"
        with open(sweep_path, newline="") as f:
            rows = list(csv.reader(f))
        assert len(rows) == sweep["n_steps"] + 1, (
            f"Expected {sweep['n_steps'] + 1} rows in {sweep['path']}, got {len(rows)}"
        )

    # 4. input_sha256 matches the actual file
    expected = hashlib.sha256(bench_json.read_bytes()).hexdigest()
    assert manifest["input_sha256"] == expected, "input_sha256 mismatch"


@pytest.mark.skipif(not HAS_DOLFINX, reason="dolfinx not available")
def test_write_artifact_cfg_sha_when_no_input_path(tmp_path):
    """When input_json_path is omitted, input_sha256 is derived from result.cfg."""
    from semi.io.artifact import write_artifact
    from semi.run import run
    from semi.schema import load as schema_load

    bench_json = BENCHMARKS_DIR / "pn_1d" / "pn_junction.json"
    cfg = schema_load(str(bench_json))
    result = run(cfg)
    run_dir = write_artifact(result, tmp_path / "runs", input_json_path=None)

    manifest = json.loads((run_dir / "manifest.json").read_text())
    expected = hashlib.sha256(json.dumps(result.cfg, sort_keys=True).encode()).hexdigest()
    assert manifest["input_sha256"] == expected
