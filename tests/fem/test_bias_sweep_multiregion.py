"""
Tests for the multi-region dispatch in ``semi.runners.bias_sweep``.

The bias_sweep runner chooses between two code paths at config-inspection
time: single-region (pn junction, resistor, one-material body) and
multi-region (semiconductor + insulator mix, e.g. MOSFET / HBT with
isolation oxide). These tests pin the dispatch decision and exercise the
multi-region machinery against a minimal 2D p-body + oxide config that
runs quickly but still touches every MR-specific path:

  - region-role detection
  - ``build_submesh_by_role`` for carrier spaces
  - cellwise ``eps_r`` DG0 assembly
  - ``build_dd_block_residual_mr`` + ``entity_maps`` block solve
  - parent-facet -> submesh-facet BC translation

The tests intentionally use a shared Si/SiO2 interface (built-in
``regions_by_box``) rather than a gmsh geometry, so they do not require
an external mesh subprocess and run under every CI matrix.
"""
from __future__ import annotations

import copy

import numpy as np
import pytest

_BASE_CFG = {
    "schema_version": "1.1.0",
    "name": "mr_dispatch_test",
    "description": "Minimal p-body + oxide 2D cap for bias_sweep MR dispatch",
    "dimension": 2,
    "mesh": {
        "source": "builtin",
        "extents": [[0.0, 1.0e-6], [0.0, 1.05e-6]],
        "resolution": [8, 21],
        "regions_by_box": [
            {"name": "silicon", "tag": 1, "bounds": [[0.0, 1.0e-6], [0.0, 1.0e-6]]},
            {"name": "oxide",   "tag": 2, "bounds": [[0.0, 1.0e-6], [1.0e-6, 1.05e-6]]},
        ],
        "facets_by_plane": [
            {"name": "body", "tag": 1, "axis": 1, "value": 0.0},
            {"name": "gate", "tag": 2, "axis": 1, "value": 1.05e-6},
        ],
    },
    "regions": {
        "silicon": {"material": "Si",   "tag": 1, "role": "semiconductor"},
        "oxide":   {"material": "SiO2", "tag": 2, "role": "insulator"},
    },
    "doping": [
        {"region": "silicon",
         "profile": {"type": "uniform", "N_D": 0.0, "N_A": 1.0e17}},
    ],
    "contacts": [
        {"name": "body", "facet": "body", "type": "ohmic", "voltage": 0.0},
        {"name": "gate", "facet": "gate", "type": "gate", "voltage": 0.0,
         "workfunction": 0.0},
    ],
    "physics": {
        "temperature": 300.0,
        "statistics": "boltzmann",
        "mobility": {"mu_n": 1400.0, "mu_p": 450.0},
        "recombination": {"srh": False, "tau_n": 1.0e-7, "tau_p": 1.0e-7, "E_t": 0.0},
    },
    "solver": {
        "type": "bias_sweep",
        "continuation": {
            "min_step": 1.0e-4, "max_halvings": 6, "max_step": 0.1,
            "easy_iter_threshold": 4, "grow_factor": 1.5,
        },
    },
    "output": {"directory": "./results/mr_dispatch_test",
               "fields": ["potential", "n", "p"]},
}


def _make_cfg(**overrides):
    cfg = copy.deepcopy(_BASE_CFG)
    cfg.update(overrides)
    return cfg


def test_multi_region_role_detection_for_mosfet_shape_config():
    """A config with both 'semiconductor' and 'insulator' roles triggers MR."""
    cfg = _make_cfg()
    regions = cfg["regions"]
    roles = {r.get("role", "semiconductor") for r in regions.values()}
    assert "insulator" in roles
    assert "semiconductor" in roles


def test_single_region_role_detection_for_resistor_shape_config():
    """A config with only 'semiconductor' roles stays on the SR path."""
    cfg = _make_cfg(regions={
        "bulk": {"material": "Si", "tag": 1, "role": "semiconductor"}
    })
    regions = cfg["regions"]
    roles = {r.get("role", "semiconductor") for r in regions.values()}
    assert "insulator" not in roles


def test_map_parent_facets_to_submesh_drops_oxide_only_facets():
    """Parent facets whose vertices are not in the semiconductor submesh get dropped."""
    from dolfinx import mesh as dmesh
    from mpi4py import MPI

    from semi.mesh import build_mesh, build_submesh_by_role, map_parent_facets_to_submesh
    from semi import schema

    cfg = _make_cfg()
    cfg = schema.validate(cfg)
    msh, cell_tags, facet_tags = build_mesh(cfg)
    submesh, _em_cell, em_vert, _ = build_submesh_by_role(
        msh, cell_tags, cfg["regions"], role="semiconductor",
    )

    # Body facet (silicon bottom) should map fully; gate facet (oxide top)
    # should map to zero submesh facets since it's not part of silicon.
    body_parent = facet_tags.find(1)
    body_sub = map_parent_facets_to_submesh(msh, submesh, em_vert, body_parent)
    assert len(body_sub) == len(body_parent), (
        f"body facets must transfer completely to submesh: "
        f"parent={len(body_parent)}, submesh={len(body_sub)}"
    )

    gate_parent = facet_tags.find(2)
    gate_sub = map_parent_facets_to_submesh(msh, submesh, em_vert, gate_parent)
    assert len(gate_sub) == 0, (
        f"gate (oxide-top) facets must NOT appear in the silicon submesh: "
        f"got {len(gate_sub)}"
    )


def test_bias_sweep_mr_dispatch_runs_to_completion():
    """End-to-end: MR config runs through bias_sweep and returns a non-empty IV."""
    from semi import schema
    from semi.runners.bias_sweep import run_bias_sweep

    cfg = _make_cfg(name="mr_bias_sweep_end_to_end")
    cfg = schema.validate(cfg)

    result = run_bias_sweep(cfg)

    # With no ohmic voltage_sweep declared, the run produces a single seed
    # solve at V=0 (see _resolve_sweep). The MR dispatch should still
    # produce a converged solve and a one-row IV.
    assert result.iv, "MR bias_sweep must produce at least one IV row"
    assert result.solver_info.get("converged", False), (
        f"MR bias_sweep seed solve did not converge: {result.solver_info}"
    )
    # phi_n / phi_p live on the silicon submesh, so their DOF counts are
    # strictly smaller than V_psi's parent-mesh DOF count.
    n_psi = int(result.psi.x.array.size)
    n_phi_n = int(result.phi_n.x.array.size)
    assert n_phi_n < n_psi, (
        f"MR path must put phi_n on a submaller submesh: "
        f"n_psi={n_psi}, n_phi_n={n_phi_n}"
    )


def test_bias_sweep_single_region_path_unchanged_on_pure_semi_config():
    """A pure-semiconductor config keeps the SR path (phi_n on full mesh)."""
    from semi import schema
    from semi.runners.bias_sweep import run_bias_sweep

    cfg = _make_cfg(
        name="sr_bias_sweep_pure_semi",
        regions={"bulk": {"material": "Si", "tag": 1, "role": "semiconductor"}},
        mesh={
            "source": "builtin",
            "extents": [[0.0, 1.0e-6], [0.0, 1.0e-6]],
            "resolution": [8, 8],
            "regions_by_box": [
                {"name": "bulk", "tag": 1,
                 "bounds": [[0.0, 1.0e-6], [0.0, 1.0e-6]]},
            ],
            "facets_by_plane": [
                {"name": "left",  "tag": 1, "axis": 0, "value": 0.0},
                {"name": "right", "tag": 2, "axis": 0, "value": 1.0e-6},
            ],
        },
        contacts=[
            {"name": "left",  "facet": "left",  "type": "ohmic", "voltage": 0.0},
            {"name": "right", "facet": "right", "type": "ohmic", "voltage": 0.0},
        ],
    )
    cfg = schema.validate(cfg)

    result = run_bias_sweep(cfg)
    assert result.iv
    # In the SR path, phi_n lives on the same mesh as psi -> equal DOF count.
    assert result.psi.x.array.size == result.phi_n.x.array.size, (
        "SR path must put phi_n on the same mesh as psi"
    )
