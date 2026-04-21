"""
Direct unit tests for `semi.mesh.build_mesh`.

The benchmarks exercise `build_mesh` transitively, but the 2D and 3D
builtin paths and the file-source NotImplementedError stub are not
hit by any pytest collection without these tests. They are cheap.
"""
from __future__ import annotations

import pytest

from semi.mesh import build_mesh


def _interval_cfg(N: int = 8):
    return {
        "dimension": 1,
        "mesh": {
            "source": "builtin",
            "extents": [[0.0, 1.0]],
            "resolution": [N],
            "facets_by_plane": [
                {"name": "left",  "tag": 1, "axis": 0, "value": 0.0},
                {"name": "right", "tag": 2, "axis": 0, "value": 1.0},
            ],
        },
    }


def _rectangle_cfg(Nx: int = 4, Ny: int = 3):
    return {
        "dimension": 2,
        "mesh": {
            "source": "builtin",
            "extents": [[0.0, 1.0], [0.0, 0.5]],
            "resolution": [Nx, Ny],
            "regions_by_box": [
                {"name": "low_y",  "tag": 1, "bounds": [[0.0, 1.0], [0.0, 0.25]]},
                {"name": "high_y", "tag": 2, "bounds": [[0.0, 1.0], [0.25, 0.5]]},
            ],
            "facets_by_plane": [
                {"name": "left",  "tag": 1, "axis": 0, "value": 0.0},
                {"name": "right", "tag": 2, "axis": 0, "value": 1.0},
            ],
        },
    }


def _box_cfg(Nx: int = 2, Ny: int = 2, Nz: int = 2):
    return {
        "dimension": 3,
        "mesh": {
            "source": "builtin",
            "extents": [[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]],
            "resolution": [Nx, Ny, Nz],
        },
    }


def test_build_mesh_interval_facets_tagged_at_endpoints():
    msh, cell_tags, facet_tags = build_mesh(_interval_cfg(N=10))
    assert msh.topology.dim == 1
    assert cell_tags is None
    # Both endpoint facets should be tagged.
    assert len(facet_tags.find(1)) == 1
    assert len(facet_tags.find(2)) == 1


def test_build_mesh_rectangle_creates_2d_mesh_with_region_and_facet_tags():
    cfg = _rectangle_cfg(Nx=4, Ny=4)
    msh, cell_tags, facet_tags = build_mesh(cfg)
    assert msh.topology.dim == 2
    # Both region tags should appear.
    assert len(cell_tags.find(1)) > 0
    assert len(cell_tags.find(2)) > 0
    # Left/right facets are non-empty (one cell column wide each).
    assert len(facet_tags.find(1)) > 0
    assert len(facet_tags.find(2)) > 0


def test_build_mesh_box_creates_3d_mesh():
    msh, cell_tags, facet_tags = build_mesh(_box_cfg())
    assert msh.topology.dim == 3
    # No regions / facets requested -> tags are None (early-return paths).
    assert cell_tags is None
    assert facet_tags is None


def test_build_mesh_unknown_source_raises():
    cfg = _interval_cfg()
    cfg["mesh"]["source"] = "magic"
    with pytest.raises(ValueError, match="Unknown mesh source"):
        build_mesh(cfg)


def test_build_mesh_file_source_raises_notimplemented():
    cfg = _interval_cfg()
    cfg["mesh"]["source"] = "file"
    with pytest.raises(NotImplementedError):
        build_mesh(cfg)


def test_build_mesh_dimension_mismatch_raises():
    cfg = _interval_cfg()
    # Declare 1D but pass 2D extents -> _build_builtin should reject.
    cfg["mesh"]["extents"] = [[0.0, 1.0], [0.0, 1.0]]
    with pytest.raises(ValueError, match="wrong dimension"):
        build_mesh(cfg)


def test_build_mesh_unsupported_dimension_raises():
    # 4D is impossible; the dispatcher should reject before dolfinx.
    cfg = {
        "dimension": 4,
        "mesh": {
            "source": "builtin",
            "extents": [[0.0, 1.0]] * 4,
            "resolution": [2] * 4,
        },
    }
    with pytest.raises(ValueError, match="Unsupported dimension"):
        build_mesh(cfg)
