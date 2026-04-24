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


def test_build_mesh_file_source_xdmf_missing_file_raises():
    """M12 wired XDMF through `dolfinx.io.XDMFFile`; a missing file
    surfaces as a dolfinx HDF5 error rather than a NotImplementedError.
    The XDMF happy-path is covered in `tests/fem/test_mesh_xdmf.py`."""
    cfg = _interval_cfg()
    cfg["mesh"] = {"source": "file", "path": "dummy.xdmf", "format": "xdmf"}
    with pytest.raises(RuntimeError):
        build_mesh(cfg)


def test_build_mesh_file_source_unknown_format_raises():
    cfg = _interval_cfg()
    cfg["mesh"] = {"source": "file", "path": "whatever.stl", "format": "stl"}
    with pytest.raises(ValueError, match="Unknown mesh file format"):
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


# ---------------------------------------------------------------------------
# Submesh + cellwise eps_r helpers (M6)
# ---------------------------------------------------------------------------

def _mos_like_cfg(Nx: int = 8, Ny: int = 10):
    """2D two-region mesh: 'silicon' (role semiconductor) + 'oxide' (insulator)."""
    return {
        "dimension": 2,
        "mesh": {
            "source": "builtin",
            "extents": [[0.0, 1.0e-7], [0.0, 1.0e-7]],
            "resolution": [Nx, Ny],
            "regions_by_box": [
                {"name": "silicon", "tag": 1,
                 "bounds": [[0.0, 1.0e-7], [0.0, 0.7e-7]]},
                {"name": "oxide", "tag": 2,
                 "bounds": [[0.0, 1.0e-7], [0.7e-7, 1.0e-7]]},
            ],
            "facets_by_plane": [
                {"name": "body", "tag": 1, "axis": 1, "value": 0.0},
                {"name": "gate", "tag": 2, "axis": 1, "value": 1.0e-7},
            ],
        },
        "regions": {
            "silicon": {"material": "Si", "tag": 1, "role": "semiconductor"},
            "oxide":   {"material": "SiO2", "tag": 2, "role": "insulator"},
        },
    }


def test_is_single_region_semiconductor_true_for_none():
    from semi.mesh import is_single_region_semiconductor
    assert is_single_region_semiconductor(None, {"silicon": {"material": "Si", "tag": 1, "role": "semiconductor"}})


def test_is_single_region_semiconductor_true_for_1d_pn():
    from semi.mesh import build_mesh, is_single_region_semiconductor

    cfg = {
        "dimension": 1,
        "mesh": {
            "source": "builtin",
            "extents": [[0.0, 2.0e-6]],
            "resolution": [8],
            "regions_by_box": [
                {"name": "silicon", "tag": 1, "bounds": [[0.0, 2.0e-6]]},
            ],
        },
    }
    msh, cell_tags, _facet_tags = build_mesh(cfg)
    regions_cfg = {"silicon": {"material": "Si", "tag": 1, "role": "semiconductor"}}
    assert is_single_region_semiconductor(cell_tags, regions_cfg)


def test_is_single_region_semiconductor_false_for_mos():
    from semi.mesh import build_mesh, is_single_region_semiconductor

    cfg = _mos_like_cfg()
    msh, cell_tags, _facet_tags = build_mesh(cfg)
    assert not is_single_region_semiconductor(cell_tags, cfg["regions"])


def test_build_eps_r_function_assigns_per_region_values():
    from semi.mesh import build_eps_r_function, build_mesh

    cfg = _mos_like_cfg()
    msh, cell_tags, _facet_tags = build_mesh(cfg)
    eps_r_fn = build_eps_r_function(msh, cell_tags, cfg["regions"])

    from dolfinx import fem

    V_DG0 = fem.functionspace(msh, ("DG", 0))
    # Pull each region's DOF values and check they match the material.
    si_cells = cell_tags.find(1)
    ox_cells = cell_tags.find(2)
    for c in si_cells:
        dof = V_DG0.dofmap.cell_dofs(int(c))[0]
        assert abs(eps_r_fn.x.array[dof] - 11.7) < 1.0e-12
    for c in ox_cells:
        dof = V_DG0.dofmap.cell_dofs(int(c))[0]
        assert abs(eps_r_fn.x.array[dof] - 3.9) < 1.0e-12


def test_build_submesh_by_role_returns_semiconductor_cells_only():
    from semi.mesh import build_mesh, build_submesh_by_role

    cfg = _mos_like_cfg()
    msh, cell_tags, _facet_tags = build_mesh(cfg)
    num_si_cells = len(cell_tags.find(1))

    import numpy as _np

    submesh, entity_map, _vmap, _gmap = build_submesh_by_role(
        msh, cell_tags, cfg["regions"], role="semiconductor",
    )
    tdim = submesh.topology.dim
    submesh.topology.create_entities(tdim)
    n_local = submesh.topology.index_map(tdim).size_local
    assert n_local == num_si_cells
    # entity_map maps submesh cells to their parent cells; every mapped
    # parent cell must carry the 'semiconductor' tag.
    sub_ids = _np.arange(n_local, dtype=_np.int32)
    parent_ids = entity_map.sub_topology_to_topology(sub_ids, inverse=False)
    parent_si_cells = set(int(c) for c in cell_tags.find(1))
    for pid in parent_ids:
        assert int(pid) in parent_si_cells


def test_build_submesh_by_role_raises_without_cell_tags():
    import pytest as _pytest

    from semi.mesh import build_mesh, build_submesh_by_role

    cfg = {
        "dimension": 1,
        "mesh": {
            "source": "builtin",
            "extents": [[0.0, 1.0]],
            "resolution": [4],
        },
    }
    msh, cell_tags, _facet_tags = build_mesh(cfg)
    assert cell_tags is None
    with _pytest.raises(ValueError, match="requires cell_tags"):
        build_submesh_by_role(msh, cell_tags, {"silicon": {"material": "Si"}})


def test_build_submesh_by_role_raises_when_no_matching_role():
    import pytest as _pytest

    from semi.mesh import build_mesh, build_submesh_by_role

    cfg = _mos_like_cfg()
    msh, cell_tags, _facet_tags = build_mesh(cfg)
    with _pytest.raises(ValueError, match="No region in regions_cfg has role"):
        build_submesh_by_role(msh, cell_tags, cfg["regions"], role="conductor")
