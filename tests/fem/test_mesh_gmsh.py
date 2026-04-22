"""
Gmsh `.msh` loader (Day 7, Phase 3).

Exercises `semi.mesh.build_mesh` on the committed
`benchmarks/resistor_3d/fixtures/box.msh` fixture. The fixture is a
1 um x 200 nm x 200 nm tetrahedral bar with physical groups:

    "silicon"       (volume, tag 1)
    "contact_left"  (surface, tag 10) at x = 0
    "contact_right" (surface, tag 20) at x = L

The builtin-mesh equivalent uses `create_box` with resolution
[64, 16, 16] and tags the two x-faces via `facets_by_plane`. These
tests check that the gmsh path (a) loads, (b) returns a 3D mesh with
bounding box agreeing to 1e-12 m, and (c) produces the same two
contact-facet tag IDs as the JSON box/plane tagger convention so
downstream BC resolution does not care which mesh source was used.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from semi.mesh import build_mesh

FIXTURE = (
    Path(__file__).resolve().parents[2]
    / "benchmarks" / "resistor_3d" / "fixtures" / "box.msh"
)

L_expected = 1.0e-6
W_expected = 2.0e-7


def _file_cfg(path: Path) -> dict:
    return {
        "dimension": 3,
        "mesh": {
            "source": "file",
            "path": str(path),
            "format": "gmsh",
        },
    }


@pytest.mark.skipif(not FIXTURE.exists(), reason="gmsh fixture missing")
def test_gmsh_loader_returns_3d_mesh_with_correct_bbox():
    cfg = _file_cfg(FIXTURE)
    msh, cell_tags, facet_tags = build_mesh(cfg)
    assert msh.topology.dim == 3

    x = msh.geometry.x
    assert abs(x[:, 0].min() - 0.0) < 1e-12
    assert abs(x[:, 0].max() - L_expected) < 1e-12
    assert abs(x[:, 1].min() - 0.0) < 1e-12
    assert abs(x[:, 1].max() - W_expected) < 1e-12
    assert abs(x[:, 2].min() - 0.0) < 1e-12
    assert abs(x[:, 2].max() - W_expected) < 1e-12

    # Physical groups survive the round trip.
    assert cell_tags is not None
    assert facet_tags is not None
    assert set(np.unique(cell_tags.values).tolist()) == {1}
    assert set(np.unique(facet_tags.values).tolist()) == {10, 20}


@pytest.mark.skipif(not FIXTURE.exists(), reason="gmsh fixture missing")
def test_gmsh_loader_contact_facets_symmetric_and_on_the_right_planes():
    cfg = _file_cfg(FIXTURE)
    msh, _, facet_tags = build_mesh(cfg)

    left = facet_tags.find(10)
    right = facet_tags.find(20)
    assert left.size > 0
    assert right.size > 0
    # The two square faces have equal area, so gmsh-generated counts match.
    assert left.size == right.size

    # Spot-check the facet centroid x-coordinate for each contact.
    fdim = msh.topology.dim - 1
    msh.topology.create_connectivity(fdim, 0)
    f2v = msh.topology.connectivity(fdim, 0)
    coords = msh.geometry.x
    for facet in left:
        verts = f2v.links(int(facet))
        cx = coords[verts, 0].mean()
        assert abs(cx - 0.0) < 1e-12
    for facet in right:
        verts = f2v.links(int(facet))
        cx = coords[verts, 0].mean()
        assert abs(cx - L_expected) < 1e-12


@pytest.mark.skipif(not FIXTURE.exists(), reason="gmsh fixture missing")
def test_gmsh_loader_resolves_relative_path_against_source_dir():
    # Emulate the schema.load convention: cfg carries `_source_dir`
    # and `mesh.path` is a bare filename relative to that directory.
    cfg = {
        "dimension": 3,
        "mesh": {
            "source": "file",
            "path": "box.msh",
            "format": "gmsh",
        },
        "_source_dir": str(FIXTURE.parent),
    }
    msh, cell_tags, facet_tags = build_mesh(cfg)
    assert msh.topology.dim == 3
    assert cell_tags is not None and facet_tags is not None
    assert set(np.unique(facet_tags.values).tolist()) == {10, 20}
