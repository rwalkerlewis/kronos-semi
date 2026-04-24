"""
XDMF mesh loader (M12).

Exercises the XDMF branch of `semi.mesh._build_from_file`, which
previously raised `NotImplementedError`. The round-trip goes through
the committed gmsh fixture used by M7: load the `.msh`, write it back
as XDMF with `cell_tags` and `facet_tags`, then re-load through
`build_mesh(format="xdmf")` and compare tag counts.
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


@pytest.mark.skipif(not FIXTURE.exists(), reason="gmsh fixture missing")
def test_xdmf_roundtrip_preserves_cell_and_facet_tags(tmp_path):
    # Load the reference gmsh mesh and its tags.
    cfg_gmsh = {
        "dimension": 3,
        "mesh": {
            "source": "file",
            "path": str(FIXTURE),
            "format": "gmsh",
        },
    }
    msh_in, cell_tags_in, facet_tags_in = build_mesh(cfg_gmsh)

    # Write it back out as XDMF with named meshtags.
    from dolfinx.io import XDMFFile

    xdmf_path = tmp_path / "mesh.xdmf"
    with XDMFFile(msh_in.comm, str(xdmf_path), "w") as f:
        f.write_mesh(msh_in)
        if cell_tags_in is not None:
            cell_tags_in.name = "cell_tags"
            f.write_meshtags(cell_tags_in, msh_in.geometry)
        if facet_tags_in is not None:
            # The XDMF writer requires the parent topology connectivity
            # to exist before writing the facet tags.
            msh_in.topology.create_connectivity(msh_in.topology.dim - 1, msh_in.topology.dim)
            facet_tags_in.name = "facet_tags"
            f.write_meshtags(facet_tags_in, msh_in.geometry)

    # Re-load through the XDMF branch.
    cfg_xdmf = {
        "dimension": 3,
        "mesh": {
            "source": "file",
            "path": str(xdmf_path),
            "format": "xdmf",
        },
    }
    msh_out, cell_tags_out, facet_tags_out = build_mesh(cfg_xdmf)

    assert msh_out.topology.dim == 3
    assert cell_tags_out is not None, "XDMF cell_tags not re-read"
    assert facet_tags_out is not None, "XDMF facet_tags not re-read"

    # Tag-value sets survive the round trip.
    assert set(np.unique(cell_tags_in.values).tolist()) == set(
        np.unique(cell_tags_out.values).tolist()
    )
    assert set(np.unique(facet_tags_in.values).tolist()) == set(
        np.unique(facet_tags_out.values).tolist()
    )

    # Per-tag facet counts (the downstream BC resolver relies on these).
    for tag in np.unique(facet_tags_in.values):
        n_in = int((facet_tags_in.values == tag).sum())
        n_out = int((facet_tags_out.values == tag).sum())
        assert n_in == n_out, (
            f"facet tag {tag}: {n_in} in fixture vs {n_out} after XDMF round-trip"
        )
