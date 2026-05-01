"""
XDMF mesh ingest round-trip test (M14.3, Phase C).

The committed gmsh fixture `benchmarks/resistor_3d/fixtures/box.msh`
is loaded via the gmsh path, then re-written to XDMF in a temporary
directory using `dolfinx.io.XDMFFile.write_mesh` plus `write_meshtags`,
then loaded back via the new XDMF branch in `semi.mesh._build_from_file`.
The two load paths must:

1. Yield meshes with bounding-boxes that agree to 1e-12 m.
2. Reproduce the same physical-group facet tags (so downstream BCs
   resolve identically).
3. Produce identical analytical resistance R for the resistor benchmark
   within 1e-12 relative when both meshes are fed through `build_mesh`
   and the V-I sweep runner.

The dolfinx round-trip (gmsh -> XDMFFile.write -> XDMFFile.read) is
the reference path; meshio is not required.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from semi.mesh import build_mesh

REPO_ROOT = Path(__file__).resolve().parents[2]
GMSH_FIXTURE = (
    REPO_ROOT / "benchmarks" / "resistor_3d" / "fixtures" / "box.msh"
)
RESISTOR_GMSH_JSON = REPO_ROOT / "benchmarks" / "resistor_3d" / "resistor_gmsh.json"


def _file_cfg(path: Path, fmt: str = "gmsh") -> dict:
    return {
        "dimension": 3,
        "mesh": {
            "source": "file",
            "path": str(path),
            "format": fmt,
        },
    }


@pytest.fixture(scope="module")
def xdmf_fixture(tmp_path_factory):
    """Convert the committed box.msh to a single-file XDMF in a tmp dir.

    Uses dolfinx XDMFFile to write the mesh + cell tags + facet tags
    into one .xdmf (which produces the .xdmf XML descriptor plus a
    paired .h5 datafile).
    """
    if not GMSH_FIXTURE.exists():
        pytest.skip("gmsh fixture missing")

    from dolfinx.io import XDMFFile

    tmp = tmp_path_factory.mktemp("xdmf_fixture")
    out_xdmf = tmp / "box.xdmf"

    cfg = _file_cfg(GMSH_FIXTURE, fmt="gmsh")
    msh, ct, ft = build_mesh(cfg)

    with XDMFFile(msh.comm, str(out_xdmf), "w") as xdmf:
        xdmf.write_mesh(msh)
        if ct is not None:
            ct.name = "cell_tags"
            xdmf.write_meshtags(ct, msh.geometry)
        if ft is not None:
            ft.name = "facet_tags"
            xdmf.write_meshtags(ft, msh.geometry)
    return out_xdmf


def test_xdmf_loader_returns_3d_mesh_with_correct_bbox(xdmf_fixture):
    cfg = _file_cfg(xdmf_fixture, fmt="xdmf")
    msh, _ct, _ft = build_mesh(cfg)
    assert msh.topology.dim == 3
    coords = msh.geometry.x
    L_x = float(coords[:, 0].max() - coords[:, 0].min())
    W_y = float(coords[:, 1].max() - coords[:, 1].min())
    W_z = float(coords[:, 2].max() - coords[:, 2].min())
    assert L_x == pytest.approx(1.0e-6, abs=1e-12)
    assert W_y == pytest.approx(2.0e-7, abs=1e-12)
    assert W_z == pytest.approx(2.0e-7, abs=1e-12)


def test_xdmf_loader_preserves_facet_tags(xdmf_fixture):
    msh_g, _ctg, ft_g = build_mesh(_file_cfg(GMSH_FIXTURE, fmt="gmsh"))
    msh_x, _ctx, ft_x = build_mesh(_file_cfg(xdmf_fixture, fmt="xdmf"))

    assert ft_g is not None
    assert ft_x is not None
    tags_g = sorted({int(t) for t in ft_g.values})
    tags_x = sorted({int(t) for t in ft_x.values})
    assert tags_g == tags_x
    # The fixture declares physical surfaces 10 and 20 for the two
    # contact faces (see benchmarks/resistor_3d/fixtures/box.geo).
    assert tags_g == [10, 20]


@pytest.mark.skipif(
    not RESISTOR_GMSH_JSON.exists(),
    reason="resistor_gmsh.json not present",
)
def test_xdmf_resistor_R_matches_msh_within_1e_minus_12(xdmf_fixture):
    """Acceptance test 3: run the resistor sweep on both load paths and
    compare R. The XDMF path must agree with the gmsh path to 1e-12
    relative.
    """
    from semi.runners.bias_sweep import run_bias_sweep
    from semi.schema import load as schema_load

    cfg_msh = schema_load(str(RESISTOR_GMSH_JSON))
    cfg_xdmf = json.loads(json.dumps(cfg_msh))
    cfg_xdmf["mesh"]["path"] = str(xdmf_fixture)
    cfg_xdmf["mesh"]["format"] = "xdmf"

    res_msh = run_bias_sweep(cfg_msh)
    res_xdmf = run_bias_sweep(cfg_xdmf)

    iv_msh = sorted(res_msh.iv, key=lambda r: r["V"])
    iv_xdmf = sorted(res_xdmf.iv, key=lambda r: r["V"])
    assert len(iv_msh) == len(iv_xdmf)

    # Use the highest-magnitude bias point for the relative check; near
    # zero bias the current is near zero and the relative norm is
    # numerically meaningless.
    V_msh = np.array([float(r["V"]) for r in iv_msh])
    J_msh = np.array([float(r["J"]) for r in iv_msh])
    J_xdmf = np.array([float(r["J"]) for r in iv_xdmf])
    idx = int(np.argmax(np.abs(V_msh)))
    R_msh = V_msh[idx] / J_msh[idx]
    R_xdmf = V_msh[idx] / J_xdmf[idx]
    rel_err = abs(R_xdmf - R_msh) / abs(R_msh)
    assert rel_err < 1.0e-12, (
        f"XDMF R does not match gmsh R within 1e-12 relative: "
        f"R_msh={R_msh:.6e}, R_xdmf={R_xdmf:.6e}, rel_err={rel_err:.3e}"
    )


def test_xdmf_loader_returns_none_tags_when_meshtags_absent(tmp_path):
    """XDMF file without meshtags yields (msh, None, None).

    Exercises the ``except RuntimeError: ... = None`` branches in
    ``semi.mesh._build_from_file`` for the XDMF path.
    """
    from dolfinx.io import XDMFFile
    from mpi4py import MPI

    # Write a mesh-only XDMF (no meshtags) from scratch.
    msh_src, _ct, _ft = build_mesh({
        "dimension": 1,
        "mesh": {
            "source": "builtin",
            "extents": [[0.0, 1.0e-6]],
            "resolution": [4],
        },
    })
    xdmf_path = tmp_path / "mesh_only.xdmf"
    with XDMFFile(MPI.COMM_WORLD, str(xdmf_path), "w") as xf:
        xf.write_mesh(msh_src)
    # Read back: cell_tags and facet_tags should be None.
    msh2, cell_tags, facet_tags = build_mesh({
        "dimension": 1,
        "mesh": {
            "source": "file",
            "path": str(xdmf_path),
            "format": "xdmf",
        },
    })
    assert cell_tags is None
    assert facet_tags is None
    assert msh2 is not None
