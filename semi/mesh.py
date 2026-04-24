"""
Mesh construction and tagging.

Supports:
    - Builtin meshes: interval (1D), rectangle (2D), box (3D)
      with region tagging via axis-aligned boxes and facet tagging via planes
    - File-based meshes: gmsh .msh via `dolfinx.io.gmsh.read_from_msh`
      and XDMF via `dolfinx.io.XDMFFile`. Physical groups stored in the
      file are returned directly as cell_tags and facet_tags; the JSON
      box/plane tagger is bypassed.
    - Parametric `geometry` input: a `.geo` file plus physical-group map
      (M12). The `.geo` is handed to gmsh in a subprocess via
      `semi.geometry.realize`, the resulting `.msh` is cached by content
      hash, then the cached file is loaded via the gmsh path above.

The output of `build_mesh` is a tuple (mesh, cell_tags, facet_tags) where
the tag meshtag objects follow the convention used throughout dolfinx.
Integer tags are assigned from the JSON config's 'tag' fields for builtin
meshes, and from the .msh physical-group IDs for file-source meshes.

This module requires dolfinx and is imported lazily. The non-FEM code in
this package (schema, materials, doping, scaling) does NOT import this.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np


def build_mesh(cfg: dict[str, Any]):
    """
    Build a dolfinx Mesh with cell tags and facet tags from a validated config.

    Returns
    -------
    msh : dolfinx.mesh.Mesh
    cell_tags : dolfinx.mesh.MeshTags (on cells), or None if no regions_by_box
    facet_tags : dolfinx.mesh.MeshTags (on facets), or None if no facets_by_plane
    """
    # Parametric `geometry` input (M12): realize via gmsh subprocess and
    # synthesize a file-source mesh block so the rest of the pipeline
    # reuses the existing gmsh-loader code path unchanged.
    if "geometry" in cfg and "mesh" not in cfg:
        from . import geometry as _geom

        geo_cfg = cfg["geometry"]
        source_dir = cfg.get("_source_dir")
        msh_path = _geom.realize(
            geo_cfg,
            source_dir=Path(source_dir) if source_dir is not None else Path("."),
            dimension=int(cfg["dimension"]),
        )
        synthetic_mesh_cfg = {
            "source": "file",
            "format": "gmsh",
            "path": str(msh_path),
        }
        return _build_from_file(
            synthetic_mesh_cfg,
            dim=int(cfg["dimension"]),
            source_dir=None,
        )

    mesh_cfg = cfg["mesh"]
    source = mesh_cfg["source"]

    if source == "builtin":
        msh = _build_builtin(mesh_cfg, cfg["dimension"])
        cell_tags = _tag_regions(msh, mesh_cfg.get("regions_by_box", []))
        facet_tags = _tag_facets(msh, mesh_cfg.get("facets_by_plane", []))
        return msh, cell_tags, facet_tags
    if source == "file":
        return _build_from_file(
            mesh_cfg,
            dim=int(cfg["dimension"]),
            source_dir=cfg.get("_source_dir"),
        )
    raise ValueError(f"Unknown mesh source {source!r}")


def _build_builtin(mesh_cfg: dict, dim: int):
    """Create an interval / rectangle / box mesh from extents and resolution."""
    from dolfinx import mesh as dmesh
    from mpi4py import MPI

    extents = mesh_cfg["extents"]
    res = mesh_cfg["resolution"]
    if len(extents) != dim or len(res) != dim:
        raise ValueError(
            f"Mesh extents/resolution have wrong dimension: "
            f"expected {dim}, got extents={len(extents)}, resolution={len(res)}"
        )

    comm = MPI.COMM_WORLD
    if dim == 1:
        return dmesh.create_interval(comm, res[0], [extents[0][0], extents[0][1]])
    if dim == 2:
        return dmesh.create_rectangle(
            comm,
            [np.array([extents[0][0], extents[1][0]]),
             np.array([extents[0][1], extents[1][1]])],
            [res[0], res[1]],
        )
    if dim == 3:
        return dmesh.create_box(
            comm,
            [np.array([extents[0][0], extents[1][0], extents[2][0]]),
             np.array([extents[0][1], extents[1][1], extents[2][1]])],
            [res[0], res[1], res[2]],
        )
    raise ValueError(f"Unsupported dimension {dim}")


def _build_from_file(mesh_cfg: dict, dim: int, source_dir: str | None = None):
    """
    Load a mesh from disk.

    Supports gmsh `.msh` via `dolfinx.io.gmsh.read_from_msh`; physical
    groups stored in the file are returned verbatim as `cell_tags` and
    `facet_tags`, so `build_mesh` skips the box/plane tagger for
    file-source meshes. XDMF is reserved for a future PR.

    The `path` field is resolved relative to the JSON source directory
    when it is not absolute, matching the convention used by other
    file-relative references loaded via `semi.schema.load`.
    """
    fmt = mesh_cfg.get("format", "gmsh")
    raw_path = Path(mesh_cfg["path"])
    if not raw_path.is_absolute() and source_dir is not None:
        raw_path = Path(source_dir) / raw_path

    if fmt == "gmsh":
        from dolfinx.io import gmsh as _gmsh_io
        from mpi4py import MPI

        meshdata = _gmsh_io.read_from_msh(
            str(raw_path), MPI.COMM_WORLD, gdim=dim,
        )
        return meshdata.mesh, meshdata.cell_tags, meshdata.facet_tags
    if fmt == "xdmf":
        from dolfinx.io import XDMFFile
        from mpi4py import MPI

        with XDMFFile(MPI.COMM_WORLD, str(raw_path), "r") as f:
            msh = f.read_mesh(name="mesh")
            try:
                cell_tags = f.read_meshtags(msh, name="cell_tags")
            except RuntimeError:
                cell_tags = None
            try:
                msh.topology.create_connectivity(msh.topology.dim - 1, msh.topology.dim)
                facet_tags = f.read_meshtags(msh, name="facet_tags")
            except RuntimeError:
                facet_tags = None
        return msh, cell_tags, facet_tags
    raise ValueError(f"Unknown mesh file format {fmt!r}")


def _tag_regions(msh, boxes: list[dict]):
    """Tag cells whose centroids lie inside one of the given axis-aligned boxes."""
    if not boxes:
        return None
    from dolfinx import mesh as dmesh

    tdim = msh.topology.dim
    msh.topology.create_entities(tdim)
    num_cells = msh.topology.index_map(tdim).size_local + msh.topology.index_map(tdim).num_ghosts
    all_cells = np.arange(num_cells, dtype=np.int32)

    # Centroid-based classification
    centroids = _cell_centroids(msh, all_cells)

    values = np.zeros(num_cells, dtype=np.int32)
    # Later boxes override earlier ones (intentional: allow nesting)
    for box in boxes:
        tag = int(box["tag"])
        bounds = np.asarray(box["bounds"], dtype=float)  # shape (dim, 2)
        inside = np.ones(num_cells, dtype=bool)
        for d in range(tdim):
            inside &= (centroids[:, d] >= bounds[d, 0]) & (centroids[:, d] <= bounds[d, 1])
        values[inside] = tag

    # Build meshtags for tagged cells only
    tagged = np.where(values != 0)[0].astype(np.int32)
    vals = values[tagged]
    # Sort by entity index (required by meshtags)
    order = np.argsort(tagged)
    tagged = tagged[order]
    vals = vals[order]
    return dmesh.meshtags(msh, tdim, tagged, vals)


def _tag_facets(msh, planes: list[dict]):
    """Tag boundary facets lying on the given axis-aligned planes."""
    if not planes:
        return None
    from dolfinx import mesh as dmesh

    tdim = msh.topology.dim
    fdim = tdim - 1

    all_tags = []
    all_facets = []

    for plane in planes:
        tag = int(plane["tag"])
        axis = int(plane["axis"])
        value = float(plane["value"])
        tol = float(plane.get("tol", 1.0e-12))

        facets = dmesh.locate_entities_boundary(
            msh, fdim,
            lambda x, a=axis, v=value, t=tol: np.isclose(x[a], v, atol=t)
        )
        all_facets.append(facets)
        all_tags.append(np.full(len(facets), tag, dtype=np.int32))

    if not all_facets:
        return None
    facets = np.concatenate(all_facets).astype(np.int32)
    tags   = np.concatenate(all_tags).astype(np.int32)

    # meshtags wants sorted, unique entities
    order = np.argsort(facets)
    facets = facets[order]
    tags   = tags[order]

    # Drop duplicates (keep first occurrence)
    _, first = np.unique(facets, return_index=True)
    facets = facets[first]
    tags   = tags[first]

    return dmesh.meshtags(msh, fdim, facets, tags)


def _cell_centroids(msh, cells: np.ndarray) -> np.ndarray:
    """Compute centroids of given cells by averaging their vertex coordinates."""
    tdim = msh.topology.dim
    # Ensure cell->vertex connectivity exists
    msh.topology.create_connectivity(tdim, 0)
    c2v = msh.topology.connectivity(tdim, 0)
    coords = msh.geometry.x  # (num_geom_nodes, 3)
    # For simplex / hex meshes, cell-to-vertex gives us what we need
    centroids = np.zeros((len(cells), tdim))
    for i, c in enumerate(cells):
        verts = c2v.links(int(c))
        centroids[i] = coords[verts, :tdim].mean(axis=0)
    return centroids


def _semiconductor_tags(regions_cfg: dict) -> set[int]:
    """Set of cell tags whose region role is 'semiconductor'."""
    out: set[int] = set()
    for r in regions_cfg.values():
        if "tag" not in r:
            continue
        role = r.get("role", "semiconductor")
        if role == "semiconductor":
            out.add(int(r["tag"]))
    return out


def is_single_region_semiconductor(cell_tags, regions_cfg: dict) -> bool:
    """
    Detect the single-region-semiconductor fast path.

    Returns True if `cell_tags` is None (no regions declared, treated as
    bare semiconductor) or if every tagged cell has role 'semiconductor'
    and there is at most one distinct tag. This is the condition under
    which the Poisson LHS coefficient can collapse to a scalar Constant
    and the continuity equations assemble directly on the full mesh,
    giving byte-identical behaviour against the M2-M5 1D path.
    """
    if cell_tags is None:
        return True
    tag_role: dict[int, str] = {}
    for r in regions_cfg.values():
        if "tag" not in r:
            continue
        tag_role[int(r["tag"])] = r.get("role", "semiconductor")
    unique_tags = set(int(t) for t in np.unique(cell_tags.values))
    if len(unique_tags) > 1:
        return False
    for t in unique_tags:
        if tag_role.get(t, "semiconductor") != "semiconductor":
            return False
    return True


def build_submesh_by_role(msh, cell_tags, regions_cfg: dict, role: str = "semiconductor"):
    """
    Construct a submesh consisting of all cells whose region tag maps to
    the given `role` in `regions_cfg`.

    Returns
    -------
    (submesh, entity_map, vertex_map, geom_map)
        `entity_map` has length equal to the number of cells in
        `submesh`; `entity_map[i]` is the parent-mesh cell index of
        submesh cell `i`. This is the object needed as
        `entity_maps={msh: inv_entity_map}` when assembling a form on
        the parent mesh that reads Functions defined on the submesh.

    Raises
    ------
    ValueError if `cell_tags` is None or no cells match the role.
    """
    from dolfinx import mesh as dmesh

    if cell_tags is None:
        raise ValueError(
            "build_submesh_by_role requires cell_tags; no regions declared"
        )

    tag_role: dict[int, str] = {}
    for r in regions_cfg.values():
        if "tag" not in r:
            continue
        tag_role[int(r["tag"])] = r.get("role", "semiconductor")

    wanted_tags = sorted(t for t, rr in tag_role.items() if rr == role)
    if not wanted_tags:
        raise ValueError(
            f"No region in regions_cfg has role={role!r}; cannot build submesh."
        )

    tdim = msh.topology.dim
    index_lists = [cell_tags.find(t) for t in wanted_tags]
    cells = np.unique(np.concatenate(index_lists)).astype(np.int32) if index_lists else np.array([], dtype=np.int32)
    if cells.size == 0:
        raise ValueError(
            f"No cells carry a tag matching role={role!r}; check regions_by_box."
        )

    submesh, entity_map, vertex_map, geom_map = dmesh.create_submesh(
        msh, tdim, cells,
    )
    return submesh, entity_map, vertex_map, geom_map


def map_parent_facets_to_submesh(parent_msh, submesh, vertex_map, parent_facets):
    """Translate parent-mesh facet indices to submesh facet indices.

    For each parent facet in `parent_facets`, look up its vertex set via
    the parent mesh's facet-to-vertex connectivity, map those vertices
    into the submesh via `vertex_map` (a `dolfinx.mesh.EntityMap` of
    dimension 0 returned by `create_submesh`), and return the submesh
    facet carrying the same vertex set. Parent facets whose vertices do
    not all lie in the submesh are silently dropped; those facets are
    not boundaries of the submesh (e.g. an oxide-only facet when the
    submesh is silicon-only).

    Used by the multi-region bias-sweep path to translate parent
    `facet_tags` into submesh-facet Dirichlet BCs for phi_n / phi_p on
    the semiconductor while psi continues to be solved on the parent.
    """
    import numpy as np

    tdim = parent_msh.topology.dim
    fdim = tdim - 1
    parent_msh.topology.create_connectivity(fdim, 0)
    submesh.topology.create_connectivity(fdim, 0)
    parent_f2v = parent_msh.topology.connectivity(fdim, 0)
    sub_f2v = submesh.topology.connectivity(fdim, 0)

    parent_facets = np.asarray(list(parent_facets), dtype=np.int32)
    if parent_facets.size == 0:
        return np.array([], dtype=np.int32)

    unique_parent_verts = set()
    for pf in parent_facets:
        for v in parent_f2v.links(int(pf)):
            unique_parent_verts.add(int(v))
    parent_vert_arr = np.asarray(
        sorted(unique_parent_verts), dtype=np.int32,
    )
    sub_vert_arr = np.asarray(
        vertex_map.sub_topology_to_topology(parent_vert_arr, True),
        dtype=np.int64,
    )
    parent_to_sub = {
        int(p): int(s) for p, s in zip(parent_vert_arr, sub_vert_arr)
    }

    sub_facet_index = submesh.topology.index_map(fdim)
    num_sub_facets = sub_facet_index.size_local + sub_facet_index.num_ghosts
    sub_by_vertices: dict[tuple, int] = {}
    for f_sub in range(num_sub_facets):
        verts = tuple(sorted(int(v) for v in sub_f2v.links(f_sub)))
        sub_by_vertices[verts] = f_sub

    out: list[int] = []
    for pf in parent_facets:
        mapped = [parent_to_sub.get(int(v), -1) for v in parent_f2v.links(int(pf))]
        if any(m < 0 for m in mapped):
            continue
        sub_f = sub_by_vertices.get(tuple(sorted(mapped)))
        if sub_f is not None:
            out.append(sub_f)
    return np.asarray(out, dtype=np.int32)


def build_eps_r_function(msh, cell_tags, regions_cfg: dict):
    """
    Build a DG0 cellwise `eps_r(x)` Function on the parent mesh.

    Maps each cell's region tag to the material's relative permittivity.
    Cells whose tag is not represented in `regions_cfg` fall back to
    eps_r = 1.0 (vacuum), which is never exercised in practice because
    `regions_by_box` defines a full partition of the mesh for every
    multi-region benchmark.
    """
    from dolfinx import fem
    from petsc4py import PETSc

    from .materials import get_material

    tag_eps_r: dict[int, float] = {}
    for r in regions_cfg.values():
        if "tag" not in r:
            continue
        mat = get_material(r["material"])
        tag_eps_r[int(r["tag"])] = float(mat.epsilon_r)

    V_DG0 = fem.functionspace(msh, ("DG", 0))
    eps_r_fn = fem.Function(V_DG0, name="eps_r")
    eps_r_fn.x.array[:] = PETSc.ScalarType(1.0)

    if cell_tags is None:
        return eps_r_fn

    for tag, eps in tag_eps_r.items():
        cells = cell_tags.find(tag)
        for c in cells:
            dof = V_DG0.dofmap.cell_dofs(int(c))[0]
            eps_r_fn.x.array[dof] = PETSc.ScalarType(eps)
    eps_r_fn.x.scatter_forward()
    return eps_r_fn
