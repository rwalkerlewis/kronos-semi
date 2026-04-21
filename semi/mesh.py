"""
Mesh construction and tagging.

Supports:
    - Builtin meshes: interval (1D), rectangle (2D), box (3D)
      with region tagging via axis-aligned boxes and facet tagging via planes
    - File-based meshes: gmsh .msh via dolfinx.io.gmshio (deferred; not needed for 1D)

The output of `build_mesh` is a tuple (mesh, cell_tags, facet_tags) where
the tag meshtag objects follow the convention used throughout dolfinx.
Integer tags are assigned from the JSON config's 'tag' fields.

This module requires dolfinx and is imported lazily. The non-FEM code in
this package (schema, materials, doping, scaling) does NOT import this.
"""
from __future__ import annotations

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
    # Lazy imports — these fail loudly if dolfinx isn't installed

    mesh_cfg = cfg["mesh"]
    source = mesh_cfg["source"]

    if source == "builtin":
        msh = _build_builtin(mesh_cfg, cfg["dimension"])
    elif source == "file":
        msh = _build_from_file(mesh_cfg)
    else:
        raise ValueError(f"Unknown mesh source {source!r}")

    cell_tags  = _tag_regions(msh, mesh_cfg.get("regions_by_box", []))
    facet_tags = _tag_facets(msh, mesh_cfg.get("facets_by_plane", []))
    return msh, cell_tags, facet_tags


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


def _build_from_file(mesh_cfg: dict):
    """Read gmsh .msh or XDMF mesh. Returns only the mesh; tags come from file."""
    # Deferred; will implement for 2D MOS and 3D cases
    raise NotImplementedError("File-based mesh loading not yet wired up. Use builtin for now.")


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
