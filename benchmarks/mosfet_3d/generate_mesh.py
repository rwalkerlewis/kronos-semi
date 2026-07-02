#!/usr/bin/env python3
"""Regenerate the M19 3D MOSFET tetrahedral meshes from mosfet3d.geo.

Usage (inside the dolfinx/gmsh Docker image):

    python benchmarks/mosfet_3d/generate_mesh.py                # 200k-DOF default
    python benchmarks/mosfet_3d/generate_mesh.py --cl 5e-8 \
        --out benchmarks/mosfet_3d/mosfet3d_coarse.msh          # smoke mesh
    python benchmarks/mosfet_3d/generate_mesh.py --cl 1.4e-8 \
        --out benchmarks/mosfet_3d/mosfet3d_fine.msh            # ~500k GPU mesh

`--cl` is the characteristic mesh length in METERS; it is converted to
the .geo's micrometer units and injected into the `cl` DefineConstant so
the geometry / physical-group definitions stay in one place. Meshes are
written in gmsh msh2 (ASCII 2.2) format, which `dolfinx.io.gmsh` reads
and which keeps the committed files diff-friendly.
"""
from __future__ import annotations

import argparse
from pathlib import Path

HERE = Path(__file__).resolve().parent
GEO = HERE / "mosfet3d.geo"


def generate(cl_m: float, out: Path, binary: bool = False) -> None:
    import gmsh

    cl_um = float(cl_m) * 1.0e6  # .geo works in micrometers

    gmsh.initialize()
    try:
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        if binary:
            gmsh.option.setNumber("Mesh.Binary", 1)
        gmsh.open(str(GEO))
        # Force a uniform target size of `cl`, overriding the .geo's
        # per-point MeshSize default. MeshSizeFromPoints=0 disables the
        # geometry point sizes so MeshSizeMax/Min bind globally.
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeMax", cl_um)
        gmsh.option.setNumber("Mesh.MeshSizeMin", cl_um)
        gmsh.model.mesh.generate(3)
        gmsh.write(str(out))

        # Report node / physical-group counts for the commit log.
        node_tags, _, _ = gmsh.model.mesh.getNodes()
        n_nodes = len(node_tags)
        print(f"[generate_mesh] cl={cl_m:.2e} m -> {out.name}: {n_nodes} nodes")
        for dim, tag in gmsh.model.getPhysicalGroups():
            name = gmsh.model.getPhysicalName(dim, tag)
            ents = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
            print(f"    physical dim={dim} tag={tag:2d} '{name}' "
                  f"({len(ents)} entities)")
    finally:
        gmsh.finalize()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cl", type=float, default=20.0e-9,
                    help="characteristic mesh length in meters (default 20 nm)")
    ap.add_argument("--out", type=Path, default=HERE / "fixtures" / "mosfet3d.msh",
                    help="output .msh path")
    ap.add_argument("--binary", action="store_true",
                    help="write a binary msh2 file (smaller; still read by "
                         "dolfinx.io.gmsh)")
    args = ap.parse_args()
    generate(args.cl, args.out, binary=args.binary)


if __name__ == "__main__":
    main()
