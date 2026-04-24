"""
Parametric geometry: `.geo` to `.msh` orchestration (M12).

The M12 `geometry` input block points at a gmsh `.geo` file and asks the
engine to produce the `.msh` at solve time. This module owns the
subprocess call: it hashes the `.geo` contents together with the
characteristic-length option, looks up the result in a filesystem
cache, and only re-runs gmsh on a cache miss.

The contract is "gmsh as a subprocess, never as a Python import" so the
engine process stays free of a gmsh Python binding at module scope and
the base-image gmsh is authoritative about meshing. The cache key is
deliberately content-based (not path-based), so two configs that share
the same bytes for their `.geo` hit the same cached `.msh` even when
they live in different directories.

Cache layout::

    <cache_dir>/<first-2-of-hash>/<full-hash>.msh

Keyed on ``sha256(geo_bytes || b"\\x00" || clmax_bytes || b"\\x00" || dim_bytes)``.
"""
from __future__ import annotations

import hashlib
import os
import subprocess
from pathlib import Path
from typing import Any


class GeometryError(RuntimeError):
    """Raised when gmsh fails to mesh a .geo file."""


def _default_cache_dir() -> Path:
    env = os.environ.get("KRONOS_MESH_CACHE")
    if env:
        return Path(env)
    return Path.home() / ".cache" / "kronos-semi" / "mesh"


def _cache_key(geo_bytes: bytes, clmax: float | None, dimension: int) -> str:
    h = hashlib.sha256()
    h.update(geo_bytes)
    h.update(b"\x00")
    if clmax is not None:
        h.update(f"{clmax!r}".encode())
    h.update(b"\x00")
    h.update(str(int(dimension)).encode())
    return h.hexdigest()


def realize(
    geo_cfg: dict[str, Any],
    source_dir: Path | str,
    cache_dir: Path | str | None = None,
    dimension: int = 2,
) -> Path:
    """
    Realize a gmsh `.geo` into a cached `.msh` file.

    Parameters
    ----------
    geo_cfg
        The ``geometry`` block from a validated input JSON. Must declare
        ``source == "gmsh_geo"`` and carry a ``path`` to the `.geo`
        file; ``characteristic_length`` is optional.
    source_dir
        Directory to resolve a relative ``geo_cfg['path']`` against. The
        schema loader populates ``cfg['_source_dir']`` with the input
        JSON's parent directory; pass that here.
    cache_dir
        Override the cache directory. Defaults to
        ``$KRONOS_MESH_CACHE`` if set, else ``~/.cache/kronos-semi/mesh``.
    dimension
        Mesh dimension. Passed to gmsh as ``-1`` / ``-2`` / ``-3``.

    Returns
    -------
    Path
        The cached `.msh` file on disk. Safe to pass to
        ``dolfinx.io.gmsh.read_from_msh``.
    """
    source = geo_cfg.get("source")
    if source != "gmsh_geo":
        raise GeometryError(
            f"Unsupported geometry.source {source!r}; only 'gmsh_geo' is implemented."
        )

    raw_path = Path(geo_cfg["path"])
    if not raw_path.is_absolute():
        raw_path = Path(source_dir) / raw_path
    geo_path = raw_path.resolve()
    if not geo_path.is_file():
        raise GeometryError(f"Geometry .geo file not found: {geo_path}")

    clmax = geo_cfg.get("characteristic_length")
    clmax_f = float(clmax) if clmax is not None else None
    if clmax_f is not None and not (clmax_f > 0.0):
        raise GeometryError(
            f"characteristic_length must be positive, got {clmax_f!r}"
        )

    if dimension not in (1, 2, 3):
        raise GeometryError(f"dimension must be 1, 2, or 3, got {dimension!r}")

    cache_root = Path(cache_dir) if cache_dir is not None else _default_cache_dir()
    geo_bytes = geo_path.read_bytes()
    key = _cache_key(geo_bytes, clmax_f, dimension)
    msh_path = cache_root / key[:2] / f"{key}.msh"

    if msh_path.is_file():
        return msh_path

    msh_path.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "gmsh",
        f"-{int(dimension)}",
        "-format", "msh22",
        str(geo_path),
        "-o", str(msh_path),
    ]
    if clmax_f is not None:
        cmd += ["-clmax", f"{clmax_f:g}"]

    try:
        proc = subprocess.run(
            cmd,
            check=False,
            capture_output=True,
            text=True,
            timeout=300,
        )
    except FileNotFoundError as exc:
        raise GeometryError(
            "gmsh executable not found on PATH. Install gmsh >= 4.11 or "
            "set the image's gmsh binary before using geometry inputs."
        ) from exc
    except subprocess.TimeoutExpired as exc:
        if msh_path.is_file():
            try:
                msh_path.unlink()
            except OSError:
                pass
        raise GeometryError(
            f"gmsh timed out after 300 s on {geo_path.name}"
        ) from exc

    if proc.returncode != 0 or not msh_path.is_file():
        if msh_path.is_file():
            try:
                msh_path.unlink()
            except OSError:
                pass
        stderr = (proc.stderr or "").strip()
        stdout = (proc.stdout or "").strip()
        tail = stderr if stderr else stdout
        raise GeometryError(
            f"gmsh failed on {geo_path.name} (exit {proc.returncode}):\n{tail}"
        )

    return msh_path
