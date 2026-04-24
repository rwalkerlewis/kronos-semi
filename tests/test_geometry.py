"""
Tests for `semi.geometry.realize` (M12).

Pure-Python: these tests exercise the gmsh subprocess + content-hashed
cache without importing dolfinx. They do require the `gmsh` binary on
PATH and are skipped when it is absent.
"""
from __future__ import annotations

import shutil
import time
from pathlib import Path

import pytest

from semi import geometry as geom

pytestmark = pytest.mark.skipif(
    shutil.which("gmsh") is None, reason="gmsh not on PATH"
)


_GEO_A = """
Point(1) = {0, 0, 0, 1};
Point(2) = {1, 0, 0, 1};
Point(3) = {1, 1, 0, 1};
Point(4) = {0, 1, 0, 1};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Physical Surface("body", 1) = {1};
Physical Line("bottom", 2) = {1};
Mesh.MshFileVersion = 2.2;
"""

_GEO_B = """
Point(1) = {0, 0, 0, 1};
Point(2) = {2, 0, 0, 1};
Point(3) = {2, 1, 0, 1};
Point(4) = {0, 1, 0, 1};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Physical Surface("body", 1) = {1};
Physical Line("bottom", 2) = {1};
Mesh.MshFileVersion = 2.2;
"""


def _write_geo(path: Path, text: str) -> Path:
    path.write_text(text)
    return path


def test_realize_hits_cache_on_identical_inputs(tmp_path):
    """Second call with the same bytes + clmax + dim reuses the cached .msh."""
    src = tmp_path / "src"
    src.mkdir()
    cache = tmp_path / "cache"

    geo_path = _write_geo(src / "dev.geo", _GEO_A)
    cfg = {"source": "gmsh_geo", "path": "dev.geo", "characteristic_length": 0.25}

    t0 = time.perf_counter()
    msh1 = geom.realize(cfg, source_dir=src, cache_dir=cache, dimension=2)
    cold_dt = time.perf_counter() - t0

    assert msh1.is_file(), f"first call did not produce a .msh: {msh1}"
    mtime_1 = msh1.stat().st_mtime
    size_1 = msh1.stat().st_size

    t1 = time.perf_counter()
    msh2 = geom.realize(cfg, source_dir=src, cache_dir=cache, dimension=2)
    warm_dt = time.perf_counter() - t1

    assert msh2 == msh1
    assert msh2.stat().st_mtime == mtime_1, "cache hit rewrote the file"
    assert msh2.stat().st_size == size_1

    # Cold call runs gmsh (subprocess spawn + mesher; >= ~0.05 s typical).
    # Warm call should be essentially free (< 0.1 s); allow a generous
    # margin for slow test runners.
    assert warm_dt < cold_dt, (
        f"cache warm call ({warm_dt:.3f} s) not faster than cold "
        f"({cold_dt:.3f} s); caching may be broken"
    )
    assert warm_dt < 0.2, (
        f"cache warm call took {warm_dt:.3f} s; expected << 1 s -- the "
        f"subprocess likely ran again despite the cache"
    )


def test_realize_cache_distinguishes_characteristic_length(tmp_path):
    """Changing characteristic_length produces a new cached .msh."""
    src = tmp_path / "src"
    src.mkdir()
    cache = tmp_path / "cache"

    _write_geo(src / "dev.geo", _GEO_A)
    cfg_a = {"source": "gmsh_geo", "path": "dev.geo", "characteristic_length": 0.25}
    cfg_b = {"source": "gmsh_geo", "path": "dev.geo", "characteristic_length": 0.125}

    msh_a = geom.realize(cfg_a, source_dir=src, cache_dir=cache, dimension=2)
    msh_b = geom.realize(cfg_b, source_dir=src, cache_dir=cache, dimension=2)
    assert msh_a != msh_b
    assert msh_a.is_file() and msh_b.is_file()
    assert msh_a.stat().st_size != msh_b.stat().st_size or msh_a.read_bytes() != msh_b.read_bytes()


def test_realize_cache_distinguishes_geo_bytes(tmp_path):
    """Different .geo content produces a different cache entry."""
    src = tmp_path / "src"
    src.mkdir()
    cache = tmp_path / "cache"

    geo_a = _write_geo(src / "a.geo", _GEO_A)
    geo_b = _write_geo(src / "b.geo", _GEO_B)
    cfg_a = {"source": "gmsh_geo", "path": str(geo_a), "characteristic_length": 0.25}
    cfg_b = {"source": "gmsh_geo", "path": str(geo_b), "characteristic_length": 0.25}

    msh_a = geom.realize(cfg_a, source_dir=src, cache_dir=cache, dimension=2)
    msh_b = geom.realize(cfg_b, source_dir=src, cache_dir=cache, dimension=2)
    assert msh_a != msh_b
    assert msh_a.is_file() and msh_b.is_file()


def test_realize_raises_on_missing_geo(tmp_path):
    """A missing .geo path raises GeometryError before any subprocess call."""
    src = tmp_path / "src"
    src.mkdir()
    cfg = {"source": "gmsh_geo", "path": "ghost.geo", "characteristic_length": 0.25}

    with pytest.raises(geom.GeometryError, match="not found"):
        geom.realize(cfg, source_dir=src, cache_dir=tmp_path / "cache", dimension=2)


def test_realize_raises_on_bad_geo(tmp_path):
    """A malformed .geo makes gmsh fail; GeometryError carries stderr."""
    src = tmp_path / "src"
    src.mkdir()
    (src / "bad.geo").write_text("this is not a valid gmsh script\n")
    cfg = {"source": "gmsh_geo", "path": "bad.geo"}

    with pytest.raises(geom.GeometryError, match=r"(gmsh|exit)"):
        geom.realize(cfg, source_dir=src, cache_dir=tmp_path / "cache", dimension=2)


def test_realize_rejects_unknown_source(tmp_path):
    cfg = {"source": "opencascade_step", "path": "x.step"}
    with pytest.raises(geom.GeometryError, match="Unsupported"):
        geom.realize(cfg, source_dir=tmp_path, dimension=2)


def test_realize_rejects_invalid_dimension(tmp_path):
    src = tmp_path / "src"
    src.mkdir()
    _write_geo(src / "dev.geo", _GEO_A)
    cfg = {"source": "gmsh_geo", "path": "dev.geo"}
    with pytest.raises(geom.GeometryError, match="dimension"):
        geom.realize(cfg, source_dir=src, cache_dir=tmp_path / "cache", dimension=5)


def test_realize_rejects_nonpositive_clmax(tmp_path):
    src = tmp_path / "src"
    src.mkdir()
    _write_geo(src / "dev.geo", _GEO_A)
    cfg = {"source": "gmsh_geo", "path": "dev.geo", "characteristic_length": -1.0}
    with pytest.raises(geom.GeometryError, match="positive"):
        geom.realize(cfg, source_dir=src, cache_dir=tmp_path / "cache", dimension=2)


def test_realize_default_cache_uses_env(tmp_path, monkeypatch):
    """KRONOS_MESH_CACHE env var overrides the ~/.cache default."""
    src = tmp_path / "src"
    src.mkdir()
    cache = tmp_path / "envcache"

    _write_geo(src / "dev.geo", _GEO_A)
    monkeypatch.setenv("KRONOS_MESH_CACHE", str(cache))
    cfg = {"source": "gmsh_geo", "path": "dev.geo"}

    msh = geom.realize(cfg, source_dir=src, dimension=2)
    assert str(cache) in str(msh)
    assert msh.is_file()
