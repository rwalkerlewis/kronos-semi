"""
MMS test fixtures and collection guard.

Tests in this directory require dolfinx and run only inside the
Dockerized FEM environment. The pure-Python CI matrix does not have
dolfinx installed, so we tell pytest to skip the entire subdirectory
when the import fails.
"""
from __future__ import annotations

try:
    import dolfinx  # noqa: F401

    HAS_DOLFINX = True
except ImportError:
    HAS_DOLFINX = False

if not HAS_DOLFINX:
    collect_ignore_glob = ["*"]
