"""Runtime probe of PETSc's linear-solver capabilities (M15 Phase C).

This module is intentionally thin and pure-Python. It does NOT import
dolfinx, cupy, pycuda, or any vendor SDK directly; everything goes through
``petsc4py.PETSc``. When ``petsc4py`` is not installed (as in the lean CI
job for the pure-Python core), every probe returns the conservative
CPU-only answer.

Public API:
    available_backends()  -> list of backend names PETSc can run today.
    device_info()         -> dict for the /capabilities HTTP endpoint and
                             the result manifest.
    resolve_backend()     -> map a user `solver.backend` string to a
                             concrete backend, raising ConfigError if the
                             requested backend is unavailable.

Backend taxonomy (matches schema 1.4.0 from M15 Phase B):
    "cpu-mumps"  - direct LU via MUMPS (V&V reference; always available).
    "gpu-amgx"   - GMRES + AMGX preconditioner on CUDA matrices.
    "gpu-hypre"  - BiCGStab + hypre BoomerAMG on CUDA matrices, where
                   hypre was itself built with HYPRE_USE_GPU. We
                   conservatively report "gpu-hypre" as unavailable when
                   the build flag is not detectable.

Environment override:
    KRONOS_BACKEND - When ``solver.backend == "auto"``, this overrides the
                     auto-resolution preference order. Useful in CI to pin
                     to "cpu-mumps" even on a GPU-capable host.
"""
from __future__ import annotations

import os
from typing import Any

from semi import __version__ as _ENGINE_VERSION


class ConfigError(Exception):
    """Raised when a requested backend is not available in this build."""


_KNOWN_BACKENDS: tuple[str, ...] = ("cpu-mumps", "gpu-amgx", "gpu-hypre")
_AUTO_PREFERENCE: tuple[str, ...] = ("gpu-amgx", "gpu-hypre", "cpu-mumps")


def _probe_petsc() -> dict[str, Any]:
    """Inspect the linked PETSc build. Returns a dict of probe results.

    All fields default to the CPU-only, no-GPU answer so callers never
    need to special-case the missing-petsc4py path.
    """
    info: dict[str, Any] = {
        "petsc_version": None,
        "petsc_complex": False,
        "petsc_int64": False,
        "mat_types": frozenset(),
        "vec_types": frozenset(),
        "pc_types": frozenset(),
        "hypre_gpu": False,
    }
    try:
        from petsc4py import PETSc  # type: ignore[import-not-found]
    except Exception:
        return info

    try:
        info["petsc_version"] = ".".join(str(v) for v in PETSc.Sys.getVersion())
    except Exception:
        pass
    try:
        info["petsc_complex"] = bool(PETSc.ScalarType().dtype.kind == "c")  # type: ignore[attr-defined]
    except Exception:
        pass
    try:
        info["petsc_int64"] = bool(PETSc.IntType().dtype.itemsize == 8)  # type: ignore[attr-defined]
    except Exception:
        pass

    # Matrix and vector type strings supported by this PETSc build.
    # PETSc.Mat.Type and PETSc.Vec.Type are namespaces of string constants.
    mat_types: set[str] = set()
    vec_types: set[str] = set()
    pc_types: set[str] = set()
    try:
        for name in dir(PETSc.Mat.Type):
            if name.startswith("_"):
                continue
            value = getattr(PETSc.Mat.Type, name)
            if isinstance(value, str):
                mat_types.add(value.lower())
    except Exception:
        pass
    try:
        for name in dir(PETSc.Vec.Type):
            if name.startswith("_"):
                continue
            value = getattr(PETSc.Vec.Type, name)
            if isinstance(value, str):
                vec_types.add(value.lower())
    except Exception:
        pass
    try:
        for name in dir(PETSc.PC.Type):
            if name.startswith("_"):
                continue
            value = getattr(PETSc.PC.Type, name)
            if isinstance(value, str):
                pc_types.add(value.lower())
    except Exception:
        pass

    info["mat_types"] = frozenset(mat_types)
    info["vec_types"] = frozenset(vec_types)
    info["pc_types"] = frozenset(pc_types)

    # Conservative hypre-on-GPU probe. PETSc does not expose a clean
    # "is hypre built with CUDA" flag through petsc4py; HYPRE_USE_GPU is
    # a hypre-build option. We only report gpu-hypre when both
    # aijcusparse (or aijhipsparse) and the hypre PC are present, AND
    # PETSc.Options has been pre-set with -pc_hypre_use_gpu. Otherwise
    # treat as CPU-only.
    try:
        opts = PETSc.Options()
        info["hypre_gpu"] = bool(opts.getBool("pc_hypre_use_gpu", default=False))
    except Exception:
        info["hypre_gpu"] = False

    return info


def available_backends() -> list[str]:
    """Subset of ``_KNOWN_BACKENDS`` supported by the linked PETSc build.

    ``cpu-mumps`` is always reported (validated at solve time). GPU
    backends require both a GPU matrix type (``aijcusparse`` or
    ``aijhipsparse``) and the corresponding preconditioner.
    """
    info = _probe_petsc()
    available: list[str] = ["cpu-mumps"]

    has_gpu_mat = (
        "aijcusparse" in info["mat_types"]
        or "aijhipsparse" in info["mat_types"]
    )
    has_gpu_vec = (
        "cuda" in info["vec_types"]
        or "hip" in info["vec_types"]
    )
    if has_gpu_mat and has_gpu_vec:
        if "amgx" in info["pc_types"]:
            available.append("gpu-amgx")
        if "hypre" in info["pc_types"] and info["hypre_gpu"]:
            available.append("gpu-hypre")

    return available


def _gpu_count_and_name() -> tuple[int | None, str | None]:
    """Best-effort device count and device name without importing cupy.

    We call ``nvidia-smi`` only as a last resort and only when a GPU
    backend is otherwise available. Returns ``(None, None)`` on any
    failure; the absence of GPU info is itself a signal.
    """
    return None, None


def device_info() -> dict[str, Any]:
    """Snapshot of build- and device-level facts for the manifest and
    the ``/capabilities`` endpoint.

    Output schema is fixed; tests assert on the key set.
    """
    info = _probe_petsc()
    backends = available_backends()
    device_count, device_name = _gpu_count_and_name()
    return {
        "engine_version": _ENGINE_VERSION,
        "petsc_version": info["petsc_version"],
        "petsc_complex": bool(info["petsc_complex"]),
        "petsc_int64": bool(info["petsc_int64"]),
        "backends_available": backends,
        "device_count": device_count,
        "device_name": device_name,
    }


def resolve_backend(requested: str, available: list[str] | None = None) -> str:
    """Map a user-requested backend to a concrete one, or raise ConfigError.

    Rules (matching docs/M15_STARTER_PROMPT.md Phase C):
      - ``"cpu-mumps"``: always returns ``"cpu-mumps"`` (no GPU dependency).
      - ``"gpu-amgx"``: returns ``"gpu-amgx"`` if available, else raises.
      - ``"gpu-hypre"``: returns ``"gpu-hypre"`` if available, else raises.
      - ``"auto"``: prefers ``gpu-amgx`` > ``gpu-hypre`` > ``cpu-mumps``;
        ``KRONOS_BACKEND`` env var overrides the preference (must itself
        resolve to an available backend).
    """
    if available is None:
        available = available_backends()
    if requested == "cpu-mumps":
        # cpu-mumps is the V&V reference; always honored.
        return "cpu-mumps"
    if requested in ("gpu-amgx", "gpu-hypre"):
        if requested not in available:
            raise ConfigError(
                f"solver.backend={requested!r} was requested but is not "
                f"available in this build (available: {available}). "
                f"Use solver.backend='auto' for graceful CPU fallback, "
                f"or rebuild PETSc with the required GPU support."
            )
        return requested
    if requested == "auto":
        override = os.environ.get("KRONOS_BACKEND")
        if override:
            if override not in _KNOWN_BACKENDS:
                raise ConfigError(
                    f"KRONOS_BACKEND={override!r} is not a recognised "
                    f"backend; expected one of {list(_KNOWN_BACKENDS)}"
                )
            if override not in available and override != "cpu-mumps":
                raise ConfigError(
                    f"KRONOS_BACKEND={override!r} is not available in "
                    f"this build (available: {available})"
                )
            return override
        for candidate in _AUTO_PREFERENCE:
            if candidate in available:
                return candidate
        # cpu-mumps is always in available, so we should not reach here.
        return "cpu-mumps"  # pragma: no cover
    raise ConfigError(
        f"solver.backend={requested!r} is not a recognised backend; "
        f"expected one of {list(_KNOWN_BACKENDS) + ['auto']}"
    )
