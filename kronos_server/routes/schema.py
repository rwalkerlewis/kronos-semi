"""Static metadata endpoints: /schema, /materials, /capabilities."""
from __future__ import annotations

from typing import Any

from fastapi import APIRouter

from ..models import (
    CapabilitiesBackends,
    CapabilitiesEngine,
    CapabilitiesPhysics,
    CapabilitiesResponse,
)

router = APIRouter()


@router.get("/schema")
def get_schema() -> dict[str, Any]:
    """The input JSON schema plus version metadata. UI form-builders consume this.

    Returns ``{"schema": <Draft-07 schema>, "version": "1.0.0",
    "supported_major": 1}``. The embedded ``schema`` matches
    ``semi.schema.SCHEMA`` verbatim; ``version`` is the advertised schema
    version (taken from the ``schema_version`` property's ``examples``);
    ``supported_major`` is the engine's accepted major version.
    """
    from semi.schema import ENGINE_SUPPORTED_SCHEMA_MAJOR, SCHEMA
    version = (
        SCHEMA.get("properties", {})
        .get("schema_version", {})
        .get("examples", ["1.0.0"])[0]
    )
    return {
        "schema": SCHEMA,
        "version": version,
        "supported_major": ENGINE_SUPPORTED_SCHEMA_MAJOR,
    }


@router.get("/materials")
def get_materials() -> dict[str, dict[str, Any]]:
    """Dump the material database: {name: parameters}."""
    from dataclasses import asdict

    from semi.materials import MATERIALS
    return {name: asdict(mat) for name, mat in MATERIALS.items()}


@router.get("/capabilities", response_model=CapabilitiesResponse)
def get_capabilities() -> CapabilitiesResponse:
    """Feature flags describing what this engine build supports.

    Values are hard-coded from the current state of the codebase as shipped
    at M10. Future milestones flip these flags as they land (see
    docs/IMPROVEMENT_GUIDE.md §4).
    """
    try:
        from semi import __version__
        version = str(__version__)
    except Exception:  # pragma: no cover
        version = "unknown"

    return CapabilitiesResponse(
        engine=CapabilitiesEngine(name="kronos-semi", version=version),
        solver_types=["equilibrium", "drift_diffusion", "bias_sweep", "mos_cv"],
        dimensions=[1, 2, 3],
        physics=CapabilitiesPhysics(
            boltzmann=True,
            fermi_dirac=False,
            srh_recombination=True,
            auger_recombination=False,
            field_dependent_mobility=False,
            transient=False,
            ac_small_signal=False,
            schottky_contacts=False,
            heterojunctions=False,
            tunneling=False,
        ),
        mesh_sources=["builtin", "gmsh"],
        backends=CapabilitiesBackends(
            cpu_mumps=True,
            gpu_amgx=False,
            gpu_hypre=False,
        ),
    )
