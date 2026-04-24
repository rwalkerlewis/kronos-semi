"""Pydantic request/response models surfaced in the OpenAPI spec."""
from __future__ import annotations

from typing import Any

from pydantic import BaseModel, Field


class SolveResponse(BaseModel):
    run_id: str = Field(..., description="Identifier for the newly-enqueued run.")
    status: str = Field(..., description="Queued/running/completed/failed.")
    status_url: str = Field(..., description="Relative URL for polling run status.")


class RunSummary(BaseModel):
    run_id: str
    status: str
    created_at: str = Field(..., description="ISO-8601 UTC timestamp parsed from run_id.")


class RunListResponse(BaseModel):
    runs: list[RunSummary]


class HealthResponse(BaseModel):
    status: str


class ReadyResponse(BaseModel):
    ready: bool


class CapabilitiesEngine(BaseModel):
    name: str
    version: str


class CapabilitiesPhysics(BaseModel):
    boltzmann: bool
    fermi_dirac: bool
    srh_recombination: bool
    auger_recombination: bool
    field_dependent_mobility: bool
    transient: bool
    ac_small_signal: bool
    schottky_contacts: bool
    heterojunctions: bool
    tunneling: bool


class CapabilitiesBackends(BaseModel):
    cpu_mumps: bool
    gpu_amgx: bool
    gpu_hypre: bool


class CapabilitiesResponse(BaseModel):
    engine: CapabilitiesEngine
    solver_types: list[str]
    dimensions: list[int]
    physics: CapabilitiesPhysics
    mesh_sources: list[str]
    backends: CapabilitiesBackends


class ErrorResponse(BaseModel):
    detail: Any
