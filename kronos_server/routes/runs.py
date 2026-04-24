"""GET /runs, GET /runs/{id}, plus field / iv / log / input / manifest downloads."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from fastapi import APIRouter, HTTPException, Request
from fastapi.responses import FileResponse, JSONResponse

from ..jobs import read_status
from ..models import RunListResponse, RunSummary

router = APIRouter()


def _created_at_from_run_id(run_id: str) -> str:
    """Parse the ISO-8601 timestamp prefix from `<ts>_<name>_<sha>` IDs.

    The run_id encodes the timestamp as `YYYY-MM-DDTHH-MM-SSZ` (hyphens in
    the time portion so the string is filename-safe). Convert back to the
    standard ISO-8601 `YYYY-MM-DDTHH:MM:SSZ`.
    """
    stamp = run_id.split("_", 1)[0]
    if "T" not in stamp:
        return run_id
    date_part, time_part = stamp.split("T", 1)
    return f"{date_part}T{time_part.replace('-', ':')}"


def _status_for(run_dir: Path) -> str:
    manifest_path = run_dir / "manifest.json"
    if manifest_path.exists():
        try:
            m = json.loads(manifest_path.read_text())
            status = m.get("status")
            if status in ("completed", "failed", "running"):
                return status
        except Exception:
            pass
    status_rec = read_status(run_dir)
    if status_rec and isinstance(status_rec.get("status"), str):
        return status_rec["status"]
    return "unknown"


def _run_summary_payload(run_id: str, run_dir: Path) -> dict[str, Any]:
    return {
        "run_id": run_id,
        "status": _status_for(run_dir),
        "created_at": _created_at_from_run_id(run_id),
    }


@router.get("/runs", response_model=RunListResponse)
def list_runs(request: Request) -> RunListResponse:
    storage = request.app.state.storage
    entries: list[RunSummary] = []
    for rid in storage.list_run_ids():
        rd = storage.run_dir(rid)
        entries.append(RunSummary(**_run_summary_payload(rid, rd)))
    return RunListResponse(runs=entries)


def _resolve_run_dir(request: Request, run_id: str) -> Path:
    storage = request.app.state.storage
    run_dir = storage.run_dir(run_id)
    if not run_dir.exists() or not run_dir.is_dir():
        raise HTTPException(status_code=404, detail=f"Run {run_id!r} not found.")
    return run_dir


@router.get("/runs/{run_id}")
def get_run(run_id: str, request: Request) -> JSONResponse:
    run_dir = _resolve_run_dir(request, run_id)
    manifest_path = run_dir / "manifest.json"
    status = _status_for(run_dir)
    if manifest_path.exists():
        try:
            manifest = json.loads(manifest_path.read_text())
        except Exception as exc:  # pragma: no cover
            raise HTTPException(
                status_code=500, detail=f"Corrupt manifest: {exc}",
            ) from exc
        manifest["status"] = status
        return JSONResponse(content=manifest)
    payload = _run_summary_payload(run_id, run_dir)
    status_rec = read_status(run_dir) or {}
    for key in ("input_sha256", "wall_time_s", "error"):
        if key in status_rec:
            payload[key] = status_rec[key]
    return JSONResponse(content=payload)


@router.get("/runs/{run_id}/manifest")
def get_manifest(run_id: str, request: Request) -> FileResponse:
    run_dir = _resolve_run_dir(request, run_id)
    path = run_dir / "manifest.json"
    if not path.exists():
        raise HTTPException(status_code=404, detail="Manifest not yet written.")
    return FileResponse(str(path), media_type="application/json")


@router.get("/runs/{run_id}/input")
def get_input(run_id: str, request: Request) -> FileResponse:
    run_dir = _resolve_run_dir(request, run_id)
    path = run_dir / "input.json"
    if not path.exists():
        raise HTTPException(status_code=404, detail="Input JSON not found.")
    return FileResponse(str(path), media_type="application/json")


@router.get("/runs/{run_id}/fields/{name}")
def get_field(run_id: str, name: str, request: Request) -> FileResponse:
    run_dir = _resolve_run_dir(request, run_id)
    fields_dir = run_dir / "fields"
    bp_dir = fields_dir / f"{name}.bp"
    xdmf_path = fields_dir / f"{name}.xdmf"
    if bp_dir.exists() and bp_dir.is_dir():
        # VTX BP directories are multi-file; stream a tar.
        import tarfile
        import tempfile
        fd, tmp = tempfile.mkstemp(prefix=f"{name}_", suffix=".tar")
        import os
        os.close(fd)
        with tarfile.open(tmp, "w") as tar:
            tar.add(str(bp_dir), arcname=f"{name}.bp")
        return FileResponse(
            tmp,
            media_type="application/x-tar",
            filename=f"{name}.bp.tar",
        )
    if xdmf_path.exists():
        return FileResponse(str(xdmf_path), media_type="application/xml")
    raise HTTPException(status_code=404, detail=f"Field {name!r} not found.")


@router.get("/runs/{run_id}/iv/{contact}")
def get_iv(run_id: str, contact: str, request: Request) -> FileResponse:
    run_dir = _resolve_run_dir(request, run_id)
    path = run_dir / "iv" / f"{contact}.csv"
    if not path.exists():
        raise HTTPException(status_code=404, detail=f"IV file for {contact!r} not found.")
    return FileResponse(str(path), media_type="text/csv")


@router.get("/runs/{run_id}/logs")
def get_logs(run_id: str, request: Request) -> FileResponse:
    run_dir = _resolve_run_dir(request, run_id)
    path = run_dir / "logs" / "engine.log"
    if not path.exists():
        raise HTTPException(status_code=404, detail="Engine log not found.")
    return FileResponse(str(path), media_type="text/plain")
