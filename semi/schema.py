"""
JSON schema for semiconductor simulation input.

The schema itself lives in ``schemas/input.v1.json`` (Draft-07); this module
loads it once at import time and exposes it as ``SCHEMA`` for engine and
server code.

Design principles:
    - Flat where possible, nested where it aids clarity
    - Units: lengths in meters, potentials in volts, densities in cm^-3
      (device-physics tradition; converted to m^-3 internally)
    - Mesh and external files are referenced by relative path from the JSON
    - Schema validates before any physics runs; missing required fields
      or wrong types fail fast with a clear error

Top-level structure:
    {
      "schema_version": "1.0.0",     MAJOR.MINOR.PATCH; major gated by engine
      "name": "...",                 case name, used in output paths
      "dimension": 1 | 2 | 3,
      "mesh": { ... },               how to get the mesh (builtin or file)
      "regions": { name: {...} },    material per region
      "doping": [ ... ],             list of doping profile specs
      "contacts": [ ... ],           contact BCs (ohmic, gate, insulating)
      "physics": { ... },            model selections
      "solver": { ... },             Newton/continuation parameters
      "sweep": { ... },              optional bias sweep
      "output": { ... }              output controls
    }
"""
from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any

_SCHEMA_PATH = Path(__file__).parent.parent / "schemas" / "input.v1.json"

# Major version of the input schema this engine build accepts. Inputs whose
# `schema_version` major differs from this number are rejected by `validate`.
# Minor/patch skew is accepted silently.
# TODO(M11+): publish schemas/input.v1.json to a stable URL on every tag
# (e.g. https://schemas.kronos-semi.org/input/v1.0.0.json).
ENGINE_SUPPORTED_SCHEMA_MAJOR = 1

# Minor version tracking. Bumped when a backward-compatible schema addition
# is made.
#   M13 (1.1.0): added the `transient` solver type.
#   M14 (1.2.0): added the `ac_sweep` solver type plus solver.dc_bias and
#                solver.ac sub-objects (frequency sweep specification).
#   M15 (1.3.0): added the top-level `coordinate_system` field with the
#                `axisymmetric` option for cylindrical 2D MOSCAP and similar
#                rotationally-symmetric devices.
#   M15 (1.4.0): added solver.backend (cpu-mumps | gpu-amgx | gpu-hypre |
#                auto) and solver.compute (device, precision, preconditioner,
#                linear_solver). Defaults preserve byte-equivalent CPU-MUMPS
#                behavior; resolution of `auto` happens at solve time.
SCHEMA_SUPPORTED_MINOR = 4


@lru_cache(maxsize=1)
def _load_schema() -> dict[str, Any]:
    with open(_SCHEMA_PATH) as f:
        return json.load(f)


SCHEMA: dict[str, Any] = _load_schema()


class SchemaError(Exception):
    """Raised when a config fails JSON schema validation."""
    pass


def validate(cfg: dict[str, Any]) -> dict[str, Any]:
    """
    Validate a config dict against SCHEMA and fill defaults.

    Uses jsonschema if available; otherwise does minimal required-key checks.
    After successful structural validation the major component of
    ``schema_version`` is checked against ``ENGINE_SUPPORTED_SCHEMA_MAJOR``;
    mismatched majors raise ``SchemaError``.

    Returns the config (possibly modified in-place) with defaults filled.

    Raises SchemaError with all errors concatenated on failure.
    """
    try:
        import jsonschema  # noqa: F401  (presence check for fallback below)
        from jsonschema import Draft7Validator
        validator = Draft7Validator(SCHEMA)
        errors = sorted(validator.iter_errors(cfg), key=lambda e: list(e.path))
        if errors:
            messages = []
            for e in errors:
                path = ".".join(str(p) for p in e.absolute_path) or "<root>"
                messages.append(f"  at {path}: {e.message}")
            raise SchemaError("JSON schema validation failed:\n" + "\n".join(messages))
    except ImportError:
        missing = [k for k in SCHEMA["required"] if k not in cfg]
        if missing:
            raise SchemaError(f"Missing required fields: {missing}") from None

    requested_version = cfg["schema_version"]
    try:
        requested_major = int(str(requested_version).split(".")[0])
    except (ValueError, AttributeError):
        raise SchemaError(
            f"schema_version {requested_version!r} is not a valid semver string"
        ) from None
    if requested_major != ENGINE_SUPPORTED_SCHEMA_MAJOR:
        raise SchemaError(
            f"schema_version {requested_version!r} has major {requested_major}; "
            f"engine supports major {ENGINE_SUPPORTED_SCHEMA_MAJOR}"
        )

    return _fill_defaults(cfg)


def _validate_coordinate_system(cfg: dict[str, Any]) -> None:
    """
    Cross-field validation for the `coordinate_system` setting.

    `axisymmetric` requires:
      - dimension == 2 (the meridian half-plane is a 2D domain),
      - all radial coordinates non-negative (mesh extents r >= 0 if builtin),
      - no Dirichlet (ohmic or gate) contact pinned at r = 0; the symmetry
        axis is a natural no-flux boundary.

    These checks happen after JSON-schema structural validation so the
    error messages can reference resolved values. Failures raise
    SchemaError with a clear message.
    """
    cs = cfg.get("coordinate_system", "cartesian")
    if cs == "cartesian":
        return
    if cs != "axisymmetric":
        raise SchemaError(
            f"coordinate_system={cs!r} is not supported; "
            f"expected 'cartesian' or 'axisymmetric'"
        )

    dim = int(cfg["dimension"])
    if dim != 2:
        raise SchemaError(
            f"coordinate_system='axisymmetric' requires dimension=2 "
            f"(meridian half-plane); got dimension={dim}"
        )

    mesh = cfg.get("mesh", {})
    if mesh.get("source") == "builtin":
        extents = mesh.get("extents", [])
        if extents and float(extents[0][0]) < 0.0:
            raise SchemaError(
                "coordinate_system='axisymmetric': radial extent extents[0] "
                f"must satisfy r_min >= 0; got r_min={extents[0][0]}"
            )

    # Forbid Dirichlet contacts on the symmetry axis. We detect axis facets
    # by looking for facets_by_plane entries with axis=0 and value=0 in
    # builtin meshes; for file meshes the user is responsible for not
    # tagging the axis as a contact (documented in the schema description).
    axis_facet_names: set[str] = set()
    for plane in mesh.get("facets_by_plane", []) if mesh.get("source") == "builtin" else []:
        if int(plane.get("axis", -1)) == 0 and abs(float(plane.get("value", 1.0))) < 1.0e-15:
            axis_facet_names.add(plane["name"])
    for c in cfg.get("contacts", []):
        if c.get("type") in ("ohmic", "gate") and c.get("facet") in axis_facet_names:
            raise SchemaError(
                f"coordinate_system='axisymmetric': contact {c.get('name')!r} "
                f"is pinned at r=0 on facet {c.get('facet')!r}; the symmetry "
                f"axis must remain a natural (Neumann) boundary"
            )


def _validate_compute(cfg: dict[str, Any]) -> None:
    """
    Cross-field validation for `solver.backend` against `solver.compute.device`.

    Mirrors the cross-field pattern used for `coordinate_system`. Rules:
      - `backend == "cpu-mumps"` requires `compute.device in {"auto", "cpu"}`.
      - `backend in {"gpu-amgx", "gpu-hypre"}` requires
        `compute.device in {"auto", "cuda", "hip"}`.
      - `backend == "auto"` is always accepted at validation time; the
        actual device is resolved at solve time by `semi.compute`
        (introduced in Phase C of M15).

    No-op when neither `backend` nor `compute` is present (default-fill
    sets `backend="cpu-mumps"` and leaves `compute` absent, which is the
    byte-equivalent legacy path).
    """
    solver = cfg.get("solver", {})
    backend = solver.get("backend")
    compute = solver.get("compute")
    if backend is None and compute is None:
        return
    backend = backend or "cpu-mumps"
    device = (compute or {}).get("device", "cpu" if backend == "cpu-mumps" else "auto")

    if backend == "cpu-mumps":
        if device not in ("auto", "cpu"):
            raise SchemaError(
                f"solver.backend='cpu-mumps' requires solver.compute.device "
                f"in {{'auto', 'cpu'}}; got {device!r}"
            )
    elif backend in ("gpu-amgx", "gpu-hypre"):
        if device not in ("auto", "cuda", "hip"):
            raise SchemaError(
                f"solver.backend={backend!r} requires solver.compute.device "
                f"in {{'auto', 'cuda', 'hip'}}; got {device!r}"
            )
    elif backend == "auto":
        return
    else:
        raise SchemaError(
            f"solver.backend={backend!r} is not a recognised backend; "
            f"expected one of {{'cpu-mumps', 'gpu-amgx', 'gpu-hypre', 'auto'}}"
        )


def _fill_defaults(cfg: dict[str, Any]) -> dict[str, Any]:
    """Fill in sensible defaults for optional sections."""
    cfg.setdefault("coordinate_system", "cartesian")
    _validate_coordinate_system(cfg)
    phys = cfg.setdefault("physics", {})
    phys.setdefault("temperature", 300.0)
    rec = phys.setdefault("recombination", {})
    rec.setdefault("srh", True)
    rec.setdefault("tau_n", 1.0e-7)
    rec.setdefault("tau_p", 1.0e-7)
    rec.setdefault("E_t", 0.0)
    rec.setdefault("auger", False)
    mob = phys.setdefault("mobility", {})
    mob.setdefault("model", "constant")
    mob.setdefault("mu_n", 1400.0)
    mob.setdefault("mu_p", 450.0)
    phys.setdefault("statistics", "boltzmann")

    solver = cfg.setdefault("solver", {})
    solver.setdefault("type", "equilibrium")
    solver.setdefault("backend", "cpu-mumps")
    solver.setdefault("max_iterations", 50)
    solver.setdefault("atol", 1.0e-10)
    solver.setdefault("rtol", 1.0e-8)
    solver.setdefault("damping", 1.0)
    ls = solver.setdefault("linear_solver", {})
    ls.setdefault("ksp_type", "preonly")
    ls.setdefault("pc_type", "lu")
    ls.setdefault("factor_mat_solver_type", "mumps")
    cont = solver.setdefault("continuation", {})
    cont.setdefault("steps", 10)
    cont.setdefault("adaptive", True)
    cont.setdefault("min_step", 1.0e-4)
    cont.setdefault("max_halvings", 6)
    cont.setdefault("easy_iter_threshold", 4)
    cont.setdefault("grow_factor", 1.5)

    _validate_compute(cfg)

    out = cfg.setdefault("output", {})
    out.setdefault("directory", "./results")
    out.setdefault("fields", ["potential", "n", "p"])
    out.setdefault("write_xdmf", True)
    out.setdefault("write_iv", True)

    return cfg


def load(path: str | Path) -> dict[str, Any]:
    """Load and validate a JSON config file; record its source path for relative lookups."""
    path = Path(path)
    with path.open() as f:
        cfg = json.load(f)
    cfg = validate(cfg)
    cfg["_source_path"] = str(path.resolve())
    cfg["_source_dir"] = str(path.parent.resolve())
    return cfg


def dumps(cfg: dict[str, Any], indent: int = 2) -> str:
    """Serialize a config back to JSON, stripping internal bookkeeping fields."""
    clean = {k: v for k, v in cfg.items() if not k.startswith("_")}
    return json.dumps(clean, indent=indent)
