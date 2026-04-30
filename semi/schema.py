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
#   M15 (1.3.0): added top-level `coordinate_system` ("cartesian" |
#                "axisymmetric"); added `oxide_thickness` and
#                `gate_workfunction_eV` to `gate` contacts; added top-level
#                `cv_sweep` block (Vg_min, Vg_max, n_points, frequency_mode).
SCHEMA_SUPPORTED_MINOR = 3


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

    cfg = _fill_defaults(cfg)
    # Cross-field invariants run AFTER defaults so the checks see the
    # final cfg state (e.g. coordinate_system populated to "cartesian").
    _check_cross_field_invariants(cfg)

    return cfg


def _check_cross_field_invariants(cfg: dict[str, Any]) -> None:
    """Validate constraints that the JSON schema cannot express directly.

    These checks live here because Draft-07 has no clean way to express
    "field A is required *only when* field B equals X". They run after
    the structural JSON-schema validation, so by this point we know the
    types are right; we are only enforcing semantic constraints.

    Raises ``SchemaError`` on the first violation found.
    """
    coord = cfg.get("coordinate_system", "cartesian")
    dim = cfg.get("dimension")
    if coord == "axisymmetric":
        if dim != 2:
            raise SchemaError(
                "coordinate_system='axisymmetric' requires dimension == 2 "
                f"(got dimension={dim}); the (r, z) half-plane is two-dimensional."
            )
        mesh = cfg.get("mesh", {})
        if mesh.get("source") == "builtin":
            extents = mesh.get("extents", [])
            # Axis 0 is r; r must start at exactly 0.0 so that the symmetry
            # axis is a mesh facet and the r-weighted weak form behaves
            # correctly (no Dirichlet on r=0; the r factor handles it).
            if extents and extents[0][0] != 0.0:
                raise SchemaError(
                    "coordinate_system='axisymmetric' requires "
                    f"mesh.extents[0][0] == 0.0 (the symmetry axis), "
                    f"got {extents[0][0]!r}."
                )

    # `gate` contacts must declare an oxide thickness; the C-V orchestrator
    # and the verifier need C_ox = eps_ox / t_ox.
    for c in cfg.get("contacts", []):
        if c.get("type") == "gate" and "oxide_thickness" not in c:
            raise SchemaError(
                f"contact {c.get('name', '<unnamed>')!r} has type='gate' but "
                f"is missing the required 'oxide_thickness' field."
            )

    # cv_sweep requires exactly one gate contact.
    if "cv_sweep" in cfg:
        n_gates = sum(1 for c in cfg.get("contacts", []) if c.get("type") == "gate")
        if n_gates != 1:
            raise SchemaError(
                f"cv_sweep requires exactly one contact with type='gate'; "
                f"found {n_gates}."
            )
        sw = cfg["cv_sweep"]
        if "Vg_min" in sw and "Vg_max" in sw and sw["Vg_min"] >= sw["Vg_max"]:
            raise SchemaError(
                f"cv_sweep.Vg_min ({sw['Vg_min']}) must be strictly less than "
                f"cv_sweep.Vg_max ({sw['Vg_max']}); a single-point sweep is "
                f"not meaningful for dQ/dV finite differences."
            )


def _fill_defaults(cfg: dict[str, Any]) -> dict[str, Any]:
    """Fill in sensible defaults for optional sections."""
    cfg.setdefault("coordinate_system", "cartesian")

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

    out = cfg.setdefault("output", {})
    out.setdefault("directory", "./results")
    out.setdefault("fields", ["potential", "n", "p"])
    out.setdefault("write_xdmf", True)
    out.setdefault("write_iv", True)

    if "cv_sweep" in cfg:
        cfg["cv_sweep"].setdefault("frequency_mode", "both")

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
