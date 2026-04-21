"""
JSON schema for semiconductor simulation input.

Design principles:
    - Flat where possible, nested where it aids clarity
    - Units: lengths in meters, potentials in volts, densities in cm^-3
      (device-physics tradition; converted to m^-3 internally)
    - Mesh and external files are referenced by relative path from the JSON
    - Schema validates before any physics runs; missing required fields
      or wrong types fail fast with a clear error

Top-level structure:
    {
      "name": "...",                case name, used in output paths
      "dimension": 1 | 2 | 3,
      "mesh": { ... },              how to get the mesh (builtin or file)
      "regions": { name: {...} },   material per region
      "doping": [ ... ],            list of doping profile specs
      "contacts": [ ... ],          contact BCs (ohmic, gate, insulating)
      "physics": { ... },           model selections
      "solver": { ... },            Newton/continuation parameters
      "sweep": { ... },             optional bias sweep
      "output": { ... }             output controls
    }
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

SCHEMA: dict[str, Any] = {
    "$schema": "http://json-schema.org/draft-07/schema#",
    "type": "object",
    "required": ["name", "dimension", "mesh", "regions", "doping", "contacts"],
    "properties": {
        "name": {"type": "string"},
        "description": {"type": "string"},
        "dimension": {"type": "integer", "enum": [1, 2, 3]},

        "mesh": {
            "type": "object",
            "oneOf": [
                {
                    "description": "Load mesh from file (gmsh .msh, xdmf)",
                    "required": ["source", "path"],
                    "properties": {
                        "source": {"const": "file"},
                        "path": {"type": "string"},
                        "format": {"type": "string", "enum": ["gmsh", "xdmf"]},
                    },
                },
                {
                    "description": "Generate a built-in mesh (interval/rectangle/box)",
                    "required": ["source", "extents", "resolution"],
                    "properties": {
                        "source": {"const": "builtin"},
                        "extents": {
                            "type": "array",
                            "description": "Bounds per dimension: [[xmin,xmax],[ymin,ymax],...]",
                            "items": {
                                "type": "array",
                                "items": {"type": "number"},
                                "minItems": 2, "maxItems": 2,
                            },
                        },
                        "resolution": {
                            "type": "array",
                            "items": {"type": "integer", "minimum": 1},
                        },
                        "regions_by_box": {
                            "type": "array",
                            "items": {
                                "type": "object",
                                "required": ["name", "tag", "bounds"],
                            },
                        },
                        "facets_by_plane": {
                            "type": "array",
                            "items": {
                                "type": "object",
                                "required": ["name", "tag", "axis", "value"],
                            },
                        },
                    },
                },
            ],
        },

        "regions": {
            "type": "object",
            "description": "Region name -> properties. Names match mesh cell tags.",
            "additionalProperties": {
                "type": "object",
                "required": ["material"],
                "properties": {
                    "material": {"type": "string"},
                    "tag": {"type": "integer"},
                    "role": {"type": "string", "enum": ["semiconductor", "insulator"]},
                },
            },
        },

        "doping": {
            "type": "array",
            "items": {
                "type": "object",
                "required": ["region", "profile"],
                "properties": {
                    "region": {"type": "string"},
                    "profile": {
                        "oneOf": [
                            {
                                "type": "object",
                                "required": ["type", "N_D", "N_A"],
                                "properties": {
                                    "type": {"const": "uniform"},
                                    "N_D": {"type": "number"},
                                    "N_A": {"type": "number"},
                                },
                            },
                            {
                                "type": "object",
                                "required": ["type", "axis", "location",
                                             "N_D_left", "N_A_left",
                                             "N_D_right", "N_A_right"],
                                "properties": {
                                    "type": {"const": "step"},
                                    "axis": {"type": "integer", "enum": [0, 1, 2]},
                                    "location": {"type": "number"},
                                    "N_D_left": {"type": "number"},
                                    "N_A_left": {"type": "number"},
                                    "N_D_right": {"type": "number"},
                                    "N_A_right": {"type": "number"},
                                },
                            },
                            {
                                "type": "object",
                                "required": ["type", "center", "sigma", "peak", "dopant"],
                                "properties": {
                                    "type": {"const": "gaussian"},
                                    "center": {"type": "array", "items": {"type": "number"}},
                                    "sigma": {"type": "array", "items": {"type": "number"}},
                                    "peak": {"type": "number"},
                                    "dopant": {"type": "string", "enum": ["donor", "acceptor"]},
                                    "background_N_D": {"type": "number"},
                                    "background_N_A": {"type": "number"},
                                },
                            },
                        ],
                    },
                },
            },
        },

        "contacts": {
            "type": "array",
            "items": {
                "type": "object",
                "required": ["name", "facet", "type"],
                "properties": {
                    "name": {"type": "string"},
                    "facet": {
                        "oneOf": [{"type": "string"}, {"type": "integer"}],
                    },
                    "type": {"type": "string", "enum": ["ohmic", "gate", "insulating"]},
                    "voltage": {"type": "number"},
                    "voltage_sweep": {
                        "type": "object",
                        "description": "Forward bias ramp on this contact.",
                        "required": ["start", "stop", "step"],
                        "properties": {
                            "start": {"type": "number"},
                            "stop": {"type": "number"},
                            "step": {"type": "number", "exclusiveMinimum": 0.0},
                        },
                    },
                    "workfunction": {"type": "number"},
                    "oxide_region": {"type": "string"},
                },
            },
        },

        "physics": {
            "type": "object",
            "properties": {
                "temperature": {"type": "number"},
                "recombination": {
                    "type": "object",
                    "properties": {
                        "srh": {"type": "boolean"},
                        "tau_n": {"type": "number"},
                        "tau_p": {"type": "number"},
                        "E_t": {
                            "type": "number",
                            "description": "Trap level measured from the intrinsic level, eV. Zero for a mid-gap trap.",
                        },
                        "auger": {"type": "boolean"},
                    },
                },
                "mobility": {
                    "type": "object",
                    "properties": {
                        "model": {"type": "string", "enum": ["constant"]},
                        "mu_n": {"type": "number"},
                        "mu_p": {"type": "number"},
                    },
                },
                "statistics": {"type": "string", "enum": ["boltzmann"]},
            },
        },

        "solver": {
            "type": "object",
            "properties": {
                "type": {
                    "type": "string",
                    "enum": [
                        "newton_block",
                        "equilibrium",
                        "drift_diffusion",
                        "bias_sweep",
                        "mos_cv",
                    ],
                },
                "max_iterations": {"type": "integer"},
                "atol": {"type": "number"},
                "rtol": {"type": "number"},
                "damping": {"type": "number"},
                "linear_solver": {"type": "object"},
                "continuation": {
                    "type": "object",
                    "properties": {
                        "steps": {"type": "integer"},
                        "adaptive": {"type": "boolean"},
                        "min_step": {"type": "number", "exclusiveMinimum": 0.0},
                        "max_halvings": {"type": "integer", "minimum": 0},
                        "max_step": {"type": "number", "exclusiveMinimum": 0.0},
                        "easy_iter_threshold": {
                            "type": "integer",
                            "minimum": 1,
                            "description": "SNES converging in strictly fewer than this many iterations counts as easy.",
                        },
                        "grow_factor": {
                            "type": "number",
                            "exclusiveMinimum": 1.0,
                            "description": "Multiplier on the continuation step after easy_iter_threshold consecutive easy solves.",
                        },
                    },
                },
            },
        },

        "sweep": {
            "type": "object",
            "properties": {
                "contact": {"type": "string"},
                "values": {"type": "array", "items": {"type": "number"}},
            },
        },

        "output": {
            "type": "object",
            "properties": {
                "directory": {"type": "string"},
                "fields": {"type": "array", "items": {"type": "string"}},
                "write_xdmf": {"type": "boolean"},
                "write_iv": {"type": "boolean"},
            },
        },
    },
}


class SchemaError(Exception):
    """Raised when a config fails JSON schema validation."""
    pass


def validate(cfg: dict[str, Any]) -> dict[str, Any]:
    """
    Validate a config dict against SCHEMA and fill defaults.

    Uses jsonschema if available; otherwise does minimal required-key checks.
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

    return _fill_defaults(cfg)


def _fill_defaults(cfg: dict[str, Any]) -> dict[str, Any]:
    """Fill in sensible defaults for optional sections."""
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
