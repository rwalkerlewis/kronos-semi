# JSON input schema reference

The input JSON contract is defined by
[`schemas/input.v1.json`](../../schemas/input.v1.json) (Draft-07).
Every object node carries a UI-facing `description`. The engine
loader is [`semi/schema.py`](../../semi/schema.py).

This page summarises the schema. For the authoritative definition,
the full annotated source of truth is the JSON file.

## Versioning

- `schema_version` is required at the top level. Format `MAJOR.MINOR.PATCH`.
- The engine refuses inputs whose major version differs from
  `ENGINE_SUPPORTED_SCHEMA_MAJOR` (currently `1` in
  [`semi/schema.py`](../../semi/schema.py)).
- Minor/patch skew is accepted silently. Current minor: `3`
  (`SCHEMA_SUPPORTED_MINOR = 3`).
- Current schema version: **1.3.0** (added `coordinate_system`).

History:
- **1.0.0** (M11) — extracted from `semi.schema.SCHEMA`.
- **1.1.0** (M13) — `solver.type = "transient"` and the BDF fields.
- **1.2.0** (M14) — `solver.type = "ac_sweep"` and the `dc_bias` / `ac` blocks.
- **1.3.0** (M14.2) — `coordinate_system` (`"cartesian"` / `"axisymmetric"`).

## Top-level fields

| Field | Required | Type | Notes |
|---|---|---|---|
| `schema_version` | ✓ | string | Semver. Currently `"1.3.0"`. |
| `coordinate_system` | – | enum | `"cartesian"` (default) or `"axisymmetric"`. |
| `name` | ✓ | string | Case name; used in output paths. |
| `description` | – | string | Free-form. |
| `dimension` | ✓ | int | 1, 2, or 3. Must equal 2 when `coordinate_system="axisymmetric"`. |
| `mesh` | ✓ | object | Builtin box or external gmsh `.msh`. |
| `regions` | ✓ | object | name → material, role. |
| `doping` | ✓ | array | Per-region doping profiles. |
| `contacts` | ✓ | array | Boundary conditions. |
| `physics` | – | object | Recombination model parameters. |
| `solver` | – | object | Solver type and tolerances. |

## `mesh` block

Two variants chosen by `source`:

- `source: "builtin"` — structured box mesh.
  - `extents`: `[[xmin, xmax], ...]` per dimension, in meters.
  - `resolution`: cells per dimension.
  - `regions_by_box`: optional list of axis-aligned sub-boxes that
    tag cells by region (required for multi-region devices).
  - `facets_by_plane`: list of axis-aligned planes that tag boundary
    facets, referenced by name in `contacts[*].facet`.

- `source: "file"` — external gmsh `.msh`.
  - `path`: relative to the JSON file's directory.
  - `format`: `"gmsh"` (XDMF reserved, not implemented).
  - Physical groups in the `.msh` supply region and facet tags.
  - For axisymmetric benchmarks, the meridian mesh is built from a
    `.geo` template (e.g.
    [`benchmarks/moscap_axisym_2d/moscap_axisym.geo`](../../benchmarks/moscap_axisym_2d/moscap_axisym.geo))
    and a generated `.msh` is committed alongside it.

## `coordinate_system` and axisymmetric constraints

When `coordinate_system == "axisymmetric"`, cross-field validation
in [`semi/schema.py`](../../semi/schema.py) enforces:

1. `dimension == 2`.
2. The radial extent is non-negative; the first mesh coordinate is `r ≥ 0`.
3. No Dirichlet contact lives on the symmetry axis `r = 0`. The axis is a natural (no-flux) boundary.
4. The outer radial wall `r = R` is treated as homogeneous Neumann; choose `R` large enough that the solution is insensitive to the cutoff (`R ≳ 5 W_dmax` for MOSCAP).

See [`docs/theory/axisymmetric.md`](../theory/axisymmetric.md) for the derivation.

## `doping` profiles

Each entry targets one region with one profile. Densities are in
cm⁻³ (device-physics convention) and are converted to m⁻³ internally.

- `uniform`: `N_D`, `N_A`.
- `step`: `axis`, `location`, `N_{D,A}_{left,right}`.
- `gaussian`: `axis`, `center`, `sigma`, `peak`, `donor` (bool).

## `contacts`

| Field | Notes |
|---|---|
| `name` | Symbolic contact name, used in output IV files. |
| `facet` | Either an int tag or a string referencing `mesh.facets_by_plane[*].name`. |
| `type` | `"ohmic"`, `"gate"`, `"insulating"`. |
| `voltage` | DC value in V (default 0). |
| `voltage_sweep` | `{start, stop, step}` for bias-sweep runners. |
| `workfunction` | Required for `type: "gate"` (V). |
| `oxide_region` | Required for gate; the insulator region key. |

## `solver`

| `solver.type` | Runner | Notes |
|---|---|---|
| `"equilibrium"` | `run_equilibrium` | Equilibrium Poisson only. |
| `"drift_diffusion"` | `run` | Coupled Slotboom DD at fixed bias. |
| `"bias_sweep"` | `run_bias_sweep` | Adaptive continuation along a contact's `voltage_sweep`. |
| `"mos_cv"` | `run_mos_cv` | 2D MOS C–V via `numpy.gradient(Q, V)`. |
| `"mos_cap_ac"` | `run_mos_cap_ac` | C–V via analytic dQ/dV (M14.1). |
| `"transient"` | `run_transient` | BDF1/BDF2 (M13). Fields: `t_end`, `dt`, `order`, `max_steps`, `output_every`. |
| `"ac_sweep"` | `run_ac_sweep` | Small-signal AC at a DC operating point (M14). |

Common `solver` sub-keys: `tolerances`, `continuation`, `jacobian_shift`.

## Minimal 1D pn-junction example

```json
{
  "schema_version": "1.0.0",
  "name": "pn_junction_1d",
  "dimension": 1,
  "mesh": {
    "source": "builtin",
    "extents": [[0.0, 2.0e-6]],
    "resolution": [400],
    "facets_by_plane": [
      {"name": "anode",   "tag": 1, "axis": 0, "value": 0.0},
      {"name": "cathode", "tag": 2, "axis": 0, "value": 2.0e-6}
    ]
  },
  "regions": {
    "silicon": {"material": "Si", "tag": 1, "role": "semiconductor"}
  },
  "doping": [
    {
      "region": "silicon",
      "profile": {
        "type": "step", "axis": 0, "location": 1.0e-6,
        "N_D_left": 0.0,  "N_A_left": 1.0e17,
        "N_D_right": 1.0e17, "N_A_right": 0.0
      }
    }
  ],
  "contacts": [
    {"name": "anode",   "facet": "anode",   "type": "ohmic", "voltage": 0.0},
    {"name": "cathode", "facet": "cathode", "type": "ohmic", "voltage": 0.0}
  ]
}
```

## Axisymmetric example

[`benchmarks/moscap_axisym_2d/moscap_axisym.json`](../../benchmarks/moscap_axisym_2d/moscap_axisym.json)
is the live reference for the axisymmetric path.

## Known schema caveats (carried forward from M11)

- **Permissive about unknown keys.** Top-level and most nested
  objects lack `additionalProperties: false`, so a UI typo in a
  field name validates and is silently ignored. The manifest
  schema is strict; flipping the input schema is a v2.0.0 change.
- **No per-block `schema_version`.** A bump in the `physics` block
  forces every caller to upgrade atomically.

## See also

- [`schemas/manifest.v1.json`](../../schemas/manifest.v1.json) — output manifest schema (M9).
- [`semi/schema.py`](../../semi/schema.py) — loader and version gate.
- [`semi/io/`](../../semi/io/) — manifest writer and reader.
