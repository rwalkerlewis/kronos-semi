# kronos-semi

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/rwalkerlewis/kronos-semi/blob/main/notebooks/01_pn_junction_1d.ipynb)

A JSON-driven finite-element semiconductor device simulator built on FEniCSx
(dolfinx 0.10). Solves Poisson-coupled drift-diffusion with SRH recombination
on 1D, 2D, and 3D meshes, with an HTTP API for programmatic access.
Originally delivered as an evaluation task; the project is now a working
production engine suitable for web UI integration.

If you're a contributor or a coding agent, start at the [Orientation](#orientation)
section.

## Status (v0.10.0, end of M10)

Milestones M1 through M10 are merged on `main`. The capability matrix is
what the engine actually does today, verified in CI:

| Capability                         | Dimensions    | Status  | Verifier                                                     |
|------------------------------------|---------------|---------|--------------------------------------------------------------|
| Equilibrium Poisson                | 1D / 2D / 3D  | shipped | MMS finest-pair rates L2 = 2.0, H1 = 1.0                     |
| Coupled Slotboom drift-diffusion   | 1D / 2D       | shipped | MMS finest-pair rates L2 >= 1.99 across variants             |
| SRH recombination                  | 1D / 2D       | shipped | verified against SNS analytical at reverse bias              |
| Ohmic contact BCs                  | 1D / 2D / 3D  | shipped | Shockley diode within 10% at forward bias                    |
| Gate contact BCs with phi_ms       | 2D            | shipped | MOS C-V within 10% in depletion window                       |
| Multi-region Poisson (Si/SiO2)     | 2D            | shipped | multi-region MMS L2 = 2.0                                    |
| File-sourced gmsh .msh meshes      | 3D            | shipped | builtin vs gmsh R-match within 1%                            |
| Adaptive bias continuation         | uni + bipolar | shipped | pn junction forward + reverse, 3D resistor                   |
| 3D ohmic V-I linearity             | 3D            | shipped | V-I linearity within 1%                                      |
| Benchmarks                         | 5             | shipped | pn_1d, pn_1d_bias, pn_1d_bias_reverse, mos_2d, resistor_3d   |
| Conservation / mesh convergence    | 1D            | shipped | charge neutrality, Cauchy rates >= 1.8/doubling              |
| V&V                                | 10 studies    | shipped | 62/62 PASS                                                   |
| On-disk result artifact (M9)       | all           | shipped | round-trip tests on all 5 benchmarks, manifest Draft-07 validated |
| `semi-run` CLI (M9)                | all           | shipped | end-to-end artifact writer exercised in CI                   |
| HTTP server (M10)                  | all           | shipped | 15 server tests, end-to-end POST /solve + WebSocket progress |
| Test suite                         | pure + FEM    | shipped | 230 tests, >=95% coverage                                    |
| CI                                 | lint+test+FEM | shipped | green on dev and main                                        |

See [CHANGELOG.md](CHANGELOG.md) for per-milestone deliverables.

## Where this is going (M11+)

M1 through M10 produced a numerically sound engine with an HTTP API. The
remaining work is mostly additional physics and the linear-solver scaling
needed for real 3D devices. The roadmap lives in
[docs/IMPROVEMENT_GUIDE.md](docs/IMPROVEMENT_GUIDE.md).

| Milestone | Summary                                                      | Blocking for                     |
|-----------|--------------------------------------------------------------|----------------------------------|
| M11       | Schema versioning + UI-facing schema companion               | Form-builder-driven UI           |
| M12       | Mesh input beyond axis-aligned boxes                         | Real (non-rectangular) devices   |
| M13       | Transient solver                                             | C-V transients, diode turn-on    |
| M14       | AC small-signal analysis                                     | True C-V, admittance spectroscopy|
| M15       | GPU linear solver path                                       | Anything above ~200k DOFs        |
| M16       | Physics completeness pass (mobility, Auger, FD, Schottky, tunneling) | Real device models       |
| M17       | Heterojunctions / position-dependent band structure          | HEMTs, HBTs                      |
| M18       | UI: first cut (separate repo)                                | Shipping a product               |

Each milestone in the guide has explicit deliverables and acceptance tests.

## Orientation

Four entry points, depending on what you want to do.

### I want to understand the physics without a TCAD background

Read [docs/PHYSICS_INTRO.md](docs/PHYSICS_INTRO.md). Tutorial-style, assumes
programming literacy and calculus but not device physics. About 30 minutes.

### I want to understand the code architecture

Read [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) for the five-layer design
and import rules, then [docs/WALKTHROUGH.md](docs/WALKTHROUGH.md) for an
end-to-end trace of what happens when you call `semi.run.run(cfg)`, with
file and line anchors.

### I want to extend the physics or solver

Read [docs/PHYSICS.md](docs/PHYSICS.md) for the full equation and scaling
reference, then the existing derivation documents under `docs/` (one per
non-trivial device class), then [docs/adr/](docs/adr/) for locked
decisions. Every new physics model needs an MMS verifier — read ADR 0006.

### I want to pick up work from the roadmap

Read [PLAN.md](PLAN.md) for current state and next task, then the specific
milestone section of [docs/IMPROVEMENT_GUIDE.md](docs/IMPROVEMENT_GUIDE.md).
Each milestone has an acceptance test that gates merge — do not negotiate
around it.

## Quick start on Colab (zero setup)

Click the badge above, or use the direct link:

```
https://colab.research.google.com/github/rwalkerlewis/kronos-semi/blob/main/notebooks/01_pn_junction_1d.ipynb
```

The first cell installs FEniCSx on Colab via [FEM on Colab](https://fem-on-colab.github.io/)
(~30 s). Subsequent cells clone this repo, load a JSON benchmark, run the
solve, and plot against analytical curves.

### Notebook catalog

| Notebook | Benchmark | Colab |
|---|---|---|
| 01 pn junction 1D | Equilibrium Poisson, depletion approximation | [![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/rwalkerlewis/kronos-semi/blob/main/notebooks/01_pn_junction_1d.ipynb) |
| 02 pn junction bias | Forward Shockley + reverse SNS sweep | [![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/rwalkerlewis/kronos-semi/blob/main/notebooks/02_pn_junction_bias.ipynb) |
| 03 MOS C-V | 2D MOS capacitor C-V sweep | [![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/rwalkerlewis/kronos-semi/blob/main/notebooks/03_mos_cv.ipynb) |
| 04 resistor 3D | 3D doped bar, builtin and gmsh meshes | [![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/rwalkerlewis/kronos-semi/blob/main/notebooks/04_resistor_3d.ipynb) |

## Local install

dolfinx installs most reliably via conda:

```bash
conda create -n kronos-semi -c conda-forge python=3.12 fenics-dolfinx mpich pyvista
conda activate kronos-semi
git clone https://github.com/rwalkerlewis/kronos-semi.git
cd kronos-semi
pip install -e ".[dev]"
```

For HTTP server support, include the `server` extra:

```bash
pip install -e ".[dev,server]"
```

The pure-Python layer (`schema`, `materials`, `scaling`, `doping`,
`constants`) does not need dolfinx and installs standalone via
`pip install -e .` — useful for schema validation and unit tests in
CI environments without a FEM toolchain. The `server` process also
imports no dolfinx at module scope; FEM work happens in spawned worker
subprocesses.

## Running the HTTP server (M10)

The `kronos-server` console entry point starts a FastAPI server that
exposes the engine over HTTP:

```bash
kronos-server                           # 127.0.0.1:8000 by default
KRONOS_SERVER_PORT=8080 kronos-server   # or via env vars
```

Under Docker Compose, the `server` service runs the same entry point:

```bash
docker compose up -d server
curl http://localhost:8000/health
curl http://localhost:8000/capabilities | python -m json.tool
```

### Example: submit a solve over HTTP

```bash
# Submit a benchmark as a JSON body
RUN_ID=$(curl -sf -X POST http://localhost:8000/solve \
  -H "Content-Type: application/json" \
  -d @benchmarks/pn_1d/pn_junction.json \
  | python -c 'import sys,json; print(json.load(sys.stdin)["run_id"])')

# Poll until complete
while true; do
  STATUS=$(curl -sf http://localhost:8000/runs/$RUN_ID \
    | python -c 'import sys,json; print(json.load(sys.stdin)["status"])')
  [ "$STATUS" = "completed" ] && break
  sleep 2
done

# Fetch the manifest and an IV CSV
curl -sf http://localhost:8000/runs/$RUN_ID/manifest
curl -sf http://localhost:8000/runs/$RUN_ID/iv/cathode
```

OpenAPI spec at `http://localhost:8000/openapi.json`, Swagger UI at
`/docs`. Full endpoint reference in
[CHANGELOG.md](CHANGELOG.md) under the `[0.10.0]` entry.

## Docker

Reproducible dev environment on top of `ghcr.io/fenics/dolfinx/dolfinx:stable`:

```bash
docker compose build
docker compose run --rm test                       # pytest (230 tests)
docker compose run --rm benchmark pn_1d            # run a benchmark
docker compose up -d dev && docker compose exec dev bash
docker compose up jupyter                          # JupyterLab on :8888
docker compose up -d server                        # HTTP server on :8000
```

## Running benchmarks

```bash
docker compose run --rm benchmark pn_1d
docker compose run --rm benchmark pn_1d_bias
docker compose run --rm benchmark pn_1d_bias_reverse
docker compose run --rm benchmark mos_2d
docker compose run --rm benchmark resistor_3d
```

Each benchmark is a JSON file under `benchmarks/<n>/`. The CLI loads it,
runs the solver, checks a registered verifier, and writes plots to
`results/<n>/`.

From Python:

```python
from semi import schema, run

cfg = schema.load("benchmarks/pn_1d/pn_junction.json")
result = run.run(cfg)

# result.psi_phys, result.n_phys, result.p_phys are numpy arrays
# result.x_dof gives the mesh coordinates
# result.iv is a list of {"V": ..., "J": ...} dicts (for bias sweeps)
```

`result` is a `SimulationResult` dataclass; see
[docs/WALKTHROUGH.md §11](docs/WALKTHROUGH.md) for the full schema and
caveats about dolfinx-backed fields not being serializable.

To persist a result to disk for later reading (including from an
environment without dolfinx), use the M9 artifact writer:

```bash
semi-run benchmarks/pn_1d/pn_junction.json --out runs/
```

The resulting `runs/<run_id>/` directory contains a schema-validated
`manifest.json`, the input JSON, mesh and field files, IV CSVs, and
convergence logs.

## JSON input example

Minimal 1D pn junction:

```json
{
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

Doping densities in JSON are in cm⁻³ (device-physics tradition); everything
else is SI. The full schema is defined in `semi/schema.py`.

## Scope vs COMSOL Semiconductor Module

kronos-semi covers the **quasi-static, steady-state** subset.

**In scope (shipped, M1–M10):**

- Poisson with multi-region dielectric (Si/SiO2)
- Drift-diffusion in Slotboom (quasi-Fermi potential) form
- SRH recombination with configurable mid-gap trap energy
- Ohmic contacts and ideal gate contacts with work-function offset
- 1D, 2D, and 3D structured meshes (builtin); 3D unstructured meshes via gmsh
- Bias sweeps with adaptive step-size continuation (unipolar and bipolar)
- Method-of-Manufactured-Solutions and conservation V&V suite
- Machine-readable on-disk result artifacts (schema-validated manifests)
- HTTP API with solve submission, progress streaming, and artifact fetch

**Out of scope today (planned for M13–M17, see [docs/IMPROVEMENT_GUIDE.md](docs/IMPROVEMENT_GUIDE.md)):**

- Caughey-Thomas / Lombardi field-dependent mobility (M16)
- Auger and radiative recombination (M16)
- Fermi-Dirac statistics (M16; Boltzmann throughout today, valid below ~10¹⁹ cm⁻³)
- Impact ionization and avalanche generation (future)
- AC small-signal analysis (M14)
- Transient (time-dependent) solver (M13)
- Band-to-band or trap-assisted tunneling (M16)
- Heterojunctions and position-dependent band structure (M17)
- Schottky contacts and contact resistance models (M16)
- Thermal coupling / self-heating (future)
- Optical generation (future)
- Full MOSFET, FinFET, or 3D transistor geometries (depends on M12 + M16)

## Design notes

### Why Slotboom variables

Galerkin FEM on raw drift-diffusion is unstable when drift dominates
diffusion, which is almost everywhere in a real device. The Slotboom
transform rewrites the continuity equations in quasi-Fermi potentials
Φ_n, Φ_p; currents become pure gradients of Φ, the weak form is
well-posed without upwind stabilization. Implemented in
`semi/physics/drift_diffusion.py`, verified by the MMS-DD suite. See
[ADR 0004](docs/adr/0004-slotboom-variables-for-dd.md) for the
alternatives considered.

### Why nondimensional scaling

A raw 1 µm device at 10¹⁷ cm⁻³ doping has permittivity ~10⁻¹¹, elementary
charge ~10⁻¹⁹, density ~10²³. Newton on that directly diverges — Jacobian
condition number ~10³⁰. All fields are scaled so ψ̂ is O(1) and carrier
ratios n/C₀ are O(1). The scaled Poisson equation has a small parameter
λ² = ε V_t / (q C₀ L₀²) ~ 10⁻⁴, which is the squared Debye-length-to-device
ratio — a singular perturbation that *is* the depletion-region physics.

### Why dolfinx 0.10

The `NonlinearProblem` class in 0.10 wraps PETSc SNES directly (the old
`NewtonSolver` is deprecated), giving us line search, convergence
diagnostics, and block-system support for the coupled (ψ, Φ_n, Φ_p)
solve. See [ADR 0003](docs/adr/0003-dolfinx-0-10-api.md).

### Why no dolfinx in the server process

`kronos_server/` never imports dolfinx, UFL, or PETSc at module scope.
FEM work happens only in worker subprocesses spawned via
`ProcessPoolExecutor(mp_context="spawn")`. A server process can start
and serve `/health`, `/ready`, `/capabilities`, and `/schema` in a
pure-Python environment without the FEM stack, which is useful for
smoke-tested deployments and for UI developers who want a static mock
server.

## Verification

Analytical results are regenerated at import time by
`tests/check_day1_math.py`, which runs without dolfinx. Nine sanity
checks including thermal voltage, Debye length, built-in potential
formulas, mass-action, charge neutrality, and peak |E| for a symmetric
junction.

```bash
python tests/check_day1_math.py    # runs offline, no dolfinx required
pytest tests/                       # 230 tests under Docker
python scripts/run_verification.py  # V&V suite, 62 studies
```

## Planning documents (read order for contributors)

1. [PLAN.md](PLAN.md) — current state, next task, invariants. **Single source of truth** for what to work on.
2. [docs/IMPROVEMENT_GUIDE.md](docs/IMPROVEMENT_GUIDE.md) — M11+ roadmap, acceptance tests, anti-goals.
3. [docs/PHYSICS_INTRO.md](docs/PHYSICS_INTRO.md) — physics tutorial for programmers.
4. [docs/PHYSICS.md](docs/PHYSICS.md) — governing equations, scaling, BCs (reference).
5. [docs/WALKTHROUGH.md](docs/WALKTHROUGH.md) — end-to-end code trace with file:line anchors.
6. [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) — five-layer component design and import rules.
7. [docs/adr/](docs/adr/) — locked decisions. Open a new ADR before changing any invariant.
8. [docs/ROADMAP.md](docs/ROADMAP.md) — per-milestone delivery history (M1–M10).

## License

MIT. See [LICENSE](LICENSE).
