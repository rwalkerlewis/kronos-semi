# kronos-semi

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/rwalkerlewis/kronos-semi/blob/main/notebooks/01_pn_junction_1d.ipynb)

A JSON-driven finite-element semiconductor device simulator built on FEniCSx
(dolfinx 0.10). Solves Poisson-coupled drift-diffusion with SRH recombination
on 1D, 2D, and 3D meshes, with an HTTP API for programmatic access.
Originally delivered as an evaluation task; the project is now a working
production engine suitable for web UI integration.

If you're a contributor or a coding agent, start at the [Orientation](#orientation)
section.

## Status (v0.14.0, end of M14.1)

Milestones M1 through M14.1 are merged on `main`. M13.1 is a tracked
follow-up (Scharfetter-Gummel rewrite of the transient continuity
discretisation, see [issue #34](https://github.com/rwalkerlewis/kronos-semi/issues/34));
the transient steady-state-limit test is currently `xfail` pending that
work.

What the engine does today, in plain terms:

- **Equilibrium Poisson** in 1D, 2D, 3D, multi-region (Si/SiO2 with
  cellwise eps_r and the natural flux-continuity interface condition).
- **Coupled Slotboom drift-diffusion bias sweeps** (M2/M3) with adaptive
  step-size continuation in `AdaptiveStepController`; SRH recombination
  with configurable trap energy; supports unipolar and zero-spanning
  bipolar sweeps.
- **2D MOSFET** with Gaussian n+ source/drain implants (M12), built on
  the multi-region Poisson + Slotboom DD stack.
- **Transient solver** with BDF1 (backward Euler) and BDF2 time
  integration (M13), in (psi, n, p) primary-density form; ADRs 0009 /
  0010. Note: the (n,p) Galerkin discretisation diverges at tight SNES
  atol; M13.1 (Scharfetter-Gummel edge-flux assembly, issue #34) is the
  fix, currently in flight.
- **MOS capacitor C-V via analytic AC** (M14.1, the just-landed cycle):
  the `mos_cap_ac` runner solves the linearised Poisson sensitivity at
  each gate bias to return dQ/dV directly in F/m^2, replacing the
  noisier `numpy.gradient(Q, V)` of `mos_cv`. Worst error 6.79% vs the
  depletion-approximation reference in the verifier window.
- **Small-signal AC sweep** for two-terminal pn diodes (M14):
  `(J + jωM) δu = -dF/dV δV` solved via real 2x2 block reformulation;
  the `rc_ac_sweep` benchmark matches analytical depletion C within
  0.4% over [1 Hz, 1 MHz]; ADR 0011.
- **3D doped resistor** with V-I linearity within 1%, builtin or
  gmsh-sourced unstructured meshes (M7).
- **On-disk result artifacts** (M9): manifest.json + fields + IV CSVs +
  convergence logs, schema-validated.
- **HTTP API** (M10): `POST /solve`, progress over WebSocket, `GET
  /runs/{id}/manifest`, etc. Server process imports no dolfinx at
  module scope.
- **Schema versioning + UI form-builder annotations** (M11):
  `schemas/input.v1.json` Draft-07, every object node carries a
  description.
- **V&V**: MMS for Poisson (1D linear/nonlinear, 2D triangles, 2D
  multi-region), mesh convergence, charge / current conservation, MMS
  for coupled DD across three variants and three grids; finest-pair
  rates at theoretical L2 = 2.000.
- **CI**: lint + pure-Python on Python 3.10/3.11/3.12 + dolfinx
  benchmarks + V&V on every push, coverage gate at 95%.

The full capability matrix is captured in [docs/ROADMAP.md](docs/ROADMAP.md);
the per-milestone deliverables are in [CHANGELOG.md](CHANGELOG.md).

## Where this is going (post-M14.1)

M1 through M14.1 produced a numerically sound engine with a versioned,
UI-facing JSON input contract, an HTTP API, transient and AC solvers,
and a 2D MOSFET benchmark. The remaining work is one transient
discretisation rewrite and the usual physics-completeness and
linear-solver scaling axes. Detailed deliverables and acceptance tests
live in [docs/IMPROVEMENT_GUIDE.md](docs/IMPROVEMENT_GUIDE.md).

| Milestone | Summary                                                              | Status              |
|-----------|----------------------------------------------------------------------|---------------------|
| M13.1     | Scharfetter-Gummel edge-flux discretisation for transient continuity | In flight (issue #34) |
| M15       | 2D-axisymmetric MOSCAP C-V (LF + HF, Hu Fig 5-18)                    | In flight           |
| M15       | GPU linear solver path                                               | Planned             |
| M16       | Physics completeness (Caughey-Thomas, Auger, FD, Schottky, tunneling) | Planned             |
| M17       | Heterojunctions / position-dependent band structure                  | Planned             |
| M18       | UI: first cut (separate repo)                                        | Out of scope here   |

The original "Day 5 -- 2D MOS capacitor" line item is **superseded by the
2D-axisymmetric MOSCAP** above. Axisymmetric exploits the body-of-revolution
symmetry of a circular gate: instead of meshing a full 3D disc, we mesh
the (r, z) half-plane and integrate against `2 pi r dr dz`. For the same
fringing-field fidelity this is `O(N^2)` DOFs instead of `O(N^3)`, and the
geometry is naturally a disc-shaped device rather than the rectangular
infinite-stripe of a 2D Cartesian model.

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
docker compose run --rm test                       # pytest (256 tests + 1 xfail)
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
else is SI. The full schema lives in `schemas/input.v1.json` (JSON
Schema Draft-07, every object node annotated with a UI-facing
description); `semi/schema.py` loads it and exposes it as
`semi.schema.SCHEMA`. Every input JSON must declare a top-level
`"schema_version": "1.0.0"`; the engine refuses inputs whose major
version does not match `ENGINE_SUPPORTED_SCHEMA_MAJOR` in
`semi/schema.py` (currently `1`).

## Scope vs COMSOL Semiconductor Module

kronos-semi covers the **quasi-static, steady-state** subset.

**In scope (shipped, M1 through M14.1):**

- Poisson with multi-region dielectric (Si/SiO2)
- Drift-diffusion in Slotboom (quasi-Fermi potential) form
- SRH recombination with configurable mid-gap trap energy
- Ohmic contacts and ideal gate contacts with work-function offset
- 1D, 2D, and 3D structured meshes (builtin); 3D unstructured meshes via gmsh
- Bias sweeps with adaptive step-size continuation (unipolar and bipolar)
- Method-of-Manufactured-Solutions and conservation V&V suite
- Machine-readable on-disk result artifacts (schema-validated manifests)
- HTTP API with solve submission, progress streaming, and artifact fetch
- 2D MOSFET (M12) with Gaussian n+ source/drain implants
- Transient BDF1/BDF2 (M13); the (n,p) Galerkin discretisation has a
  known atol-driven divergence (issue #34), tracked by M13.1
- Small-signal AC sweep (M14) on two-terminal pn devices, plus the
  M14.1 `mos_cap_ac` runner for analytic differential C(V_gate)

**Out of scope today (planned for M15 through M17, see [docs/IMPROVEMENT_GUIDE.md](docs/IMPROVEMENT_GUIDE.md)):**

- GPU linear solver path; CPU-LU only above ~200k DOFs is a hard wall (M15)
- Caughey-Thomas / Lombardi field-dependent mobility (M16)
- Auger and radiative recombination (M16)
- Fermi-Dirac statistics (M16; Boltzmann throughout today, valid below ~10¹⁹ cm⁻³)
- Impact ionization and avalanche generation (future)
- Band-to-band or trap-assisted tunneling (M16)
- Heterojunctions and position-dependent band structure (M17)
- Schottky contacts and contact resistance models (M16)
- Thermal coupling / self-heating (future)
- Optical generation (future)
- Full FinFET or 3D transistor geometries (depends on M16 + meshing work)

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
pytest tests/                       # 256 tests + 1 xfail under Docker
python scripts/run_verification.py  # V&V suite
```

### Axisymmetric MOSCAP (M15, in flight)

`tests/check_moscap_axi_math.py` runs the closed-form Hu-Ch.5 reference
formulas for the n+poly / SiO2 / p-Si stack at the benchmark defaults
(`t_ox = 5 nm, N_A = 1e17 cm^-3, T = 300 K`):

- `phi_F = V_t ln(N_A / n_i)` and `2 phi_F` (~ 0.83 V at N_A = 1e17)
- `V_FB = phi_m - (chi_Si + Eg/2 + phi_F)` for n+ poly (~ -0.95 V)
- `W_dep,max = sqrt(2 eps_Si 2 phi_F / (q N_A))` (~ 104 nm)
- `C_ox = eps_ox / t_ox` and `C_dep,min = eps_Si / W_dep,max`
- `C_min/C_ox = 1 / (1 + eps_ox W_dep,max / (eps_Si t_ox))`
- `V_T = V_FB + 2 phi_F + sqrt(2 eps_Si q N_A 2 phi_F) / C_ox`

Wired into `pytest` at `tests/test_moscap_axi_math.py` and exercised by
the CV sweep in the notebook (Slice 3+). FEM acceptance criteria from
the problem statement: LF and HF curves coincide for V_g < V_T to within
1%, the HF curve floor matches `C_min/C_ox` to within 5%, and the LF
curve climbs back to >= 0.95 C_ox at `V_g = V_T + 1 V`.

```bash
python tests/check_moscap_axi_math.py            # standalone sanity
pytest tests/test_moscap_axi_math.py -v          # 13 closed-form assertions
```

## Planning documents (read order for contributors)

1. [PLAN.md](PLAN.md) — current state, next task, invariants. **Single source of truth** for what to work on.
2. [docs/IMPROVEMENT_GUIDE.md](docs/IMPROVEMENT_GUIDE.md) — M15+ roadmap, acceptance tests, anti-goals.
3. [docs/PHYSICS_INTRO.md](docs/PHYSICS_INTRO.md) — physics tutorial for programmers.
4. [docs/PHYSICS.md](docs/PHYSICS.md) — governing equations, scaling, BCs (reference).
5. [docs/WALKTHROUGH.md](docs/WALKTHROUGH.md) — end-to-end code trace with file:line anchors.
6. [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) — five-layer component design and import rules.
7. [docs/adr/](docs/adr/) — locked decisions. Open a new ADR before changing any invariant.
8. [docs/ROADMAP.md](docs/ROADMAP.md) — per-milestone delivery history (M1 through M14.1).

## License

MIT. See [LICENSE](LICENSE).
