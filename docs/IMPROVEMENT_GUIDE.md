# kronos-semi: Improvement Guide

**Purpose.** This document is the ground-truth roadmap for transforming kronos-semi
from a KronosAI submission artifact into a production engine that sits behind a
COMSOL-Semiconductor-style web GUI. It is written to be consumed by both human
contributors and coding agents (Claude, Codex, Cursor, etc.); read it alongside
`PLAN.md` and `docs/ROADMAP.md`.

Nothing in here is aspirational hand-waving. Every milestone below states its
input, its acceptance test, and its dependency on prior milestones.

---

## 1. Honest current state (as of v0.14.1)

What exists and works:

- `semi/` package with a five-layer architecture; Layer 3 (pure-Python
  core) has no dolfinx dependency and is independently testable.
- JSON schema 1.3.0 (Draft-07, enforced via `jsonschema`) covering
  mesh, regions, doping, contacts, physics, solver, output, plus the
  top-level `coordinate_system` field added in M14.2.
- Builtin meshes in 1D/2D/3D plus gmsh `.msh` loader with physical
  groups.
- Equilibrium Poisson, coupled drift-diffusion in Slotboom form, SRH
  recombination with configurable trap energy.
- Ohmic and ideal-gate contacts; multi-region Poisson with Si/SiO2
  via `create_submesh`.
- Adaptive bias-ramp continuation with halving and growth, bipolar
  sweep support.
- Transient solver in Slotboom primary unknowns: backward-Euler and
  BDF2, with the deep-steady-state runner matching `bias_sweep`
  within 1e-4 relative error and the `pn_1d_turnon` benchmark within
  5% (M13, M13.1; ADRs 0010, 0014).
- Small-signal AC sweep for two-terminal devices via the linearised
  `(J + jωM) δu = -dF/dV δV` system in real 2x2 block form; the
  `rc_ac_sweep` benchmark matches analytical depletion C within
  0.4% over [1 Hz, 1 MHz] (M14; ADR 0011).
- MOSCAP differential capacitance via AC admittance: the
  `mos_cap_ac` runner returns dQ/dV directly from `Im(Y) / (2πf)`,
  audit case 03 confirms byte-identity with `mos_cv` Q_gate (M14.1).
- Axisymmetric (cylindrical) 2D path: r-weighted Poisson and Slotboom
  drift-diffusion on the meridian half-plane; `benchmarks/moscap_axisym_2d/`
  reproduces Hu Fig. 5-18; runner dispatch in `mos_cap_ac.py` r-weights
  charge and sensitivity forms (M14.2; PRs #64, #65).
- Eight shipped benchmarks (`pn_1d`, `pn_1d_bias`, `pn_1d_bias_reverse`,
  `pn_1d_turnon`, `mos_2d`, `resistor_3d`, `rc_ac_sweep`,
  `moscap_axisym_2d`), each with a verifier asserting closed-form
  physical correctness.
- MMS verification for Poisson and DD across multiple variants and
  grids, discrete-conservation checks, mesh-convergence ladder, plus
  AC consistency and audit suites; 237+ pure-Python tests passing.
- Five Colab notebooks that install FEniCSx via fem-on-colab and run
  each benchmark end-to-end.
- On-disk result artifact contract: `manifest.json`, fields, IV CSVs,
  and convergence logs under `runs/<run_id>/`, schema-validated; a
  pure-Python reader round-trips without dolfinx (M9).
- HTTP server surface: `kronos_server/` with `POST /solve`,
  `GET /runs/{id}`, WebSocket `/runs/{id}/stream`, `GET /materials`,
  `GET /schema`; the server process imports no dolfinx at module
  scope (M10).
- Schema versioning: `schemas/input.v1.json` Draft-07 with full
  UI-facing annotations, every benchmark validates, the major-version
  gate refuses mismatched majors (M11).
- 2D MOSFET with Gaussian n+ source/drain implants on the multi-region
  Poisson + Slotboom DD stack (M12).

What does *not* exist:

- **No GPU linear solve.** Every solve is MUMPS LU on CPU via PETSc.
  Above ~200k DOFs this is the wall. M15 closes the gap.
- **No field-dependent mobility, Auger, Fermi-Dirac, Schottky
  contacts, or tunneling.** COMSOL Semiconductor has all of these.
  M16, one PR per model.
- **No heterojunctions.** Position-dependent χ and Eg are not yet
  supported. M17 (depends on M16.4 Fermi-Dirac).
- **No Cartesian-2D MOSCAP variant** and no rigorous gate-driven HF
  C-V method. Tracked as M14.2.x in
  [`PLAN.md`](../PLAN.md) §Backlog and `docs/ROADMAP.md`.
- **No `GET /capabilities` endpoint.** UI integration item; M15 will
  close it because GPU availability is itself a capability.

Bottom line: the numerics are sound, the engine is addressable from
outside Python, results are consumable without a dolfinx install, and
the schema is versioned and richly annotated. The remaining work is
scaling the linear solver (M15) and adding the physics models
engineers actually use daily (M16+).

---

## 2. Design principles (non-negotiable)

These are the constraints every subsequent milestone must respect. They are
the invariants that let a GUI be built without constantly re-ploughing the
engine.

**P1. The JSON is the contract.** The UI never calls Python functions
directly. The UI emits JSON, the engine consumes JSON, the engine emits a
result artifact. No shared memory, no pickled objects, no language coupling.

**P2. Results are plain files.** A completed simulation is a directory on
disk (or object-store prefix) containing a `manifest.json`, a `mesh.vtu` (or
`.bp` / `.xdmf+h5`), one `field_<name>.vtu` per scalar/vector field, and one
`iv.csv` per swept quantity. A UI reads these with stock libraries (vtk.js,
numpy, pandas); no dolfinx import ever happens on the UI side.

**P3. The engine is stateless between jobs.** Each `POST /solve` gets a fresh
MPI process (or worker pool slot). State that must persist across jobs lives
in the object store, never in RAM.

**P4. The schema versions.** Every JSON input carries
`"schema_version": "1.2.0"`. Engine refuses to run mismatched major versions
and warns on minor skew.

**P5. GPU is opt-in, not load-bearing.** The CPU path stays as the reference
and V&V truth. GPU accelerates the linear solver and the assembly of large
forms; it never changes what the result *is*, only how fast.

**P6. Nothing in the engine API leaks PETSc or UFL types.** The boundary types
are: JSON in, numpy arrays + mesh files + manifest dict out. Anything else is
an implementation detail.

---

## 3. Result artifact contract (ship this first, everything else depends on it)

Every successful solve writes a directory with this layout:

```
<run_id>/
├── manifest.json           # schema-versioned metadata
├── input.json              # exact JSON that was solved (copy)
├── mesh/
│   ├── mesh.xdmf           # topology + geometry
│   └── mesh.h5
├── fields/
│   ├── psi.bp              # potential (ADIOS2 BP5 or vtu)
│   ├── n.bp                # electron density
│   ├── p.bp                # hole density
│   ├── N_net.bp            # net doping
│   ├── E.bp                # electric field (vector)
│   └── J.bp                # current density (vector, per bias step)
├── iv/
│   └── <contact>.csv       # V, J_n, J_p, J_tot per sweep step
├── convergence/
│   └── snes.csv            # per-step Newton iter count, residual norms
└── logs/
    └── engine.log
```

### `manifest.json` schema (v1)

```json
{
  "schema_version": "1.0.0",
  "engine": {"name": "kronos-semi", "version": "0.9.0", "commit": "..."},
  "run_id": "2026-04-23T16-32-11Z_pn_bias_a1b2c3",
  "status": "completed",
  "wall_time_s": 12.4,
  "input_sha256": "...",
  "solver": {
    "type": "bias_sweep",
    "backend": "petsc-cpu-mumps",
    "n_dofs": 1200,
    "n_steps": 21,
    "converged": true
  },
  "fields": [
    {"name": "psi",   "units": "V",    "path": "fields/psi.bp",   "rank": 0},
    {"name": "n",     "units": "m^-3", "path": "fields/n.bp",     "rank": 0},
    {"name": "p",     "units": "m^-3", "path": "fields/p.bp",     "rank": 0},
    {"name": "E",     "units": "V/m",  "path": "fields/E.bp",     "rank": 1},
    {"name": "J",     "units": "A/m^2","path": "fields/J.bp",     "rank": 1}
  ],
  "sweeps": [
    {"kind": "voltage", "contact": "cathode",
     "path": "iv/cathode.csv", "n_steps": 21, "bipolar": false}
  ],
  "mesh": {
    "path": "mesh/mesh.xdmf", "dimension": 1,
    "n_cells": 400, "n_vertices": 401,
    "regions": [{"name": "silicon", "tag": 1, "material": "Si"}]
  },
  "warnings": []
}
```

The UI reads the manifest first; everything else is lazy-loaded.

---

## 4. Milestones (M9 onward; M1–M8 already shipped)

Each milestone is self-contained, has a single deliverable, and has an
acceptance test that can be demanded of a coding agent. Do them roughly in
order, but M9/M10 are load-bearing — everything after assumes them.

### M9 — Result artifact writer + manifest

**Status: Done (2026-04-23).** `semi/io/artifact.py`, `semi-run` CLI,
`schemas/manifest.v1.json`, pure-Python reader. Full deliverable,
acceptance tests, and scope preserved in §10.

---

### M10 — HTTP server: `POST /solve`, `GET /runs/{id}`

**Status: Done (M10).** `kronos_server/` package, FastAPI surface,
WebSocket progress channel, `ProcessPoolExecutor`-backed workers.
Full deliverable, acceptance tests, and scope preserved in §10.

---

### M11 — Schema versioning + UI-facing schema companion

**Status: Done (M11).** `schemas/input.v1.json` Draft-07 with full
UI annotations; `schema_version` enforced via the major-version gate
in `semi/schema.py`. Full deliverable, acceptance tests, and scope
preserved in §10.

---

### M12 — Mesh input beyond axis-aligned boxes

**Status: Done (M12).** Gmsh `.geo` / `.msh` ingest with content-hash
caching; 2D MOSFET benchmark with Gaussian n+ source/drain implants.
Full deliverable, acceptance tests, and scope preserved in §10.

---

### M13 — Transient solver

**Status: Done (2026-04-25).** Backward-Euler and BDF2 time stepping
on the coupled (psi, phi_n, phi_p) system; `pn_1d_turnon` benchmark.
Full deliverable, acceptance tests, and scope preserved in §10.

---

### M13.1 — Slotboom transient close-out

**Status: Done (2026-04-27).** Slotboom primary unknowns (ADR 0014,
supersedes ADR 0009) replace the M13 (n, p) Galerkin discretisation;
both `test_transient_steady_state_limit` and the BDF rate tests are
active and passing as of v0.14.1; the deep-steady-state runner
matches `bias_sweep` within 1e-4 relative error and the
`pn_1d_turnon` benchmark is within 5%.

---

### M14 — AC small-signal analysis

**Status: Done (2026-04-26).** Linearised
`(J + jωM) δu = -dF/dV δV` in real 2x2 block form (PETSc-real build);
`rc_ac_sweep` benchmark matches analytical depletion C within 0.4%
over [1 Hz, 1 MHz]; ADR 0011 with sign-convention errata. Full
deliverable, acceptance tests, and scope preserved in §10.

---

### M14.1 — AC differential capacitance for MOSCAP

**Status: Done (2026-04-26).** `semi/runners/mos_cap_ac.py` returns
dQ/dV directly from `Im(Y) / (2πf)`, replacing the noisier
`numpy.gradient(Q, V)` of `mos_cv` for the `mos_2d` C-V verifier;
audit case 03 confirms byte-identity with `mos_cv` Q_gate at all 42
gate voltages.

---

### M14.2 — Axisymmetric (cylindrical) 2D MOSCAP

**Status: Done (2026-04-30).** Top-level `coordinate_system` field
(schema 1.3.0, `cartesian` default or `axisymmetric`) with cross-
field validation; r-weighted Poisson and Slotboom drift-diffusion on
the meridian half-plane (`semi/physics/axisymmetric.py`); pure-Python
MOSCAP analytical helpers (`semi/cv.py`); `benchmarks/moscap_axisym_2d/`
reproduces Hu Fig. 5-18; runner dispatch in `mos_cap_ac.py` r-weights
the charge and sensitivity forms (PRs #64, #65).

---

### M15 — GPU linear solver path

**Why:** Above ~200k DOFs the MUMPS LU factorization dominates wall time. 3D
FinFET-class problems will be in the millions of DOFs and CPU-LU is hopeless.
GPU-accelerated iterative solvers change the scaling law.

**What can actually run on GPU** (be realistic about this):

| Stage                          | GPU-suitable?                | Library                       |
|--------------------------------|------------------------------|-------------------------------|
| Mesh build, tagging            | No (serial graph operations) | n/a                           |
| UFL → assembly                 | Yes, with dolfinx-cuda/KOKKOS| dolfinx-cuda (experimental)   |
| Jacobian assembly              | Partial; the dominant O(N)   | dolfinx-cuda                  |
| Linear solve (Krylov)          | **Yes — biggest win**        | PETSc + HIP/CUDA, AMGX, hypre |
| Preconditioner setup (AMG)     | Yes                          | AMGX, hypre-BoomerAMG         |
| Nonlinear SNES loop overhead   | No (negligible cost)         | n/a                           |
| I/O (writing fields)           | No                           | n/a                           |

**Realistic plan:** the 80/20 is "PETSc built against CUDA/HIP, Jacobian
solved on GPU via AMGX or hypre on device." Do not try to port assembly to
GPU yet — dolfinx-cuda is experimental and moving fast.

**Deliverable:**

- `semi/solver.py` grows a `backend: str = "cpu-mumps" | "gpu-amgx" | "gpu-hypre"`
  parameter.
- CPU path untouched, bit-identical to current output.
- GPU path: Jacobian assembled on CPU, copied to GPU, solved with AMGX (or
  PETSc-GAMG-on-device); solution copied back.
- Benchmark: 3D FinFET-scale mesh (≥500k DOFs) showing ≥5× speedup over
  CPU-MUMPS on an A100 or MI300X.
- CI covers CPU path only; GPU path has a separate self-hosted runner or is
  left to nightly.

**Acceptance tests:**

1. On the existing benchmarks, GPU backend produces results within 1e-8
   relative L2 of the CPU-MUMPS reference.
2. On the new large-3D benchmark, wall-time speedup ≥5× for the linear-solve
   portion alone.
3. Fallback: GPU backend gracefully degrades to CPU if no device is found,
   with a warning in the manifest.

**Estimated scope:** 5–7 days, dominated by environment setup.

**Dependencies:** M9 (needs artifact writer for timing instrumentation).

---

### M16 — Physics completeness pass

**Why:** To mimic COMSOL Semiconductor properly we need the models engineers
actually use daily.

**Deliverable (do these in order, each a separate PR):**

1. **Caughey-Thomas field-dependent mobility.** Closed-form velocity-saturation
   model. ~200 LOC.
2. **Lombardi surface-mobility model.** Needed for realistic MOSFET I-V in
   the inversion regime. ~300 LOC.
3. **Auger recombination.** `R_Auger = (C_n n + C_p p)(np - n_i^2)`.
   ~100 LOC.
4. **Fermi-Dirac statistics.** Replace `exp(±psi)` with Fermi integrals
   (use Blakemore approximation for speed, full FDI for validation). Gate
   with a `statistics: "fermi-dirac"` schema option; Boltzmann remains the
   default. ~400 LOC.
5. **Schottky contacts.** Thermionic-emission boundary condition. ~300 LOC.
6. **Band-to-band (Kane) and trap-assisted (Hurkx) tunneling.** Required for
   any useful diode or floating-gate simulation. ~600 LOC.

Each ships with its own benchmark and verifier.

**Acceptance tests:** one analytical-or-SPICE comparison per model.

**Estimated scope:** 2–4 days per item. Don't batch; one PR per model.

**Dependencies:** M9, M11.

---

### M17 — Heterojunction / position-dependent band structure

**Why:** AlGaAs/GaAs, InGaAs/InP, SiGe/Si — none of this works without
position-varying electron affinity `chi(x)` and bandgap `Eg(x)`.

**Deliverable:** extend the material model and Poisson/continuity kernels
to read `chi` and `Eg` as cellwise DG0 fields. Benchmark: AlGaAs/GaAs HEMT
or SiGe HBT.

**Estimated scope:** 4 days.

**Dependencies:** M9, M16.4 (Fermi-Dirac, because heterojunctions break the
nondegenerate approximation at the barrier).

---

### M18 — UI: first cut

**Why:** At this point the engine has everything; the UI is what turns this
into a product.

**Not in scope of this repo.** Mentioned here so the engine team knows what
invariants the UI relies on.

Recommended stack: React + TypeScript, three.js for 3D mesh visualization,
vtk.js or plotly for field rendering, JSONForms for schema-driven device
setup. The UI is a separate repo that consumes this engine's JSON schema
and HTTP API.

---

## 5. GPU section (expanded)

### What is worth accelerating

Profile first, then accelerate. But a general ordering of the biggest wins:

1. **Linear solve (Krylov + AMG preconditioner).** 60–85% of wall time on
   problems above ~100k DOFs. PETSc has first-class CUDA and HIP backends;
   AMGX (NVIDIA) and rocALUTION (AMD) are drop-in options. **Biggest single
   lever.**
2. **Jacobian assembly.** 10–25% of wall time. dolfinx-cuda (still
   experimental as of 2026) handles this via batched cellwise kernels. The
   UFL-to-GPU-kernel path is not yet production-grade; revisit in a year.
3. **Field interpolation + post-processing.** Minor (<5%) but trivially GPU
   with CuPy for things like computing `E = -grad psi` on a dense grid.
4. **Multi-job parallelism.** Less "GPU" and more "throughput": run 100
   parameter-sweep jobs, each on its own small GPU slice. MIG partitions on
   A100/H100 or MPS on older cards. This is actually the most impactful
   pattern for a parameter-exploration GUI.

### What is *not* worth GPU-ifying (yet)

- **Meshing.** Gmsh is not GPU-friendly; time spent here is negligible
  anyway.
- **SNES nonlinear outer loop.** Host-side Python; totally negligible.
- **I/O.** Writing XDMF/BP is CPU-bound and fine that way.
- **Schema validation, JSON parsing, job dispatch.** Obviously host.

### Concrete engineering path

Phase 1 (M15): PETSc built with `--with-cuda` or `--with-hip`, add
`mat_type: aijcusparse` and `vec_type: cuda` as solver options. Use AMGX
(nvidia/AMGX) as the preconditioner via `-pc_type amgx`. Engine exposes a
`solver.backend` field. All assembly stays CPU; matrix and vectors are
moved to device for the Krylov solve.

Phase 2 (post-M15, optional): experiment with dolfinx-cuda for assembly.
Only if Phase 1 profiling shows assembly is now the bottleneck.

Phase 3 (post-M18, throughput): MIG partitioning for parameter sweeps.
Relevant when the UI ships and users run design-of-experiments studies.

### Honest caveat

PETSc-on-GPU is mature for structured problems (Poisson, elliptic). For
coupled drift-diffusion with exponential nonlinearities, the Jacobian
conditioning is terrible and AMG can stall. Expect to spend real time tuning
preconditioner parameters per problem class. A 5× speedup is realistic; a
50× speedup is marketing.

---

## 6. UI integration checklist (what the UI team will ask for)

The engine is ready for a UI when all of these are green:

- [x] M9 complete: runs produce a stable on-disk artifact.
- [x] M10 complete: HTTP API exposes solve, status, stream, artifact fetch.
- [x] M11 complete: schema is versioned, richly annotated, published.
- [x] `GET /materials` returns the material DB.
- [x] `GET /schema` returns the input schema.
- [ ] `GET /capabilities` returns which physics models this engine build
      supports (e.g., `{"transient": true, "ac": false, "gpu": true, ...}`).
      (M15 will close this; GPU availability is a capability.)
- [ ] WebSocket stream emits at least one message per SNES iteration with
      `{"residual_norm": ..., "iteration": ..., "bias_step": ...}`.
- [ ] CORS configured on the HTTP server for the UI's dev origin.
- [ ] OpenAPI spec auto-generated (FastAPI does this for free).
- [ ] At least one end-to-end integration test where a browser-like client
      submits a JSON and visualizes the result manifest.

---

## 7. What coding agents should do when picking up work here

1. **Read PLAN.md and ROADMAP.md before touching code.** "Current state" and
   "Next task" are the source of truth for what is done.
2. **Pick exactly one milestone.** Do not mix M9 and M10 in a single PR. Do
   not start M11 before M9 is merged.
3. **Write the acceptance test first.** Every milestone above states its
   acceptance test explicitly. Make it fail, then make it pass.
4. **Respect the layer boundaries.** Layer 3 (pure-Python core) never imports
   dolfinx. The new `kronos_server/` package is Layer 6 and never imports
   PETSc or UFL types into its public API.
5. **Do not rewrite the physics kernels.** The UFL forms are correct and
   MMS-verified. If you think you see a bug, write a failing MMS test first;
   90% of the time the "bug" is a scaling-convention confusion that the
   existing tests already cover.
6. **Never skip `make_scaling_from_config`.** The raw equations have a
   Jacobian condition number >1e30. This is not negotiable.
7. **Use `snes_rtol: 1.0e-14` only in tests.** Production defaults are looser
   and documented in `DEFAULT_PETSC_OPTIONS`.

---

## 8. Anti-goals (things to actively *not* do)

- **Do not build a bespoke mesh format.** Use XDMF/ADIOS2-BP and gmsh. They
  work, they're standard, the VTK/ParaView/vtk.js ecosystem reads them for
  free.
- **Do not embed a physics DSL.** The JSON schema is the DSL. If it gets
  inadequate, extend the schema; do not invent a second configuration layer.
- **Do not pickle `SimulationResult`.** It contains live PETSc handles; pickles
  break across processes. Always go through `write_artifact`.
- **Do not add `requests`-based auth to the HTTP server before a real
  deployment need.** Rate limiting and API keys can wait until there is a
  second user.
- **Do not port the whole codebase to Rust or C++ for performance.** The
  bottleneck is the linear solve, not Python. M15 addresses it in the right
  place.
- **Do not claim COMSOL parity.** COMSOL has 30 years of model variations,
  convergence hacks, and user-reported-edge-cases. kronos-semi will have a
  useful, correct subset. Be honest about that in marketing.

---

## 9. Change log for this document

- **2026-04-23**, initial version, written post-M8.
- **2026-04-30**, refresh of §1, §4, §6 to reflect v0.14.x reality
  post-M14.2; introduce §10 shipped-milestone appendix.

---

## 10. Shipped milestone detail

This appendix preserves the original §4 deliverable, acceptance, and
scope blocks for milestones that have shipped. The §4 entries above
carry only a one-line status badge; full text lives here so that
acceptance criteria remain available as regression references.

### M9 — Result artifact writer + manifest

**Why:** Nothing downstream works without a stable, machine-readable output.
This is the single highest-value piece of work on the list.

**Deliverable:** New module `semi/io/artifact.py` exposing:

```python
def write_artifact(result: SimulationResult, out_dir: Path, run_id: str) -> Path:
    """Write the full result tree under out_dir/run_id/. Returns the run dir."""
```

Plus a new CLI entry point `semi-run <input.json> [--out runs/]` that wraps
`semi.run.run` and calls `write_artifact`.

**Acceptance tests (gate for merge):**

1. For each of the five shipped benchmarks, `semi-run` produces a run dir
   validated against a new JSON-schema file `schemas/manifest.v1.json`.
2. A pure-Python reader `semi.io.artifact.read_manifest(run_dir)` round-trips.
3. A Python-only (no dolfinx) consumer script can read `fields/psi.bp` via
   `adios2` and `iv/cathode.csv` via `numpy.loadtxt` and plot them.
4. New test file `tests/test_artifact.py` with at least 10 assertions.

**Estimated scope:** 2 days. ~400 LOC + schema + tests.

**Dependencies:** none.

---

### M10 — HTTP server: `POST /solve`, `GET /runs/{id}`

**Why:** This is what lets a browser UI talk to the engine.

**Deliverable:** New top-level package `kronos_server/` (kept out of `semi/` to
preserve the pure-engine boundary):

```
kronos_server/
├── __init__.py
├── app.py              # FastAPI app
├── jobs.py             # job queue + worker spawning
├── storage.py          # LocalFS + S3 backends for run dirs
├── models.py           # pydantic models mirroring schemas/
└── routes/
    ├── solve.py        # POST /solve  -> {run_id, status_url}
    ├── runs.py         # GET  /runs/{id}, GET /runs/{id}/manifest
    └── stream.py       # WebSocket /runs/{id}/stream for SNES progress
```

**Endpoints:**

| Method | Path                       | Purpose                                   |
|-------:|----------------------------|-------------------------------------------|
| POST   | `/solve`                   | Submit JSON, get back `{run_id, status}`  |
| GET    | `/runs/{id}`               | Status + manifest                         |
| GET    | `/runs/{id}/fields/{name}` | Stream a field file                       |
| GET    | `/runs/{id}/iv/{contact}`  | Stream CSV                                |
| WS     | `/runs/{id}/stream`        | Live SNES residual norms during solve     |
| GET    | `/materials`               | Material DB dump for UI autocomplete      |
| GET    | `/schema`                  | The JSON schema itself, for UI validation |

**Worker model:** start with a `ProcessPoolExecutor` of size `N_CPU / 4` (each
worker can use 4-way MPI). Redis + Celery is a later upgrade, not a
requirement at first ship.

**Acceptance tests:**

1. End-to-end: `curl -X POST /solve -d @pn_junction.json` then polling
   `/runs/{id}` yields `status: "completed"` in under 30 s.
2. WebSocket client sees at least `n_steps` progress messages during a bias
   sweep.
3. A second client hitting `/runs/{id}/manifest` while the job runs gets
   `status: "running"`, not a 500.
4. Integration test matrix covers all five benchmarks via the HTTP API.

**Estimated scope:** 3 days. ~700 LOC + tests.

**Dependencies:** M9.

---

### M11 — Schema versioning + UI-facing schema companion

**Why:** A UI form-builder needs introspectable, versioned schemas, not just
a validator.

**Deliverable:**

- Add `schema_version` as a required top-level JSON field; validator refuses
  mismatched majors.
- Split `semi/schema.py` into `schemas/input.v1.json` (Draft-07) plus a
  Python loader that imports it. The JSON-file schema is what the UI pulls
  from `GET /schema` and feeds to a form generator (JSONForms, rjsf, etc.).
- Add `title`, `description`, `default`, `examples`, and `enum` annotations
  to every schema node so form-builders render sensible labels and help text
  without a separate UI config.
- Publish the schema on every git tag to a stable URL (a static GitHub Pages
  site is sufficient).

**Acceptance tests:**

1. `jsonschema-cli validate schemas/input.v1.json` against every benchmark
   JSON passes.
2. All existing tests still pass; the schema refactor is a lift-and-shift.
3. New test asserts every node with `type: object` has a `description`.

**Estimated scope:** 1 day.

**Dependencies:** M9 (so the manifest schema and input schema can share a
layout).

---

### M12 — Mesh input beyond axis-aligned boxes

**Why:** Real devices are not boxes. UI will export gmsh `.geo` or `.msh`
files; current engine support is partial.

**Deliverable:**

- Extend `semi/mesh.py::_build_from_file` to cover XDMF (previously
  `NotImplementedError`).
- Add a `geometry` input alternative to `mesh` in the schema:

```json
"geometry": {
  "source": "gmsh_geo",
  "path": "device.geo",
  "characteristic_length": 5.0e-8,
  "physical_groups": {
    "silicon": {"dim": 2, "tag": 1},
    "oxide":   {"dim": 2, "tag": 2},
    "gate":    {"dim": 1, "tag": 10}
  }
}
```

The engine calls gmsh in a subprocess to produce the `.msh`, caches it by
content hash, then loads it.

- Add a 2D MOSFET benchmark (`benchmarks/mosfet_2d/`) that exercises the new
  geometry input: source/drain/gate/body contacts, non-axis-aligned doping.

**Acceptance tests:**

1. `.geo` to `.msh` to solve to IV curve matches SPICE-level analytical
   approximation within 20% at V_DS = 0.1 V.
2. Cache hit avoids re-meshing on identical `(geo_hash, char_length)`.
3. Mesh-convergence test on the MOSFET benchmark.

**Estimated scope:** 3 days.

**Dependencies:** M9, M11.

---

### M13 — Transient solver

**Why:** Every C-V, capacitance-transient-spectroscopy, and realistic device
characterization involves time. Missing transient is the single biggest
physics gap relative to COMSOL.

**Deliverable:**

- `semi/runners/transient.py` implementing backward-Euler and BDF2 time
  stepping on the coupled (psi, phi_n, phi_p) system.
- Schema extension: `"solver": {"type": "transient", "t_end": 1e-9, "dt": 1e-12}`
  plus adaptive time-stepping options mirroring the bias-continuation options.
- Benchmark `benchmarks/pn_1d_turnon/`: forward-bias step response, verify
  charge storage and RC time constant.

**Acceptance tests:**

1. Steady-state limit of the transient solver matches `run_bias_sweep` within
   numerical noise.
2. Turn-on benchmark: minority-carrier lifetime extracted from the transient
   matches the SRH `tau_n` used as input within 5%.
3. BDF2 convergence rate verified against MMS.

**Estimated scope:** 4 days.

**Dependencies:** M9.

**M13.1 follow-up (closed in v0.14.1).** Slotboom primary unknowns
(ADR 0014) replace the M13 (n, p) Galerkin discretisation; both
`test_transient_steady_state_limit` and the BDF rate tests are active
and passing as of v0.14.1; deep-steady-state runner matches
`bias_sweep` within 1e-4 relative error and `pn_1d_turnon` within 5%.

---

### M14 — AC small-signal analysis

**Why:** Capacitance-voltage, admittance-spectroscopy, and high-frequency
characterization all need AC. The math reuses the transient Jacobian.

**Deliverable:**

- `semi/runners/ac_sweep.py`: at each DC operating point, assemble
  `(J + jωM)` and solve the linearised system for the AC response. The
  PETSc-real build forces a real 2x2 block reformulation
  (`[[J, -ωM], [ωM, J]] [[Re δu], [Im δu]] = [[Re b], [Im b]]`); ADR
  0011 documents the decision and lists complex PETSc as deferred
  M16+ work.
- Terminal admittance includes the displacement current at the
  contact in addition to the linearised conduction current.
- MOS C-V benchmark upgraded to compute true C(V) via AC admittance
  (M14.1 follow-up: `mos_cap_ac` returns dQ/dV directly from
  `Im(Y) / (2πf)`).

**Acceptance tests:**

1. Low-frequency limit of AC capacitance matches `dQ/dV` from the existing
   MOS benchmark within 2%.
2. Frequency sweep produces the expected MOS C-V frequency dispersion for
   interface-trap models (flat until traps are added).
3. `rc_ac_sweep` benchmark matches analytical depletion C within 0.4%
   over [1 Hz, 1 MHz].

**Estimated scope:** 4 days.

**Dependencies:** M9, M13 (for the mass matrix infrastructure).

**M14 sign-convention errata (2026-04-28).** Resolved Phase 1 audit
Cases 02 and 05. `run_ac_sweep` now reports terminal current INTO the
device (matching `bias_sweep` and
`semi/postprocess.evaluate_current_at_contact`); `result.C =
+Im(Y) / (2πf)` is bit-identical to the pre-fix value, and ADR 0011
gained an Errata section formalising the convention. Audit Cases 02
and 05 assertions tightened from NaN-guards to sign-equality plus
1% / 5% relative-error gates.
