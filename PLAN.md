# PLAN

Authoritative planning document for kronos-semi. Human and AI contributors
should read this file before writing any code.

## How to use this document

- **Read before writing code.** The "Current state" and "Next task" sections
  are the ground truth for what is done and what to work on next. Do not
  duplicate finished work.
- **Check invariants before architectural changes.** Anything in "Invariants"
  is locked. To change an invariant, open an ADR under `docs/adr/` first,
  get it accepted, and then update this file.
- **Update on completion.** When a task in "Next task" is finished, move its
  summary line to "Completed work log" (append, never delete), replace
  "Current state" with the new reality, and set "Next task" to the following
  item from the roadmap.
- **Never delete history.** The completed work log is append-only.
- **Keep it short.** Target under 300 lines. If a section grows, move the
  detail to `docs/ROADMAP.md` or an ADR and leave a one-line summary here.

## Project one-liner

A JSON-driven finite-element semiconductor device simulator built on
FEniCSx (dolfinx 0.10), solving Poisson-coupled drift-diffusion with SRH
recombination for 1D/2D/3D devices.

## Repository

- URL: https://github.com/rwalkerlewis/kronos-semi
- License: MIT
- Primary branch: `main`
- Active dev branch: `dev/final-housekeeping` (M9: Result artifact writer
  in flight; M8 merged via PRs #10/#11)

## Current state

M1 through M10 are merged into `main`. M9 delivered the on-disk artifact
contract; M10 now exposes the engine over HTTP via a new top-level
`kronos_server/` package. `CHANGELOG.md` records the `[0.10.0] - M10:
HTTP server` entry. The project has delivered the KronosAI evaluation
milestones plus the first two pieces of the production-engine track
(artifact writer and HTTP API).

The capability matrix (verified in CI) is authoritative: see `README.md`
§Status or `docs/ROADMAP.md`.

### What works (verified in Docker on current `main`)

- Everything from M1–M8: Docker dev environment, pure-Python core, JSON
  schema, builtin meshes, gmsh `.msh` loader, equilibrium Poisson,
  coupled Slotboom drift-diffusion, SRH recombination, adaptive bias
  continuation (including bipolar sweeps), five shipped benchmarks
  (pn_1d, pn_1d_bias, pn_1d_bias_reverse, mos_2d, resistor_3d), MMS
  and conservation V&V suite (62/62 PASS), four Colab notebooks,
  GitHub Actions CI.
- **M9 deliverables (v0.9.0):**
  - `schemas/manifest.v1.json` — JSON Schema Draft-07 for run manifests,
    `additionalProperties: false` throughout.
  - `semi/io/artifact.py` — `write_artifact(result, out_dir, run_id)`
    produces a validated on-disk run tree.
  - `semi/io/reader.py` — `read_manifest(run_dir)` runs stdlib-only (no
    numpy, no dolfinx), suitable for M10 subprocess use.
  - `semi/io/cli.py` + `semi-run` console entry point.
  - `tests/test_artifact.py` — 3 pure-Python + 5 FEM round-trip tests
    across all five benchmarks.
  - Verified in Docker: gates 3, 4, 5 and sanity checks S1–S5 all pass
    (see `scripts/m9_verify.sh`).
- **M10 deliverables (new in v0.10.0):**
  - `kronos_server/` top-level package (FastAPI app, pydantic DTOs,
    LocalFS storage, in-process `ProcessPoolExecutor` worker pool,
    file-backed progress event stream).
  - Endpoints: `POST /solve`, `GET /runs`, `GET /runs/{id}`, plus
    `/runs/{id}/{manifest,input,fields/{name},iv/{contact},logs}`,
    `WS /runs/{id}/stream`, `GET /schema`, `GET /materials`,
    `GET /capabilities`, `GET /health`, `GET /ready`. CORS configurable
    via `KRONOS_SERVER_CORS_ORIGINS`. OpenAPI at `/openapi.json` and
    Swagger UI at `/docs` (auto-generated).
  - `kronos-server` console entry point and `server` docker-compose
    service on port 8000.
  - Optional `progress_callback=None` threaded through
    `semi/runners/{equilibrium,bias_sweep,mos_cv}.py` (only change to
    `semi/` in M10).
  - `tests/test_kronos_server.py` — 15 tests covering HTTP API,
    WebSocket progress, and concurrent solves; all pass alongside the
    existing suite (230 total).

### What does not work / not yet built

The gaps between the current state and a production UI-backed engine
are enumerated and sequenced in
[`docs/IMPROVEMENT_GUIDE.md`](docs/IMPROVEMENT_GUIDE.md), milestones
M11 through M18. Summary:

- **No schema versioning.** Input JSON still lacks `schema_version`.
  M11.
- **Mesh input is axis-aligned-only** except for imported gmsh files
  that came pre-meshed. M12.
- **No transient, no AC small-signal.** Steady-state only. M13, M14.
- **Linear solver is CPU-LU only.** Unusable above ~200k DOFs. M15.
- **Physics gaps:** no field-dependent mobility, Auger, Fermi-Dirac,
  Schottky contacts, tunneling, heterojunctions. M16, M17.

## Next task

**M11: Schema versioning + UI-facing schema companion.** Scoped in full
in [`docs/IMPROVEMENT_GUIDE.md`](docs/IMPROVEMENT_GUIDE.md) §4 (M11).

- **Goal:** make the input JSON schema introspectable and versioned so
  a UI form-builder can render labels, defaults, and help text without
  a separate config layer, and so the engine can refuse inputs from an
  incompatible major version.
- **Scope, in:** a required top-level `schema_version` field in every
  input; a standalone `schemas/input.v1.json` file loaded by
  `semi/schema.py`; `title`/`description`/`default`/`examples` added to
  every schema node; an assertion test that every `type: object` node
  has a `description`.
- **Scope, out:** any physics change, M10 server surface changes,
  manifest schema changes.
- **Preconditions:** M10 merged on `main`.

## Roadmap

| Milestone | Summary | Status |
|:-----------------------------|--------------------------------------------------------------|------------|
| M1: Equilibrium Poisson | 1D pn junction, Docker env | Done |
| M2: Coupled drift-diffusion | Slotboom, coupled Newton, bias sweep | Done |
| M3: Adaptive continuation | Bias ramping, Shockley IV verifier | Done |
| M4: V&V suite | MMS, mesh convergence, conservation, CI | Done |
| M5: Refactor and test pass | run.py split, bcs.py extracted, coverage 96.25% | Done |
| M6: 2D MOS capacitor | Oxide + silicon multi-region, C-V sweep | Done |
| M7: 3D doped resistor | gmsh loader, bipolar sweep, V-I 1% | Done |
| M8: Submission polish | Notebooks, catalog, CHANGELOG, v0.8.0 tag | Done |
| M9: Result artifact writer | manifest.json, on-disk field/IV files, semi-run CLI | Done |
| M10: HTTP server | POST /solve, GET /runs/{id}, WebSocket progress | Done |
| M11: Schema versioning | UI-facing schema companion, form-builder annotations | Next |
| M12: Mesh beyond boxes | XDMF loader, gmsh .geo ingress, 2D MOSFET benchmark | Planned |
| M13: Transient solver | Backward-Euler + BDF2, diode turn-on benchmark | Planned |
| M14: AC small-signal | Complex-frequency admittance, true C-V | Planned |
| M15: GPU linear solver | PETSc CUDA/HIP, AMGX preconditioner, 500k-DOF 3D benchmark | Planned |
| M16: Physics completeness | Caughey-Thomas, Lombardi, Auger, FD, Schottky, tunneling | Planned |
| M17: Heterojunctions | Position-dependent chi, Eg; HEMT or HBT benchmark | Planned |
| M18: UI (separate repo) | React + vtk.js + JSONForms, consumes M9+M10+M11 contracts | Out of scope (this repo) |

Per-milestone acceptance tests and scope definitions live in
[`docs/IMPROVEMENT_GUIDE.md`](docs/IMPROVEMENT_GUIDE.md). The M1–M8
delivery history is in [`docs/ROADMAP.md`](docs/ROADMAP.md).

## Invariants

These decisions are locked. Changing any of them requires an accepted ADR
under `docs/adr/`.

1. **JSON is the only supported input format.** No Python DSL, no TOML,
   no YAML. See `docs/adr/0001-json-as-input-format.md`.
2. **Nondimensionalization is mandatory.** The raw equations have a
   Jacobian condition number exceeding 10^30 and Newton fails. See
   `docs/adr/0002-nondimensionalization-mandatory.md`.
3. **Mesh coordinates stay in meters.** Only psi, densities, and
   quasi-Fermi potentials are scaled. The scaled Poisson coefficient is
   therefore `L_D^2 * eps_r`, not `lambda2 * eps_r`. See `docs/PHYSICS.md`
   and ADR 0002.
4. **Pure-Python core must not depend on dolfinx.** The modules
   `constants`, `materials`, `scaling`, `doping`, `schema`,
   `continuation`, `diode_analytical`, and `bcs` must import without
   dolfinx so CI and users can run them standalone. See
   `docs/ARCHITECTURE.md` and ADR 0007.
5. **Target dolfinx 0.10 API only.** Use `dolfinx.fem.petsc.NonlinearProblem`
   with `petsc_options_prefix`. The `dolfinx.nls.petsc.NewtonSolver` class
   is deprecated and must not be introduced. See
   `docs/adr/0003-dolfinx-0-10-api.md`.
6. **Slotboom variables for drift-diffusion.** No SUPG or streamline
   diffusion stabilization of raw continuity equations. See
   `docs/adr/0004-slotboom-variables-for-dd.md`.
7. **Physics-style variable names are allowed and expected.** `N_A`,
   `V_t`, `psi_L`, `eps_Si`, etc. Ruff is configured to accept them
   (N rules disabled). Do not rename to PEP 8.
8. **No em dashes in prose anywhere in the repo.** Use commas, periods,
   or parentheses.
9. **Docker for dev, conda/FEM-on-Colab for users.** The supported dev
   path is `docker compose`. See
   `docs/adr/0005-docker-for-dev-conda-fem-on-colab-for-users.md`.

## Non-goals

The following are explicitly out of scope for the evaluation submission.
They may be added after submission as stretch goals (see
`docs/ROADMAP.md`, "Post-submission" section).

- Field-dependent mobility (Caughey-Thomas, Canali, saturation velocity).
- Auger and radiative recombination.
- Fermi-Dirac statistics.
- Heterojunctions (multiple semiconductor materials with different band
  alignments in the same device).
- Incomplete ionization of dopants.
- AC small-signal analysis.
- Band-to-band or trap-assisted tunneling.
- Transient (time-dependent) solver.
- Full MOSFET with source/drain/gate/body contacts. **Post-submission stretch goal only.**
- FinFET or any 3D transistor geometry. **Post-submission stretch goal only.**
- GUI or web frontend.

## Post-submission cleanups

Pre-existing doc/code inconsistencies discovered during the M8: Submission polish
polish pass and deferred per the M8 anti-pattern rule against
expanding scope on the final submission PR. None of these affect any
verifier or benchmark result; they are all surface-level alignment
items.

- **`docs/mos_derivation.md` §6 ψ-reference convention.** Section 6
  derives `psi_s = psi(x, y_int) - psi(x, y_bulk)` (psi referenced to
  the bulk Fermi level), but the shipped MOS code uses the project-
  wide `psi = 0` at the intrinsic Fermi level convention enforced by
  the ohmic-contact equilibrium BC. Both conventions yield the same
  C(V) curve once V_FB is computed against the matching reference;
  the M6 verifier uses the shipped code's convention and passes
  at 9.25 percent worst-case in [V_FB + 0.2, V_T - 0.1] V. Action:
  rewrite §6 in the intrinsic-reference convention so the derivation
  reads against the same baseline as the code, and add a short
  appendix mapping between the two conventions for readers who learn
  MOS from textbooks (Sze, Pierret) that use the bulk reference.
- **`semi/__init__.py::__version__` still reads `"0.1.0"`.** This was
  the M1 placeholder and never bumped through M2 through M7.
  Action: bump to `"0.8.0"` as part of Phase 9 (CHANGELOG `[0.8.0] -
  M8: Submission polish` entry) so the package version, the changelog header, and
  the eventual `v0.2.0` tag (post-merge per the M8 prompt's
  Phase 11) are internally consistent. Note the discrepancy between
  package SemVer `0.8.0` and release tag `v0.2.0`: the tag tracks
  the *submission* version (v0.1.0 = M1 baseline, v0.2.0 = M8
  submission), the package SemVer tracks days-of-work shipped.

## Completed work log

Append-only. Newest entries on top.

- **M10 (2026-04-23):** HTTP server; new `kronos_server/` top-level package
  (FastAPI app, `ProcessPoolExecutor` worker pool, file-backed progress
  stream), `kronos-server` console entry, `server` docker-compose
  service, `[server]` optional extra in `pyproject.toml`, 15 new server
  tests (`tests/test_kronos_server.py`). Added optional
  `progress_callback=None` to `semi/runners/{equilibrium,bias_sweep,mos_cv}.py`
  (only permitted `semi/` change).

- **M9 (2026-04-23):** Result artifact writer; `semi/io/artifact.py`, `schemas/manifest.v1.json`, `semi-run` CLI.

- **M7 (2026-04-21):** 3D doped resistor benchmark, first
  dimension extension; delivered on `dev/day7-resistor-3d` (PR #9):
  - **`docs/resistor_derivation.md`** (pre-approved derivation-lite
    gate): device geometry, analytical ohmic resistance, V-I
    linearity metric with 1% tolerance rationale, 3D slice-plot
    strategy, gmsh loader test strategy (466a820).
  - **`semi/mesh.py::_build_from_file`** wired via
    `dolfinx.io.gmsh.read_from_msh` (e395830, docs in 3abd170);
    physical groups in the file are returned verbatim as
    `cell_tags` and `facet_tags`; the builtin JSON box-tagger is
    bypassed for file-source meshes. XDMF remains a clear
    `NotImplementedError`.
  - **`semi/runners/bias_sweep.py`** two-leg walk for sign-spanning
    bias sweeps (7c947b2, extracted into
    `compute_bipolar_legs` in Phase 6 for unit-testability).
    Unipolar sweeps fall through to the original single-endpoint
    ramp unchanged.
  - **`benchmarks/resistor_3d/resistor.json`**: 1 um x 200 nm x
    200 nm bar, uniform N_D = 1e18 cm^-3, ohmic contacts on the
    x = 0 and x = L faces, 5-point sweep in [-0.01, +0.01] V.
    `resistor_gmsh.json` plus `fixtures/box.geo` / `box.msh`
    committed fixture for the unstructured path (62f52d2).
  - **`scripts/run_benchmark.py`**: `verify_resistor_3d` V-I
    linearity verifier (1% tolerance, zero-bias-noise + sign
    sanity checks); 3D slice plotters for psi and |J_n| at the
    y = W/2 midplane plus an I-V scatter.
  - **Tests**: `tests/fem/test_mesh_gmsh.py` (round-trip,
    physical-group preservation, builtin-vs-gmsh V-I equivalence);
    `tests/fem/test_resistor_3d.py` (coarsened smoke + builtin
    and gmsh production-JSON theory match); pure-Python
    `tests/test_bipolar_sweep.py` pinning the leg computation.
  - **`docs/PHYSICS.md` Section 7** (3D extension notes, ~1 page):
    cites Poisson and drift-diffusion as dimension-independent,
    documents the bipolar-sweep driver path, references
    `docs/resistor_derivation.md` for the full derivation.
  - **CI / verification**: On this branch, commits e395830 and
    3abd170 landed a first pass of the gmsh loader that failed
    the Dockerized FEM CI job (wrong dolfinx namespace); the fix
    in 62f52d2 re-wired `_build_from_file` against
    `dolfinx.io.gmsh.read_from_msh`. 466a820 (Phase 2 derivation)
    and 62f52d2 (Phase 4 benchmark) both greened the full matrix.
    No `main`-branch artifact ever saw a failing build.
  - **No regressions**: 1D and 2D benchmarks (`pn_1d`,
    `pn_1d_bias`, `pn_1d_bias_reverse`, `mos_2d`) byte-identical
    to the M6 baseline; V&V suite green; coverage stays >= 95%.

- **M6 (2026-04-21):** 2D MOS capacitor (first 2D benchmark, first
  multi-region device). Merged via PR #8 from `dev/day6-mos-2d`:
  - **`docs/mos_derivation.md`** (pre-approved derivation gate):
    device geometry, per-region equations, Si/SiO2 interface
    conditions, submesh formulation, gate BC, MOS C-V theory with
    10% verifier rationale, multi-region Poisson MMS construction.
  - **`semi/mesh.py`**: `build_submesh_by_role` via
    `dolfinx.mesh.create_submesh`; DG0 cellwise `eps_r(x)` Function
    on the parent mesh (719b9e8).
  - **`semi/bcs.py`**: gate branch in `build_psi_dirichlet_bcs`
    (Dirichlet `(V_gate - phi_ms) / V_t` on psi, no Slotboom BC);
    `build_dd_dirichlet_bcs` explicitly skips gate contacts
    (b24c650).
  - **`semi/physics/drift_diffusion.py`**: `DDBlockSpacesMR`,
    `make_dd_block_spaces_mr`, `build_dd_block_residual_mr` for
    `V_phi_n`, `V_phi_p` on the submesh with `entity_maps` threading
    the parent<->submesh mapping through `fem.form` and
    `NonlinearProblem` (455196a). Covered by `test_dd_submesh.py`.
  - **`semi/physics/poisson.py`**: `build_equilibrium_poisson_form_mr`
    for multi-region equilibrium Poisson (stiffness on full mesh
    with cellwise eps_r, space-charge restricted to silicon via
    `dx(subdomain_id=semi_tag)`). The scalar single-region path is
    preserved byte-identically.
  - **`semi/runners/mos_cv.py`**: new `run_mos_cv` runner driven
    by `solver.type == "mos_cv"`. Sweeps the gate contact, solves
    at each V_gate, integrates silicon space charge to produce
    `(V_gate, Q_gate)` rows.
  - **`benchmarks/mos_2d/mos_cap.json`**: 500 nm p-type Si
    (N_A = 1e17 cm^-3) / 5 nm SiO2; uniform 1 nm vertical mesh
    (505 cells) respects the Si/SiO2 interface as a grid line;
    V_gate sweep [-0.9, +1.2] V in 0.05 V steps (43 points).
    Ideal gate (phi_ms = 0), V_FB = -0.417 V, V_T = +0.658 V.
  - **`scripts/run_benchmark.py`**: 2D `plot_mos_2d` (tricontourf
    of psi, central-column psi(y), C-V and Q-V curves);
    `verify_mos_2d` (C_sim via centered FD, depletion-approx theory
    via closed-form psi_s inversion, 10% tolerance in
    [V_FB + 0.2, V_T - 0.1] V, monotone non-increasing check). 1D
    plotters untouched.
  - **`semi/verification/mms_poisson.py`**:
    `run_mms_poisson_2d_multiregion` / `run_mr_convergence_study`
    on the Si/SiO2 coefficient jump with eps-weighted flux
    continuity (99fcd2b). Finest-pair rate_L^2 >= 1.99, rate_H^1
    >= 0.95 enforced in `scripts/run_verification.py mms_poisson`.
    Pytest `test_mms_poisson_2d_multiregion_convergence` gates the
    looser >= 1.85 / >= 0.85 thresholds.
  - **`docs/PHYSICS.md`** Section 6 (condensed MOS reference:
    device, equilibrium model, BC-convention shift for V_FB,
    capacitance extraction, verifier result). Section 5.5 adds the
    multi-region Poisson MMS entry.
  - **`tests/fem/test_mos_cv.py`** (4 tests): gate sweep Q_gate
    behaviour, flatband bulk equilibrium, multi-region form
    byte-identity on single-region mesh, malformed-config error
    path.
  - **CI**: `.github/workflows/ci.yml` adds the `mos_2d` benchmark
    step after the three `pn_1d*` benchmarks.
  - **Verifier result**: worst |C_sim - C_theory|/C_theory = 9.25%
    at V_gate = -0.20 V (window edge), inside the fixed 10%
    tolerance. The initial V_FB+0.1 window hit 10.06% at
    V_gate = -0.25 V; per the reviewer's guidance the fix was to
    shrink the window to V_FB+0.2, not loosen the tolerance.
  - **Regressions**: 1D benchmarks (`pn_1d`, `pn_1d_bias`,
    `pn_1d_bias_reverse`) are byte-identical to the M5 baseline; V&V suite
    clears all gates including the new multi-region MMS; pytest
    195/195 passes; coverage 95.43% (95% CI gate passes).

- **M5 (2026-04-21):** Refactor pass, coverage to 95%+, completed
  `docs/PHYSICS.md` Section 2.5. Merged via PR #7 from
  `dev/day5-refactor`:
  - **`semi/bcs.py`** extracted (pure-Python core, no dolfinx at
    module scope). Public API: `ContactBC` dataclass,
    `resolve_contacts(cfg, facet_tags=None, voltages=None)`,
    `build_psi_dirichlet_bcs`, `build_dd_dirichlet_bcs`. Bias-sweep
    drivers pass `voltages={name: V_step}` rather than mutating
    config. Validated against an embedded byte-for-byte copy of the
    legacy inline implementation in `tests/fem/test_bcs.py`. Design
    rationale: ADR 0007.
  - **`semi/run.py` split** from 580 lines to a 74-line dispatcher.
    Solver paths moved to `semi/runners/{equilibrium,bias_sweep}.py`;
    facet/current/IV helpers moved to `semi/postprocess.py`. No module
    in `semi/` exceeds 300 lines (max is `bias_sweep.py` at 294).
    `SimulationResult` and `run(cfg)` stay in `semi.run`;
    `run_equilibrium`, `run_bias_sweep`, `_fmt_tag`, `_resolve_sweep`
    remain importable from `semi.run` via a `__getattr__` shim.
  - **Coverage to 96.25%** across 1598 statements (60 missed). 38 new
    tests across pure-Python and FEM matrices; total 177/177 pass.
    CI `docker-fem` adds `pytest --cov=semi --cov-fail-under=95`.
    `pyproject.toml` excludes pragma markers, gmsh
    `NotImplementedError`, `__main__` guards, and TYPE_CHECKING; the
    V&V CLI driver functions are `# pragma: no cover` because they
    are exercised end-to-end by `scripts/run_verification.py all`.
  - **`docs/PHYSICS.md` Section 2.5** completed with the full scaled
    drift-diffusion derivation. ADR 0007 records the BC interface
    design.
  - **No behavioral regression.** `pn_1d`, `pn_1d_bias`,
    `pn_1d_bias_reverse` byte-identical to the M4 baseline; V&V 53
    PASS / 0 FAIL with rates byte-identical. PR #7.

- **M4 (2026-04-21):** Verification & Validation suite. Replaces
  the originally-planned refactor (pushed to M5) with four
  verification activities, each gated in CI via
  `scripts/run_verification.py`:
  - **MMS for equilibrium Poisson** (`semi/verification/mms_poisson.py`):
    1D linear, 1D nonlinear, 2D triangles sweep with UFL weak-form
    forcing. Finest-pair L^2 rates at theoretical 2.000, H^1 at
    0.999-1.000. `2d_quad_smoke` sanity ratio 0.394.
  - **Mesh convergence on `pn_1d`** (`semi/verification/mesh_convergence.py`):
    N in [50..1600] sweep; Cauchy self-convergence >= 1.99x per
    doubling on E_peak and W over the first four levels; honest-flag
    plateau documented because the depletion-approximation reference
    is itself approximate.
  - **Conservation** (`semi/verification/conservation.py`): charge
    conservation on `pn_1d` equilibrium at rel 1.5e-17; current
    continuity on `pn_1d_bias` (V in {0.30, 0.45, 0.60} V) and
    `pn_1d_bias_reverse` (V in {-0.50, -1.00, -2.00} V) with worst
    max_rel 1.91% forward / 0.020% reverse, inside the 5% / 15% tols.
  - **MMS for coupled drift-diffusion** (`semi/verification/mms_dd.py`):
    three variants (psi-only / full coupling no R / full coupling with
    SRH) x three grids (1D default, 1D nonlinear, 2D); every gated
    block rate >= 1.99 L^2 and >= 0.996 H^1. Derivation artifact:
    `docs/mms_dd_derivation.md`. Implementation-time amendment to the
    derivation: SNES `atol = 0.0`, `stol = 1e-12` (the originally-
    proposed `atol = 1e-16` was below the 2D continuity-block initial
    residual and terminated Newton prematurely).
  - **CLI driver** `scripts/run_verification.py` with subcommands
    `mms_poisson`, `mesh_convergence`, `conservation`, `mms_dd`, `all`.
  - **CI integration:** `docker-fem` job adds a V&V step running
    `run_verification.py all` inside the dolfinx image; timeout
    tightened from 30 to 15 minutes; `benchmark-plots` artifact
    renamed `fem-results` and keeps uploading the full `results/`
    tree.
  - **Documentation:** new "Verification & Validation" section in
    `docs/PHYSICS.md` with current finest-pair rates and honest flags
    (depletion plateau, Variant A scope, block residual scale
    disparity); new ADR
    `docs/adr/0006-verification-and-validation-strategy.md`. PR:
    `dev/day4-vnv`.

- **M3 (2026-04-21):** adaptive bias continuation, Sah-Noyce-
  Shockley recombination in the forward verifier, and a dedicated
  reverse-bias benchmark. Added `semi/continuation.py`
  (`AdaptiveStepController`: grow on easy_iter_threshold easy solves,
  halve on failure, clamp to sweep endpoint) and rewrote
  `run_bias_sweep` to drive it, cutting `pn_1d_bias` forward (0 to
  0.6 V) from 42 SNES iterations to 31 (26.2% reduction). Extracted
  Shockley, depletion-width, SNS, and SRH-generation reference curves
  into `semi/diode_analytical.py` (pure-Python, testable). The
  `pn_1d_bias` verifier now matches J_sim to the SNS total
  (J_diff + J_rec with Sze f = 2 V_t/(V_bi - V) correction) within
  15% on [0.15, 0.55] V and to Shockley diffusion within 10% at
  V = 0.6 V. New `pn_1d_bias_reverse` benchmark sweeps the anode 0 to
  -2 V and demands |J| within 20% of
  (q n_i / 2 tau_eff)(W(V) - W(0)) on [-2, -0.5] V; this replaces the
  originally-proposed "|J| saturates to J_0" acceptance criterion,
  which does not hold for tau = 1e-8 s devices where SRH thermal
  generation dominates reverse current over Shockley diffusion by
  ~5 orders of magnitude. Schema extended with
  `continuation.{max_step, easy_iter_threshold, grow_factor}`; docs
  add a "Bias continuation strategy" subsection to `PHYSICS.md`; CI
  now also runs `pn_1d_bias_reverse`. 26 new tests (continuation
  controller, diode analytical helpers, schema fields). PR:
  `dev/day3-bias-hardening`.
- **M2 (2026-04-20):** coupled Slotboom drift-diffusion with SRH
  recombination and forward-bias sweep. Added `semi/physics/slotboom.py`
  (n/p from (psi, phi_n, phi_p), UFL and NumPy),
  `semi/physics/recombination.py` (SRH with E_t trap level),
  `semi/physics/drift_diffusion.py` (three-block P1 Galerkin residual
  with `L_D^2 * eps_r` on Poisson and `L_0^2 * mu_hat` on continuity);
  extended `semi/solver.py` with `solve_nonlinear_block` (blocked
  NonlinearProblem wrapper) and `semi/run.py` with `run_bias_sweep`
  (adaptive-halving voltage continuation, snapshot/restore on SNES
  failure, UFL facet-integral current evaluation, IV table recording).
  Extended the JSON schema with `recombination.E_t`, per-contact
  `voltage_sweep`, and `solver.type in {drift_diffusion, bias_sweep}`
  plus `solver.continuation.{min_step, max_halvings}`. Added the
  `pn_1d_bias` benchmark (20 um symmetric 1e17 junction, anode swept
  0->0.6 V) and Shockley long-diode verifier. At V=0.6 V the simulated
  current matches Shockley within 10%; at lower bias SRH depletion-
  region current dominates, so the verifier only requires qualitative
  properties there. 34 new tests (schema, Slotboom, SRH, bias BC
  helpers) pass alongside the M1 suite (70 total). PR:
  `dev/day2-drift-diffusion`.
- **M1 (2026-04-20):** equilibrium Poisson, 1D pn junction benchmark,
  Docker dev environment, benchmark runner CLI, physics coefficient fix
  (scaled Poisson LHS uses `L_D^2 = lambda2 * L_0^2`, not `lambda2`, since
  the mesh is in physical meters). All 6 `pn_1d` verifier checks pass.
  PR: `dev/docker-day1-fix`.
