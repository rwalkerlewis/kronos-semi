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
- Active dev branch: `dev/day4-vnv` (Day 4 V&V suite in flight)

## Current state

Day 1, Day 2, and Day 3 are merged into `main`. CI hardening (branch
glob plus Dockerized FEM job) merged on `ci/docker-benchmark-matrix`.
Day 3 (adaptive bias continuation, Sah-Noyce-Shockley verifier,
reverse-bias generation check) merged via PR #5 from
`dev/day3-bias-hardening`. Day 4 has been re-scoped from a refactor pass
to a Verification & Validation suite (MMS, mesh convergence,
conservation checks); the original refactor is pushed to Day 5.

### What works (verified in Docker on current `main`)

- Docker dev environment on `ghcr.io/fenics/dolfinx/dolfinx:stable`
  (dolfinx 0.10). `docker compose run --rm test` runs 70/70 pytest.
- Pure-Python core: constants, materials (Si, Ge, GaAs, SiO2, HfO2, Si3N4),
  nondimensional scaling, doping profiles (uniform/step/gaussian), JSON
  schema with validate/load.
- JSON input schema (jsonschema Draft-07) with defaults, including
  `recombination.E_t`, per-contact `voltage_sweep`, and
  `solver.continuation.{min_step, max_halvings}`.
- Builtin mesh generation (interval/rectangle/box) with region and facet
  tagging from axis-aligned boxes and planes.
- Equilibrium Poisson under Boltzmann statistics, solved via PETSc SNES
  through `dolfinx.fem.petsc.NonlinearProblem`. Quadratic Newton
  convergence (5 iterations, residual 1.5e-6 to 8.1e-15).
- Coupled Slotboom drift-diffusion with SRH recombination. Three-block
  (psi, phi_n, phi_p) residual with `L_D^2 * eps_r` on Poisson and
  `L_0^2 * mu_hat` on continuity. Forward-bias sweep with snapshot/
  restore adaptive halving continuation; UFL facet-integral current
  evaluation.
- Benchmark runner CLI (`scripts/run_benchmark.py`) producing plots and
  running registered verifiers.
- `pn_1d` benchmark (equilibrium): V_bi, peak |E|, bulk densities, mass
  action, all within 5-10% of depletion-approximation theory.
- `pn_1d_bias` benchmark (forward bias): J at V = 0.6 V within 10% of
  Shockley long-diode theory; J_sim within 15% of the Sah-Noyce-Shockley
  total (J_diff + J_rec, with the Sze f = 2 V_t/(V_bi - V) correction)
  on V in [0.15, 0.55] V. The original Day 2 "qualitative below 0.5 V"
  caveat was retired by the Day 3 SNS verifier.
- `pn_1d_bias_reverse` benchmark (reverse bias): |J| within 20% of the
  net SRH generation current (q n_i / 2 tau_eff)(W(V) - W(0)) on V in
  [-2, -0.5] V.
- GitHub Actions runs pure-Python tests on 3.10/3.11/3.12, ruff, and a
  Dockerized FEM job that runs pytest plus both benchmark verifiers on
  every push to `main`, `dev/**`, `ci/**`, `docs/**`.

### What does not work / not yet built

- Verification & Validation suite (Day 4, in flight on
  `dev/day4-vnv`): MMS for Poisson and DD, mesh convergence study,
  current/charge conservation checks, CI integration.
- Refactor of `semi/run.py` and BC construction (Day 5).
- 2D MOS capacitor benchmark (Day 6).
- 3D doped resistor benchmark (Day 7).
- Gmsh `.msh` mesh loader (stubbed, raises `NotImplementedError`).
- Gate contacts (`type: "gate"`), Schottky contacts.
- Field-dependent mobility, Auger/radiative recombination, Fermi-Dirac
  statistics (see Non-goals).

## Next task

**Day 4: Verification & Validation suite.** In flight on
`dev/day4-vnv`. Replaces the originally-planned refactor pass; the
refactor is pushed to Day 5.

- **Branch:** `dev/day4-vnv` (cut from `main` after PR #5 merged).
- **Rationale:** the existing `pn_1d` / `pn_1d_bias` /
  `pn_1d_bias_reverse` "verifiers" are single-point physical sanity
  tests, not verification in the Roache / Oberkampf-Roy sense. We have
  no mesh convergence evidence, no manufactured solutions, and no
  conservation checks. Closing this gap is a higher-value use of Day 4
  than refactoring.
- **Scope, in (five phases):**
  - **Phase 1: MMS for Poisson.** New module
    `semi/verification/mms_poisson.py` with 1D and 2D smooth
    manufactured solutions, observed L^2 / H^1 convergence rate
    measurement, CSV+PNG output under `results/mms_poisson/`. Tests in
    `tests/fem/test_mms_poisson.py` (skipped when dolfinx missing).
    Acceptance: L^2 rate >= 1.85 on the finest mesh pair.
  - **Phase 2: Mesh convergence on `pn_1d`.** New module
    `semi/verification/mesh_convergence.py`. Sweep N in
    [50, 100, 200, 400, 800, 1600], record V_bi / peak |E| /
    depletion-width error, Newton iterations, solve time. Acceptance:
    monotone error reduction; documented convergence order.
  - **Phase 3: Conservation checks.** Current continuity check
    (J_n + J_p constant in space, +-5% forward / +-15% reverse) added
    to `pn_1d_bias` and `pn_1d_bias_reverse` verifiers; charge
    conservation check (integrated rho == 0 to solver tolerance) added
    to `pn_1d`.
  - **Phase 4: MMS for coupled drift-diffusion.** New module
    `semi/verification/mms_dd.py`. Three tests: psi-only with frozen
    quasi-Fermis, full three-block coupling, full coupling with SRH.
    Acceptance: L^2 rate >= 1.75 on each field on the finest pair.
  - **Phase 5: CI integration.** Add `run_verification.py all` to the
    `docker-fem` job; upload `results/mms_*` and
    `results/mesh_convergence/` as workflow artifacts.
  - **Phase 6 (Final): Documentation.** New "Verification & Validation"
    section in `docs/PHYSICS.md`, new ADR
    `0006-verification-and-validation-strategy.md`, CHANGELOG update.
- **Scope, out:**
  - Refactor of `run.py` / BC extraction (Day 5).
  - 2D MOS capacitor (Day 6) and 3D resistor (Day 7).
  - Colab notebook updates.
  - `devsim` code-to-code comparison (post-submission).
- **Preconditions:** Day 3 (PR #5) merged into `main`. Done.
- **Hard invariants for this PR:** mesh stays in meters with
  L_D^2 = lambda2 * L_0^2 on Poisson LHS; pure-Python core stays
  dolfinx-free (`semi/verification/` may import dolfinx); dolfinx 0.10
  API only; Slotboom variables for DD; no em dashes in new prose;
  no Colab work.

## Roadmap

| Day | Milestone                                                    | Status     | Notes                                                                 |
|----:|--------------------------------------------------------------|------------|-----------------------------------------------------------------------|
| 1   | Equilibrium Poisson, 1D pn junction, Docker env              | Done       | 6/6 verifier checks pass; PR `dev/docker-day1-fix`                    |
| 2   | Slotboom drift-diffusion, coupled Newton, bias sweep         | Done       | 6/6 `pn_1d_bias` checks pass; PR `dev/day2-drift-diffusion`           |
| 3   | Bias ramping continuation, Shockley IV verifier hardening    | Done       | Adaptive ramp (-26.2% iters), SNS verifier, reverse-bias gen check    |
| 4   | Verification & Validation suite (MMS, conv, conservation)    | In flight  | `dev/day4-vnv`; replaces refactor; MMS-Poisson, MMS-DD, mesh, CI      |
| 5   | Refactor pass, expanded test coverage, physics docs updates  | Planned    | Ruff-clean, coverage target, ADR review (was Day 4)                   |
| 6   | 2D MOS capacitor (oxide + silicon multi-region)              | Planned    | Uses submesh for carriers; verify C-V curve (was Day 5)               |
| 7   | 3D doped resistor                                            | Planned    | Framework extension; verify Ohmic V-I linearity (was Day 6)           |
| 8   | Final polish, submission packaging                           | Planned    | Regenerate notebooks, tag release (was Day 7)                         |

See `docs/ROADMAP.md` for the full per-day breakdown.

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
   `constants`, `materials`, `scaling`, `doping`, `schema` must import
   without dolfinx so CI and users can run them standalone. See
   `docs/ARCHITECTURE.md`.
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

## Completed work log

Append-only. Newest entries on top.

- **Day 3 (2026-04-21):** adaptive bias continuation, Sah-Noyce-
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
- **Day 2 (2026-04-20):** coupled Slotboom drift-diffusion with SRH
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
  helpers) pass alongside the Day 1 suite (70 total). PR:
  `dev/day2-drift-diffusion`.
- **Day 1 (2026-04-20):** equilibrium Poisson, 1D pn junction benchmark,
  Docker dev environment, benchmark runner CLI, physics coefficient fix
  (scaled Poisson LHS uses `L_D^2 = lambda2 * L_0^2`, not `lambda2`, since
  the mesh is in physical meters). All 6 `pn_1d` verifier checks pass.
  PR: `dev/docker-day1-fix`.
