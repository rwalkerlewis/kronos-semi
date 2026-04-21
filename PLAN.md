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
- Active dev branch: `dev/day3-bias-hardening` (Day 3 in flight)

## Current state

Day 1 and Day 2 are merged into `main`. CI hardening (branch glob plus
Dockerized FEM job) merged on `ci/docker-benchmark-matrix`. Day 3
(adaptive bias continuation, Sah-Noyce-Shockley verifier, reverse-bias
check) is in flight on `dev/day3-bias-hardening`.

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
  Shockley long-diode theory. Low-bias behaviour is qualitative in Day 2
  because depletion-region SRH recombination raises ideality toward 2;
  Day 3 will make that regime quantitative by adding the Sah-Noyce-
  Shockley term.
- GitHub Actions runs pure-Python tests on 3.10/3.11/3.12, ruff, and a
  Dockerized FEM job that runs pytest plus both benchmark verifiers on
  every push to `main`, `dev/**`, `ci/**`, `docs/**`.

### What does not work / not yet built

- Reverse-bias saturation check (Day 3).
- Adaptive step growth (halving works; growth after easy solves does
  not, Day 3).
- Sah-Noyce-Shockley reference curve in the bias verifier (Day 3).
- 2D MOS capacitor benchmark (Day 5).
- 3D doped resistor benchmark (Day 6).
- Gmsh `.msh` mesh loader (stubbed, raises `NotImplementedError`).
- Gate contacts (`type: "gate"`), Schottky contacts.
- Field-dependent mobility, Auger/radiative recombination, Fermi-Dirac
  statistics (see Non-goals).

## Next task

**Day 3: Bias ramping continuation, IV verifier hardening.** In flight.

- **Branch:** `dev/day3-bias-hardening` (cut from `main` at the
  post-CI-hardening commit).
- **Scope, in:**
  - Adaptive continuation step growth: after a configurable number of
    consecutive easy SNES solves, grow the step by a configurable
    factor, bounded above by the continuation max step. Preserve the
    halving-on-failure behavior from Day 2.
  - Extend the Sah-Noyce-Shockley (depletion-region SRH) current term
    into the `pn_1d_bias` verifier so the [0.15, 0.6] V range becomes
    quantitative (15% tolerance on J_diff + J_rec).
  - Reverse-bias saturation check, -2 V to -0.5 V, within 20% of J_0.
  - Document the bias continuation strategy in `docs/PHYSICS.md`.
  - Consider Scharfetter-Gummel box-scheme discretization as an ADR if
    Galerkin Slotboom residuals show accuracy issues at high doping or
    wider bias ranges (not expected in Day 3 scope).
- **Scope, out:**
  - Field-dependent mobility, Auger, transient solver (post-submission).
  - 2D or 3D benchmarks (Days 5-6).
  - Colab notebook updates (separate PR, deferred).
- **Acceptance criteria:**
  - `docker compose run --rm benchmark pn_1d` still green.
  - `docker compose run --rm benchmark pn_1d_bias` green with the
    tightened forward verifier over [0.15, 0.6] V and the new reverse-
    bias saturation check.
  - Adaptive ramp reduces total SNES iterations on the Day 2 forward
    sweep by at least 25% relative to the halving-only baseline.
  - At least 4 new tests (step-growth logic, SNS reference curve,
    schema additions). Total pytest >= 74.
  - `docs/PHYSICS.md` has a "Bias continuation strategy" subsection.
  - Ruff clean; no em dashes in new prose or comments.
- **Preconditions (satisfied):**
  - Day 2 merged to `main`.
  - `ci/docker-benchmark-matrix` merged to `main` (PR #4).
  - 70/70 pytest green on `main` in Docker.
  - Both `pn_1d` and `pn_1d_bias` benchmarks green on `main`.

## Roadmap

| Day | Milestone                                                    | Status  | Notes                                                                 |
|----:|--------------------------------------------------------------|---------|-----------------------------------------------------------------------|
| 1   | Equilibrium Poisson, 1D pn junction, Docker env              | Done    | 6/6 verifier checks pass; PR `dev/docker-day1-fix`                    |
| 2   | Slotboom drift-diffusion, coupled Newton, bias sweep         | Done    | 6/6 `pn_1d_bias` checks pass; PR `dev/day2-drift-diffusion`           |
| 3   | Bias ramping continuation, Shockley IV verifier hardening    | Planned | Tighten tolerance, add reverse bias (saturation region)               |
| 4   | Refactor pass, expanded test coverage, physics docs updates  | Planned | Ruff-clean, coverage target, ADR review                               |
| 5   | 2D MOS capacitor (oxide + silicon multi-region)              | Planned | Uses submesh for carriers; verify C-V curve                           |
| 6   | 3D doped resistor                                            | Planned | Framework extension; verify Ohmic V-I linearity                       |
| 7   | Final polish, submission packaging                           | Planned | Regenerate notebooks, tag release                                     |

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
