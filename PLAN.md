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
- Active dev branch: `dev/docker-day1-fix` (Day 1 complete, pending merge)
- Docs branch: `docs/planning-scaffolding` (this change)

## Current state

### What works (verified in Docker)

- Docker dev environment on `ghcr.io/fenics/dolfinx/dolfinx:stable`
  (dolfinx 0.10). `docker compose run --rm test` runs 36/36 pytest.
- Pure-Python core: constants, materials (Si, Ge, GaAs, SiO2, HfO2, Si3N4),
  nondimensional scaling, doping profiles (uniform/step/gaussian), JSON
  schema with validate/load. 36/36 pytest passing.
- JSON input schema (jsonschema Draft-07) with defaults.
- Builtin mesh generation (interval/rectangle/box) with region and facet
  tagging from axis-aligned boxes and planes.
- Equilibrium Poisson under Boltzmann statistics, solved via PETSc SNES
  through `dolfinx.fem.petsc.NonlinearProblem`. Quadratic Newton
  convergence observed (5 iterations, residual 1.5e-6 to 8.1e-15).
- Benchmark runner CLI (`scripts/run_benchmark.py`) producing plots and
  running registered physical verifiers. `pn_1d` benchmark passes all
  six checks:

  | Check                                      | Sim                 | Theory             | Rel err |
  |--------------------------------------------|--------------------:|-------------------:|--------:|
  | Built-in voltage V_bi                      | 0.8334 V            | 0.8334 V           | 0.00%   |
  | Peak \|E\|                                 | 104.84 kV/cm        | 113.53 kV/cm       | 7.65%   |
  | p-side bulk hole density                   | 1.00e23 m^-3        | N_A = 1.00e23      | ratio 1.00 |
  | n-side bulk electron density               | 1.00e23 m^-3        | N_D = 1.00e23      | ratio 1.00 |
  | Mass action n*p in p-side bulk             | 1.000e32            | n_i^2 = 1.000e32   | 0.00%   |
  | Mass action n*p in n-side bulk             | 1.000e32            | n_i^2 = 1.000e32   | 0.00%   |

### What does not work / not yet built

- Drift-diffusion under applied bias (Day 2 target).
- Slotboom variable formulation (Day 2).
- SRH recombination kernel (Day 2).
- Coupled block Newton for (psi, Phi_n, Phi_p) (Day 2).
- Bias ramping continuation (Day 3).
- Forward-bias IV curve verification against Shockley diode equation
  (Day 3).
- 2D MOS capacitor benchmark (Day 5).
- 3D doped resistor benchmark (Day 6).
- Gmsh `.msh` mesh loader (stubbed, raises `NotImplementedError`).
- Gate contacts (`type: "gate"`), Schottky contacts.
- Field-dependent mobility, Auger/radiative recombination, Fermi-Dirac
  statistics (see Non-goals).

## Next task

**Day 2: Slotboom drift-diffusion and bias support.**

- **Branch:** `dev/day2-drift-diffusion` (to be created off `main` after
  `dev/docker-day1-fix` merges).
- **Preconditions:**
  - `dev/docker-day1-fix` is merged to `main`.
  - `docker compose run --rm benchmark pn_1d` exits 0 on the fresh `main`.
  - `docker compose run --rm test` is 36/36 green on the fresh `main`.
- **Scope, in:**
  - Implement Slotboom quasi-Fermi variables (Phi_n, Phi_p) with
    Boltzmann carrier expressions `n = n_i exp((psi - Phi_n) / V_t)`,
    `p = n_i exp((Phi_p - psi) / V_t)`.
  - Add SRH recombination term with tau_n, tau_p from JSON.
  - Build a coupled (psi, Phi_n, Phi_p) block residual using dolfinx
    0.10 blocked function space support
    (`NonlinearProblem(..., kind="nest")` or equivalent).
  - Extend the schema with optional `recombination` (SRH parameters) and
    per-contact applied bias sweeps.
  - Add a bias ramping driver that reuses the previous solution as the
    initial guess.
  - New benchmark `benchmarks/pn_1d_bias/` with a forward-bias IV sweep.
  - New verifier: Shockley diode equation `J = J_0 (exp(V/V_t) - 1)`
    matched within 10% over V in [0.2, 0.6] V at 0.05 V steps.
  - Unit tests for Slotboom math and the SRH kernel.
- **Scope, out:**
  - Reverse bias / breakdown.
  - Field-dependent mobility or Caughey-Thomas.
  - AC small-signal.
  - 2D and higher (stays Day 5+).
  - Notebook updates (deferred to after the code runs clean in Docker).
- **Acceptance criteria:**
  - `docker compose run --rm benchmark pn_1d` still green (no
    regression on Day 1).
  - `docker compose run --rm benchmark pn_1d_bias` exits 0 with IV
    curve within 10% of Shockley across the specified bias range.
  - pytest green with new Slotboom and SRH tests (target: at least 10
    new tests).
  - No em dashes in any new prose.
  - No new non-dolfinx imports in the pure-Python core.
- **Estimated effort:** one focused work day in Docker.

## Roadmap

| Day | Milestone                                                    | Status  | Notes                                                                 |
|----:|--------------------------------------------------------------|---------|-----------------------------------------------------------------------|
| 1   | Equilibrium Poisson, 1D pn junction, Docker env              | Done    | 6/6 verifier checks pass; PR `dev/docker-day1-fix`                    |
| 2   | Slotboom drift-diffusion, coupled Newton, bias sweep         | Planned | Target benchmark `pn_1d_bias` vs Shockley                             |
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

- **Day 1 (2026-04-20):** equilibrium Poisson, 1D pn junction benchmark,
  Docker dev environment, benchmark runner CLI, physics coefficient fix
  (scaled Poisson LHS uses `L_D^2 = lambda2 * L_0^2`, not `lambda2`, since
  the mesh is in physical meters). All 6 `pn_1d` verifier checks pass.
  PR: `dev/docker-day1-fix`.
