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
- Active dev branch: `dev/day8-polish` (Day 8 final polish and
  submission packaging in flight; Day 7 merged via PR #9, `a604b12`)

## Current state

Days 1 through 7 are merged into `main`. Day 8 (final polish and
submission packaging) is in flight on `dev/day8-polish`.
CI hardening (branch glob plus Dockerized FEM job) merged on
`ci/docker-benchmark-matrix`. Day 3 (adaptive bias continuation,
Sah-Noyce-Shockley verifier, reverse-bias generation check) merged
via PR #5 from `dev/day3-bias-hardening`. Day 4 (Verification &
Validation suite: MMS-Poisson, mesh convergence, discrete
conservation, MMS for coupled drift-diffusion, CI integration,
documentation) merged via PR #6 from `dev/day4-vnv`. Day 5 (refactor
pass: `semi/bcs.py` extraction, `semi/run.py` split into a 74-line
dispatcher plus `runners/` and `postprocess.py`, coverage to 96.25%
with a 95% CI gate, completed `docs/PHYSICS.md` Section 2.5, ADR 0007)
merged via PR #7 from `dev/day5-refactor`. Day 6 (2D MOS capacitor:
multi-region Poisson over oxide plus silicon, gate contact, continuity
on a semiconductor submesh, C-V verifier matching depletion-approximation
MOS theory within 10% in [V_FB + 0.2, V_T - 0.1] V, multi-region MMS)
merged via PR #8 from `dev/day6-mos-2d`. Day 7 (3D doped resistor with
gmsh `.msh` loader, bipolar-sweep driver path, V-I linearity verifier
at 1%, and 3D slice plots) merged via PR #9 from
`dev/day7-resistor-3d` at `a604b12`.

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
- `mos_2d` benchmark (C-V sweep, Day 6): 500 nm p-type Si / 5 nm SiO2
  capacitor; |C_sim - C_theory|/C_theory < 10% in the depletion-regime
  window [V_FB + 0.2, V_T - 0.1] V (worst 9.25% at V_gate = -0.20 V
  under ideal phi_ms = 0 with psi = 0 at intrinsic level). Sweep
  covers [-0.9, +1.2] V; accumulation and strong-inversion regimes
  are rendered on the plot but excluded from the verifier per
  depletion-approximation scope.
- GitHub Actions runs pure-Python tests on 3.10/3.11/3.12, ruff, and a
  Dockerized FEM job that runs pytest, three benchmark verifiers, and
  the full V&V suite on every push to `main`, `dev/**`, `ci/**`,
  `docs/**`. Job timeout is 15 minutes once the dolfinx image is cached.
- Verification & Validation suite (Day 4): MMS for Poisson (1D linear,
  1D nonlinear, 2D triangles, 2D quad smoke) with finest-pair rates at
  theoretical 2.0/1.0; mesh convergence on `pn_1d` with Cauchy ratios
  >= 1.99x per doubling; discrete conservation (charge on `pn_1d`
  equilibrium at relative 1.5e-17; current continuity on forward and
  reverse bias sweeps); MMS for coupled drift-diffusion (three variants
  x three grids) with all gated block rates >= 1.99.

- `resistor_3d` benchmark (Day 7): 3D rectangular silicon bar (1 um x
  200 nm x 200 nm), uniform n-type N_D = 1e18 cm^-3, two ohmic
  contacts on the x = 0 and x = L faces, symmetric bias sweep in
  [-0.01, +0.01] V. V-I linearity within 1% of
  `R_theory = L / (q N_D mu_n A) = 1115 Ohm` on both the builtin
  `create_box` mesh and the committed gmsh fixture. 3D slice plots
  of psi and |J_n| at the y = W/2 midplane are written per run.
- Gmsh `.msh` loader (Day 7): `semi/mesh.py::_build_from_file` via
  `dolfinx.io.gmsh.read_from_msh`; physical groups in the file are
  returned verbatim as `cell_tags` and `facet_tags`, so `build_mesh`
  bypasses the JSON box-tagger for file-source meshes. XDMF remains
  a clear `NotImplementedError` for a future PR.
- Bipolar (sign-spanning) bias sweep (Day 7): the driver walks
  `V = 0 -> min(V) -> max(V)` with a fresh
  `AdaptiveStepController` on each leg. Unipolar pn-junction and
  MOS sweeps fall through to the original single-endpoint ramp
  unchanged.

### What does not work / not yet built

- Schottky contacts: deferred (Non-goals).
- Field-dependent mobility, Auger/radiative recombination, Fermi-Dirac
  statistics (see Non-goals).

## Next task

**Day 8: final polish and submission packaging.** In flight on
`dev/day8-polish`, targeting PR #10 against `main`.

- **Goal:** make the submission presentation-ready. No new physics,
  no new benchmarks; this is a documentation, notebook, and release
  pass only. Scope creep on the final PR is the submission-day
  failure mode to avoid.
- **Scope, in:**
  - Sync `PLAN.md` "Current state" and `docs/ROADMAP.md` statuses to
    reflect the Day 7 merge and Day 8 in flight.
  - Rewrite the `README.md` status section as an end-of-Day-7
    capability matrix (Days 1-7 shipped across PRs #2-#9).
  - Regenerate `notebooks/01_pn_junction_1d.ipynb` using
    `scripts/build_notebook_01.py` (current-state framing, not
    "Day 1").
  - Author `notebooks/02_pn_junction_bias.ipynb` (Day 2-3 content:
    forward Shockley + reverse SNS sweeps).
  - Author `notebooks/03_mos_cv.ipynb` (Day 6 content: MOS C-V with
    the V_FB + 0.2 verifier-window disclosure in the narrative).
  - Author `notebooks/04_resistor_3d.ipynb` (Day 7 content: bipolar
    sweep, builtin-vs-gmsh comparison).
  - Verify every notebook by actually executing it on Colab; record
    wall time and final plot observation in the PR body. Local
    `jupyter nbconvert --execute` is a convenience check, not a
    verification, because FEM-on-Colab's dolfinx pin can drift from
    the local Docker pin.
  - Add a README "Notebooks" catalog with Colab badges.
  - Append `[0.8.0] - Day 8` to `CHANGELOG.md`.
  - Open PR #10, wait for CI, do not self-merge.
  - Post-merge and only on explicit human prompt, tag `v0.2.0`.
- **Scope, out:** any new physics code, new verifier, new benchmark
  JSON, or schema change. If a notebook reveals a bug, note it in
  the PR body and defer the fix to a Day 9 PR.
- **Preconditions:** Day 7 (PR #9) merged into `main` at `a604b12`.

## Day 8 execution

The full Day 8 narrative, per-phase commit history, reviewer-caught
decisions (MOS ψ-reference convention, verifier window shift, Option-A
Colab gmsh install plus the `libGLU` apt dependency), verification
status (CI, local Docker, Colab QA), and deferred Day 9+ cleanups live
in [`docs/day8-submission-log.md`](docs/day8-submission-log.md). This
PLAN section carries only the phase index; the log is the reference
document.

| Phase | Scope                                                                 | Commit SHA   | Status                          |
|------:|-----------------------------------------------------------------------|--------------|---------------------------------|
| 0     | Verify `main` baseline (pytest 206, V&V 62/62, ruff clean)            | N/A          | N/A (verification-only)         |
| 1     | Sync `PLAN.md` "Current state" and `docs/ROADMAP.md` statuses         | `1d86504`    | Done                            |
| 2     | Rewrite `README.md` status section as end-of-Day-7 capability matrix  | `c7343c9`    | Done                            |
| 3     | Regenerate `notebooks/01_pn_junction_1d.ipynb` (end-of-Day-7 framing) | `8ce6b10`    | Done                            |
| 4     | Author `notebooks/02_pn_junction_bias.ipynb` (Day 2-3 content)        | `f7c83b0`    | Done                            |
| 5     | Author `notebooks/03_mos_cv.ipynb` (Day 6 C-V content)                | `bc3409b`    | Done                            |
| 6     | Author `notebooks/04_resistor_3d.ipynb` (Day 7 content)               | `365da06`    | Done                            |
| 7     | Cross-notebook consistency pass (headings, install cells, artifacts)  | `d52ca20`    | Done                            |
| 8     | Colab QA on all four notebooks, record wall times                     | pending      | In flight (NB04 outstanding)    |
| 9     | CHANGELOG `[0.8.0] - Day 8`, `semi/__init__.py` version bump          | pending      | Pending                         |
| 10    | Open PR #10, wait for CI, do not self-merge                           | pending      | Pending                         |
| 11    | Post-merge and only on explicit prompt, tag `v0.2.0`                  | pending      | Pending                         |

## Roadmap

| Day | Milestone                                                    | Status     | Notes                                                                 |
|----:|--------------------------------------------------------------|------------|-----------------------------------------------------------------------|
| 1   | Equilibrium Poisson, 1D pn junction, Docker env              | Done       | 6/6 verifier checks pass; PR `dev/docker-day1-fix`                    |
| 2   | Slotboom drift-diffusion, coupled Newton, bias sweep         | Done       | 6/6 `pn_1d_bias` checks pass; PR `dev/day2-drift-diffusion`           |
| 3   | Bias ramping continuation, Shockley IV verifier hardening    | Done       | Adaptive ramp (-26.2% iters), SNS verifier, reverse-bias gen check    |
| 4   | Verification & Validation suite (MMS, conv, conservation)    | Done       | `dev/day4-vnv`; all four phases green, CI V&V step within 15 min      |
| 5   | Refactor pass, expanded test coverage, physics docs updates  | Done       | `dev/day5-refactor`; run.py 580->74, bcs.py extracted, coverage 96.25% |
| 6   | 2D MOS capacitor (oxide + silicon multi-region)              | Done       | `dev/day6-mos-2d`; mos_cv runner, 4/4 C-V checks green, coverage 95.43% |
| 7   | 3D doped resistor                                            | Done       | Merged via PR #9 (`a604b12`): gmsh loader, bipolar sweep, V-I 1%      |
| 8   | Final polish, submission packaging                           | In flight  | `dev/day8-polish` (PR #10): README, 4 notebooks, CHANGELOG, tag prep  |

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

## Day 9+ cleanups

Pre-existing doc/code inconsistencies discovered during the Day 8
polish pass and deferred per the Day 8 anti-pattern rule against
expanding scope on the final submission PR. None of these affect any
verifier or benchmark result; they are all surface-level alignment
items.

- **`docs/mos_derivation.md` §6 ψ-reference convention.** Section 6
  derives `psi_s = psi(x, y_int) - psi(x, y_bulk)` (psi referenced to
  the bulk Fermi level), but the shipped MOS code uses the project-
  wide `psi = 0` at the intrinsic Fermi level convention enforced by
  the ohmic-contact equilibrium BC. Both conventions yield the same
  C(V) curve once V_FB is computed against the matching reference;
  the Day 6 verifier uses the shipped code's convention and passes
  at 9.25 percent worst-case in [V_FB + 0.2, V_T - 0.1] V. Action:
  rewrite §6 in the intrinsic-reference convention so the derivation
  reads against the same baseline as the code, and add a short
  appendix mapping between the two conventions for readers who learn
  MOS from textbooks (Sze, Pierret) that use the bulk reference.
- **`semi/__init__.py::__version__` still reads `"0.1.0"`.** This was
  the Day 1 placeholder and never bumped through Days 2 through 7.
  Action: bump to `"0.8.0"` as part of Phase 9 (CHANGELOG `[0.8.0] -
  Day 8` entry) so the package version, the changelog header, and
  the eventual `v0.2.0` tag (post-merge per the Day 8 prompt's
  Phase 11) are internally consistent. Note the discrepancy between
  package SemVer `0.8.0` and release tag `v0.2.0`: the tag tracks
  the *submission* version (v0.1.0 = Day 1 baseline, v0.2.0 = Day 8
  submission), the package SemVer tracks days-of-work shipped.

## Completed work log

Append-only. Newest entries on top.

- **Day 7 (2026-04-21):** 3D doped resistor benchmark, first
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
    to Day 6; V&V suite green; coverage stays >= 95%.

- **Day 6 (2026-04-21):** 2D MOS capacitor (first 2D benchmark, first
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
    `pn_1d_bias_reverse`) are byte-identical to Day 5; V&V suite
    clears all gates including the new multi-region MMS; pytest
    195/195 passes; coverage 95.43% (95% CI gate passes).

- **Day 5 (2026-04-21):** Refactor pass, coverage to 95%+, completed
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
    `pn_1d_bias_reverse` byte-identical to Day 4 baseline; V&V 53
    PASS / 0 FAIL with rates byte-identical. PR #7.

- **Day 4 (2026-04-21):** Verification & Validation suite. Replaces
  the originally-planned refactor (pushed to Day 5) with four
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
