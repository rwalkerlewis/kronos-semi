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
- Active branch: `main` (PR #65 merged via `799934d` on 2026-04-30;
  v0.14.2 tag still pending alongside the package version bump)

## Current state

M1 through M14.2 are merged into `main`. Current package version is
`0.14.1`; M14.2 (axisymmetric MOSCAP, schema 1.3.0) merged via PR #64
(`a4649be`), and the post-merge documentation cleanup, scipy base-
dependency promotion, and axisymmetric runner dispatch in
`semi/runners/mos_cap_ac.py` followed in PR #65 (`799934d`). The
v0.14.2 tag and the corresponding `pyproject.toml` version bump are
outstanding administrative items, deferred to land alongside the M15
close-out (which itself bumps the minor to 0.15.0). M13.1 closed in v0.14.1: the 1D transient runner
uses Slotboom primary unknowns (ADR 0014, supersedes ADR 0009) and
matches bias_sweep at deep steady state. SG flux primitives ship in
`semi/fem/` but are not the active CD path. ADRs 0012 (SG flux) and
0013 (BC-ramp continuation) remain in force. The xfailed
`test_transient_steady_state_limit` and BDF rate tests are now active
and passing.

The capability matrix (verified in CI) is authoritative: see `README.md`
§Status or `docs/ROADMAP.md`.

### What works (verified in Docker on current `main`)

- Everything from M1–M13.1 (see `CHANGELOG.md` for per-version detail).
- **M14 — small-signal AC sweep** (v0.14.0): `semi/runners/ac_sweep.py`
  in Slotboom primary form, real 2x2 block reformulation, displacement
  + conduction current at the contact, schema 1.2.0 (`solver.type =
  "ac_sweep"`, `solver.dc_bias`, `solver.ac`); `benchmarks/rc_ac_sweep`
  matches analytical depletion C within 0.4% over [1 Hz, 1 MHz];
  `tests/fem/test_ac_dc_limit.py` and `tests/mms/test_ac_consistency.py`
  green; ADR 0011 (with errata for the sign-convention fix).
- **M14.1 — differential capacitance via AC admittance** (PR #38):
  `semi/runners/mos_cap_ac.py` returns `dQ/dV` from `Im(Y)/(2πf)` so
  the `mos_2d` C-V benchmark is verified by both the dQ/dV (`mos_cv`)
  and AC paths; byte-identical Q_gate confirmed in audit case 03.
- **M14.2 — axisymmetric (cylindrical) 2D MOSCAP** (PR #64, schema
  1.3.0): top-level `coordinate_system` field accepts `"cartesian"`
  (default) or `"axisymmetric"`; cross-field validation enforces
  `dimension == 2`, non-negative radial extent, and rejects Dirichlet
  contacts on r = 0. r-weighted Poisson and Slotboom forms in
  `semi/physics/axisymmetric.py`. `semi/cv.py` provides pure-Python
  MOSCAP analytical helpers (V_fb, V_t, |phi_B|, W_dmax, C_ox, C_min)
  and LF/HF C-V curves. `benchmarks/moscap_axisym_2d/` reproduces Hu
  Fig. 5-18; `notebooks/05_moscap_axisym_cv.ipynb` runs end-to-end on
  Colab. Axisymmetric dispatch in `mos_cap_ac.py` (PR #65) replaces
  `W_lat` with `L_gate = ∫_gate r ds` and r-weights the charge and
  sensitivity forms.

### What does not work / not yet built

The gaps between the current state and a production UI-backed engine
are enumerated and sequenced in
[`docs/IMPROVEMENT_GUIDE.md`](docs/IMPROVEMENT_GUIDE.md), milestones
M15 through M18. Summary:

- **Linear solver is CPU-LU only.** Unusable above ~200k DOFs. M15.
- **Physics gaps:** no field-dependent mobility, Auger, Fermi-Dirac,
  Schottky contacts, tunneling. M16.
- **Heterojunctions:** position-dependent χ and Eg not yet supported.
  M17.
- **Cartesian-2D MOSCAP variant** and a rigorous AC small-signal HF
  C-V method (driving the gate at high f) are tracked as M14.2.x open
  items in `docs/ROADMAP.md`.

## Next task

**M15: GPU linear solver path.** Add a PETSc-CUDA / PETSc-HIP linear
solver backend behind `solver.backend`, keep the CPU-MUMPS path bit-
identical, ship a 3D Poisson benchmark at >=500k DOFs, and prove a 5x
speedup of the linear-solve portion alone. Full deliverable, phased
plan, and acceptance tests in
[`docs/IMPROVEMENT_GUIDE.md`](docs/IMPROVEMENT_GUIDE.md) §4 M15 and
§5; ordered work plan in
[`docs/M15_STARTER_PROMPT.md`](docs/M15_STARTER_PROMPT.md).

## Backlog

Items deferred behind the next task. Order is informational, not
binding; the owner picks one explicitly when the next task ships.

- **Physics validation suite, Phase 2.** External validation against
  Sze and Nicollian-Brews. Phase 1 is complete: cases 01-04 pass
  internal-consistency checks (audit case 03 confirms `mos_cv` and
  `mos_cap_ac` byte-identity); case 02/05 sign-convention findings
  were resolved in PR #62 (M14 sign fix); case 06 deferred. Audit
  suite is CI-gated via the `docker-fem-audit` job.
- **v0.14.2 tag and `pyproject.toml` version bump.** Outstanding
  administrative items from the PR #65 merge; will fold into the M15
  close-out which bumps to 0.15.0.

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
| M11: Schema versioning | UI-facing schema companion, form-builder annotations | Done |
| M12: MOSFET n+ + SNES amendment | Gaussian implants, relaxed SNES tols, ADR 0008 | Done |
| M13: Transient solver | Backward-Euler + BDF2, diode turn-on benchmark | Done |
| M13.1: Slotboom transient | Slotboom primary unknowns (ADR 0014); xfails closed; v0.14.1 | Done |
| M14: AC small-signal | Linearised (J + jωM) δu = -dF/dV δV; rc_ac_sweep verifier within 0.4% of analytical C_dep; ADR 0011 | Done |
| M14.1: AC differential C-V | `mos_cap_ac` runner returns dQ/dV via Im(Y)/(2πf); audit case 03 byte-identity | Done |
| M14.2: Axisymmetric MOSCAP | Cylindrical 2D path; schema 1.3.0 `coordinate_system`; Hu Fig. 5-18 benchmark | Done |
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

The following are explicitly out of scope for the current release. They
are tracked as stretch milestones (M16/M17) in
[`docs/IMPROVEMENT_GUIDE.md`](docs/IMPROVEMENT_GUIDE.md):

- Field-dependent mobility (Caughey-Thomas, Canali, saturation velocity).
- Auger and radiative recombination.
- Fermi-Dirac statistics.
- Heterojunctions (multiple semiconductor materials with different band
  alignments in the same device).
- Incomplete ionization of dopants.
- Band-to-band or trap-assisted tunneling.
- Schottky contacts.
- Full MOSFET with source/drain/gate/body contacts beyond the M12
  benchmark (M16+).
- FinFET or any 3D transistor geometry (M16+).
- GUI or web frontend (M18, separate repo).

Capabilities previously listed here that have since shipped:
**transient solver** (M13/M13.1), **AC small-signal analysis** (M14),
**axisymmetric (cylindrical) 2D devices** (M14.2).

## Post-merge follow-ups

Small alignment items deferred from the M8 / M14.2 polish passes; none
affect any verifier or benchmark result.

- **`docs/mos_derivation.md` §6 ψ-reference convention.** Section 6
  derives `psi_s = psi(x, y_int) - psi(x, y_bulk)` (bulk-Fermi
  reference), but the shipped MOS code uses the project-wide
  intrinsic-Fermi convention `psi = 0` enforced by the ohmic-contact
  equilibrium BC. Both yield the same C(V) once V_FB matches the
  reference; the M6 verifier passes at 9.25% worst-case. Action:
  rewrite §6 in the intrinsic-reference convention with a short
  appendix mapping the two conventions for textbook (Sze, Pierret)
  readers.

## Completed work log

Append-only. Newest entries on top.

- **M14.2 axisymmetric (cylindrical) 2D MOSCAP (2026-04-30):** Merged
  via PR #64 (`a4649be`). Schema 1.3.0 introduces a top-level
  `coordinate_system` field (`"cartesian"` default, `"axisymmetric"`)
  with cross-field validation: dimension must equal 2, radial extent
  is non-negative, Dirichlet contacts on the symmetry axis r = 0 are
  rejected. New `semi/physics/axisymmetric.py` provides r-weighted
  Poisson and Slotboom drift-diffusion forms on the meridian half-
  plane. New `semi/cv.py` provides pure-Python (no dolfinx) MOSCAP
  analytical helpers (`analytical_moscap_params`, `lf_cv_quasistatic`,
  `hf_cv_depletion_approximation`, FEM postprocessors
  `compute_lf_cv_fem`, `compute_hf_cv_depletion_clamp`). New
  `benchmarks/moscap_axisym_2d/` (gmsh `.geo` for the meridian mesh,
  `moscap_axisym.json` config, analytical `reference_cv.csv`)
  reproduces Hu Fig. 5-18. New `notebooks/05_moscap_axisym_cv.ipynb`
  runs end-to-end on Colab. New tests:
  `tests/check_axisym_moscap_math.py`,
  `tests/test_coordinate_system.py`,
  `tests/test_moscap_axisym_cv.py`,
  `tests/test_moscap_axisym_cv_fem.py`. PR #65 followed up with the
  axisymmetric runner dispatch in `semi/runners/mos_cap_ac.py`
  (replaces `W_lat = extents[0][1] - extents[0][0]` with
  `L_gate = ∫_gate r ds_meridian`, r-weights `charge_form` and
  `sensitivity_form`), promoted scipy from `[test]` to base
  dependencies (used by `semi.cv` and `semi.runners.ac_sweep` at
  runtime), and made the Colab path self-contained.

- **M14.1 differential capacitance via AC admittance (2026-04-26):**
  PR #38. New `semi/runners/mos_cap_ac.py` returns dQ/dV directly
  from `Im(Y) / (2π f)` at low frequency, replacing the noisier
  `numpy.gradient(Q, V)` of `mos_cv` for the `mos_2d` C-V verifier.
  Worst error 6.79% vs the depletion-approximation reference in the
  verifier window. Audit case 03 (`tests/audit/`) confirms `mos_cv`
  and `mos_cap_ac` agree on Q_gate to machine precision (rel_err =
  0.000 at all 42 gate voltages).

- **M14 sign-convention fix (2026-04-28):** Resolved Phase 1 audit
  Cases 02 and 05 (Class C). `run_ac_sweep` was reporting terminal
  admittance Y in the OUT-of-device convention, while `bias_sweep`
  / `semi/postprocess.evaluate_current_at_contact` reports current
  INTO the device. The two linearisations therefore disagreed in
  sign at the same DC operating point. Fixed by negating the
  assembled total terminal current at the assembly site in
  `run_ac_sweep` and flipping the sign of the C read-out
  (`C = +Im(Y)/(2πf)`) so `result.C` is bit-identical and
  `benchmarks/rc_ac_sweep` plus `tests/fem/test_ac_dc_limit.py`
  pass without modification. ADR 0011 grew an Errata section
  documenting the convention as "current INTO device" with
  `Y = +jωC` for an ideal capacitor. Audit Cases 02 and 05
  assertions tightened from NaN-guards to sign-equality + 1% / 5%
  relative-error gates. `docs/PHYSICS_AUDIT.md` Resolution section
  added.

- **M13.1 (2026-04-27):** 1D transient close-out via Slotboom primary
  unknowns (ADR 0014), BC-ramp continuation (ADR 0013), SG primitives
  (ADR 0012), and a Jacobian shift plus MUMPS workspace bump
  (`mat_mumps_icntl_14=200`, shift 1e-14) for the rank-deficient
  phi_n row in the deep p-bulk of the Slotboom lumped-mass time
  term. PR #52's pivot-threshold approach was superseded by PR #54
  (factor-mat options never reached MUMPS in dolfinx 0.10). Both
  previously-xfailed tests (test_transient_steady_state_limit and
  the BDF rate tests) now pass without xfail. PRs #44, #45, #46,
  #47, #48, #49, #50, #51, #52, #54. v0.14.1.

- **M14 (2026-04-26):** Small-signal AC sweep runner. Added
  `semi/runners/ac_sweep.py` solving the linearised system
  (J + jωM) δu = -dF/dV δV around a converged DC operating point. The
  Jacobian J is the steady-state DD operator at u_0 in (ψ, n_hat,
  p_hat) primary-density form; M is the lumped diagonal carrier-density
  mass shipped with M13 (`semi/fem/mass.py`); -dF/dV is obtained via
  finite-difference DC sensitivity (eps_V = 1e-3 V) under the hood.
  PETSc-real build forces a real 2x2 block reformulation
  ([[J, -ωM], [ωM, J]] [[Re δu], [Im δu]] = [[Re b], [Im b]]); ADR
  0011 documents the decision and lists complex PETSc as deferred
  M16+ work. Terminal admittance includes the displacement current
  jω·ε·grad(δψ)·n at the contact in addition to the linearised
  conduction current. Added `AcSweepResult` to `semi/results.py`,
  bumped `schemas/input.v1.json` to 1.2.0 and `SCHEMA_SUPPORTED_MINOR`
  to 2 (new `solver.dc_bias` and `solver.ac` sub-objects). Acceptance
  benchmark `benchmarks/rc_ac_sweep/` matches analytical pn depletion
  C within 0.41 % worst-case at V_DC = -1 V across [1 Hz, 1 MHz];
  `tests/fem/test_ac_dc_limit.py` matches within 5 % at V_DC ∈
  {-2.0, -1.0, -0.5} V; `tests/mms/test_ac_consistency.py` checks
  three internal AC-formulation invariants. Version bumped to
  0.14.0.

- **M12 (2026-04-25):** MOSFET n+ doping + SNES tolerance amendment.
  Added `benchmarks/mosfet_2d/mosfet_2d.json` with three doping entries:
  uniform p-type body (N_A = 1 × 10¹⁶ cm⁻³) plus Gaussian n+ source and
  drain implants (peak 5 × 10¹⁹ cm⁻³, σ = (0.4 µm, 0.15 µm)). Amended
  SNES tolerances in `semi/runners/bias_sweep.py`: `rtol` 1e-14 → 1e-10,
  `atol` 1e-14 → 1e-7, `max_it` 60 → 100 (stol unchanged). Documented
  the decision in `docs/adr/0008-snes-tolerances.md`. Added
  `tests/fem/test_bias_sweep_multiregion.py` with a 3-step forward ramp
  test asserting ≥ 3 IV points and non-decreasing J_n. Version bumped to
  0.12.0. PR closes #25.

- **M11 (2026-04-23):** Schema versioning + UI-facing schema companion.
  Extracted the input schema verbatim from `semi/schema.py` into
  `schemas/input.v1.json` (Draft-07); `semi/schema.py` now loads that
  file through an `@lru_cache`'d loader and exposes the same public
  `SCHEMA`/`validate`/`load` surface. Added required top-level
  `schema_version` (semver) with a major-version gate in
  `semi.schema.validate` (`ENGINE_SUPPORTED_SCHEMA_MAJOR = 1`); minor/
  patch skew is accepted silently. Annotated every `type: object`
  node with a UI-facing description (22/22) and hand-copied
  `_fill_defaults` defaults into the schema via `default`/`examples`.
  `GET /schema` now returns `{schema, version, supported_major}`
  instead of the bare schema. Added a symmetric manifest major-version
  gate in `semi/io/reader.py`
  (`ENGINE_SUPPORTED_MANIFEST_MAJOR = 1`). Added
  `tests/test_schema_versioning.py` (8 pure-Python tests).
  `pyproject.toml` version bumped to `0.11.0`; the wheel now
  `force-include`s `schemas/` so pip-installed kronos-semi ships the
  schema files. 238 tests pass, coverage 95.5%, ruff clean.

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
