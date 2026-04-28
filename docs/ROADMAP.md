# Roadmap

kronos-semi is a FEniCSx-based finite-element semiconductor device simulator that mimics the capabilities of the COMSOL Semiconductor Module. Simulations are driven by a single JSON file (optionally referencing an external geometry/mesh artifact) and are otherwise plain text. The project ships 1D, 2D, and 3D benchmark problems and a zero-setup Colab notebook (cloud-hosted, no local install) so reviewers can run everything from a browser.

## Capability matrix

Milestones M1 through M14.1 are shipped (v0.14.1, 2026-04-27). The
table below reflects the state on `main` after the M13.1 close-out
and is regenerated whenever a milestone merges. PLAN.md tracks the
in-flight task; this table is the "what does this project do"
snapshot.

| Capability | Dimensions | Status | Verifier |
|---|---|---|---|
| Equilibrium Poisson | 1D / 2D / 3D | shipped (M1, M6, M7) | MMS finest-pair rates L2 = 2.0, H1 = 1.0 |
| Coupled drift-diffusion (Slotboom) | 1D / 2D | shipped (M2) | MMS L2 rates >= 1.99 across variants |
| SRH recombination | 1D / 2D | shipped (M2, M3) | SNS analytical at reverse bias |
| Ohmic contact BCs | 1D / 2D / 3D | shipped (M2, M7) | Shockley diode within 10% at forward bias |
| Gate contact BCs with phi_ms | 2D | shipped (M6) | MOS C-V within 10% in depletion window |
| Multi-region Poisson (Si/SiO2) | 2D | shipped (M6) | multi-region MMS L2 >= 1.99 |
| File-sourced gmsh meshes | 3D | shipped (M7) | builtin vs gmsh V-I match within 1% |
| Adaptive bias continuation | uni + bipolar | shipped (M3, M7) | pn forward + reverse, 3D resistor |
| 3D ohmic V-I linearity | 3D | shipped (M7) | V-I linearity within 1% |
| Result artifact writer + manifest | n/a | shipped (M9) | manifest schema validation, semi-run CLI round-trip |
| HTTP server (POST /solve, GET /runs/{id}) | n/a | shipped (M10) | 15-test FastAPI suite incl. WebSocket progress |
| Schema versioning + UI-facing companion | n/a | shipped (M11) | Draft-07 validation, major-version gate |
| MOSFET 2D (n+ implants, SNES amendment) | 2D | shipped (M12) | mosfet_2d verifier within +-20%; ADR 0008 |
| Transient solver (BDF1/BDF2, Slotboom) | 1D | shipped (M13, M13.1) | steady-state-limit < 1e-4; BDF rate tests; pn_1d_turnon within 5% |
| AC small-signal analysis | 1D / 2D | shipped (M14) | rc_ac_sweep C(f) within 0.41% over [1 Hz, 1 MHz] |
| Differential capacitance via AC admittance | 2D | shipped (M14.1) | mos_cap_ac vs mos_cv on Q_gate (PR #38) |
| Test suite | pure + FEM | shipped | ~280 test functions across pure-Python and FEM; coverage gate >= 95% |
| V&V | 10 studies | shipped | 62/62 PASS |
| CI | pure-python + lint + docker-fem (parallelized) | shipped (M5, PR #56) | green on `main` |

The "shipped (Mx)" column lists the milestone(s) in which a
capability landed. ADR cross-references appear in the per-milestone
delivery history below. The exact test count and coverage figure on
`main` after a given milestone merge are recorded in the
corresponding `CHANGELOG.md` entry.

## Scope vs. COMSOL Semiconductor Module

kronos-semi covers the quasi-static, steady-state, and small-signal
subset of the COMSOL Semiconductor Module, plus a 1D transient
turn-on solver:

**In scope (shipped):**
- Poisson equation with multi-region dielectric (Si/SiO2)
- Drift-diffusion in Slotboom (quasi-Fermi potential) form
- SRH recombination with configurable trap energy
- Ohmic contacts and ideal gate contacts (with phi_ms)
- 1D interval, 2D rectangle, and 3D box meshes (builtin)
- 3D unstructured tetrahedral meshes via gmsh .msh files
- Bias sweeps with adaptive step-size continuation (uni- and bipolar)
- Method-of-Manufactured-Solutions and conservation V&V suite
- Result artifact writer with manifest schema (M9), HTTP server
  (M10), versioned input/manifest schemas (M11)
- 2D MOSFET with Gaussian source/drain implants (M12)
- Transient (time-dependent) solver, BDF1 and BDF2, Slotboom
  primary unknowns with BC-ramp continuation (M13, M13.1)
- AC small-signal analysis with frequency sweep (linearised
  `(J + j*omega*M) du = -dF/dV dV`, real 2x2 block reformulation)
  (M14)
- Differential capacitance via AC admittance on the MOS capacitor
  (M14.1)

**Explicitly out of scope (post-submission stretch goals, sequenced
as M15-M17 in `docs/IMPROVEMENT_GUIDE.md`):**
- Caughey-Thomas / Lombardi field-dependent mobility (M16.1, M16.2)
- Auger and radiative recombination (M16.3)
- Fermi-Dirac statistics; Boltzmann is used throughout (M16.4)
- Schottky contacts (M16.5)
- Band-to-band or trap-assisted tunneling (M16.6)
- Heterojunctions and position-dependent band structure (M17)
- GPU linear solver (CPU-LU only today; unusable above ~200k DOFs).
  M15.
- Full FinFET or 3D transistor geometries (depends on M15 GPU
  scaling)
- Incomplete dopant ionization (freeze-out)
- GUI or web frontend (kronos-semi engine ships M9-M11 contracts;
  the consumer UI is a separate repo, M18)

## JSON input contract

Every simulation is driven by a JSON file validated against the jsonschema Draft-07 schema in `semi/schema.py`. The JSON specifies: mesh source (builtin geometry or path to an external .msh file), regions with materials and tags, doping profiles, contacts, physics options, and solver options. A JSON file plus the package version fully determines the simulation. See `docs/adr/0001-json-as-input-format.md` for the rationale.

## Colab onboarding

Each benchmark notebook opens from a Colab badge. The first cell installs FEniCSx via FEM-on-Colab (~30 s), clones the repo, and loads a JSON benchmark file. No local installation is required. See the README for Colab links.

## Delivery history

The milestone entries below record the goal, deliverables,
verification, and dependencies for each shipped milestone. ADRs and
PR numbers are cross-referenced; ADR contents are not duplicated
here. The capability matrix above is the concise entry point for
reviewers; this section is the historical record of how each
capability was built.

Read `PLAN.md` for the current in-flight task and the next milestone.

## M1: Equilibrium Poisson, 1D pn junction, Docker env

- **Status:** Done (2026-04-20). PR `dev/docker-day1-fix`.
- **Goal:** prove the end-to-end pipeline works for a known analytical
  test case.
- **Deliverables:**
  - `Dockerfile`, `docker-compose.yml`, `.devcontainer/devcontainer.json`,
    `.dockerignore`.
  - `semi/constants.py`, `semi/materials.py`, `semi/scaling.py`,
    `semi/doping.py`, `semi/schema.py`.
  - `semi/mesh.py`, `semi/physics/poisson.py`, `semi/solver.py`,
    `semi/run.py`.
  - `benchmarks/pn_1d/pn_junction.json`.
  - `scripts/run_benchmark.py` with a `pn_1d` verifier.
  - 36 pytest tests covering the pure-Python core.
- **Verification:**
  - `docker compose run --rm test` exits 0 with 36/36 pass.
  - `docker compose run --rm benchmark pn_1d` exits 0 with all six
    verifier checks green.
- **Dependencies:** none.

## M2: Coupled drift-diffusion (Slotboom, coupled Newton, bias sweep)

- **Status:** Done (2026-04-20). Delivered on
  `dev/day2-drift-diffusion`; all deliverables landed, the
  `pn_1d_bias` verifier is green (V=0.6 V within 10% of Shockley),
  and M1 did not regress. The low-bias ideal-Shockley mismatch is
  physical (SRH depletion-region recombination raises the ideality
  factor) and is handled by the M2 verifier with qualitative
  checks below 0.5 V. M3 will replace those qualitative checks
  with a Sah-Noyce-Shockley term.
- **Goal (as delivered):** solved the full coupled (psi, Phi_n, Phi_p)
  system with SRH recombination under forward bias, verified against
  Shockley at high bias.
- **Deliverables:**
  - `semi/physics/drift_diffusion.py` with Slotboom continuity forms
    (electron and hole continuity as separate form builders, plus a
    combined block residual).
  - `semi/physics/recombination.py` with the SRH kernel as a UFL
    expression builder.
  - Extension of `semi/solver.py` to dispatch between scalar and
    blocked `NonlinearProblem` via the `kind` argument.
  - Update `semi/run.py` to build the coupled system when the config's
    `solver.type` is `"drift_diffusion"` or `"bias_sweep"`.
  - Schema extension: `recombination` block (tau_n, tau_p, E_t),
    `contacts[*].voltage_sweep` for bias ramping.
  - `benchmarks/pn_1d_bias/pn_junction_bias.json` with a forward-bias
    IV sweep from 0.0 to 0.6 V in 0.05 V increments.
  - `pn_1d_bias` verifier in `scripts/run_benchmark.py` checking IV vs
    Shockley `J = J_0 (exp(V/V_t) - 1)` with the saturation current
    $J_0$ computed from $D_n, D_p, \tau_n, \tau_p, N_A, N_D, n_i$.
  - At least 10 new unit tests: Slotboom math round-trip (n, p
    recovery), SRH kernel limits (R = 0 at equilibrium, R > 0 under
    bias, correct saturation behavior).
- **Verification:**
  - `docker compose run --rm benchmark pn_1d` still green (M1
    regression).
  - `docker compose run --rm benchmark pn_1d_bias` exits 0 with IV
    curve within 10% of Shockley over V in [0.2, 0.6] V. Below 0.2 V
    the current is dominated by numerical noise and Shockley is not a
    meaningful target.
  - `docker compose run --rm test` green with the new tests.
  - ruff clean.
- **Dependencies:** M1 merged to `main`.

## M3: Adaptive continuation (bias ramping hardening, Shockley IV polish)

- **Status:** Done (2026-04-21). PR `dev/day3-bias-hardening`.
- **Goal:** make the bias sweep robust, add reverse-bias check,
  document the ramp-continuation logic.
- **Delivered:**
  - Adaptive step-size control: `AdaptiveStepController` grows the
    step by `grow_factor` (default 1.5) after `easy_iter_threshold`
    (default 4) consecutive easy solves, halves on SNES failure, and
    clamps to the sweep endpoint. Cuts total SNES iterations on the
    forward `pn_1d_bias` sweep from 42 to 31 (26.2% reduction).
  - Sah-Noyce-Shockley recombination term (with the leading-order
    Sze f = 2 V_t/(V_bi - V) correction) in the forward verifier:
    J_sim matches J_diff + J_rec within 15% on [0.15, 0.55] V and
    Shockley diffusion within 10% at V = 0.6 V.
  - New `benchmarks/pn_1d_bias_reverse/` (0 to -2 V) with a
    verifier demanding |J| within 20% of the net SRH generation
    current (q n_i / 2 tau_eff)(W(V) - W(0)) on [-2, -0.5] V. For
    tau = 1e-8 s this replaces the originally-proposed "J saturates
    to J_0" bar, which does not apply in the SRH-generation regime.
  - Extracted analytical reference curves into
    `semi/diode_analytical.py`; 26 new pytest tests.
  - New "Bias continuation strategy" subsection in `docs/PHYSICS.md`.
  - CI `docker-fem` job runs the reverse benchmark on every push.
- **Dependencies:** M2.

## M4: V&V suite (Verification & Validation)

- **Status:** Done (2026-04-21). PR `dev/day4-vnv`.
- **Goal:** replace single-point physical sanity tests with proper
  verification (Roache / Oberkampf-Roy sense): manufactured solutions,
  mesh-refinement convergence studies, and discrete conservation
  checks. The originally-planned refactor moves to M5.
- **Rationale:** the existing benchmark "verifiers" prove the code
  reproduces analytical reference values at one mesh resolution. They
  do not prove convergence, do not exercise the discretization away
  from operating points, and do not catch sign or coefficient errors
  that happen to cancel at a single resolution. MMS is the standard
  remedy.
- **Deliverables:**
  - **Phase 1: MMS for Poisson.**
    - `semi/verification/mms_poisson.py`: 1D and 2D smooth
      manufactured solutions, UFL-based forcing, mesh sweep,
      L^2 / H^1 error and observed-rate computation, CSV+PNG output.
    - `tests/fem/test_mms_poisson.py` with `tests/fem/conftest.py`
      that skips when dolfinx is missing.
  - **Phase 2: Mesh convergence on `pn_1d`.**
    - `semi/verification/mesh_convergence.py`: sweep N in
      [50, 100, 200, 400, 800, 1600], record V_bi error, peak |E|
      error, depletion-width error, Newton iterations, solve time.
    - `tests/fem/test_mesh_convergence.py` (3 levels for speed).
    - Documented convergence order, including any plateau caused by
      the depletion-approximation reference itself being
      approximate.
  - **Phase 3: Conservation checks.**
    - Current continuity check (max |J_total(x) - mean| / mean)
      added to the `pn_1d_bias` and `pn_1d_bias_reverse` verifiers.
      Threshold +-5% forward, +-15% reverse.
    - Charge conservation check (integrated rho over the device
      bounded by solver tolerance) added to `pn_1d`.
    - Unit tests on the check functions in isolation, with mock
      `SimulationResult` objects.
  - **Phase 4: MMS for coupled drift-diffusion.**
    - `semi/verification/mms_dd.py`: manufactured psi, Phi_n, Phi_p
      with three forcing terms in UFL.
    - Three pytest tests: psi-only with frozen quasi-Fermis;
      full three-block coupling; full coupling with nonzero SRH.
    - Acceptance: L^2 rate >= 1.75 on each field on the finest pair.
  - **Phase 5: CI integration.**
    - `scripts/run_verification.py` CLI parallel to
      `scripts/run_benchmark.py` (subcommands `mms_poisson`,
      `mms_dd`, `mesh_convergence`, `conservation`, `all`).
    - `.github/workflows/ci.yml`: run V&V step inside `docker-fem`,
      upload `results/mms_*` and `results/mesh_convergence/`
      artifacts. Total CI runtime under 15 minutes.
  - **Final phase: Documentation.**
    - New "Verification & Validation" section in `docs/PHYSICS.md`.
    - New ADR `docs/adr/0006-verification-and-validation-strategy.md`.
    - PLAN.md, ROADMAP.md, CHANGELOG.md final update.
- **Delivered:**
  - MMS-Poisson 1D linear / 1D nonlinear / 2D triangles with
    finest-pair L^2 at theoretical 2.000 and H^1 at 0.999-1.000;
    `2d_quad_smoke` sanity ratio 0.394.
  - Mesh convergence on `pn_1d` over N in [50..1600] with
    self-convergence Cauchy ratios at >= 1.99x per doubling for
    E_peak and W over the first four levels (before the
    physics-model plateau from the depletion-approximation
    reference, which is honest-flagged in both the runner output
    and in `docs/PHYSICS.md`). V_bi is set by the Ohmic BCs and
    so is reported but not gated.
  - Conservation: charge on `pn_1d` equilibrium at rel 1.5e-17
    (threshold 1e-10); current continuity on `pn_1d_bias` V in
    {0.30, 0.45, 0.60} V and `pn_1d_bias_reverse` V in {-0.50,
    -1.00, -2.00} V, worst max_rel 1.91% forward / 0.020%
    reverse inside the 5% / 15% tolerances.
  - MMS-DD nine studies (three variants x three grids) with
    every gated block rate >= 1.99 L^2 and >= 0.996 H^1.
    Variant A gates only the psi block because the continuity
    rows collapse to machine roundoff when the quasi-Fermis are
    identically zero (documented in `mms_dd_derivation.md`
    Amendment section).
  - `scripts/run_verification.py` with subcommands `mms_poisson`,
    `mesh_convergence`, `conservation`, `mms_dd`, `all`.
  - CI `docker-fem` job now runs V&V inside the dolfinx image,
    timeout tightened to 15 minutes, artifacts uploaded as
    `fem-results` (renamed from `benchmark-plots`).
  - New ADR `0006-verification-and-validation-strategy.md`, new
    "Verification & Validation" section in `docs/PHYSICS.md`,
    derivation artifact `docs/mms_dd_derivation.md` amended with
    the SNES atol-0 fix rationale.
- **Dependencies:** M3 (PR #5) merged into `main`.

## M5: Refactor and test pass

- **Status:** Done (2026-04-21). Delivered on `dev/day5-refactor`.
- **Goal:** paid down technical debt accumulated across Days 1-3 and
  expanded coverage before moving to higher dimensions.
- **Deliverables (as delivered):**
  - `semi/bcs.py` extracted (pure-Python core tier, no dolfinx import
    at module scope) with `ContactBC` dataclass,
    `resolve_contacts(cfg, facet_tags=None, voltages=None)`,
    `build_psi_dirichlet_bcs(...)`, and `build_dd_dirichlet_bcs(...)`.
    Inline `_build_ohmic_bcs_psi` and `_build_dd_ohmic_bcs` removed
    from `semi/run.py`.
  - `semi/run.py` split into a thin dispatcher (74 lines) plus
    `semi/runners/` (`equilibrium.py`, `bias_sweep.py`, `_common.py`)
    and `semi/postprocess.py`. `bias_sweep.py` 294 lines; no module
    in `semi/` exceeds 300 lines (was 580 in `run.py`).
  - Coverage 96.25% (1598 statements, 60 missed) across 177 tests
    (139 prior + 38 new). CI gate `--cov-fail-under=95` enforced
    in the `docker-fem` job.
  - `docs/PHYSICS.md` Section 2.5 completed: full scaled
    drift-diffusion derivation (continuity rows pick up `L_0^2`
    explicitly because the mesh stays in meters per Invariant 3),
    scaled SRH kernel, residual-sign block summary, J_0 numerical
    check.
  - `docs/adr/0007-contact-bc-interface.md` records the dataclass-
    over-dict / voltages-override / insulating-skip / optional-facet-
    verification design choices.
- **Verification:**
  - `pytest --cov=semi --cov-fail-under=95` exits 0 at 96.25%.
  - `pn_1d`, `pn_1d_bias`, `pn_1d_bias_reverse` benchmarks all pass
    with byte-identical numerics to the M4 baseline.
  - `python scripts/run_verification.py all` reports 53 PASS / 0 FAIL
    with finest-pair MMS rates byte-identical to M4.
- **Dependencies:** M4 V&V suite merged (PR #6).

## M6: 2D MOS capacitor

- **Status:** Done (2026-04-21). Delivered on `dev/day6-mos-2d`.
  All eight deliverables landed; the `mos_2d` benchmark exits 0
  with 4/4 verifier checks green; M4-M5 V&V suite stayed green
  with the new multi-region Poisson MMS study clearing
  rate_L^2 >= 1.99; 1D benchmarks (`pn_1d`, `pn_1d_bias`,
  `pn_1d_bias_reverse`) are byte-identical to the M5 baseline.
- **Goal (as delivered):** first 2D benchmark and first multi-region
  device. Equilibrium Poisson assembles over the full mesh with
  cellwise DG0 eps_r; space charge is restricted to silicon via
  `dx(subdomain_id=semi_tag)`; gate contact applies a Dirichlet BC on
  psi at the oxide top face; a C-V verifier differentiates gate charge
  with respect to V_gate and matches depletion-approximation MOS
  theory within 10% in the verifier window.
- **Deliverables (as landed):**
  - `docs/mos_derivation.md`: derivation-first gate, approved M6
    prerequisite (9da5fa6).
  - `semi/mesh.py` submesh + cellwise eps_r helpers (719b9e8).
  - `semi/bcs.py` gate contact wiring (b24c650).
  - `semi/physics/poisson.py:build_equilibrium_poisson_form_mr`
    (this commit).
  - `semi/physics/drift_diffusion.py` submesh block residual via
    `entity_maps` (455196a).
  - `benchmarks/mos_2d/mos_cap.json` device spec: 500 nm p-type Si
    (N_A = 1e17 cm^-3) / 5 nm SiO2, uniform 1 nm vertical mesh,
    V_gate sweep [-0.9, +1.2] V.
  - `scripts/run_benchmark.py`: 2D `plot_mos_2d` (tricontourf of
    psi, central-column psi(y), C-V and Q-V) and `verify_mos_2d`
    (10% tolerance in [V_FB + 0.2, V_T - 0.1] V).
  - `semi/verification/mms_poisson.py`:
    `run_mms_poisson_2d_multiregion` with Si/SiO2 coefficient jump
    (99fcd2b).
  - `docs/PHYSICS.md` Section 6 (MOS reference).
- **Verification (observed):**
  - `benchmark mos_2d` exits 0; worst |C_sim - C_theory|/C_theory in
    the verifier window is 9.25% at V_gate = -0.20 V.
  - 4 plots written: `psi_2d.png`, `potentials_1d.png`, `cv.png`,
    `qv.png`.
  - V&V suite rates: see `PHYSICS.md` sections 5.1-5.5.
  - `pytest --cov=semi --cov-fail-under=95` exits 0 at 95.43%
    (195 tests pass).
- **Dependencies:** M5 refactor (PR #7).

## M7: 3D doped resistor

- **Status:** Completed (2026-04-21). Merged via PR #9 from
  `dev/day7-resistor-3d` at `a604b12`.
  All deliverables landed; the `resistor_3d` benchmark exits 0 on
  both the builtin `create_box` mesh and the committed gmsh fixture
  with worst |R_sim - R_theory|/R_theory well inside the fixed 1%
  tolerance on both paths; M1-M6 benchmarks (`pn_1d`,
  `pn_1d_bias`, `pn_1d_bias_reverse`, `mos_2d`) stay byte-identical;
  the M4-M6 V&V suite stays green (every MMS rate within 0.01 of
  the post-M6 values, including `2d_multiregion`).
- **Goal (as delivered):** confirmed the simulator framework extends
  to 3D with no physics changes. Equilibrium Poisson and the
  Slotboom drift-diffusion forms (`semi/physics/poisson.py`,
  `semi/physics/drift_diffusion.py`) are dimension agnostic and
  were left untouched. The only new machinery is the gmsh `.msh`
  loader, a 3D benchmark device, a bipolar-sweep driver path, a
  V-I linearity verifier, and 3D slice plots.
- **Delivered:**
  - `docs/resistor_derivation.md` (derivation-lite gate; 466a820):
    device geometry, analytical ohmic resistance, 1% V-I linearity
    rationale, 3D slice-plot strategy, gmsh loader test strategy.
  - `semi/mesh.py::_build_from_file` wired via
    `dolfinx.io.gmsh.read_from_msh` (e395830, docs alignment in
    3abd170). Physical groups stored in the `.msh` file are
    returned verbatim as `cell_tags` and `facet_tags`; the JSON
    box-tagger is bypassed for file-source meshes.
  - `semi/runners/bias_sweep.py` two-leg walk for bipolar sweeps
    (7c947b2): when `v_sweep_list` spans zero, the ramp runs from
    `V = 0 -> min(V) -> max(V)` with a fresh
    `AdaptiveStepController` on each leg; unipolar sweeps behave
    exactly as before.
  - `benchmarks/resistor_3d/resistor.json`: 3D rectangular bar
    (1 um x 200 nm x 200 nm), uniform n-type N_D = 1e18 cm^-3,
    ohmic contacts on x=0 and x=L, 5-point sweep in
    [-0.01, +0.01] V (62f52d2).
  - `benchmarks/resistor_3d/resistor_gmsh.json` plus a committed
    fixture `fixtures/box.geo`/`box.msh` (reproducible via
    `gmsh -3 box.geo -o box.msh`). Same device on an unstructured
    tetrahedral mesh.
  - `scripts/run_benchmark.py`: `verify_resistor_3d` V-I linearity
    verifier (1% tol, zero-bias-noise and sign-consistency sanity
    checks); 3D slice plotters for psi and |J_n| at the y = W/2
    midplane plus an I-V scatter.
  - Tests: `tests/fem/test_mesh_gmsh.py` (round-trip + physical
    group preservation) and `tests/fem/test_resistor_3d.py`
    (coarse smoke plus builtin/gmsh production-JSON theory match).
  - `docs/PHYSICS.md` Section 7 (3D extension notes, ~1 page)
    citing `docs/resistor_derivation.md` for the full derivation.
- **Verification:**
  - `docker compose run --rm benchmark resistor_3d` exits 0 with
    V-I linearity within 1% on both builtin and gmsh variants;
    three plots written (`psi_slice_y_midplane.png`,
    `jn_slice_y_midplane.png`, IV scatter).
  - All prior benchmarks byte-identical; all prior V&V gates
    green; coverage stays >= 95%.
- **Dependencies:** M6 MOS merged (PR #8).

## M8: Submission polish

- **Status:** Done (2026-04-22). Merged via PRs #10/#11 from
  `dev/submission-polish`. Originally tagged v0.8.0.
- **Goal:** make the submission presentation-ready. No new physics,
  no new benchmarks; this was a documentation, notebook, and release
  pass only.
- **Deliverables:**
  - Sync `PLAN.md` "Current state" and this roadmap to reflect the
    M7 merge and M8 in flight.
  Rewrite the `README.md` status section as an end-of-M7
    capability matrix (M1-M7 shipped across PRs #2-#9), with a
    short scope-out prose block enumerating COMSOL Semiconductor-Module
    features that are deliberately out of scope.
  - Regenerate `notebooks/01_pn_junction_1d.ipynb` via
    `scripts/build_notebook_01.py` with current-state framing (no
    framing (no "M1" heading).
  - Author `notebooks/02_pn_junction_bias.ipynb` (M2-M3 content:
    forward Shockley + reverse SNS sweeps).
  - Author `notebooks/03_mos_cv.ipynb` (M6 content: MOS C-V, with
    the `V_FB + 0.2` verifier-window disclosure surfaced in the
    narrative).
  - Author `notebooks/04_resistor_3d.ipynb` (M7 content: bipolar
    sweep, builtin-vs-gmsh comparison).
  - Add a README "Notebooks" catalog with Colab badges.
  - Append `[0.8.0] - M8: Submission polish` to `CHANGELOG.md`.
  - Open PR #10, wait for CI, do not self-merge.
  - Post-merge and only on explicit human prompt, create annotated
    tag `v0.2.0`.
- **Verification:**
  - Every notebook executed on a real Colab runtime (not just
    `jupyter nbconvert --execute` locally) with wall time and
    final-cell observation recorded in the PR body.
  README and CHANGELOG accurately describe the end-of-M7
    capability set.
  - CI green on every commit on the branch
    (`gh run list --branch dev/submission-polish` verbatim).
  - All five benchmarks (`pn_1d`, `pn_1d_bias`, `pn_1d_bias_reverse`,
    `mos_2d`, `resistor_3d`) stay green; coverage stays >= 95%;
    V&V suite stays green; no physics, schema, or verifier changes.
- **Dependencies:** M6 and M7.

## M9: Result artifact writer

- **Status:** Done (2026-04-23). Tagged v0.9.0.
- **Goal:** define a versioned on-disk run artifact contract so the
  M10 HTTP server (and any future UI) can read engine output without
  re-running the simulation, and so artifacts are reproducible across
  engine versions.
- **Deliverables:**
  - `schemas/manifest.v1.json`: Draft-07 manifest schema covering
    engine metadata, solver summary, field list, mesh topology,
    optional bias-sweep IV entries, and warnings. `additionalProperties:
    false` throughout.
  - `semi/io/artifact.py::write_artifact(result, out_dir, run_id,
    input_json_path)` writes a versioned run directory with
    `manifest.json`, `input.json`, `mesh/mesh.xdmf`, `fields/<name>.bp`
    (or `.xdmf` fallback), `iv/<contact>.csv`, `convergence/snes.csv`.
  - `semi/io/reader.py::read_manifest(run_dir)` validates and loads
    a manifest using stdlib JSON plus optional jsonschema; no dolfinx
    or numpy dependency, safe for the M10 server subprocess.
  - `semi-run` CLI entry point (`semi/io/cli.py`) wired through
    `[project.scripts]`.
  - `tests/test_artifact.py` covering pure-Python schema tests plus
    FEM-heavy parametrized benchmark round-trips (skipif when
    dolfinx is unavailable).
- **Verification:** schema-validated round-trip on every shipped
  benchmark; `read_manifest` accepts artifacts written by `write_artifact`
  and rejects invalid manifests.
- **Dependencies:** M8 polish merged.

## M10: HTTP server

- **Status:** Done (2026-04-23). Tagged v0.10.0.
- **Goal:** expose the M9 engine over HTTP with a clean separation
  between the server process (no dolfinx import) and the worker
  subprocesses that do the actual FEM work.
- **Deliverables:**
  - `kronos_server/` top-level package (FastAPI app factory
    `build_app`, in-process `ProcessPoolExecutor` job manager,
    LocalFS storage backend, file-backed progress event stream
    `progress.ndjson`, pydantic request/response models, CORS).
    No dolfinx, UFL, or PETSc imports at module scope; FEM-heavy
    imports happen only inside `mp_context="spawn"` worker
    subprocesses.
  - HTTP endpoints: `POST /solve`, `GET /runs`, `GET /runs/{id}`
    (manifest + status), `GET /runs/{id}/{manifest,input,fields/{name},iv/{contact},logs}`,
    `GET /schema`, `GET /materials`, `GET /capabilities`, `GET /health`,
    `GET /ready`.
  - `WS /runs/{id}/stream` WebSocket emitting `step_done` and
    `run_done` messages.
  - `kronos-server` console entry point with env-driven configuration;
    `server` service in `docker-compose.yml`; `[server]` extra in
    `pyproject.toml`.
  - Optional `progress_callback=None` parameter added to
    `semi/runners/{equilibrium,bias_sweep,mos_cv}.py` (the only
    permitted change to `semi/` in M10).
  - `tests/test_kronos_server.py` (15 tests including end-to-end
    solve, WebSocket progress on `pn_1d_bias`, and concurrent solves).
- **Verification:** 15-test FastAPI suite green; manifest returned
  via `GET /runs/{id}` matches the on-disk manifest produced by M9.
- **Dependencies:** M9 artifact writer.

## M11: Schema versioning + UI-facing companion

- **Status:** Done (2026-04-23). Tagged v0.11.0.
- **Goal:** define a stable input-schema contract a UI form-builder
  can consume, with explicit major/minor version semantics on both
  input and manifest.
- **Deliverables:**
  - `schemas/input.v1.json`: Draft-07 input schema extracted verbatim
    from the old `semi.schema.SCHEMA` dict. Every `type: object` node
    carries a `description`; leaf properties carry descriptions and
    `default`/`examples` annotations.
  - Required top-level `schema_version` field on every input JSON.
    Engine constant `semi.schema.ENGINE_SUPPORTED_SCHEMA_MAJOR = 1`;
    major-mismatch raises `SchemaError`, minor/patch skew is accepted
    silently.
  - Symmetric `semi.io.reader.ENGINE_SUPPORTED_MANIFEST_MAJOR = 1`
    output-side gate; pre-M9 manifests without `schema_version`
    are still loadable with a warning.
  - `GET /schema` returns
    `{schema, version, supported_major}` instead of the bare schema.
  - `schema_version: "1.0.0"` as the first key of every shipped
    benchmark JSON; `pyproject.toml` `force-include`s `schemas/` so
    `pip install kronos-semi` ships the schema files.
  - `tests/test_schema_versioning.py` (8 pure-Python tests).
- **Verification:** Draft-07 validity, description presence on every
  object node, semver compliance on benchmark JSONs, major-mismatch
  rejection, minor-skew acceptance, manifest major-mismatch rejection.
- **Dependencies:** M9, M10.

## M12: MOSFET 2D + SNES tolerance amendment

- **Status:** Done (2026-04-25). Tagged v0.12.0. Closes #25.
- **Goal:** ship a first MOSFET geometry (2D, n+ source/drain
  Gaussian implants) and amend the SNES tolerances inherited from
  the M2-era 1D pn-junction settings to remain achievable on
  multi-region FEM problems with high-injection regions.
- **Deliverables:**
  - `benchmarks/mosfet_2d/mosfet_2d.json`: 2D n-channel MOSFET, p-Si
    body (5 um x 2 um, N_A = 1e16 cm^-3) with two Gaussian n+
    source/drain implants (peak 5e19 cm^-3). Gate sweeps 0 -> 1.5 V
    at V_DS = 0.05 V.
  - `semi/runners/bias_sweep.py` SNES tolerances:
    `snes_rtol` 1e-14 -> 1e-10; `snes_atol` 1e-14 -> 1e-7;
    `snes_max_it` 60 -> 100; `snes_stol` unchanged.
  - `tests/fem/test_bias_sweep_multiregion.py` 3-step ramp test on a
    minimal 1D multi-region config; asserts >= 3 IV points and
    non-decreasing electron current.
  - `docs/adr/0008-snes-tolerances.md`.
- **Verification:** mosfet_2d benchmark verifier within +-20% of the
  textbook MOS theory in the operating window; M1-M7 benchmarks
  remain byte-identical to the M11 baseline.
- **Dependencies:** M11.

## M13: Transient solver (BDF1 / BDF2)

- **Status:** Done (2026-04-26). Tagged v0.13.0.
- **Goal:** add a 1D transient (time-dependent) drift-diffusion solver
  with BDF1 (backward Euler) and BDF2 time integration and a diode
  turn-on benchmark. Steady-state runner and equilibrium runner
  unchanged.
- **Deliverables:**
  - `semi/timestepping.py` (`BDFCoefficients`; pure Python).
  - `semi/fem/__init__.py`, `semi/fem/mass.py`
    (`assemble_lumped_mass(V_n, V_p, dx)`).
  - `semi/results.py::TransientResult` dataclass.
  - `semi/runners/transient.py::run_transient` (initial design used
    `(psi, n_hat, p_hat)` primary-density form; superseded by ADR
    0014 in M13.1, see below).
  - `semi/postprocess.py::evaluate_partial_currents` returning
    `(J_n, J_p)` separately at a contact.
  - `schemas/input.v1.json` schema extension (`solver.type` enum
    extended with `"transient"`; new `t_end`, `dt`, `order`,
    `max_steps`, `output_every` fields).
  - `benchmarks/pn_1d_turnon/`: 1D pn diode transient turn-on
    benchmark; `verify.py` fits `J_anode(t) = J_ss (1 - exp(-t / tau))`
    and asserts `|tau_eff - tau_p| / tau_p < 5%`.
  - `tests/fem/test_transient_steady_state.py` (originally xfail
    pending M13.1).
  - `tests/mms/test_transient_convergence.py` (BDF1 rate >= 0.95,
    BDF2 rate >= 1.9; originally xfail pending M13.1).
  - ADRs `0009-transient-formulation.md` (superseded by 0014 in
    M13.1) and `0010-bdf-time-integration.md`.
- **Verification:** see M13.1 close-out.
- **Dependencies:** M12.

## M13.1: Transient close-out (Slotboom + MUMPS workspace + Jacobian shift)

- **Status:** Done (2026-04-27). Tagged v0.14.1. PRs
  #44, #45, #46, #47, #48, #49, #50, #51, #52, #54.
- **Goal:** close the gap between `run_transient` (deep steady state)
  and `run_bias_sweep` to within 1e-4 relative error and remove the
  xfail markers on `test_transient_steady_state_limit` and the BDF
  rate tests.
- **Deliverables:**
  - Slotboom `(psi, phi_n, phi_p)` primary unknowns in `run_transient`
    (ADR 0014, supersedes ADR 0009).
  - BC-ramp continuation in the transient runner (ADR 0013, applies
    in either formulation).
  - SG flux primitives in `semi/fem/scharfetter_gummel.py` and
    `semi/fem/sg_assembly.py` (ADR 0012; primary-unknown assumption
    updated by ADR 0014). Shipped for use by future work; not the
    active continuity-discretization path in `run_transient`.
  - Configurable `solver.jacobian_shift` (default 1e-14) applied
    after each Jacobian assembly via the SNES Jacobian callback;
    addresses the rank-deficient `phi_n` row in the deep p-bulk
    where every term carries `exp(psi - phi_n)` below floating-point
    precision.
  - MUMPS workspace bump (`mat_mumps_icntl_14=200`) on the transient
    factorization path. PR #52's pivot-threshold options
    (`pc_factor_zeropivot`, `mat_mumps_cntl_3`) were superseded by
    PR #54's direct petsc4py path because the factor Mat is created
    lazily after `NonlinearProblem.__init__` pushes the options DB
    in dolfinx 0.10.
  - Slotboom-native MMS rate tests added; the previously-xfailed
    `test_transient_steady_state_limit` and BDF rate tests now pass
    without xfail.
  - Removed `(psi, n_hat, p_hat)` primary-unknown transient path
    and the `use_sg_flux` opt-in flag.
- **Verification:**
  - `pn_1d_turnon`: `|tau_eff - tau_p| / tau_p < 5%`.
  - `test_transient_steady_state_limit`: relative error < 1e-4 at
    V = 0.3 V vs `run_bias_sweep`.
  - BDF rate tests: BDF1 >= 0.95, BDF2 >= 1.9 on `||error||_inf` for
    n and p.
- **Dependencies:** M13. References: M13.1 blocker docs
  `docs/m13.1-followup-5-blocker.md`,
  `docs/m13.1-ci-followup-blocker.md`.

## M14: AC small-signal analysis

- **Status:** Done (2026-04-26). Tagged v0.14.0.
- **Goal:** linearised frequency-response analysis around a converged
  DC operating point: `(J + j*omega*M) du = -dF/dV dV`.
- **Deliverables:**
  - `semi/runners/ac_sweep.py::run_ac_sweep`. At each frequency:
    converts the Slotboom DC solution to `(psi, n_hat, p_hat)`
    primary-density form via Boltzmann statistics; assembles the
    steady-state DD Jacobian J in primary-density form and the lumped
    mass M (carrier rows only; psi has no time derivative);
    builds and solves the **real 2x2 block reformulation** of the
    complex linear system (forced by the dolfinx-real PETSc build
    used in CI); evaluates terminal current including the displacement
    term `j*omega*epsilon*grad(d_psi).n` at the contact.
  - `semi/results.py::AcSweepResult` (frequencies, Y, Z, C, G,
    dc_bias, meta). Complex numbers serialise as
    `{"re": ..., "im": ...}` per the JSON-as-contract invariant.
  - `schemas/input.v1.json` `solver.type` enum extended with
    `"ac_sweep"`; new `solver.dc_bias` and `solver.ac` sub-objects.
    `frequencies` accepts an explicit `list`, a `logspace` spec, or a
    `linspace` spec. `SCHEMA_SUPPORTED_MINOR` bumped 1 -> 2;
    `schema_version` 1.2.0.
  - `benchmarks/rc_ac_sweep/`: 1D pn diode at V_DC = -1 V swept 1 Hz
    to 1 GHz logspace 41 points. C(f) within 5% of analytical
    depletion C over [1 Hz, 1 MHz]; current measurement matches to
    0.41% worst-case.
  - `tests/mms/test_ac_consistency.py` (Y(0) real; Re(Y) stable at
    low omega; C(f) omega-independent in the depletion regime).
  - `tests/fem/test_ac_dc_limit.py` (depletion-C agreement within 5%
    at V_DC in {-2.0, -1.0, -0.5} V).
  - `docs/adr/0011-ac-small-signal.md`.
- **Verification:** rc_ac_sweep verifier within 0.41%; three internal
  AC consistency invariants; depletion-limit agreement within 5%.
- **Dependencies:** M13 (lumped mass, primary-density form).

## M14.1: Differential capacitance via AC admittance (mos_cap_ac)

- **Status:** Done (2026-04-27). PR #38.
- **Goal:** ship an AC-admittance-based MOS capacitance runner
  (`mos_cap_ac`) and validate it against the existing `mos_cv`
  numerical-dQ/dV runner on the `mos_2d` benchmark.
- **Deliverables:**
  - `semi/runners/mos_cap_ac.py`: applies the M14 AC admittance
    machinery to the MOS capacitor C-V curve and reports
    `(V_gate, Q_gate)` from the imaginary part of the admittance at
    a low frequency.
  - PR #38 commit message records byte-identity to `mos_cv`'s
    numerical dQ/dV on `Q_gate` at matching gate voltages on the
    `mos_2d` benchmark.
- **Verification:** mos_cap_ac vs mos_cv on Q_gate (Phase 1 audit
  case 03 formalises this comparison; see `docs/PHYSICS_AUDIT.md`).
- **Dependencies:** M14. References: PR #38 commit body.

## Forward-looking

This document (`docs/ROADMAP.md`) is the **"what does this project
do"** snapshot, updated when milestones merge.

For the current in-flight task and the immediate next milestone, see
`PLAN.md` (`Current state` and `Next task` sections, plus the
append-only `Completed work log`).

For per-milestone acceptance specifications of the remaining
milestones (M15 GPU linear solver, M16 physics completeness, M17
heterojunctions, M18 UI in a separate repo) see
`docs/IMPROVEMENT_GUIDE.md`. ROADMAP.md will be extended with delivery
history entries for each as they merge.

The Phase 1 cross-runner consistency audit (`tests/audit/`,
`pytest -m audit`) ships its findings in `docs/PHYSICS_AUDIT.md`.
Phases 2 (external validation against Sze and Nicollian-Brews) and
3 (adversarial robustness) follow in subsequent PRs.
