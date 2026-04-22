# Roadmap

Full day-by-day plan for the kronos-semi evaluation submission. This
expands on the roadmap table in `PLAN.md`. Each day lists a goal,
concrete deliverables, verification criteria (the commands you run to
prove it is done), and dependencies on prior days.

Read `PLAN.md` for the short version and the current in-flight task.

## Day 1. Equilibrium Poisson, 1D pn junction, Docker env

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

## Day 2. Slotboom drift-diffusion, coupled Newton, bias sweep

- **Status:** Done (2026-04-20). Delivered on
  `dev/day2-drift-diffusion`; all deliverables landed, the
  `pn_1d_bias` verifier is green (V=0.6 V within 10% of Shockley),
  and Day 1 did not regress. The low-bias ideal-Shockley mismatch is
  physical (SRH depletion-region recombination raises the ideality
  factor) and is handled by the Day 2 verifier with qualitative
  checks below 0.5 V. Day 3 will replace those qualitative checks
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
  - `docker compose run --rm benchmark pn_1d` still green (Day 1
    regression).
  - `docker compose run --rm benchmark pn_1d_bias` exits 0 with IV
    curve within 10% of Shockley over V in [0.2, 0.6] V. Below 0.2 V
    the current is dominated by numerical noise and Shockley is not a
    meaningful target.
  - `docker compose run --rm test` green with the new tests.
  - ruff clean.
- **Dependencies:** Day 1 merged to `main`.

## Day 3. Bias ramping hardening, Shockley IV polish

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
- **Dependencies:** Day 2.

## Day 4. Verification & Validation suite

- **Status:** Done (2026-04-21). PR `dev/day4-vnv`.
- **Goal:** replace single-point physical sanity tests with proper
  verification (Roache / Oberkampf-Roy sense): manufactured solutions,
  mesh-refinement convergence studies, and discrete conservation
  checks. The originally-planned refactor moves to Day 5.
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
- **Dependencies:** Day 3 (PR #5) merged into `main`.

## Day 5. Refactor, expanded test coverage, docs pass

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
    with byte-identical numerics to the Day 4 baseline.
  - `python scripts/run_verification.py all` reports 53 PASS / 0 FAIL
    with finest-pair MMS rates byte-identical to Day 4.
- **Dependencies:** Day 4 V&V suite merged (PR #6).

## Day 6. 2D MOS capacitor

- **Status:** Done (2026-04-21). Delivered on `dev/day6-mos-2d`.
  All eight deliverables landed; the `mos_2d` benchmark exits 0
  with 4/4 verifier checks green; Day 4-5 V&V suite stayed green
  with the new multi-region Poisson MMS study clearing
  rate_L^2 >= 1.99; 1D benchmarks (`pn_1d`, `pn_1d_bias`,
  `pn_1d_bias_reverse`) are byte-identical to Day 5.
- **Goal (as delivered):** first 2D benchmark and first multi-region
  device. Equilibrium Poisson assembles over the full mesh with
  cellwise DG0 eps_r; space charge is restricted to silicon via
  `dx(subdomain_id=semi_tag)`; gate contact applies a Dirichlet BC on
  psi at the oxide top face; a C-V verifier differentiates gate charge
  with respect to V_gate and matches depletion-approximation MOS
  theory within 10% in the verifier window.
- **Deliverables (as landed):**
  - `docs/mos_derivation.md`: derivation-first gate, approved Day 6
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
- **Dependencies:** Day 5 refactor (PR #7).

## Day 7. 3D doped resistor

- **Status:** Completed (2026-04-21). Merged via PR #9 from
  `dev/day7-resistor-3d` at `a604b12`.
  All deliverables landed; the `resistor_3d` benchmark exits 0 on
  both the builtin `create_box` mesh and the committed gmsh fixture
  with worst |R_sim - R_theory|/R_theory well inside the fixed 1%
  tolerance on both paths; Day 1-6 benchmarks (`pn_1d`,
  `pn_1d_bias`, `pn_1d_bias_reverse`, `mos_2d`) stay byte-identical;
  the Day 4-6 V&V suite stays green (every MMS rate within 0.01 of
  the post-Day-6 values, including `2d_multiregion`).
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
- **Dependencies:** Day 6 MOS merged (PR #8).

## Day 8. Final polish and submission packaging

- **Status:** In flight. In progress on `dev/day8-polish`, targeting
  PR #10 against `main`. Was originally Day 7 in the pre-resistor
  roadmap.
- **Goal:** make the submission presentation-ready. No new physics,
  no new benchmarks; this is a documentation, notebook, and release
  pass only.
- **Deliverables:**
  - Sync `PLAN.md` "Current state" and this roadmap to reflect the
    Day 7 merge and Day 8 in flight.
  - Rewrite the `README.md` status section as an end-of-Day-7
    capability matrix (Days 1-7 shipped across PRs #2-#9), with a
    short scope-out prose block enumerating COMSOL Semiconductor-Module
    features that are deliberately out of scope.
  - Regenerate `notebooks/01_pn_junction_1d.ipynb` via
    `scripts/build_notebook_01.py` with current-state framing (no
    "Day 1" heading).
  - Author `notebooks/02_pn_junction_bias.ipynb` (Day 2-3 content:
    forward Shockley + reverse SNS sweeps).
  - Author `notebooks/03_mos_cv.ipynb` (Day 6 content: MOS C-V, with
    the `V_FB + 0.2` verifier-window disclosure surfaced in the
    narrative).
  - Author `notebooks/04_resistor_3d.ipynb` (Day 7 content: bipolar
    sweep, builtin-vs-gmsh comparison).
  - Add a README "Notebooks" catalog with Colab badges.
  - Append `[0.8.0] - Day 8` to `CHANGELOG.md`.
  - Open PR #10, wait for CI, do not self-merge.
  - Post-merge and only on explicit human prompt, create annotated
    tag `v0.2.0`.
- **Verification:**
  - Every notebook executed on a real Colab runtime (not just
    `jupyter nbconvert --execute` locally) with wall time and
    final-cell observation recorded in the PR body.
  - README and CHANGELOG accurately describe the end-of-Day-7
    capability set.
  - CI green on every commit on the branch
    (`gh run list --branch dev/day8-polish` verbatim).
  - All five benchmarks (`pn_1d`, `pn_1d_bias`, `pn_1d_bias_reverse`,
    `mos_2d`, `resistor_3d`) stay green; coverage stays >= 95%;
    V&V suite stays green; no physics, schema, or verifier changes.
- **Dependencies:** Days 6 and 7.

## Post-submission (stretch goals)

Explicitly out of scope for the evaluation submission. Listed here so
contributors and reviewers know the intended direction.

- **Field-dependent mobility:** Caughey-Thomas parameterization with
  saturation velocity, needed for any quantitatively accurate
  high-field device simulation.
- **Auger and radiative recombination:** becomes important in heavily
  doped regions and direct-gap materials.
- **Fermi-Dirac statistics:** needed for degeneracy in source/drain
  regions and for wide-gap materials at low temperature.
- **Incomplete dopant ionization:** freeze-out effects at low T.
- **AC small-signal:** linearized frequency response around a DC
  operating point; enables impedance and C-V analysis.
- **Transient solver:** time-dependent drift-diffusion for
  switching analysis.
- **Band-to-band and trap-assisted tunneling:** breakdown and gate
  leakage.
- **Heterojunctions:** multi-semiconductor stacks (GaAs/AlGaAs,
  Si/SiGe).
- **Full MOSFET:** source/drain/gate/body with submesh for carriers
  and proper work-function handling.
- **FinFET and 3D transistor geometries:** geometry generation and
  meshing, plus whatever new physics the above stretch goals bring.
- **GUI or web frontend:** a JavaScript SPA that builds the JSON
  input interactively and streams results back.
