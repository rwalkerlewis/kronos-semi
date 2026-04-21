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

- **Status:** In flight on `dev/day4-vnv`.
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
- **Verification:**
  - `pn_1d`, `pn_1d_bias`, `pn_1d_bias_reverse` still green
    (no regression).
  - MMS-Poisson 1D and 2D: L^2 rate >= 1.85 on the finest pair.
  - MMS-DD full coupling and with-SRH: L^2 rate >= 1.75 on each
    of psi, Phi_n, Phi_p on the finest pair.
  - Mesh convergence on `pn_1d`: monotone error reduction; at least
    2x error reduction per mesh doubling on V_bi and peak |E|.
  - Current continuity: <=5% non-uniformity at V in [0.3, 0.6] V.
  - Charge conservation: |Q_net| <= 10^-10 of scaled charge.
  - CI green on `dev/day4-vnv` with V&V step running.
- **Dependencies:** Day 3 (PR #5) merged into `main`.

## Day 5. Refactor, expanded test coverage, docs pass

- **Status:** Planned. Was originally Day 4; pushed back one slot to
  make room for the V&V suite.
- **Goal:** pay down technical debt accumulated across Days 1-3 and
  expand coverage before moving to higher dimensions.
- **Deliverables:**
  - Break `semi/run.py` into `run_equilibrium` and `run_bias_sweep`
    functions with a thin `run(cfg)` dispatcher.
  - Factor BC construction into `semi/bcs.py` so ohmic, gate, and
    future Schottky contacts share a common interface.
  - Raise pure-Python coverage to 95%+; add integration tests for
    `run` itself that can be parameterized over small in-memory
    configs.
  - Update `docs/PHYSICS.md` with the scaled drift-diffusion
    derivation now that it is implemented.
  - Add `docs/adr/0007-bc-construction-interface.md` if the BC
    refactor introduces a new decision (note: ADR 0006 is reserved
    for the Day 4 V&V strategy).
- **Verification:**
  - `pytest --cov=semi --cov-report=term-missing` shows core coverage
    at or above 95%.
  - No behavioral regressions in `pn_1d` or `pn_1d_bias`.
- **Dependencies:** Day 4 V&V suite merged.

## Day 6. 2D MOS capacitor

- **Status:** Planned. Was originally Day 5.
- **Goal:** extend to 2D and multi-region (oxide plus silicon).
- **Deliverables:**
  - `benchmarks/mos_2d/mos_cap.json` with a gate, oxide, silicon
    substrate, and body contact.
  - Submesh handling in `semi/mesh.py` so carriers live only on the
    semiconductor region.
  - Gate boundary condition (Dirichlet on psi, no continuity
    equations in oxide).
  - C-V curve verifier: depletion C-V should match MOS theory within
    10% across the depletion-to-inversion transition.
- **Verification:**
  - `benchmark mos_2d` exits 0.
  - 2D plots rendered (psi contour, |E| contour, n and p contours).
- **Dependencies:** Day 5 refactor.

## Day 7. 3D doped resistor

- **Status:** Planned. Was originally Day 6.
- **Goal:** confirm the framework extends to 3D unstructured meshes
  with no physics changes.
- **Deliverables:**
  - `benchmarks/resistor_3d/resistor.json` with a doped bar.
  - Gmsh `.msh` loading via `dolfinx.io.gmshio` (currently stubbed
    out in `semi/mesh.py`).
  - Ohmic V-I linearity verifier: current scales linearly with
    applied voltage to within 1% over a small voltage range.
- **Verification:**
  - `benchmark resistor_3d` exits 0.
  - At least two 3D slice plots (current density magnitude, potential).
- **Dependencies:** Day 5.

## Day 8. Final polish and submission packaging

- **Status:** Planned. Was originally Day 7.
- **Goal:** make the submission presentation-ready.
- **Deliverables:**
  - Regenerate `notebooks/01_pn_junction_1d.ipynb` using
    `scripts/build_notebook_01.py` after all Docker tests pass.
  - Add `notebooks/02_pn_junction_bias.ipynb` for the Day 2-3 content.
  - Update `README.md` status section with the final capability
    matrix.
  - Tag release `v0.2.0` (or similar) on `main` after final review.
  - Update `CHANGELOG.md`.
- **Verification:**
  - Both notebooks execute top-to-bottom on Colab with no errors.
  - README and CHANGELOG accurately describe final state.
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
