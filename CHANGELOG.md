# Changelog

## [0.4.0] - Day 4

### Added
- Verification & Validation suite under `semi/verification/`, driven
  by `scripts/run_verification.py` with subcommands `mms_poisson`,
  `mesh_convergence`, `conservation`, `mms_dd`, and `all`:
  - `mms_poisson`: 1D linear, 1D nonlinear, 2D triangles + 2D quad
    smoke. Finest-pair rates at theoretical L^2 = 2.000,
    H^1 = 0.999-1.000; quad/triangle ratio 0.394.
  - `mesh_convergence`: `pn_1d` sweep over
    N in [50, 100, 200, 400, 800, 1600] with Cauchy self-convergence
    reported as the primary convergence metric (depletion-approximation
    errors plateau because the reference is itself approximate; this
    is honest-flagged in runner output and in `docs/PHYSICS.md`).
    Cauchy ratios >= 1.99x per doubling on E_peak and W over the
    first four levels.
  - `conservation`: charge conservation on `pn_1d` equilibrium
    (|Q_net|/Q_ref = 1.5e-17 vs. 1e-10 threshold) and current
    continuity on `pn_1d_bias` forward (5% tol) and
    `pn_1d_bias_reverse` (15% tol); worst max_rel 1.91% / 0.020%.
  - `mms_dd`: three variants (psi-only, full coupling no R, full
    coupling with SRH) x three grids (1D default, 1D nonlinear, 2D).
    Every gated block rate >= 1.99 L^2 and >= 0.996 H^1.
- Derivation artifact `docs/mms_dd_derivation.md` (gate-first:
  approved before any MMS-DD code landed). Amended with the
  implementation-time SNES tolerance fix: `atol = 0.0` with
  `stol = 1e-12` (the originally-proposed `atol = 1e-16` was below
  the 2D continuity-block initial residual and terminated Newton
  prematurely).
- New ADR `docs/adr/0006-verification-and-validation-strategy.md`
  recording the four-phase V&V choice and rejected alternatives
  (devsim code-to-code, single-MMS-for-everything).
- New "Verification & Validation" section in `docs/PHYSICS.md`
  (Section 5) with current finest-pair rates and honest flags.
- CI `docker-fem` job now runs
  `python scripts/run_verification.py all` inside the dolfinx
  container; timeout tightened from 30 to 15 minutes.

### Changed
- `benchmark-plots` CI artifact renamed to `fem-results`; path
  stays `results/` and now covers benchmarks + V&V together.
- Day-4 scope re-purposed: the originally-planned refactor pass is
  pushed to Day 5 to make room for the V&V work (see
  `docs/adr/0006-*.md` context).

## [0.3.0] - Day 3

### Added
- Adaptive step-size controller for bias continuation
  (`semi.continuation.AdaptiveStepController`): grows the step by
  `grow_factor` after `easy_iter_threshold` consecutive easy SNES
  solves, halves on failure, clamps to the sweep endpoint.
- Pure-Python diode analytical helpers (`semi.diode_analytical`):
  Shockley long-diode saturation, depletion width, Sah-Noyce-Shockley
  forward total reference, SRH generation reference.
- Reverse-bias benchmark `benchmarks/pn_1d_bias_reverse/` with a
  verifier that compares |J| against the net SRH generation current
  `(q n_i / 2 tau_eff)(W(V) - W(0))` within 20% on [-2, -0.5] V.
- CI `docker-fem` job now runs the reverse benchmark on every push.
- Schema fields `solver.continuation.{max_step, easy_iter_threshold,
  grow_factor}`.
- New "Bias continuation strategy" subsection in `docs/PHYSICS.md`.
- 26 new pytest tests (continuation controller, diode analytical
  helpers, schema fields); total test count 96.

### Changed
- `pn_1d_bias` forward verifier now demands J_sim within 15% of the
  SNS reference J_diff + J_rec (with Sze f = 2 V_t/(V_bi - V)
  correction) over V in [0.15, 0.55] V, and within 10% of Shockley
  diffusion at V = 0.6 V. Day 2 "qualitative below 0.5 V" caveat is
  removed.
- `run_bias_sweep` walks the ramp adaptively instead of recording at
  every `voltage_sweep.step` point. Total SNES iterations on the
  Day 2 sweep (0 to 0.6 V) dropped from 42 to 31 (26.2% reduction).

## [0.2.0] - Day 2

### Added
- Coupled Slotboom drift-diffusion with SRH recombination
  (`semi.physics.slotboom`, `recombination`, `drift_diffusion`).
- Blocked NonlinearProblem wrapper `solve_nonlinear_block` and
  `run_bias_sweep` with halving-based voltage continuation.
- `pn_1d_bias` benchmark and Shockley long-diode verifier.
- CI hardening: `docker-fem` job runs pytest plus benchmark verifiers
  inside the dolfinx container; `branches` glob fixed to match
  `dev/**`, `ci/**`, `docs/**`.
- Schema fields `recombination.{tau_n, tau_p, E_t}`, per-contact
  `voltage_sweep`, `solver.type in {drift_diffusion, bias_sweep}`,
  and `continuation.{min_step, max_halvings}`.

## [0.1.0] - Day 1

### Added
- JSON schema with Draft-07 validation (`semi.schema`)
- Material database: Si, Ge, GaAs, SiO₂, HfO₂, Si₃N₄ (`semi.materials`)
- Physical constants module with SI/cm conversions (`semi.constants`)
- Nondimensional scaling class with auto-inference from config (`semi.scaling`)
- Doping profile evaluators: uniform, step, Gaussian (`semi.doping`)
- Builtin mesh generation with region/facet tagging (`semi.mesh`)
- Equilibrium Poisson form under Boltzmann statistics (`semi.physics.poisson`)
- PETSc SNES solver wrapper for dolfinx 0.10+ (`semi.solver`)
- Top-level `run(cfg)` entry point (`semi.run`)
- 1D pn junction benchmark with JSON input (`benchmarks/pn_1d/`)
- Colab notebook: clone repo, import package, run, plot, verify (`notebooks/01_pn_junction_1d.ipynb`)
- Test suite: 36 tests covering materials, schema, scaling, doping, constants
- Offline math sanity check script (runs without dolfinx)
- GitHub Actions CI for pure-Python tests across Python 3.10/3.11/3.12
- Ruff lint configuration

### Verified
- V_bi = 0.833 V for symmetric 10¹⁷/10¹⁷ Si junction (matches analytical)
- Peak |E| = 113.5 kV/cm (within few % of depletion approximation)
- Mass-action np = n_i² in bulk (to numerical precision)
- Charge neutrality in quasi-neutral regions (to 10⁻¹⁵)

### Not yet implemented (see Roadmap in README)
- Drift-diffusion under applied bias
- Coupled (ψ, Φₙ, Φₚ) block Newton
- Bias sweep continuation
- Multi-region (oxide + semiconductor)
- 2D MOS capacitor benchmark
- 3D resistor benchmark
- MOSFET and FinFET
