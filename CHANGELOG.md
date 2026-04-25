# Changelog

## [0.12.0] - M12: MOSFET n+ doping + SNES tolerance amendment

### Added
- `benchmarks/mosfet_2d/mosfet_2d.json`: 2D n-channel MOSFET benchmark.
  P-type Si body (5 µm × 2 µm, N_A = 1 × 10¹⁶ cm⁻³) with two Gaussian
  n+ source/drain implants (center at x = 0.5 µm and 4.5 µm, σ = (0.4 µm,
  0.15 µm), peak = 5 × 10¹⁹ cm⁻³, donor). Gate sweeps 0 → 1.5 V at
  V_DS = 0.05 V.
- `tests/fem/test_bias_sweep_multiregion.py`: new test module with
  `test_bias_sweep_multiregion_multistep`. Runs a 3-step forward bias ramp
  (0 → 0.15 V) on a minimal 1D multi-region config and asserts at least
  3 IV points and non-decreasing electron current.
- `docs/adr/0008-snes-tolerances.md`: ADR documenting the SNES tolerance
  change rationale, validation, consequences, deferred work, and
  alternatives.

### Changed
- `semi/runners/bias_sweep.py` SNES tolerances in `run_bias_sweep`:
  - `snes_rtol`: 1.0e-14 → 1.0e-10 (achievable at high injection).
  - `snes_atol`: 1.0e-14 → 1.0e-7 (absolute floor for continuity block).
  - `snes_stol`: unchanged at 1.0e-14 (displacement convergence kept tight).
  - `snes_max_it`: 60 → 100 (headroom for fine MOSFET meshes).
  Added an inline comment referencing ADR 0008.
- `README.md` Status section updated to v0.12.0.
- `pyproject.toml` version bumped `0.11.0 → 0.12.0`.
- `semi/__init__.py` `__version__` bumped `0.11.0 → 0.12.0`.

## [0.11.0] - M11: Schema versioning

### Added
- `schemas/input.v1.json`: the project's Draft-07 input schema, extracted
  verbatim from the old `semi.schema.SCHEMA` dict. Every `type: object`
  node carries a `description`; leaf properties carry descriptions and,
  where applicable, `default`/`examples` annotations that a UI form
  builder can render as labels and help text.
- Required top-level `schema_version` field on every input JSON
  (pattern `^\d+\.\d+\.\d+$`). The engine constant
  `semi.schema.ENGINE_SUPPORTED_SCHEMA_MAJOR = 1` is compared against
  the input's major version; mismatched majors raise `SchemaError`.
  Minor/patch skew is accepted silently.
- Symmetric major-version gate on the *output* side:
  `semi.io.reader.ENGINE_SUPPORTED_MANIFEST_MAJOR = 1` is compared
  against the `schema_version` of a manifest being read back;
  mismatches raise `ValueError`. Manifests without `schema_version`
  (pre-M9 artifacts) are still loadable with a warning.
- `GET /schema` now returns
  `{"schema": <Draft-07 schema>, "version": "1.0.0", "supported_major": 1}`
  instead of the bare schema, so UI clients can discover the advertised
  schema version and the engine's supported major in one round trip.
- `schema_version: "1.0.0"` as the first key of every shipped benchmark
  JSON (`pn_1d`, `pn_1d_bias`, `pn_1d_bias_reverse`, `mos_2d`,
  `resistor_3d`, `resistor_gmsh`).
- `tests/test_schema_versioning.py`: 8 pure-Python tests gating
  Draft-07 validity of the extracted schema, description presence on
  every object node, benchmark-JSON semver compliance, major-mismatch
  rejection, minor-skew acceptance, missing-`schema_version`
  rejection, loader caching (identity), and manifest major-mismatch
  rejection.

### Changed
- `semi/schema.py` is now a loader: the module-level `SCHEMA` dict is
  populated from `schemas/input.v1.json` via an `@lru_cache`'d
  `_load_schema()`. The public API (`SCHEMA`, `validate`, `load`,
  `dumps`, `SchemaError`) is unchanged; downstream imports need no
  edits.
- `tests/test_kronos_server.py::test_schema` asserts the new
  `{schema, version, supported_major}` shape rather than a bare
  schema.
- `pyproject.toml` version bumped `0.10.0 -> 0.11.0`.
  `[tool.hatch.build.targets.wheel]` now `force-include`s `schemas/`
  so `pip install kronos-semi` ships both `input.v1.json` and
  `manifest.v1.json` inside the wheel.
- `semi/__init__.py` `__version__` bumped `0.10.0 -> 0.11.0`.

### Known issues
- **Input schema is permissive about unknown keys.** The top-level
  schema and every nested object lack `"additionalProperties": false`,
  so a UI typo in a field name (e.g., `"voltag": 0.5` instead of
  `"voltage"`) validates successfully and is silently ignored. The
  manifest schema (M9) is strict via `additionalProperties: false`
  throughout; the input schema should be too, eventually. Flipping
  this is a breaking change -- every existing external JSON with
  extra keys would fail. Defer to a major schema bump (v2.0.0) in a
  later milestone.
- **No separate `schema_version` for sub-blocks.** Today, if the
  `physics` block gains a new required key in v1.1.0, every caller
  must upgrade atomically. A per-block version would allow partial
  upgrades. Not needed now; flag for when it actually bites.
- **Defaults documented in schema annotations but synthesized by
  `_fill_defaults`.** Two sources of truth that must stay in sync.
  M11 hand-copies them; a future refactor could drive defaults from
  the schema directly (e.g., using a library like
  `jsonschema-default`). Low priority unless the list of defaults
  grows significantly.

## [0.10.0] - M10: HTTP server

### Added
- `kronos_server/` top-level package exposing the M9 engine over HTTP.
  FastAPI app factory (`build_app`), in-process `ProcessPoolExecutor`
  worker pool (`JobManager`), LocalFS storage backend, file-backed
  progress event stream (`progress.ndjson`), pydantic request/response
  models, and CORS middleware. No dolfinx, UFL, or PETSc imports at
  module scope in the server process; FEM-heavy imports happen only
  inside worker subprocesses spawned with `mp_context="spawn"`.
- HTTP endpoints: `POST /solve` (202/400), `GET /runs`,
  `GET /runs/{id}` (returning the full manifest plus status),
  `GET /runs/{id}/{manifest,input,fields/{name},iv/{contact},logs}`,
  `GET /schema`, `GET /materials`, `GET /capabilities`, `GET /health`,
  `GET /ready` (200/503).
- WebSocket endpoint: `WS /runs/{id}/stream` tails the progress file and
  emits `step_done` and `run_done` messages; closes with code 4404 on
  unknown run, 1000 on completion.
- `kronos-server` console entry point (`pyproject.toml`) with env-var
  driven host/port/workers/runs-dir/cors-origins configuration.
- `server` service in `docker-compose.yml`, bound to port 8000.
- `[project.optional-dependencies]` extra `server = [fastapi, uvicorn,
  httpx, pydantic]`.
- Optional `progress_callback=None` parameter on
  `semi/runners/{equilibrium,bias_sweep,mos_cv}.py`, called per bias
  step; this is the only permitted change to `semi/` in M10.
- `tests/test_kronos_server.py` — 15 tests covering pure-Python and
  FEM-backed paths: health/ready, capabilities, schema, materials,
  OpenAPI, input validation, 404 handling, CORS preflight, end-to-end
  solve+manifest match, field download, IV download, input download,
  WebSocket progress (≥ 5 `step_done` on `pn_1d_bias`), and concurrent
  solves (3 simultaneous).

### Changed
- `Dockerfile`: `pip install -e ".[dev,server]"` so the baked image
  includes the server extra.
- `pyproject.toml` `[tool.hatch.build.targets.wheel]` packages now
  include `kronos_server`.
- `PLAN.md`: M10 moved to completed work log; next task set to M11.

## [0.9.0] - M9: Result artifact writer

### Added
- `schemas/manifest.v1.json`: JSON-schema Draft-07 describing the run-artifact
  manifest contract. Covers engine metadata, solver summary, field list,
  mesh topology, optional bias-sweep entries, and warnings.
- `semi/io/artifact.py`: `write_artifact(result, out_dir, run_id, input_json_path)`
  writes a versioned run directory containing `manifest.json`, `input.json`,
  `mesh/mesh.xdmf`, `fields/<name>.bp` (or `.xdmf` fallback), `iv/<contact>.csv`,
  and `convergence/snes.csv`.
- `semi/io/reader.py`: `read_manifest(run_dir)` loads and schema-validates
  `manifest.json` using only stdlib JSON plus optional jsonschema; no dolfinx
  or numpy dependency, safe for M10 server subprocess.
- `semi/io/__init__.py`: exposes `write_artifact` and `read_manifest`.
- `semi/io/cli.py`: `main()` function powering the `semi-run` CLI entry point.
- `tests/test_artifact.py`: assertions across pure-Python schema tests and
  FEM-heavy parametrized benchmark tests (skipif when dolfinx unavailable).
- `[project.scripts]` entry in `pyproject.toml`: `semi-run = "semi.io.cli:main"`.

### Changed
- `PLAN.md`: M9 moved to completed work log; next task set to M10.
- `CHANGELOG.md`: M9 entry added.

## [0.8.0] - M8: Submission polish

### Added
- Four Colab notebooks in `notebooks/`: `01_pn_junction_1d.ipynb`,
  `02_pn_junction_bias.ipynb`, `03_mos_cv.ipynb`,
  `04_resistor_3d.ipynb`. Each notebook installs FEniCSx via
  FEM on Colab, clones the repo, loads a benchmark JSON, runs the
  solver, and plots results against analytical references. No local
  install required.
- `scripts/build_notebook_0[1-4].py`: generator scripts that produce
  the notebooks from Python source so the cell content is reviewable
  and version-controlled outside the notebook JSON.
- `docs/submission-polish-log.md`: engineering decision log capturing
  the M8 phase plan, per-phase commit history, reviewer-caught issues
  (MOS psi-reference convention, verifier window shift, Colab gmsh
  libGLU dependency), and Colab QA results.

### Changed
- `README.md`: replaced the M1-era status section with a shipped-feature
  capability matrix covering all 14 capabilities across PRs #2-#9.
  Added a notebook catalog table with Colab badges for all four
  notebooks. Added an explicit scope/out-of-scope section.
- `PLAN.md`, `docs/ROADMAP.md`: updated to reflect M8 complete.
- `CHANGELOG.md`, `docs/`, `scripts/`, `tests/`: replaced "Day N"
  internal sprint labels with milestone names (M1-M8) throughout.
- `semi/__init__.py`: version bumped from `0.1.0` to `0.8.0`.

### Fixed
- `notebooks/04_resistor_3d.ipynb`: added `apt-get install libglu1-mesa`
  before `pip install gmsh` to resolve `OSError: libGLU.so.1` on fresh
  Colab sessions. Confirmed working with dolfinx 0.10.0.post5 and
  gmsh 4.15.2.

### Verified
- All four notebooks execute end-to-end on Colab (dolfinx 0.10.0.post5):
  NB01 (Vbi 0.834 V, W 147 nm), NB02 (J(0.6V) ~1.6e3 A/m²),
  NB03 (C-V worst error 9.25% at V_gate = -0.200 V),
  NB04 (R_sim 1.115 kΩ, both mesh variants).
- pytest 206/206 passed, V&V 62/62 PASS, ruff clean on main.

## [0.7.0] - M7: 3D doped resistor

### Added
- `benchmarks/resistor_3d/resistor.json`: 3D doped rectangular bar
  (1 um x 200 nm x 200 nm), uniform n-type N_D = 1e18 cm^-3, two
  ohmic contacts on the x=0 and x=L faces, 5-point bias sweep on
  `contact_right` in [-0.01, +0.01] V. First benchmark in the
  project to exercise the solver stack in 3D.
- `benchmarks/resistor_3d/resistor_gmsh.json` plus a committed
  `fixtures/box.geo` / `fixtures/box.msh` pair. The `.geo` is the
  reproducible source; `gmsh -3 box.geo -o box.msh` regenerates
  the `.msh` under 100 KB. Same device on an unstructured
  tetrahedral mesh, exercised by the end-to-end test matrix.
- `semi/mesh.py::_build_from_file`: `.msh` loader wired via
  `dolfinx.io.gmsh.read_from_msh`. Physical groups stored in the
  file are returned verbatim as `cell_tags` and `facet_tags`, so
  `build_mesh` skips the JSON box/plane tagger when the mesh
  brings its own tags. XDMF is left as a clear
  `NotImplementedError` for a future PR.
- `semi/runners/bias_sweep.py`: two-leg walk for bipolar sweeps.
  When the resolved `v_sweep_list` spans zero, the bias ramp runs
  from `V = 0 -> min(V) -> max(V)` with a fresh
  `AdaptiveStepController` on each leg. Unipolar sweeps (the
  pn-junction and MOS benchmarks) fall through to the original
  single-endpoint ramp unchanged, so 1D/2D numerics are
  byte-identical to the M6 baseline.
- `scripts/run_benchmark.py`: `verify_resistor_3d` V-I linearity
  verifier (max |R_sim - R_theory|/R_theory < 1% tolerance;
  sanity checks: |I(V=0)| under numerical-noise floor, and
  `sign(I(V)) == sign(V)` at every nonzero V). 3D slice plotters
  for psi and |J_n| at the `y = W/2` midplane, plus an I-V
  scatter against the theoretical line; the 1D and 2D plotters
  are untouched.
- `docs/resistor_derivation.md` (derivation-lite gate, approved
  before any M7 implementation): device geometry, analytical
  ohmic resistance `R = L / (q N_D mu_n A) = 1115 Ohm`, 1% V-I
  linearity rationale tied to the `V -> 0, uniform-everything`
  limit, 3D slice-plot strategy, gmsh loader test strategy.
- `docs/PHYSICS.md` Section 7 (3D extension notes, ~1 page):
  cites that Poisson and drift-diffusion are dimension-independent
  and that 3D is a recipe-level extension, not a physics change;
  references `docs/resistor_derivation.md` for the full
  derivation. Covers the device, analytical R, bipolar-sweep
  driver, V-I verifier, and gmsh loader.
- `tests/fem/test_mesh_gmsh.py`: round-trip load of `box.msh`,
  physical-group preservation, and equivalent-V-I result vs the
  builtin `create_box` path.
- `tests/fem/test_resistor_3d.py`: coarsened smoke test
  (16x4x4 mesh, shrunken sweep) plus builtin and gmsh
  production-JSON tests that run the verifier end-to-end and
  assert max |R_sim - R_theory|/R_theory < 1%.
- `tests/test_bipolar_sweep.py`: pure-Python test that constructs
  a `v_sweep_list` crossing zero and asserts the bias-sweep
  driver populates `bipolar_legs` at the expected endpoints; a
  companion unipolar test asserts `bipolar_legs` stays empty so
  1D/2D numerics remain on the single-endpoint ramp path.

### Changed
- `semi/mesh.py` module docstring promotes the gmsh loader from
  "stubbed / NotImplementedError" to a supported file source.
- `docs/ROADMAP.md`: M7 moves from Planned to Done; M8 is
  queued next.

### Verified
- `docker compose run --rm benchmark resistor_3d` exits 0 on both
  builtin and gmsh variants; V-I linearity well inside the fixed
  1% tolerance. Three plots written per variant.
- 1D and 2D benchmarks (`pn_1d`, `pn_1d_bias`, `pn_1d_bias_reverse`,
  `mos_2d`) byte-identical to the M6 baseline.
- `python scripts/run_verification.py all` clears all gates;
  finest-pair MMS rates (including `2d_multiregion`) within 0.01
  of the post-M6 values.
- `pytest --cov=semi --cov-fail-under=95` exits 0.

### Notes
- No changes to `semi/physics/poisson.py` or
  `semi/physics/drift_diffusion.py` were required on this branch:
  both forms are written against the abstract `nabla` operator
  and the cell measure `dx` and are dimension agnostic. The
  scaled coefficients `L_D^2 eps_r` (Poisson) and `L_0^2 mu_hat`
  (continuity) pick up no new dimensional factors in 3D; every
  `nabla` still contributes one `1 / L_0` and every `dx` still
  contributes `L_0^d`, with the `d` cancelling once the weak
  form is divided through by the reference density rate.
- On this branch two commits landed a first pass of the gmsh
  loader that failed the Dockerized FEM CI job (e395830 wired
  the loader against a wrong dolfinx namespace; 3abd170 aligned
  the docstrings but did not fix the import). Both runs failed
  in the `docker-fem` step before any V-I verifier ran; the fix
  landed in 62f52d2, which re-wired `_build_from_file` against
  `dolfinx.io.gmsh.read_from_msh` and greens the full benchmark
  matrix. No `main`-branch artifact ever saw a failing build.

## [0.6.0] - M6: 2D MOS capacitor

### Added
- `semi/physics/poisson.py:build_equilibrium_poisson_form_mr`:
  multi-region equilibrium Poisson residual. Stiffness integrates
  over the full mesh with a cellwise DG0 eps_r Function; the
  Boltzmann space-charge term is restricted to semiconductor cells
  via `dx(subdomain_id=semi_tag)`, so the oxide region carries only
  the Laplacian. The interface natural condition
  eps_r_Si grad psi . n = eps_r_ox grad psi . n is enforced
  automatically by the piecewise eps_r in the bilinear form (see
  `docs/mos_derivation.md` section 3).
- `semi/runners/mos_cv.py:run_mos_cv`: sweeps a gate contact through
  its `voltage_sweep`, solves multi-region equilibrium Poisson at
  each V_gate, and records (V_gate, Q_gate) rows by integrating the
  silicon space charge and dividing by the lateral extent. Registers
  under the new `solver.type == "mos_cv"` dispatch.
- `benchmarks/mos_2d/mos_cap.json`: 500 nm p-type Si substrate
  (N_A = 1e17 cm^-3) with 5 nm SiO2 gate oxide. Uniform 1 nm vertical
  mesh (505 cells) lands the Si/SiO2 interface on a grid line.
  V_gate sweeps [-0.9, +1.2] V in 0.05 V steps (43 points).
- `scripts/run_benchmark.py`: 2D `plot_mos_2d` renders four PNGs
  (`psi_2d.png` tricontourf, `potentials_1d.png` central-column
  psi(y), `cv.png` simulation vs depletion-approximation theory with
  verifier-window overlay, `qv.png` Q_gate(V_gate)). `verify_mos_2d`
  extracts C_sim by centered finite difference and asserts
  |C_sim - C_theory|/C_theory < 10% in the window
  [V_FB + 0.2, V_T - 0.1] V. The 1D plotters (`plot_pn_1d`,
  `plot_pn_1d_bias`, `plot_pn_1d_bias_reverse`) are untouched.
- `docs/PHYSICS.md` Section 6 (condensed MOS reference: device
  structure, equilibrium-Poisson model, V_FB / V_T formulae under
  our psi=0-at-intrinsic convention, capacitance extraction, verifier
  results). Section 5.5 adds the multi-region Poisson MMS entry.
- `tests/fem/test_mos_cv.py` (4 tests): gate sweep behaviour,
  flatband bulk equilibrium, multi-region form byte-identity vs
  scalar form on a single-region mesh, and malformed-config error
  path.
- CI `.github/workflows/ci.yml` runs the `mos_2d` benchmark after
  the three `pn_1d*` benchmarks.

### Changed
- `semi/schema.py`: adds `"mos_cv"` to the `solver.type` enum.
- `semi/run.py` dispatcher routes `solver.type == "mos_cv"` to
  `run_mos_cv`. `semi/runners/__init__.py` exports the new runner.

### Notes
- Worst C-V error in the verifier window is 9.25% at V_gate = -0.20 V
  (V_FB + 0.217 V, the window's low edge). The window low edge was
  shifted from V_FB + 0.1 to V_FB + 0.2 during development: the
  initial 0.1 V margin hit 10.06% at V_gate = -0.25 V, and per the
  reviewer's guidance the intended fix is to shrink the window, not
  to loosen the 10% tolerance. Regimes outside the window
  (accumulation, strong inversion) are physically real but are not
  captured by the depletion approximation; the simulator does render
  them on the C-V plot.
- Body ohmic BC sets psi_body = -phi_F (our psi=0-at-intrinsic
  convention), not 0. This shifts V_FB to `phi_ms - phi_F`, not the
  textbook `phi_ms`, and V_T shifts correspondingly. Documented in
  `docs/PHYSICS.md` Section 6.3.
- Coverage is 95.43% on `semi/` (195 tests pass). The new
  `mos_cv.py` contributes ~100 statements, of which 8 error-branch
  lines are uncovered; the 95% CI gate still passes.

## [0.5.0] - M5: Refactor and test pass

### Added
- `semi/bcs.py`: pure-Python boundary-condition module (joins the
  pure-Python core tier per Invariant 4). Public API:
  - `ContactBC` dataclass holding (name, kind, facet_tag, V_applied,
    work_function) for one resolved contact.
  - `resolve_contacts(cfg, facet_tags=None, voltages=None)` walks the
    JSON config, validates kinds, resolves string facet refs via
    `mesh.facets_by_plane`, applies an optional voltages override (the
    bias-sweep code path), and optionally asserts the tag has facets.
  - `build_psi_dirichlet_bcs(...)` for the equilibrium Poisson row.
  - `build_dd_dirichlet_bcs(...)` for the coupled (psi, phi_n, phi_p)
    block.
- `semi/runners/` package (`equilibrium.py`, `bias_sweep.py`,
  `_common.py`) and `semi/postprocess.py` (terminal-current evaluation,
  IV recording, facet metadata helpers). `semi/run.py` is now a 74-line
  thin dispatcher; `bias_sweep.py` is 294 lines; no module in `semi/`
  exceeds 300 lines (was 580 in `run.py`).
- `tests/test_bcs.py` (11 pure-Python tests),
  `tests/fem/test_bcs.py` (3 FEM tests comparing the extracted
  builders byte-for-byte against an embedded copy of the legacy
  inline implementation),
  `tests/fem/test_mesh.py` (7 direct unit tests for `build_mesh`),
  `tests/fem/test_solver.py` (3 SNES-wrapper tests on linear Poisson
  and a synthetic 2-block problem),
  `tests/fem/test_recombination_ufl.py` (2 tests of the UFL
  `srh_rate` helper),
  `tests/test_convergence_helpers.py` (14 pure-Python tests for
  `_convergence` CSV/plot/table helpers, previously at 0% coverage),
  `tests/test_run_shim.py` (5 pure-Python tests for the
  `semi.run` dispatcher and backward-compat `__getattr__` shim),
  plus 7 new pure-Python tests in `tests/test_scaling.py` and
  `tests/test_doping.py`.
- ADR 0007 `docs/adr/0007-contact-bc-interface.md` records the BC
  interface design (dataclass over dict, voltages override over
  config mutation, insulating-skip semantics, optional facet-tag
  verification).
- `docs/PHYSICS.md` Section 2.5 completed (was an M2 placeholder):
  full scaled-DD derivation showing the L_0^2 coefficient on the
  continuity rows, scaled SRH kernel, residual-sign block summary,
  and J_0 numerical reality check. Section 3.1 BC reference
  redirects from the now-removed inline helpers in `run.py` to
  `semi.bcs`.

### Changed
- `semi/run.py` reduced from 580 lines (M4) to 74 (M5).
  `SimulationResult` and `run(cfg)` stay here; everything else moved
  to `runners/` and `postprocess.py`.
- `semi/run.py` exposes `run_equilibrium`, `run_bias_sweep`,
  `_fmt_tag`, and `_resolve_sweep` via a `__getattr__` shim so the
  V&V suite and existing pure-Python tests keep working unchanged.
- `pyproject.toml` adds `[tool.coverage.run]` / `[tool.coverage.report]`
  with `pragma: no cover`, `raise NotImplementedError`,
  `if __name__ == "__main__":`, and `if TYPE_CHECKING` exclusions.
- `to_table_rows`, `report_table`, `write_artifacts`, `run_cli_study`
  in the verification modules and `run_bias_sweep_with_continuity`
  in `semi/verification/conservation.py` are marked
  `# pragma: no cover`: they exist only to drive the end-to-end V&V
  CLI (`scripts/run_verification.py all`) which the docker-fem CI
  step exercises on every push.
- CI `docker-fem` job adds a coverage step that fails the build if
  `pytest --cov=semi --cov-fail-under=95` drops below 95%.

### Verified
- `pytest`: 177 / 177 pass (139 prior + 38 new).
- `docker compose run --rm benchmark pn_1d`: V_bi 0.8334 V, peak |E|
  104.84 kV/cm (byte-identical to the M4 baseline).
- `docker compose run --rm benchmark pn_1d_bias`: J(V=0.6) 1.635e+03
  A/m^2, J_total continuity 1.91% at V=0.6 (byte-identical).
- `docker compose run --rm benchmark pn_1d_bias_reverse`: |J| / SRH-gen
  worst 16.6% at V=-0.55 V (byte-identical).
- `python scripts/run_verification.py all`: 53 PASS / 0 FAIL,
  finest-pair MMS rates byte-identical to the M4 baseline (no rate
  moved by more than 0.000).
- Coverage: 96.25% across 1598 statements, 60 missed; gate at
  `--cov-fail-under=95`.

## [0.4.0] - M4: V&V suite

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
- M4 scope re-purposed: the originally-planned refactor pass is
  pushed to M5 to make room for the V&V work (see
  `docs/adr/0006-*.md` context).

## [0.3.0] - M3: Adaptive continuation

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
  diffusion at V = 0.6 V. M2 "qualitative below 0.5 V" caveat is
  removed.
- `run_bias_sweep` walks the ramp adaptively instead of recording at
  every `voltage_sweep.step` point. Total SNES iterations on the
  M2 sweep (0 to 0.6 V) dropped from 42 to 31 (26.2% reduction).

## [0.2.0] - M2: Coupled drift-diffusion

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

## [0.1.0] - M1: Equilibrium Poisson

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
