# M19 starter prompt: 3D MOSFET capstone benchmark

## Context

You are working in `kronos-semi` at `v0.25.0` (post-M18, package version
`0.25.0`, schema `2.9.0`). The repo is at
`https://github.com/rwalkerlewis/kronos-semi`. `main` is at the commit
that closed M18 (adaptive-dt transient runner).

**IMPORTANT:** There may be unstaged modifications in the working tree
from prior in-progress work. Before creating the M19 branch, stash these
with `git stash` so the branch starts clean from main. Restore the stash
on main after the milestone is complete (or leave it stashed; do not
commit those changes as part of M19).

Your assignment is **M19: 3D MOSFET capstone benchmark**.

## Scope

M19 ships a real 3D n-channel MOSFET on a gmsh-sourced unstructured
tetrahedral mesh, exercising:

- The M15 GPU linear-solver path on a real device (not just on Poisson).
- The M16.1 Caughey-Thomas mobility under non-trivial 3D fields.
- The M16.4 Fermi-Dirac statistics (optional but recommended at n+ doping).
- The M14.3 XDMF / gmsh mesh ingest path (already stable).

M19 does **not** include MPI parallel orchestration (that is M19.1), HTTP
hardening (M20), or the bias-sweep SNES stabilization that would unblock
`nmos_idvgs` (separate follow-up). Do not touch `nmos_idvgs`, `mosfet_2d`,
or any existing benchmark.

## Device spec

| Parameter | Value |
|---|---|
| Channel length L | 250 nm |
| Channel width W | 1 um |
| Gate oxide t_ox | 5 nm SiO2 |
| Body | p-type Si, N_A = 1e16 cm^-3 |
| Source/drain | Gaussian n+ implants, peak N_D = 5e19 cm^-3, sigma_x = 0.4 um, sigma_y = 0.15 um |
| Mobility | Lombardi (M16.2) with Caughey-Thomas bulk (M16.1) |
| Statistics | Fermi-Dirac / Blakemore (M16.4) |
| Starting mesh size | ~200k DOFs (scalable to ~500k for GPU acceptance) |

## Branch and PR rules

- Stash uncommitted working-tree changes first: `git stash`.
- Work on a fresh branch: `git checkout -b dev/m19-mosfet-3d`.
- One milestone, one PR. Open the PR after Phase 0.
- **Run all phases consecutively without pausing for confirmation, status
  reports, or "should I continue" prompts.** Do not ask permission between
  phases. Run Phase 0 through Phase F end-to-end, push commits as they land.

## Phases

### Phase 0 - Starter prompt + schema bump

1. Write this prompt verbatim to `docs/M19_STARTER_PROMPT.md` and commit:
   `docs: ship M19 starter prompt (M19)`.
2. Bump schema in `schemas/input.v2.json` and `semi/schema.py`: schema
   `2.9.0` -> `2.10.0`. The only schema addition is allowing
   `solver.backend` to coexist with `solver.type: "bias_sweep"` in the 3D
   MOSFET JSON (it already does; just bump `SCHEMA_SUPPORTED_MINOR` to 10
   and add a `schema_version: "2.10.0"` example to the schema). Commit:
   `feat(schema): schema 2.10.0; M19 3D MOSFET config (M19)`.

### Phase A - gmsh geometry and mesh

1. Create `benchmarks/mosfet_3d/` directory.
2. Write `benchmarks/mosfet_3d/mosfet3d.geo`: a gmsh `.geo` file describing
   the 3D NMOS device:
   - Silicon body region (x in [0, L_total], y in [0, W], z in [0, H_si]
     where L_total = 1 um, W = 1 um, H_si = 200 nm).
   - SiO2 gate-oxide region (same x/y footprint, z in [H_si, H_si + t_ox]).
   - Physical volumes: `silicon` (tag 1) and `oxide` (tag 2).
   - Physical surfaces: `source` (tag 10), `drain` (tag 11), `gate`
     (tag 12), `body` (tag 13), all remaining surfaces `insulating`
     (tag 0).
   - Mesh size characteristic length cl = 20 nm (yields ~200k DOFs on a
     uniform tetrahedral mesh).
3. Generate and commit the mesh: `python -c "import gmsh; gmsh.initialize();
   gmsh.open(benchmarks/mosfet_3d/mosfet3d.geo);
   gmsh.model.mesh.generate(3); gmsh.write(benchmarks/mosfet_3d/mosfet3d.msh);
   gmsh.finalize()"`.
4. Write `benchmarks/mosfet_3d/mosfet_3d.json`: the simulation JSON for a
   V_GS sweep [0.0, 2.0] V at 0.2 V step, V_DS = 0.05 V (linear regime),
   Lombardi+CT mobility, Fermi-Dirac statistics, `solver.backend: "auto"`.
5. Commit: `feat(benchmarks): 3D MOSFET geometry, mesh, and JSON config (M19)`.

### Phase B - Verifier and analytical reference

1. In `semi/diode_analytical.py`, add `mosfet_3d_paosah_iv(V_GS_arr, V_DS,
   mu_eff, C_ox, L, W, V_T, vsat)` - the Pao-Sah linear-regime formula
   `I_D = (W/L) * mu_eff * C_ox * max(V_GS - V_T, 0) * V_DS` plus the
   velocity-saturation correction
   `mu_eff / (1 + mu_eff * V_DS / (vsat * L))`.
2. In `scripts/run_benchmark.py`, add `verify_mosfet_3d`:
   - Run `mosfet_3d.json` via `semi.run.run(cfg)`.
   - Extract I_D at each V_GS from the drain contact IV.
   - Compare to the analytical Pao-Sah reference within 25% in the linear
     regime ([V_T + 0.2, V_T + 0.8] V window).
   - Assert I_D increases monotonically with V_GS above threshold.
   - Report PASS/FAIL with observed vs expected values.
3. Commit: `feat(verifier): Pao-Sah 3D MOSFET verifier with 25% tolerance (M19)`.

### Phase C - Saturation regime benchmark

1. Add `benchmarks/mosfet_3d/mosfet_3d_sat.json`: same device but
   V_DS = 1.0 V (saturation regime), V_GS sweep [0.0, 2.0] V at 0.4 V step
   (5 points to keep runtime tolerable).
2. In `scripts/run_benchmark.py`, extend `verify_mosfet_3d` to also run the
   saturation config and check I_DSAT within 30% of the velocity-saturation
   analytical reference `I_DSAT = (W / (2*L)) * mu_eff * C_ox *
   (V_GS - V_T)^2 / (1 + (V_GS - V_T) / (2 * vsat * L / mu_eff))`.
3. Add `tests/test_mosfet_3d_verifier.py`: 4 pure-Python assertions on the
   analytical helpers (threshold, linear regime, saturation regime,
   monotonicity).
4. Commit: `feat(benchmarks): saturation-regime 3D MOSFET config and verifier (M19)`.

### Phase D - GPU acceptance test

1. Add `benchmarks/mosfet_3d/mosfet_3d_gpu.json`: same as `mosfet_3d.json`
   but `solver.backend: "gpu-amgx"` and mesh refined to ~500k DOFs (use a
   smaller cl in the `.geo`).
2. In `scripts/run_benchmark.py`, add `verify_mosfet_3d_gpu`:
   - Run the GPU config; if `gpu-amgx` is unavailable (backend resolves to
     `cpu-mumps`), skip the GPU timing gate and report `SKIP (no GPU)`.
   - If GPU is available: assert CPU/GPU linear-solve wall-clock ratio >= 5x.
   - Assert `psi` field is finite and non-NaN.
3. Add `.github/workflows/gpu-nightly.yml` matrix entry for `mosfet_3d_gpu`
   (gate on `vars.GPU_RUNNER_AVAILABLE == "true"`).
4. Commit: `feat(benchmarks): GPU acceptance config for 3D MOSFET at ~500k DOFs (M19)`.

### Phase E - CI wiring and tests

1. Add `mosfet_3d` and `mosfet_3d_sat` to `.github/workflows/ci.yml`
   docker-fem matrix (no `allow-failure`).
2. Add `tests/fem/test_mosfet_3d.py`: a coarsened (cl = 50 nm, ~20k DOFs)
   smoke test that runs the linear-regime sweep for 3 V_GS points and
   asserts finite I_D and monotone I_D(V_GS). Mark with
   `@pytest.mark.slow` if runtime > 60s at 20k DOFs.
3. Run `ruff check semi/ tests/ scripts/ benchmarks/ && pytest tests/ -x -q`
   (pure-Python suite) and fix any issues.
4. Commit: `test(mosfet_3d): coarsened FEM smoke + CI matrix wiring (M19)`.

### Phase F - Closeout

1. Update `PLAN.md`:
   - Move `Next task` pointer to "Bias-sweep SNES line-search stabilization
     (unblocks `nmos_idvgs`)" (the other candidate from PLAN.md).
   - Add M19 to the completed work log (append-only, newest on top).
   - Update `Current state` to reflect M19 shipped in v0.26.0.
2. Update `docs/IMPROVEMENT_GUIDE.md`:
   - Mark M19 row as Done in the capability matrix.
   - Update the "Honest gap" section.
3. Update `docs/ROADMAP.md` M19 row from Planned to Done.
4. Update `CHANGELOG.md` with a `[0.26.0]` entry.
5. Bump `pyproject.toml` and `semi/__init__.py` from `0.25.0` -> `0.26.0`.
6. Commit: `docs: close out M19 (PLAN, IMPROVEMENT_GUIDE, PHYSICS, ROADMAP,
   CHANGELOG, ADR placeholder) (M19)`.
7. Push all commits and open a PR titled
   "M19: 3D MOSFET capstone benchmark (#XX)".

## Acceptance tests (summary)

1. `python scripts/run_benchmark.py mosfet_3d` exits 0, linear-regime I_D
   within 25% of Pao-Sah reference, I_D monotone above V_T.
2. `python scripts/run_benchmark.py mosfet_3d_sat` exits 0, I_DSAT within
   30% of velocity-saturation reference.
3. GPU run exits 0 (or `SKIP` if no GPU hardware); if GPU available,
   linear-solve speedup >= 5x vs CPU-MUMPS at ~500k DOFs.
4. Existing benchmarks bit-identical: `pn_1d_bias` J(V=0.6 V) = 1.635e+03
   A/m^2; all other anchors unchanged.
5. `pytest tests/ -x -q` passes (pure-Python suite green, coverage gate
   >= 95).

## Invariants (do not violate)

- JSON is the only input format (ADR 0001).
- Slotboom primary variables only (ADR 0004).
- dolfinx 0.10 API only (ADR 0003).
- Pure-Python core must not import dolfinx (ADR 0007).
- Schema changes are additive minor bumps only.
- No em dashes in prose.
- Physics-style variable names (N_A, V_t, etc.) are correct; do not rename
  to PEP 8.
- No co-authored-by Claude credits in commit messages.

## Notify on completion

When completely finished (all phases committed and pushed, PR open), run:
`openclaw system event --text "M19 3D MOSFET capstone complete: PR open,
v0.26.0, schema 2.10.0" --mode now`
