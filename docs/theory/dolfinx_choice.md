# Why dolfinx 0.10

The `NonlinearProblem` class in dolfinx 0.10 wraps PETSc SNES
directly (the old `NewtonSolver` is deprecated), giving us:

- Line search (`newtonls + bt` is the engine default).
- Convergence diagnostics from SNES (residual histories, divergence
  reasons) accessible via `petsc4py`.
- Block-system support for the coupled $(\psi, \Phi_n, \Phi_p)$
  solve through `kind="nest"` or `kind=None` on a sequence of
  residual forms.

## API expectations the codebase relies on

- `NonlinearProblem` accepts `F, u` as sequences plus
  `kind=None|"nest"` for blocked problems.
- SNES options can be set via the petsc4py options DB *before*
  `NonlinearProblem.__init__`, but factor-matrix options
  (e.g. MUMPS `cntl_3`) must be set on the factor `Mat` directly,
  not via the options DB. This was the root cause of the M13.1
  PR #52 misadventure (see [adr/0014](../adr/) and PR #54).

## Default SNES settings

- `newtonls + bt`, `rtol = atol = 1e-10/1e-12` for scalar Poisson.
- For the DD block, `rtol = atol = stol = 1e-14` so convergence is
  not declared prematurely while continuity residuals are still
  large.
- Transient factorisation: 200 % MUMPS workspace
  (`mat_mumps_icntl_14=200`) and a small Jacobian shift
  (`solver.jacobian_shift = 1e-14`) applied via the SNES Jacobian
  callback to handle rank-deficient quasi-Fermi rows in the deep
  p-bulk.

## Where this lives in the code

- [`semi/solver.py`](../../semi/solver.py) — SNES wrapper.
- [`semi/runners/`](../../semi/runners/) — every runner builds its
  own `NonlinearProblem` with runner-specific tolerances.

## ADR

- [adr/0003-dolfinx-0-10-api.md](../adr/0003-dolfinx-0-10-api.md).
