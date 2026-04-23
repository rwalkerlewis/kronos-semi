# 0003. Target dolfinx 0.10 API only

- Status: Accepted
- Date: 2026-04-20

## Context

FEniCSx has a rapidly evolving API. The relevant API surface for this
project is the nonlinear solver interface:

- **dolfinx 0.9 and earlier:** `dolfinx.nls.petsc.NewtonSolver` held a
  hand-rolled Newton loop. Setting PETSc options required either
  manipulating global `PETSc.Options()` or reaching through the solver
  object. Block/nested systems required extra plumbing in user code.
- **dolfinx 0.10:** `dolfinx.fem.petsc.NonlinearProblem` wraps PETSc
  SNES directly. It takes a `petsc_options_prefix` and a dict of PETSc
  options, giving clean per-problem configuration. It also supports
  blocked and nested function spaces via `kind="blocked"` or
  `kind="nest"`, which the coupled (psi, Phi_n, Phi_p) system will
  need for M2: Coupled drift-diffusion.

The M2: Coupled drift-diffusion block Newton implementation needs the 0.10 blocked
support. Supporting both APIs simultaneously would fork the solver
driver and double the test matrix.

## Decision

Target `dolfinx >= 0.10, < 0.11` as the only supported version. Use
`dolfinx.fem.petsc.NonlinearProblem` with `petsc_options_prefix` and
`petsc_options` (see `semi/solver.py`). Do not introduce or import
`dolfinx.nls.petsc.NewtonSolver`; the class is deprecated in 0.10 and
will be removed.

## Consequences

Easier:

- Single solver driver, single set of PETSc option semantics.
- Block-Newton (M2: Coupled drift-diffusion) uses the same `NonlinearProblem` API as
  scalar Newton, with an additional `kind` argument.
- SNES line search (`snes_linesearch_type`) and advanced convergence
  criteria are available out of the box.

Harder:

- Users on older dolfinx environments (for example, stale conda
  environments, legacy FEM-on-Colab kernels) cannot run the FEM layer.
  The Docker image in the repo pins to
  `ghcr.io/fenics/dolfinx/dolfinx:stable`, which currently resolves to
  0.10, so the supported path is unambiguous.
- When dolfinx 0.11 lands, every form and driver will need to be
  retested. We accept that cost rather than hedging.

## Related

- `semi/solver.py`, the single use site for `NonlinearProblem`.
- Invariant 5 in `PLAN.md`.
