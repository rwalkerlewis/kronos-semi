# ADR 0010: BDF Time Integration and Lumped Mass Matrix

**Date:** 2026-04-25
**Status:** Accepted
**Milestone:** M13 — Transient solver

---

## Context

The transient solver needs a time-integration scheme for the coupled
(ψ, n, p) system. The choice of scheme affects:

1. **Accuracy order** — how fast the temporal error converges as dt → 0.
2. **Stability** — whether the scheme remains bounded for reasonable dt.
3. **Start-up** — whether the method requires history from multiple past
   steps (multi-step methods need seeding).
4. **Mass matrix** — whether to use the standard consistent mass matrix
   (full bilinear form ∫ n v dx) or the lumped diagonal approximation.

## Decision

### Time integration: BDF1 / BDF2

Use the **backward differentiation formula** (BDF) family at order 1 or 2.
The scheme is specified by the BDF coefficients:

| Order | α₀  | α₁  | α₂  |
|-------|-----|-----|-----|
| 1     | 1   | −1  |     |
| 2     | 3/2 | −2  | 1/2 |

The discrete time derivative at step n+1 is:

```
(1/dt) · Σ_{k=0}^{K} α_k · u^{n+1−k}
```

For BDF2, a single BDF1 step seeds the two-level history before switching
to BDF2.

**Rationale for BDF family:**
- BDF1 and BDF2 are both **A-stable** (unconditionally stable for linear
  problems), which is critical for semiconductor problems where the
  time scales span many orders of magnitude.
- The implicit formulation integrates naturally with the existing SNES
  nonlinear block solver: the residual at each timestep is assembled and
  passed to `solve_nonlinear_block` unchanged.
- BDF1 provides 1st-order temporal accuracy; BDF2 provides 2nd-order.
- Higher-order BDF (3–6) exist but are not A-stable for order ≥ 3 and
  are deferred to a future milestone.

**Adaptive time stepping** is not implemented in M13. Fixed dt is used
throughout. The `max_steps` safety cap prevents infinite loops.

### Mass matrix: consistent (Galerkin) in UFL forms

The time derivative term is assembled using the **standard consistent mass
matrix** via UFL:

```python
alpha_0 / dt * n * v_n * ufl.dx  # → ∫ n v_n dx (consistent mass)
```

The `assemble_lumped_mass` function in `semi/fem/mass.py` computes the
**lumped (row-sum) mass diagonal** via the functional ∫ v dx (by the
partition-of-unity property of P1 elements). This vector is stored in
`TransientResult.meta` and is available for M14 (AC small-signal analysis)
where the mass matrix appears explicitly in the `(J + iωM)` system.

**Rationale for consistent mass:**
- Consistent mass is the standard Galerkin approach; it requires no
  special assembly and is directly supported by UFL.
- For implicit BDF time integration (unconditionally stable), consistent
  vs lumped mass does not affect stability — only the error constant.
- The convergence rates (BDF1: 1st order, BDF2: 2nd order) are the same
  for consistent and lumped mass on uniform meshes.
- Lumped mass is beneficial mainly for (a) explicit schemes (improves
  stability threshold) and (b) conditioning of the mass matrix when
  used in an M-matrix context. Neither benefit applies here.

## Consequences

**Positive:**
- Simple implementation: no post-assembly diagonal replacement needed.
- SNES and NonlinearProblem work unchanged.
- `assemble_lumped_mass` is available for M14 without any transient
  runner changes.

**Negative:**
- Consistent mass slightly increases the bandwidth of the Jacobian
  compared to lumped mass, but the Poisson-DD Jacobian already has
  full P1×P1 coupling so this is not a new concern.
- Consistent mass with BDF2 can introduce mild dispersive errors on
  very coarse meshes; lumped mass is sometimes more accurate there.
  For the P1 element on 1D uniform meshes this effect is negligible.

**Deferred:**
- Adaptive time stepping (dt controller based on local error estimate).
- BDF3+ for 3rd-order accuracy.
- Explicit use of the lumped mass diagonal in the transient residual.

## Alternatives considered

1. **Crank-Nicolson (θ-method, θ=0.5)**: 2nd-order accurate but only
   A-stable for θ ≥ 0.5; the standard CN scheme introduces numerical
   dispersion on coarse meshes. BDF2 is preferred because it is also
   A-stable and has better damping properties.

2. **Runge-Kutta DIRK**: higher-order implicit RK is well-studied but
   requires evaluating the residual multiple times per step. Not
   compatible with the single-residual NonlinearProblem API. Deferred.

3. **Operator splitting**: decouple Poisson from continuity and advance
   each independently. Splits the physical coupling and introduces a
   splitting error; not suitable for a production solver.
