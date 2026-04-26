# ADR 0009: Transient Continuity Formulation — Density (n, p) vs Slotboom (φ_n, φ_p)

**Date:** 2026-04-25
**Status:** Accepted
**Milestone:** M13 — Transient solver

---

## Context

The steady-state drift-diffusion solver (`semi/runners/bias_sweep.py`) uses
Slotboom quasi-Fermi potentials (ψ, φ_n, φ_p) as primary unknowns. This
choice is standard in steady-state FEM semiconductor solvers because it
turns the potentially ill-conditioned carrier transport equation into a
formally coercive form (see ADR 0004).

For the transient solver we need to add a time-derivative term to the
continuity equations. The natural form of the continuity equations in
time is:

```
∂n/∂t + div(J_n) = -R
∂p/∂t + div(J_p) = +R
```

where the time derivative acts on **carrier densities** n and p, not on
the quasi-Fermi potentials φ_n and φ_p. If we keep φ_n as the primary
unknown, the time derivative becomes:

```
∂n/∂t = ∂/∂t [n_i exp(ψ - φ_n)] = n_i exp(ψ - φ_n) * (∂ψ/∂t - ∂φ_n/∂t)
```

This couples the Poisson and continuity blocks through ∂ψ/∂t and results
in a non-symmetric mass matrix that is expensive to linearise.

## Decision

For the **transient runner only**, switch the primary unknowns from
(ψ, φ_n, φ_p) to **(ψ, n_hat, p_hat)**, where n_hat = n/C₀ and
p_hat = p/C₀ are the scaled carrier densities.

The continuity residuals in n/p form are:

```
Electron:  ∂n/∂t + div(J_n) = -R
           with J_n = L₀² μ_n (n∇ψ - ∇n)   [scaled units]

Hole:      ∂p/∂t + div(J_p) = +R
           with J_p = L₀² μ_p (∇p + p∇ψ)   [scaled units]
```

The BDF weak form of the time derivative for the electron equation is:

```
α₀/dt · ∫ n · v_n + L₀² μ_n · ∫ (∇n − n∇ψ) · ∇v_n + ∫ R · v_n
+ ∫ f_hist_n · v_n = 0
```

where `f_hist_n = Σ_{k=1}^K α_k/dt · n^{n+1-k}` is updated each timestep
from the BDF history.

**Ohmic contact boundary conditions** in n/p form:
- ψ_bc: unchanged (depends on applied bias)
- n_hat_bc = n_i_hat · exp(arcsinh(N_net/(2 n_i))): **constant**, independent of V
- p_hat_bc = n_i_hat · exp(−arcsinh(N_net/(2 n_i))): **constant**, independent of V

The carrier BCs do not change with the applied bias because both ψ and
φ_n/φ_p shift by the same amount V/V_t at an ohmic contact, leaving
n = n_i · exp(ψ − φ_n) invariant.

The steady-state runner (`run_bias_sweep`) is **not modified**. It
continues to use the Slotboom (φ_n, φ_p) formulation.

## Consequences

**Positive:**
- The time derivative in (n, p) form is naturally linear in the primary
  unknowns, giving a block-diagonal mass matrix.
- The consistent mass matrix `∫ n · v_n` is symmetric positive definite,
  improving SNES convergence for the transient block.
- The carrier contact BCs in n/p form are constant (no update needed
  during time stepping), simplifying the time loop.

**Negative:**
- The n/p form has a different coercivity structure than the Slotboom
  form. For very large carrier gradients (extreme forward bias) the
  drift term n · ∇ψ may cause minor instabilities compared to the
  Slotboom form.
- Current evaluation post-processing must recover φ_n and φ_p from n,
  p via Slotboom inversion (φ_n = ψ − log(n/n_i)), adding a small
  overhead.

**Neutral:**
- The Poisson block is identical in both formulations; only the carrier
  density expressions differ (n_hat directly vs n_i_hat · exp(ψ − φ_n)).

## Alternatives considered

1. **Keep (ψ, φ_n, φ_p) for transient** and add `M · ∂n/∂(φ_n) · ∂φ_n/∂t`
   to the residual. Rejected because this requires assembling a coupled
   (ψ, φ_n) mass block and coupling the nonlinear function derivative
   into the transient term, adding significant complexity.

2. **Mixed formulation**: use φ_n for the stiffness block and n for the
   time derivative. Rejected for the same coupling reason.

3. **Explicit predictor-corrector**: avoid the implicit issue by using
   an explicit scheme for the time derivative. Rejected because explicit
   methods are conditionally stable with a CFL-like constraint on dt that
   would require very small steps for typical semiconductor problems.

## Known limitation (M13.1)

The (ψ, n_hat, p_hat) Galerkin discretization adopted above does **not**
converge to the same discrete steady state as the Slotboom (φ_n, φ_p)
discretization used by `run_bias_sweep`, even on identical meshes with
identical boundary conditions and arbitrarily tight SNES tolerance.

**Symptom.** A 1D pn junction held at forward bias for many minority-
carrier lifetimes settles to a discrete fixed point of the (n,p) form
that disagrees with `bias_sweep`'s Slotboom fixed point by ~3 × 10¹⁸
m⁻³ in carrier density at the depletion edge (orders of magnitude
larger than the minority-carrier density itself). The terminal IV from
the transient runner is the wrong magnitude and, depending on
post-processing, the wrong sign. With sufficiently tight SNES tolerance
(`atol = 1e-15`), the (n,p) Newton iterate can drive carriers to
~10²⁴ m⁻³ and ψ_hat to ~90, far outside the physical regime, before
locking at a non-physical fixed point.

**Cause.** This is a **spatial discretization** failure, not a time
integration one. The standard Galerkin discretization of the
convection-diffusion form `div(∇n − n∇ψ) = R` has no positivity guarantee
and is unstable for non-trivial cell Péclet numbers, exactly the regime
the depletion region sits in. The Slotboom form sidesteps the issue by
using exponentially-varying integrand `n_i exp(ψ − φ_n)` at quadrature
points, which absorbs the cell-Péclet variation into the basis. An
intermediate fix using Slotboom primary unknowns plus a chain-rule
expansion of ∂n/∂t was attempted on `dev/m13-mesh-and-slotboom-test`
(commit 323d30a) and hit MUMPS LU zero-pivot failures from a different
conditioning issue.

**Test status.**
* `tests/mms/test_transient_convergence.py` — **passing.** A pairwise-
  difference temporal-convergence test recovers BDF1 ≈ 1.0 and BDF2 ≈ 2.0
  on the (n,p) form, verifying that the time integrator and history
  bookkeeping are correct.
* `tests/fem/test_transient_steady_state.py` — **xfail** with `strict=False`,
  pending M13.1 fix. The test still encodes the correct acceptance
  criterion (relative IV error < 1e-4 vs `bias_sweep`); it will flip
  to passing without code change once SG lands.

**Resolution path (M13.1).** Implement Scharfetter-Gummel
exponential-fitting edge fluxes for the carrier continuity blocks. SG
absorbs the n_hat scaling into Bernoulli-weighted basis functions, is
positivity-preserving by construction, and is the standard choice in
production drift-diffusion codes. Estimated scope: a new edge-flux
assembler in `semi/fem/scharfetter_gummel.py`, ~600 LOC including
tests, plus a new ADR superseding the (n,p) Galerkin choice in this
document.

**Why not block on M13.1.** The transient infrastructure (BDF1/BDF2,
history bookkeeping, lumped mass for M14, IV recording, snapshot
output) is correct and is verified by the temporal convergence test.
Blocking the milestone on a new spatial discretizer would push
M14 (small-signal AC) by an estimated 2-3 sessions. Shipping with
the steady-state limit xfail'd is honest about the limitation while
unblocking downstream work. The default SNES `atol` for the transient
runner has been tightened from `1e-7` to `1e-10` so the (n,p) failure
is exposed as honest Newton work rather than silently masked by SNES
exiting at iteration 0.
