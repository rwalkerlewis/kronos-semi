# ADR 0009: Transient Continuity Formulation — Density (n, p) vs Slotboom (φ_n, φ_p)

**Date:** 2026-04-25
**Status:** Superseded by ADR 0011
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
