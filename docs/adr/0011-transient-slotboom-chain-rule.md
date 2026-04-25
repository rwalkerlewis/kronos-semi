# ADR 0011: Transient Slotboom Continuity with Chain-Rule Temporal Form

**Date:** 2026-04-25
**Status:** Accepted (supersedes ADR 0009)
**Milestone:** M13 — Transient solver, post-merge fix

---

## Context

ADR 0009 selected the (ψ, n, p) primary-density formulation for the
transient continuity equations. The motivation was that ∂n/∂t and
∂p/∂t are linear in (n, p) under that choice, giving a symmetric
positive-definite mass matrix.

In practice the M13 transient tests (`tests/fem/test_transient_steady_state`,
`tests/mms/test_transient_convergence`) fail because the SNES solver
drives n or p negative during a Newton iterate. The (ψ, n, p) Newton
update has no built-in positivity guarantee, and the moment any DOF
goes negative the SRH denominator τ_p (n + n₁) + τ_n (p + p₁) can
collapse and either generate a NaN or lock the iterate at a
non-physical fixed point. Mesh refinement does not fix this — the
issue is structural to the formulation, not a discretisation error.

## Decision

Switch the transient continuity equations to **Slotboom (ψ, φ_n, φ_p)
primary unknowns** matching the steady-state runner, and apply BDF
temporal discretisation through the **chain rule** on the carrier
densities:

    n = n_i exp(ψ - φ_n)        [scaled, V_t = 1]
    p = n_i exp(φ_p - ψ)

    ∂n/∂t = n (∂ψ/∂t - ∂φ_n/∂t)
    ∂p/∂t = p (∂φ_p/∂t - ∂ψ/∂t)

with each ∂/∂t approximated by BDF on the corresponding primary
unknown. Define `BDF[u] = α₀/dt · u^{n+1} + f_hist_u` where
`f_hist_u = Σ_{k≥1} α_k/dt · u^{n+1-k}` carries past values stored as
numpy arrays. The three-block residual is:

    F_ψ:    L_D² ε_r ∇ψ · ∇v_ψ - (p - n + N) v_ψ                = 0
    F_φn:   n · (BDF[ψ] - BDF[φ_n]) · v_n
            + L₀² μ_n n ∇φ_n · ∇v_n - R · v_n                    = 0
    F_φp:   p · (BDF[φ_p] - BDF[ψ]) · v_p
            + L₀² μ_p p ∇φ_p · ∇v_p + R · v_p                    = 0

Boundary conditions match the steady-state Slotboom runner: ohmic
contacts pin φ_n = φ_p = V_applied / V_t (Shockley boundary) and
ψ = arcsinh(N_net / 2 n_i) + V_applied / V_t.

History storage tracks (ψ, φ_n, φ_p) arrays directly, sized to
`order + 1` levels. Initial condition: ψ from a V=0 Poisson
equilibrium solve, φ_n = φ_p = 0.

## Consequences

**Positive:**

- n and p are exponentials of the primary unknowns and are therefore
  strictly positive at every Newton iterate. The negative-density
  failure mode of ADR 0009 is structurally eliminated.
- (ψ, φ_n, φ_p) at an ohmic contact differs by exactly the equilibrium
  arcsinh value regardless of bias, so the contact carrier densities
  are pinned to their equilibrium majority/minority values throughout
  the transient — natural Shockley boundary conditions.
- Code reuse: the same `build_dd_dirichlet_bcs`,
  `evaluate_partial_currents`, and Slotboom helpers used by
  `run_bias_sweep` apply unchanged.

**Negative:**

- The chain-rule time term carries a factor n_ufl (or p_ufl) that
  varies by ~25 orders of magnitude across the device (intrinsic
  region vs. doped region), which makes the Newton Jacobian
  ill-conditioned near BDF fixed points. MUMPS LU may fail with
  `DIVERGED_LINEAR_SOLVE` (reason -3); UMFPACK may produce NaN
  (`DIVERGED_FNORM_NAN`, reason -4).
- This is the textbook Slotboom-Jacobian-conditioning issue: standard
  remedy is Scharfetter-Gummel edge-flux discretisation (planned as
  the M13 follow-up "Fix 3"), which absorbs the n_ufl scaling into
  exponentially-weighted basis functions.

**Neutral:**

- The Poisson block is identical across the steady-state and
  transient runners.
- Recombination R is the same UFL expression as the steady-state
  Slotboom runner.

## Alternatives considered

1. **Direct BDF on n_ufl** (Selberherr / Markowich textbook form):
   replace the chain-rule time term with `α₀/dt · n_ufl + f_hist_n`
   where `f_hist_n[i] = Σ alpha_k/dt · (n_i exp(ψ^k - φ_n^k))[i]`.
   Tested empirically: Newton stagnates at the first time step
   because the time-mass term has the same ~25 OOM range as the
   chain-rule form, and the effective diagonal scaling is no better.

2. **Keep ADR 0009's (ψ, n, p) form and add positivity projection**
   (clamp n, p ≥ 0 inside the SNES residual). Rejected: introduces a
   non-smooth nonlinearity that breaks Newton's quadratic
   convergence and makes Jacobian assembly conditional.

3. **Scharfetter-Gummel discretisation.** This is the standard fix
   for the wide-range n_ufl conditioning problem and is the planned
   "Fix 3" if the chain-rule form does not meet the M13 acceptance
   tests on its own. Deferred to a follow-up because it is a
   larger refactor (per-edge flux assembly rather than UFL
   stiffness assembly) and warrants its own ADR.
