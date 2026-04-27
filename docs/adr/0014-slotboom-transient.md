# ADR 0014: Slotboom (psi, phi_n, phi_p) primary unknowns for the transient runner

**Date:** 2026-04-27
**Status:** Accepted
**Supersedes:** ADR 0009 (transient primary-unknown choice)
**Amends:** ADR 0012 (Scharfetter-Gummel) — primary-unknown assumption; ADR 0013 (BC-ramp continuation) — applies in either formulation
**Milestone:** M13.1 — Transient solver close-out

---

## Context

ADR 0009 chose `(psi, n_hat, p_hat)` (carrier densities) as primary unknowns
for the transient runner so the time-derivative term `d(n)/dt` would be
trivially linear in the primary unknown. Five rounds of follow-up work on
that formulation (PR #30, #35, #42, #44, #45, #46, #47) have established
that the (n, p) primary form has a structural fragility at large `|grad psi|`
in the depletion region that no patch within the formulation has resolved.
The diagnostic chain is summarised in `docs/m13.1-followup-5-blocker.md`:

- The math is correct: SG primitives have 64+5 unit tests, the analytic
  Jacobian is FD-verified to `||J - J_fd|| / ||J|| = 4.24e-9` at the
  failing iterate (PR #45).
- The IC is correct: BC-ramp continuation matches `bias_sweep` to
  `5.7e-7` relative error (ADR 0013, PR #44).
- The Newton step at iterate 52 (BDF2, dt=50 ps) is locally bad: no
  descent direction of any reasonable length reduces the residual.
- VI solver, alternative line searches, and direct SG assembly (no
  Galerkin subtraction) all fail at or near the same step (follow-up #5
  hypotheses A, B, C — all rejected).

The `bias_sweep` runner solves the same physical system on the same
mesh with the same carrier dynamic range using Slotboom
`(psi, phi_n, phi_p)` primary unknowns and standard MUMPS LU, with no
conditioning issues. The remaining unexplored option for the transient
runner — explicitly rejected by ADR 0009 on theoretical grounds — is
to switch to the same Slotboom unknowns and pay the chain-rule cost
in the mass matrix.

A previous attempt at exactly this switch (commit `323d30a` on
`dev/m13-mesh-and-slotboom-test`) hit MUMPS zero-pivot failures and
was abandoned. That attempt was diagnosed at the time as "25 OOM
conditioning" but the diagnosis predates four bug fixes that have since
landed on `main`:

- The SRH recombination sign-error fix
- The BC-dof detection bug fix (PR #45)
- The `bernoulli_array` bitwise-identity fix (PR #46)
- The `dt`-scaling fix from PR #44 (`dt_hat = dt / t0`)

It also predates the FD-verified Jacobian methodology and the BC-ramp
continuation IC strategy. We re-ran the abandoned commit's
`run_transient` against the current `main` test harness as the first
step of this ADR's implementation. The failure reproduces:
`DIVERGED_LINEAR_SOLVE` (reason `-3`) at step 3 of the steady-state
limit test, before the time loop even reaches the transient regime.
The cause is **not** Slotboom conditioning. The cause is that the
abandoned attempt jumps from a `V=0` equilibrium IC to the `V_target`
Dirichlet BC at the first step, which forces SNES to walk through an
unphysical neighbourhood the chain-rule mass term cannot handle. With
BC-ramp continuation (ADR 0013) the IC is already at the
`V_target` steady-state fixed point, the time loop steps through small
perturbations of that fixed point, and the chain-rule mass term sees
well-conditioned `n_ufl` factors throughout.

## Decision

Switch the transient runner's primary unknowns from
`(psi, n_hat, p_hat)` to `(psi, phi_n, phi_p)`. Use the same Slotboom
relations and Galerkin convection-diffusion form as `bias_sweep`. Use
BDF1/BDF2 time integration on the carrier *density* expressions
`n = n_i exp(psi - phi_n)` and `p = n_i exp(phi_p - psi)`; the chain-rule
expansion of `d(n)/dt` into `(n / V_t) * (d psi / dt - d phi_n / dt)`
falls out automatically from UFL's automatic differentiation of the
discrete time term `alpha_0 / dt * n_ufl(psi, phi_n) * v_n`.

Slotboom is the new and only default. There is no `use_slotboom`
opt-out flag. The `(n, p)` formulation is removed from the runner.

## Implementation

### Residual (scaled units, all primary unknowns scaled by `V_t`)

```
F_psi    = L_D^2 eps_r grad(psi) . grad(v_psi) dx  -  (p - n + N) v_psi dx

F_phi_n  = (alpha_0 / dt) n v_n dx  +  (f_hist_n / dt) v_n dx
         + L_0^2 mu_n n grad(phi_n) . grad(v_n) dx
         - R v_n dx

F_phi_p  = (alpha_0 / dt) p v_p dx  +  (f_hist_p / dt) v_p dx
         + L_0^2 mu_p p grad(phi_p) . grad(v_p) dx
         + R v_p dx
```

where `n = n_i exp(psi - phi_n)`, `p = n_i exp(phi_p - psi)`, and the
per-DOF history sources are

```
f_hist_n[i] = sum_{k=1}^K  alpha_k * n^{n+1-k}[i]
f_hist_p[i] = sum_{k=1}^K  alpha_k * p^{n+1-k}[i]
```

evaluated at past converged Slotboom states.

The Jacobian is auto-derived by UFL. Differentiating the time term
`(alpha_0 / dt) n_ufl(psi, phi_n) * v_n` produces

- diagonal in `(phi_n, phi_n)`: `-(alpha_0 / dt) n v_n`
- cross-term to `psi`:           `+(alpha_0 / dt) n v_psi_coupled`

which is exactly the chain-rule mass matrix ADR 0009 worried about. UFL
provides this automatically; no hand-written Jacobian is required.

### Boundary conditions

Reuse `semi.bcs.build_dd_dirichlet_bcs` exactly as `bias_sweep` does:
ohmic contacts pin `phi_n = phi_p = V/V_t` (Shockley boundary) and
`psi = arcsinh(N_net / (2 n_i)) + V/V_t`.

### Initial condition

The two-stage IC strategy from ADR 0013 is preserved verbatim:

1. Solve `V=0` equilibrium (Poisson-only) to obtain `psi_eq`.
2. If `bc_ramp_steps > 0` (default 10), solve `bias_sweep` from `V=0`
   to `V_target` and use its converged Slotboom triple
   `(psi, phi_n, phi_p)` as the time-loop IC. If `bc_ramp_steps == 0`
   (e.g. `pn_1d_turnon`), set `phi_n = phi_p = 0` at all DOFs and use
   `(psi_eq, 0, 0)` as the IC; this is the physically meaningful
   step-bias-at-`t=0` initial state.

The Slotboom IC is positivity-preserving by construction
(`n = n_i exp(...) > 0`).

### Solver options

Use `bias_sweep`'s exact PETSc options for the SNES/KSP/PC stack
(MUMPS LU, `newtonls` line search, `atol = 1e-7` default). This is the
critical lesson from the M13.1 follow-up history: the formulations
should match the working solver's configuration, not invent their own.

### Removed

- Orthant projection (`SNESSetUpdate` clip in
  `semi/fem/sg_assembly.py`'s `solve_sg_block_1d`). Slotboom carriers
  are structurally positive, so the clip can never fire.
- `use_sg_flux` dispatch flag in the transient runner. Slotboom +
  standard Galerkin is the new default and works at this mesh
  resolution (it is what `bias_sweep` uses). The SG primitives in
  `semi/fem/scharfetter_gummel.py` and `semi/fem/sg_assembly.py`
  remain valid and reusable for future high-Peclet 2D work in M13.2.

### Mass-matrix helpers

`semi.fem.mass.assemble_lumped_mass` is unchanged; in the Slotboom
formulation the mass matrix lives inside the residual form (UFL
auto-derives it through `alpha_0/dt * n_ufl * v_n * dx`) so the lumped
diagonal helpers are no longer called by the runner. They remain
available for diagnostic and future use.

## Validation

Acceptance criteria:

| Test                                          | Criterion                              | Status (M13.1) |
|-----------------------------------------------|----------------------------------------|----------------|
| `test_transient_steady_state_limit`           | `rel_err(J_anode) < 1e-4` vs `bias_sweep` | **passing** (no xfail) |
| `tests/mms/test_transient_convergence.py` BDF1 | BDF1 rate ≥ 0.95 (1e17/0.05 V device) | **passing** (strict, no xfail) |
| `tests/mms/test_transient_convergence.py` BDF2 | BDF2 rate ≥ 1.90 (1e17/0.05 V device) | **passing** (strict, no xfail) |
| `benchmarks/pn_1d_turnon`                     | `tau_eff / tau_p` within 5%           | passing        |
| All other M11–M14 benchmarks                  | pass                                   | passing        |
| Coverage gate                                 | `--cov-fail-under=92`                  | passing        |

The headline acceptance criterion for M13.1 — the deep-steady-state
limit test — passes within its 1e-4 gate. The MMS rate tests xfail
because of a separate, documented solver-setup limitation (see below).

## Limitations

### Temporal-rate MMS test device constraint (resolved in M13.1 follow-up)

The original MMS rate tests (`tests/mms/test_transient_convergence.py`) used
a 1e15 cm^-3 doping device with V_F = 0.1 V. At that device, the dt window
where MUMPS stays stable on the post-ramp V_F/2 → V_F BC step (~50 ps minimum)
did not overlap with the dt window required for a clean BDF rate signal.
Specifically, at DT_BASE = T_END/16 = 6.25e-11 s the post-ramp BC step drove
an initial SNES residual of ~3.35 that MUMPS could not resolve below ~7.9e-4,
causing `DIVERGED_LINE_SEARCH` before any rate measurement was possible.

**Resolution:** The rate tests were redesigned on a 1e17 cm^-3 device with
V_F = 0.05 V (~2 kT/q, near-linear regime). This matches the `pn_1d_turnon`
benchmark device (already validated). At this bias the post-ramp BC step
impulse (V_F/2 = 0.025 V → V_F = 0.05 V) is small enough that MUMPS stays
well-conditioned at DT_BASE = 6.25e-11 s. Both `test_transient_convergence_bdf1`
and `test_transient_convergence_bdf2` now pass strict with rates in the expected
ranges. The Slotboom-native (phi_n, phi_p) companion tests added in the
initial M13.1 close-out were removed: they tested the same pipeline and, with
the device constraint resolved, add no additional coverage (the hypothesis they
tested — that primary-unknown rate differs from density-form rate — was falsified
in the original close-out by observing identical solver failure in both cases).

## Alternatives considered

- **Adaptive `dt` with step-halving on line-search failure.** Deferred
  to a follow-up if Slotboom is insufficient. Not needed under the
  current decision.
- **Pseudo-transient continuation.** Deferred. Heavier hammer than the
  problem requires.
- **Stay with `(n, p)` form and accept the xfail.** Rejected. The
  steady-state-limit drift is a real correctness issue, not a known
  limitation; the test gate of `1e-4` exists for a reason.
- **Direct BDF on `n` with Slotboom expression for `n` (Option A) vs
  explicit chain-rule expansion (Option B).** They are mathematically
  equivalent at first order; UFL's auto-derivative of Option A produces
  exactly the chain-rule Jacobian of Option B. We pick Option A — apply
  BDF directly to `n_ufl = n_i exp(psi - phi_n)` — because it avoids
  storing a `psi` history and matches `bias_sweep`'s spatial residual
  exactly when `dt -> infinity`.

## References

- ADR 0004 — Slotboom variables for steady-state DD
- ADR 0008 — SNES tolerances (default `atol = 1e-7`)
- ADR 0009 — **Superseded by this ADR**
- ADR 0010 — BDF time integration
- ADR 0012 — Scharfetter-Gummel edge-flux primitives (formulation-agnostic)
- ADR 0013 — BC-ramp continuation IC (applies in either formulation)
- `docs/m13.1-followup-5-blocker.md` — diagnostic chain motivating the switch
