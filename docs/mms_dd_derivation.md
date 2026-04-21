# MMS Derivation: Coupled Drift-Diffusion (Phase 4 Gate)

This document is the mathematical specification for the Phase-4 Method of
Manufactured Solutions (MMS) verification of the coupled drift-diffusion
block residual in `semi/physics/drift_diffusion.py`. It is a gate
artifact: no Python is written until this document is approved.

The MMS target is the three-field Slotboom-form system
(`build_dd_block_residual`, lines 140-151), whose scaled residual is

```
-div( L_D^2 eps_r grad psi_hat )            - (p_hat - n_hat + N_hat) = 0
-div( L_0^2 mu_n_hat n_hat grad phi_n_hat ) - R_hat                   = 0
-div( L_0^2 mu_p_hat p_hat grad phi_p_hat ) + R_hat                   = 0
```

with Slotboom relations

```
n_hat = ni_hat * exp(psi_hat - phi_n_hat)
p_hat = ni_hat * exp(phi_p_hat - psi_hat)
```

and SRH recombination

```
R_hat = ( n_hat * p_hat - ni_hat^2 ) / ( tau_p_hat * (n_hat + n1_hat)
                                       + tau_n_hat * (p_hat + p1_hat) )
n1_hat = ni_hat * exp( E_t/V_t ),  p1_hat = ni_hat * exp(-E_t/V_t)
```

For MMS we set `N_hat = 0` and add manufactured sources `f_psi, f_n, f_p`
to each block so that a chosen smooth exact triple `(psi_e, phi_n_e,
phi_p_e)` is the solution of the modified system. Each block then has a
discretization error that decays at the FE rate, and finest-pair rates
are asserted against the thresholds in Section 5.

---

## 1. Exact solutions

Proposed triple (all scaled units; Omega = (0, L)^d with d in {1, 2}):

**1D (d = 1):**

```
psi_e(x)   = A_psi * sin(2*pi * x / L)
phi_n_e(x) = A_n   * sin(  pi * x / L)
phi_p_e(x) = A_p   * sin(  pi * x / L)
```

**2D (d = 2):**

```
psi_e(x,y)   = A_psi * sin(2*pi * x/L) * sin(3*pi * y/L)
phi_n_e(x,y) = A_n   * sin(  pi * x/L) * sin(  pi * y/L)
phi_p_e(x,y) = A_p   * sin(  pi * x/L) * sin(  pi * y/L)
```

**Amplitudes** (default, held fixed across the mesh sweep):

```
A_psi = 0.5     # potential swing ~ 13 mV at 300 K; keeps exp(+/-psi) ~ [0.6, 1.7]
A_n   =  0.3    # electron quasi-Fermi swing
A_p   = -0.3    # hole quasi-Fermi swing, sign-flipped vs A_n
```

Rationale:

- All three fields vanish identically on the boundary of Omega, so the
  homogeneous-Dirichlet BC in Section 4 is exact (no interpolation noise).
- Distinct wavenumbers in psi vs. (phi_n, phi_p) prevent accidental
  cancellation in the combination `psi - phi_n` that appears in the
  Slotboom exponent; the two continuity blocks therefore exercise
  genuinely different gradients.
- `A_p = -A_n` (equivalently, phi_p_e and phi_n_e share a spatial
  profile but opposite sign) is required for Variant C to be a
  meaningful test of the SRH code path. With `A_p = +A_n` the
  difference `phi_p_e - phi_n_e` would vanish pointwise and the SRH
  numerator `ni_hat^2 * (exp(phi_p_e - phi_n_e) - 1)` would be
  identically zero, collapsing Variant C onto Variant B. The sign
  flip also matches physical convention: under forward bias the two
  quasi-Fermi levels split in opposite directions from equilibrium.
  Section 7.4's worst-case analysis (`phi_p_e - phi_n_e \in [-0.6, 0.6]`)
  is already correct for this sign choice.
- The 2D form uses two different wavenumbers on x and y for psi
  (2*pi, 3*pi) following the existing Poisson MMS convention
  (`semi/verification/mms_poisson.py:143-164`), which guarantees the
  discretization sees a 2D function and not a tensor product of 1D
  modes hit by one axis alone.
- Amplitudes keep `|psi_e - phi_n_e| <= A_psi + |A_n| = 0.8` and
  `|phi_p_e - psi_e| <= A_psi + |A_p| = 0.8`, so `exp(...)` stays in
  `[exp(-0.8), exp(0.8)] ≈ [0.45, 2.23]`. Newton is well-conditioned
  in this range and R_hat is neither machine-zero (near-equilibrium)
  nor catastrophically cancelling.

A second `nonlinear` amplitude set `(A_psi, A_n, A_p) = (1.0, 0.5, -0.5)`
is reserved for the CLI sweep only (Section 6). The pytest gate runs the
default amplitudes. The `A_p = -A_n` sign convention carries over.

## 2. Scaling reference values (consistency with Phase 1)

MMS reuses the Phase-1 scaling constants exactly (no re-derivation). The
values come from `semi/verification/mms_poisson.py:56-60` and are
carried unchanged into the DD MMS harness:

| Symbol          | Meaning                                    | Value                              |
|-----------------|--------------------------------------------|------------------------------------|
| `T_REF`         | Reference temperature                      | 300 K                              |
| `V_T = kT/q`    | Thermal voltage (V0 scale)                 | 0.025852 V                         |
| `L_0_REF`       | Device length (matches pn_1d benchmark)    | 2.0e-6 m                           |
| `C_0_REF`       | Density scale (1e16 cm^-3 floor)           | 1.0e22 m^-3                        |
| `EPS_R`         | Relative permittivity (Si)                 | 11.7                               |
| `mu_0 = mu_n`   | Reference mobility (Si electrons)          | from `materials.get_material("Si")`|
| `n_i`           | Intrinsic density (Si)                     | from `materials.get_material("Si")`|
| `D_0 = V_t mu_0`| Einstein diffusivity                       | derived                            |
| `t_0 = L_0^2/D_0`| Diffusion time                            | derived                            |
| `lambda2`       | `eps_0 V_t / (q C_0 L_0^2)` (bare)         | ~1.6e-3 (asserted < 1)             |
| `ni_hat`        | `n_i / C_0`                                | ~1.0e-6                            |

The scaling object is built by `build_mms_scaling(L=L_0_REF)`
(`mms_poisson.py:91-112`); the DD MMS will call the same factory to
guarantee Phase-1 consistency. The assertion `lambda2 < 1` already
present in that factory protects against future edits that would
trivialize the carrier nonlinearity.

**Mobility ratios.** `mu_n_over_mu0 = 1.0` (reference is electron
mobility itself); `mu_p_over_mu0 = mu_p / mu_n = 0.45/1.4 ≈ 0.321` using
the Si values in `materials.py`. These ratios are fixed across the sweep
so the convergence rates measure discretization error only.

**SRH lifetimes.** For Variant C: `tau_n = tau_p = 1.0e-7 s` (Si
mid-gap default), yielding `tau_hat = tau / t_0 ≈ tau * D_0 / L_0^2`.
For Variants A and B: `tau_hat = 1.0e+20` (effectively infinity so
`R_hat ≈ 0`). `E_t / V_t = 0` (mid-gap traps).

## 3. Forcing derivations (weak form)

Let `v_psi, v_n, v_p` denote the three test functions on P1 Lagrange
spaces. The modified block residual is

```
F_psi_mms = F_psi_prod - f_psi_weak
F_n_mms   = F_n_prod   - f_n_weak
F_p_mms   = F_p_prod   - f_p_weak
```

Each `f_*_weak` is constructed by substituting the exact triple into the
production strong form and integrating by parts against the test
function. Boundary terms vanish because (a) the exact fields vanish on
the boundary, and (b) the test functions vanish on the Dirichlet
boundary. We follow the same weak-form construction as the Poisson MMS
(`mms_poisson.py:167-200`): we do NOT form `ufl.div(ufl.grad(...))`
directly, because constant prefactors can collapse to numerical zero on
coarse meshes. Writing the source in weak form uses the same gradient
discretization as the production residual, so MMS verifies the
discretization actually in use.

Let

```
n_e = ni_hat * exp(psi_e - phi_n_e)
p_e = ni_hat * exp(phi_p_e - psi_e)
R_e = (n_e*p_e - ni_hat^2) / ( tau_p_hat*(n_e + n1_hat)
                             + tau_n_hat*(p_e + p1_hat) )
```

### 3.1 Psi block (all variants)

Strong:  `f_psi = -div(L_D^2 eps_r grad psi_e) - (p_e - n_e)`
(recall `N_hat = 0`). Weak form applied as `-f_psi_weak`:

```
f_psi_weak =   L_D^2 * eps_r * inner(grad(psi_e), grad(v_psi)) * dx
             - (p_e - n_e) * v_psi * dx
```

This is identical in structure to the Poisson MMS weak source but with
`(p_e - n_e)` computed from the full Slotboom triple rather than
`ni_hat * (exp(-psi_e) - exp(psi_e))` (which is the `phi_n_e = phi_p_e =
0` special case).

### 3.2 Electron continuity block

Strong:  `f_n = -div(L_0^2 mu_n_hat n_e grad phi_n_e) - R_e`. Weak:

```
f_n_weak =   L_0^2 * mu_n_hat * n_e * inner(grad(phi_n_e), grad(v_n)) * dx
           - R_e * v_n * dx
```

Note the `n_e` weighting on the diffusion integrand comes from the
Slotboom current `J_n = -q D_n n grad phi_n`; it is pointwise nonzero
everywhere in Omega because `ni_hat * exp(...) > 0`, so the block is
coercive (ADR 0004).

### 3.3 Hole continuity block

Strong:  `f_p = -div(L_0^2 mu_p_hat p_e grad phi_p_e) + R_e`. Weak:

```
f_p_weak =   L_0^2 * mu_p_hat * p_e * inner(grad(phi_p_e), grad(v_p)) * dx
           + R_e * v_p * dx
```

Sign of the `R_e` term is opposite to the electron block, matching the
production residual at line 150.

### 3.4 The three variants

**Variant A -- psi-only (Poisson block sanity).**
Set `phi_n_e = phi_p_e = 0` (identically). Then
`grad(phi_n_e) = grad(phi_p_e) = 0`, the continuity diffusion integrals
drop out, `n_e = ni_hat * exp(psi_e)`, `p_e = ni_hat * exp(-psi_e)`,
and `n_e * p_e = ni_hat^2` so `R_e = 0`. The exact triple
`(psi_e, 0, 0)` is the solution of the modified coupled system with

```
f_psi_weak = (as in 3.1)
f_n_weak   = 0
f_p_weak   = 0
```

Variant A is essentially the existing Poisson MMS pulled through
`build_dd_block_residual` instead of `build_equilibrium_poisson_form`.
Purpose: detect any block-assembly error in the Poisson row (sign,
scaling, coupling to a zero Fermi level) that would be invisible to the
standalone Poisson MMS.

**Variant B -- full coupling, no recombination.**
All three fields nontrivial. Set `tau_n_hat = tau_p_hat = 1.0e+20` so
`R_e ≈ 0` to machine precision but R is not *identically* zero in the
code path (the SRH kernel is still evaluated). The forcings are the
full 3.1/3.2/3.3 weak forms with the `R_e * v` terms numerically
negligible. Purpose: verify the drift-diffusion operators (the
`n grad phi_n`, `p grad phi_p` pieces) in their coupled form without
having to isolate recombination numerics.

**Variant C -- full coupling with SRH.**
All three fields nontrivial, realistic lifetimes
`tau_n = tau_p = 1e-7 s`. Forcings are the full 3.1/3.2/3.3 weak forms
with a physically nonzero `R_e`. Purpose: verify the full residual,
including the recombination coupling between the two continuity blocks
through the shared `R_e` term.

## 4. Boundary conditions

All three fields vanish identically on `boundary(Omega)` by construction
(each exact-solution factor is a `sin(k*pi*x/L)` with `k` a positive
integer). The homogeneous-Dirichlet BC is therefore the exact BC and is
enforced identically to the Poisson MMS:

```
zero = fem.Constant(mesh, 0.0)
bdofs_psi   = locate_dofs_topological(V_psi,   fdim, boundary_facets)
bdofs_phi_n = locate_dofs_topological(V_phi_n, fdim, boundary_facets)
bdofs_phi_p = locate_dofs_topological(V_phi_p, fdim, boundary_facets)
bcs = [
    fem.dirichletbc(zero, bdofs_psi,   V_psi),
    fem.dirichletbc(zero, bdofs_phi_n, V_phi_n),
    fem.dirichletbc(zero, bdofs_phi_p, V_phi_p),
]
```

This choice (vs. interpolating the exact solution onto the boundary) is
deliberate: interpolation of a nonlinear trig function onto a coarse P1
boundary mesh introduces an O(h^2) boundary error that pollutes the
measured L^2 rate, especially on the coarsest levels. Using
`Constant(0.0)` keeps the BC exact across all refinements.

## 5. Expected convergence rates

Both spaces are P1 Lagrange on simplices (1D intervals or 2D triangles).
Standard FE theory on a smooth problem gives

```
|| u - u_h ||_{L^2}      = O(h^2)     (rate 2)
|| grad(u - u_h) ||_{L^2} = O(h^1)     (rate 1)
```

per-block. The coupled system inherits the same rates provided Newton
converges to machine precision on each mesh so that iteration error
does not floor the FE error. The derivation originally proposed

```
snes_rtol = 1.0e-14
snes_atol = 1.0e-16
snes_stol = 0.0
snes_max_it = 80
```

matching the Poisson MMS (`mms_poisson.py:259-264`) with a larger
`max_it` allowance because the coupled system needs more iterations.

**Amendment: atol = 0, stol = 1e-12 (implementation-time fix).** The
derivation's `snes_atol = 1e-16` was tuned to the Poisson block
residual scale and turned out to terminate SNES prematurely on the
continuity blocks in 2D. The cause is the block-residual scale
disparity already flagged in Section 7.3: the Poisson block carries
`L_D^2 ~ 1.6e-3 * L_0^2` while the continuity blocks carry `L_0^2`
alone, and in 2D the integration area further shrinks the continuity
residual by a factor of ~`L_0^2 ~ 4e-12`. In the 2D Variant-B and
Variant-C tests the continuity blocks' initial absolute residual
sat at roughly 1e-18 (already below the originally-proposed
`atol = 1e-16`), so SNES would quit after the very first iteration
without ever touching the Poisson block, which still had a ~1e-6
initial residual to burn down. The Poisson row then recorded
pre-Newton interpolation error and the measured rates came out
noisy.

The fix applied in `semi/verification/mms_dd.py:368-373` is:

```
snes_rtol   = 1.0e-14
snes_atol   = 0.0      # was 1.0e-16; removes the premature-termination floor
snes_stol   = 1.0e-12  # step-norm takes over as the stopping criterion
snes_max_it = 80
```

Rationale:

- Setting `atol = 0.0` defers absolute-residual termination
  entirely. Newton iterates until either (a) the relative residual
  drops by `rtol` across all blocks, or (b) the step norm drops
  below `stol`, whichever comes first.
- `stol = 1.0e-12` is the meaningful floor for this problem: the
  scaled unknowns are O(1) (all amplitudes stay in [-1, 1]), so a
  step norm of 1e-12 means every DOF has settled to within 10^-12
  of its converged value. This is well below the FE error at the
  finest mesh pair (L^2 ~ 1e-8 for N = 320 in 1D), so Newton is
  fully converged relative to the discretization error.
- `rtol = 1.0e-14` is retained from the original. In practice the
  `stol` criterion fires first on 1D and the `rtol` criterion
  fires first in 2D where the residual still has headroom.

**Sanity check.** After the amendment the 2D runs converge in 2 or
4 Newton iterations (2 for Variant A where only the Poisson block
has a non-trivial residual; 4 for B/C where all three blocks must
be relaxed), the Poisson block reaches its rate floor, and the
continuity blocks in Variants B and C record rates within 0.01 of
the theoretical 2.0 / 1.0 targets. The previous tuning produced
rates below 1.0 on the continuity blocks at N=64 because Newton
had stopped before the linear system had actually been solved.

**Variant A noise floor.** With `phi_n_e = phi_p_e = 0`, the
continuity-block manufactured sources are identically zero and
there is nothing to drive those fields off the initial zero. The
continuity-block L^2 errors therefore sit at pure roundoff
(1e-35 to 1e-41 across the sweep). Dividing one roundoff level
by another produces meaningless rates (they hop around zero),
so the Variant A continuity-block rates are reported in the CSV
but not gated. The Phase-4 CLI (`run_verification.py mms_dd`)
and the pytest gate both check only the psi block for Variant A;
all three blocks for Variants B and C. This is consistent with
the purpose of Variant A as stated in Section 3.4: it is a
Poisson-row sanity check routed through the coupled solver, not
a convergence test for the continuity rows in isolation.

**Gate thresholds** (finest-pair rate, asserted on all three blocks
independently):

| Block           | L^2 rate  | H^1 rate  |
|-----------------|-----------|-----------|
| psi             | >= 1.75   | >= 0.80   |
| phi_n           | >= 1.75   | >= 0.80   |
| phi_p           | >= 1.75   | >= 0.80   |

The 1.75 target (vs. the 1.85 used for standalone Poisson) leaves 0.25
headroom below theoretical 2 to absorb:
- block coupling: one block's iteration residual bleeds into another's
  error on any given mesh;
- SRH nonlinearity (variant C): the cross-block R term has a compound
  exponential that is more sensitive to mesh discretization than a
  bare exp;
- drift weighting in continuity: the `n_e grad phi_n_e` factor is
  mesh-variable, not constant, so the effective operator is not a
  clean Laplacian (see Section 7).

H^1 threshold is 0.80 rather than 0.85 for the same reasons. The
monotonicity check (every level's error strictly smaller than the
previous) is also asserted, as in the Poisson gate.

## 6. Mesh plan

**Pytest gate (CI):**

```
Ns_1d = [40, 80, 160]          # three levels, 2 rate estimates
Ns_2d = [16, 32, 64]           # three levels, 2D triangles only
```

Three levels is the minimum to measure a rate plus a second-finest-pair
sanity check against monotonicity. Four-level 1D `[40, 80, 160, 320]`
is reserved for the CLI. Pytest wall-clock budget: each DD MMS solve
is ~4-10x slower than Poisson MMS at the same N (three unknowns, SRH
exp stack in R), so three levels keeps the full DD MMS gate under
~60 s on a developer laptop.

**CLI artifact sweep (`scripts/run_verification.py`):**

```
Ns_1d_linear    = [40, 80, 160, 320]        # default amplitudes
Ns_1d_nonlinear = [40, 80, 160, 320]        # larger amplitudes
Ns_2d           = [16, 32, 64]              # triangles, default amplitudes
```

Writes `convergence.csv` and `convergence.png` per (dim, variant)
combination into `verification/mms_dd/<dim>_<variant>/`, mirroring the
Poisson MMS artifact layout.

## 7. Known risks

### 7.1 exp() nonlinearity and amplitude cap

The Slotboom form means every residual contains `exp(+/- psi_hat)` or
`exp(psi - phi_n)` in a way that is *pointwise* nonlinear and
*multiplicatively* coupled to `grad(phi_n)`. For amplitudes above ~1,
three things go wrong:

1. Newton ill-conditioning: the Jacobian entry d(n)/d(psi) = n itself,
   which spans 2+ orders of magnitude across Omega when `A_psi > 1`.
2. R_e denominator floor: at points where `n_e + p_e` is small (high
   |psi|), the SRH denominator can dip below machine epsilon.
3. Manufactured source magnitude asymmetry: the `p_e - n_e` term in
   f_psi grows exponentially with A_psi; the weak form still converges
   but the residual norms span 5+ orders of magnitude between blocks.

Mitigation:
- Pytest gate holds `(A_psi, A_n, A_p) = (0.5, 0.3, -0.3)`, keeping
  `exp(...) \in [0.45, 2.23]`.
- CLI nonlinear variant `(1.0, 0.5, -0.5)` still keeps arguments
  bounded by 1.5, `exp \in [0.22, 4.48]`, comfortably inside PETSc
  double precision. We do not add an `A_psi = 2` case here;
  the Poisson MMS already covers that stress regime for the scalar
  nonlinearity.

### 7.2 Galerkin drift-dominant rate degradation

The continuity block `L_0^2 mu_n_hat n_e inner(grad phi_n, grad v_n)`
is a diffusion-only operator in Slotboom form -- there is no
`grad(psi) . grad(phi_n)` term, because the drift has been absorbed
into `n_e = ni_hat exp(psi - phi_n)`. On its face, this means the
drift-dominant degradation that plagues direct (n, p) formulations
(where standard Galerkin loses 0.5 order at Peclet > 1) does NOT apply
here.

However, the *effective* operator seen by Newton at each step is the
linearization

```
J_n(phi_n_delta) = L_0^2 mu_n_hat * (
      n_e * inner(grad phi_n_delta, grad v_n)
    - n_e * phi_n_delta * inner(grad phi_n_e, grad v_n)   (from d n_e / d phi_n)
)
```

and the second term is an anti-drift-like reaction, so the Jacobian has
both a diffusion and an O(|grad phi_n|) reaction piece. For bounded
amplitudes this is harmless. For large amplitudes with rapid phi_n
variation the Peclet number equivalent `h * |grad phi_n_e| / 1` can
approach 1 at N=40 (h/L = 0.025, |grad phi_n_e| <= pi*A_n/L = pi*0.3/L
in scaled units h/L is the relevant factor so Pe ~ pi*0.3*0.025 ~ 0.024
-- well below 1). So at the amplitudes and meshes chosen we expect no
degradation.

Mitigation / monitoring:
- Keep `|A_n|, |A_p| <= 0.5` in all gated runs.
- If a future tightening to large `A_n` surfaces a rate floor at
  ~1.5 in the continuity blocks, document it here as the expected
  Galerkin signature rather than a bug, and consider SUPG or
  Scharfetter-Gummel stabilization (deferred, out of scope for Phase 4).

### 7.3 Block residual scale disparity

The Poisson block residual carries `L_D^2 ~ 1.6e-3 * L_0^2`, while the
continuity blocks carry `L_0^2` alone. Absolute SNES tolerances at
`1e-16` are needed so the Poisson block keeps iterating after the
continuity blocks have already reached their floors (and vice versa).
This is the same tolerance rationale as the Poisson MMS; we reuse it.

### 7.4 R_e near-equilibrium cancellation

In Variant C, `n_e * p_e - ni_hat^2 = ni_hat^2 * (exp(phi_p_e - phi_n_e)
- 1)`. For `A_n = -A_p = 0.3`, `phi_p_e - phi_n_e = -2*A_n * sin(pi*x/L)`
(1D) or `-2*A_n * sin(pi*x/L)*sin(pi*y/L)` (2D), so the worst-case
argument is `phi_p_e - phi_n_e \in [-0.6, 0.6]` and the numerator swings in
`[exp(-0.6)-1, exp(0.6)-1] = [-0.45, 0.82]` times `ni_hat^2 ~ 1e-12`.
The denominator is `~ 2*tau_hat * ni_hat ~ 2e-6` (at tau_hat ~ 1e-1),
so `R_e ~ 1e-6`. This is orders of magnitude below the diffusion term
in the continuity residual (`L_0^2 mu_n_hat * ni_hat * |grad phi_n|^2
~ 4e-12 * 1 * 1e-6 * (pi*A_n/L)^2`). Net: the SRH term in Variant C
contributes a genuine but sub-dominant coupling. This is desired for
gating -- it exercises the code path without dominating the error and
masking a Poisson/continuity bug.

If the measured rate on Variant C is lower than Variants A/B by more
than 0.1, that is a signal that the R_e path is introducing
discretization noise beyond the FE theory assumptions (likely through
the exp cancellation in the numerator); treat it as evidence of a bug
or poor kernel numerics, not as a legitimate rate degradation.

---

## Summary (for the gate reviewer)

- Exact solutions are sin-products that vanish on Omega's boundary;
  amplitudes `(A_psi, A_n, A_p) = (0.5, 0.3, -0.3)` keep exp(...) in
  `[0.45, 2.23]`. `A_p = -A_n` ensures `phi_p_e - phi_n_e` is non-zero
  pointwise so Variant C's SRH term is a genuine test of the code path.
- Scaling reuses `build_mms_scaling(L=2e-6)` from the Poisson MMS with
  the existing `lambda2 < 1` guard. No new constants.
- Three variants exercise Poisson-only, full coupling without R, and
  full coupling with SRH, each with closed-form symbolic forcings in
  weak form (Sections 3.1--3.4).
- BCs are homogeneous Dirichlet on all three fields, enforced exactly
  via `Constant(0.0)`.
- Expected rates are L^2 >= 1.75 and H^1 >= 0.80 per block, with
  monotonicity required across the sweep.
- Pytest uses 1D `[40, 80, 160]` and 2D `[16, 32, 64]`; CLI extends 1D
  to `[40, 80, 160, 320]` with an additional nonlinear-amplitude run.
- The main known risks are exp() blow-up (mitigated by amplitude cap)
  and drift-induced Galerkin degradation (shown to not apply at the
  chosen amplitudes by a Peclet-number argument in 7.2).

If all of the above is approved, Phase 4 implementation will create
`semi/verification/mms_dd.py` and `tests/fem/test_mms_dd.py` following
the exact structure of the Poisson MMS harness.
