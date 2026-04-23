# Physics reference

Stable reference for the governing equations, nondimensionalization, and
boundary conditions used throughout kronos-semi. This document changes
only when the physics content of the code changes. When adding a new
physics module, verify it uses the conventions here.

## 1. Governing equations (dimensional, SI)

Internal units are SI: potentials in V, lengths in m, densities in
$\mathrm{m}^{-3}$, current density in $\mathrm{A\,m^{-2}}$. JSON input
uses $\mathrm{cm}^{-3}$ for densities and $\mathrm{cm^2\,V^{-1}\,s^{-1}}$
for mobility; those are converted at ingest by helpers in
`semi/constants.py`.

### 1.1 Poisson equation

$$
-\nabla \cdot \bigl( \varepsilon_0 \varepsilon_r(\mathbf{x})\, \nabla \psi \bigr)
    = q\,\bigl( p - n + N_D^+ - N_A^- \bigr)
$$

For equilibrium problems (complete ionization, steady state) we take
$N_D^+ = N_D$ and $N_A^- = N_A$.

### 1.2 Carrier statistics (Boltzmann)

Under Boltzmann statistics, carrier densities are expressed in terms of
the electrostatic potential $\psi$ and the quasi-Fermi potentials
$\Phi_n$ (electrons) and $\Phi_p$ (holes):

$$
n = n_i \exp\!\left( \frac{\psi - \Phi_n}{V_t} \right), \qquad
p = n_i \exp\!\left( \frac{\Phi_p - \psi}{V_t} \right).
$$

At thermal equilibrium $\Phi_n = \Phi_p = 0$ and the product
$n p = n_i^2$ identically, which is the mass-action law verified in the
`pn_1d` benchmark.

### 1.3 Continuity equations (Slotboom current form)

Starting from the drift-diffusion currents
$\mathbf{J}_n = q \mu_n n \mathbf{E} + q D_n \nabla n$ (with Einstein
relation $D_n = \mu_n V_t$) and substituting the Boltzmann expression for
$n$, the electron current simplifies to

$$
\mathbf{J}_n = q\, \mu_n\, n_i \exp\!\left(\frac{\psi - \Phi_n}{V_t}\right)\, \nabla \Phi_n
    \;\equiv\; -q\, \mu_n\, n\, \nabla \Phi_n.
$$

That is, the current is a pure gradient of the electron quasi-Fermi
potential times an exponential coefficient. The hole current is
symmetric:

$$
\mathbf{J}_p = -q\, \mu_p\, p\, \nabla \Phi_p.
$$

The steady-state continuity equations are then

$$
\nabla \cdot \mathbf{J}_n = \phantom{-}q\, R(n, p), \qquad
\nabla \cdot \mathbf{J}_p = -q\, R(n, p),
$$

where $R$ is net recombination.

### 1.4 SRH recombination kernel

Shockley-Read-Hall recombination through a single mid-gap trap level:

$$
R_\mathrm{SRH}(n, p) = \frac{n p - n_i^2}{\tau_p (n + n_1) + \tau_n (p + p_1)}
$$

with $n_1 = n_i \exp(E_t / V_t)$, $p_1 = n_i \exp(-E_t / V_t)$, and
$E_t$ the trap energy referenced to the intrinsic level (zero for a
mid-gap trap).

## 2. Nondimensionalization

### 2.1 Scales

| Quantity                | Symbol        | Scale                                           | Scaled name    |
|-------------------------|---------------|-------------------------------------------------|----------------|
| Spatial coordinate      | $x$           | none (kept in meters)                           | $x$            |
| Electrostatic potential | $\psi$        | $V_0 = V_t = k_B T / q$                         | $\hat\psi$     |
| Quasi-Fermi potentials  | $\Phi_n, \Phi_p$ | $V_0 = V_t$                                  | $\hat\Phi_n$, $\hat\Phi_p$ |
| Carrier densities       | $n, p$        | $C_0 = \max\|N\|$ across profiles, floored at $10^{16}\,\mathrm{cm}^{-3}$ | $\hat n$, $\hat p$ |
| Doping                  | $N_D - N_A$   | $C_0$                                           | $\hat N$       |
| Intrinsic density       | $n_i$         | $C_0$                                           | $\hat n_i$     |

See `semi/scaling.py` for the implementation.

### 2.2 Scaled Poisson

Substituting $\psi = V_t \hat\psi$ and $n, p, N \to C_0 \hat n$, etc. into
the Poisson equation, **keeping the spatial coordinate in meters**:

$$
-\nabla \cdot \bigl( \varepsilon_0 \varepsilon_r\, V_t\, \nabla \hat\psi \bigr)
    = q\, C_0 \bigl( \hat p - \hat n + \hat N \bigr).
$$

Divide both sides by $q C_0$:

$$
-\nabla \cdot \bigl( L_D^2\, \varepsilon_r\, \nabla \hat\psi \bigr)
    = \hat p - \hat n + \hat N,
\qquad
\boxed{\; L_D^2 \;\equiv\; \frac{\varepsilon_0 V_t}{q C_0} \;}
$$

where $L_D$ is the extrinsic Debye length at the reference density
$C_0$. This is the coefficient used in `semi/physics/poisson.py`.

### 2.3 Relationship to `lambda2`

`Scaling.lambda2` returns the **dimensionless** ratio

$$
\lambda^2 = \frac{\varepsilon_0 V_t}{q C_0 L_0^2} \;=\; \frac{L_D^2}{L_0^2}.
$$

That is, $L_D^2 = \lambda^2 L_0^2$. `lambda2` is only the correct
stiffness coefficient if the spatial coordinate is itself scaled by
$L_0$.

> **Cautionary note (M1 coefficient bug).** An earlier version of
> `semi/physics/poisson.py` used `sc.lambda2` directly as the Laplacian
> coefficient. Because the mesh is in physical meters (not in units of
> $L_0$), this suppressed diffusion by a factor of $L_0^2$ (of order
> $10^{-12}$ for a micron-scale device) and Newton converged on the
> Laplace solution instead of the nonlinear Poisson one. V_bi came out
> right by coincidence (it is set entirely by the Dirichlet BCs), but
> peak $|E|$ was off by more than 20x. The fix was to use
> $L_D^2 = \lambda^2 L_0^2$. Any future scaled forms must follow the
> same convention.

### 2.4 Bias continuation strategy

The coupled (psi, phi_n, phi_p) system is stiff at applied bias. The
Slotboom carrier expression $n = n_i \exp(\hat\psi - \hat\Phi_n)$ can
underflow in deep depletion and overflow in heavy accumulation (see
`docs/adr/0004-slotboom-variables-for-dd.md`). Newton applied directly
to a high-bias problem from a zero initial guess does not converge.
The driver (`semi.run.run_bias_sweep`) solves this by ramping the
applied bias from equilibrium through a sequence of intermediate
targets, using the previous converged solution as the initial guess at
each step.

**Coordinate convention reminder.** The mesh stays in physical meters,
so the scaled Poisson stiffness coefficient is $L_D^2 \varepsilon_r =
\lambda^2 L_0^2 \varepsilon_r$, not $\lambda^2 \varepsilon_r$. Any new
scaled form added to the DD block must follow the same convention (see
Section 2.3).

**Adaptive step control.** The continuation stepper is a signed
step-size controller (`semi.continuation.AdaptiveStepController`) with
three knobs:

- On SNES convergence with iteration count strictly below
  `easy_iter_threshold` (default 4), an easy-solve counter is
  incremented. Once the counter reaches the threshold, the step is
  multiplied by `grow_factor` (default 1.5) and the counter resets.
- On SNES divergence (including line-search failure), the step is
  halved and the counter resets. Before halving, the block unknowns
  are restored to the last converged state.
- The step is bounded above by `solver.continuation.max_step`
  (default: the voltage_sweep spacing) and below by
  `solver.continuation.min_step` (default 1e-4 V). Halving below
  `min_step` aborts the run; growing above `max_step` clamps.

The controller is signed: negative `initial_step` drives a reverse
sweep. `clamp_to_endpoint` ensures the final step exactly lands on
the sweep endpoint instead of overshooting.

**Why bias ramping is necessary.** At V = 0.6 V on the M2 pn
junction, the peak electric field is roughly $2 \times 10^5$ V/cm and
the minority carrier density at the depletion edge jumps by a factor
$\exp(0.6/V_t) \approx 10^{10}$ above its equilibrium value. A Newton
step from an equilibrium initial guess to this state crosses the
domain where carrier densities underflow to zero, the Jacobian loses
rank, and line search fails to find an improvement. Marching in
increments of $V_t$ (or smaller) keeps the updates within the radius
of quadratic convergence.

**When to override defaults.**

- If a sweep halts on "Bias ramp failed after N halvings," first
  increase `continuation.max_halvings` (default 6). If halving to
  `min_step` still fails, the regime is genuinely outside what
  constant-mobility Boltzmann DD can resolve (high injection,
  breakdown); tighten the sweep range rather than forcing the solver.
- If the sweep succeeds but spends too many Newton iterations per
  step, raise `continuation.max_step` above the sweep spacing to let
  growth skip over redundant record points.
- For deep reverse bias, lower `easy_iter_threshold` or
  `grow_factor` to keep the step conservative; the Slotboom
  formulation becomes stiff as minority carrier densities underflow
  and larger jumps lose Newton convergence.

### 2.5 Scaled drift-diffusion

This section completes the nondimensionalization that was a placeholder
through M2-M4. The implementation is now in place
(`semi/physics/drift_diffusion.py`) and the V&V suite
(`semi/verification/mms_dd.py`) confirms theoretical L^2 = 2 / H^1 = 1
rates on every block, so we can fix the scaled forms here without risk
of revision.

#### Dimensional starting point

From Section 1, the steady-state Slotboom drift-diffusion system on a
semiconductor region is

$$
\begin{aligned}
-\nabla \cdot \bigl( \varepsilon_0 \varepsilon_r \nabla \psi \bigr)
    &= q\, ( p - n + N_D - N_A ), \\
+\nabla \cdot \bigl( q\, \mu_n\, n\, \nabla \Phi_n \bigr)
    &= +\, q\, R(n, p), \\
+\nabla \cdot \bigl( q\, \mu_p\, p\, \nabla \Phi_p \bigr)
    &= -\, q\, R(n, p),
\end{aligned}
$$

with $n = n_i \exp((\psi - \Phi_n)/V_t)$ and $p = n_i
\exp((\Phi_p - \psi)/V_t)$ from Boltzmann statistics. The continuity
sign convention follows from $\mathbf{J}_n = -q\, \mu_n\, n\, \nabla
\Phi_n$ and $\nabla \cdot \mathbf{J}_n = +q\, R$ (Section 1.3).

#### Scales

Reusing the Section 2.1 conventions and adding two derived scales for
the continuity rows:

| Quantity        | Scale                         | Definition                         |
|-----------------|-------------------------------|------------------------------------|
| Length          | none (mesh stays in meters)   | (Invariant 3)                      |
| Potential       | $V_0 = V_t = k_B T / q$       |                                    |
| Density         | $C_0 = \max\|N\|$ (floored)   |                                    |
| Mobility        | $\mu_0 = \mu_n^\mathrm{ref}$  | reference material electron mobility |
| Diffusivity     | $D_0 = V_t \mu_0$             | Einstein at the reference          |
| Time            | $t_0 = L_0^2 / D_0$           | sets the recombination scale       |
| Current density | $J_0 = q\, D_0\, C_0 / L_0$   |                                    |

The substitutions are $\psi = V_t \hat\psi$, $\Phi_{n,p} = V_t
\hat\Phi_{n,p}$, $n = C_0 \hat n$, $p = C_0 \hat p$, $N_{D,A} = C_0
\hat N_{D,A}$, $\mu_{n,p} = \mu_0 \hat\mu_{n,p}$, $\tau_{n,p} = t_0
\hat\tau_{n,p}$, $n_i = C_0 \hat n_i$. The mesh coordinate $x$ is
left in physical meters per Invariant 3, which is why $L_0$ shows up
explicitly in the coefficients below rather than being absorbed.

#### Scaled Poisson row

Substituting and dividing by $q C_0$ exactly as in Section 2.2:

$$
-\nabla \cdot \bigl( L_D^2\, \varepsilon_r\, \nabla \hat\psi \bigr)
    = \hat p - \hat n + \hat N,
\qquad L_D^2 = \frac{\varepsilon_0 V_t}{q C_0}.
$$

This is the form that lands in `build_dd_block_residual` for the psi
row. The coefficient is `L_D^2 * eps_r`, *not* `lambda2 * eps_r` (see
Invariant 3 and the cautionary note in Section 2.3 below).

#### Scaled continuity rows

Starting from the electron continuity equation,

$$
\nabla \cdot \bigl( q\, \mu_n\, n\, \nabla \Phi_n \bigr) = q\, R.
$$

Substituting $\Phi_n = V_t \hat\Phi_n$, $n = C_0 \hat n$, $\mu_n =
\mu_0 \hat\mu_n$ leaves the gradient operator on the left in physical
meters (Invariant 3) and yields

$$
\nabla \cdot \bigl( q\, \mu_0 V_t\, \hat\mu_n\, C_0\, \hat n\, \nabla \hat\Phi_n \bigr)
    = q\, R.
$$

Pull the constants out and divide by $q\, C_0 / t_0 = q C_0 D_0 /
L_0^2$ (the natural rate scale for $R$, since lifetimes scale as
$\tau / t_0$):

$$
\nabla \cdot \bigl( \tfrac{\mu_0 V_t}{D_0/L_0^2}\, \hat\mu_n\, \hat n\,
                    \nabla \hat\Phi_n \bigr)
    = \hat R, \qquad
\hat R \equiv R \cdot \tfrac{t_0}{C_0}.
$$

The coefficient simplifies via $D_0 = V_t \mu_0$:

$$
\frac{\mu_0 V_t}{D_0 / L_0^2}
    = \frac{\mu_0 V_t L_0^2}{V_t \mu_0}
    = L_0^2.
$$

So the scaled electron continuity row is

$$
\nabla \cdot \bigl( L_0^2\, \hat\mu_n\, \hat n\, \nabla \hat\Phi_n \bigr)
    = \hat R,
$$

and by symmetry the hole row is

$$
\nabla \cdot \bigl( L_0^2\, \hat\mu_p\, \hat p\, \nabla \hat\Phi_p \bigr)
    = -\, \hat R.
$$

The scaled SRH kernel (from Section 1.4 with the same C0/t0 division) is

$$
\hat R(\hat n, \hat p)
    = \frac{\hat n \hat p - \hat n_i^2}
           {\hat\tau_p (\hat n + \hat n_1) + \hat\tau_n (\hat p + \hat p_1)},
$$

with $\hat n_1 = \hat n_i \exp(E_t / V_t)$ and $\hat p_1 = \hat n_i
\exp(-E_t / V_t)$. This is exactly the inline expression in
`build_dd_block_residual` (the standalone UFL helper
`semi.physics.recombination.srh_rate` returns the same expression).

The matching residual signs in `build_dd_block_residual` are

```
-div( L_D^2 eps_r grad psi_hat ) - ( p_hat - n_hat + N_hat )         (psi row)
-div( L_0^2 mu_n_hat n_hat grad phi_n_hat ) + R_hat                   (phi_n row)
-div( L_0^2 mu_p_hat p_hat grad phi_p_hat ) - R_hat                   (phi_p row)
```

(Each row is rewritten as $F = 0$ for SNES.)

#### Current scale

The terminal-current evaluation in `semi.postprocess.evaluate_current_at_contact`
assembles $\mathbf{J}_n = q\, \mu_n\, n\, \nabla \Phi_n$ and
$\mathbf{J}_p = q\, \mu_p\, p\, \nabla \Phi_p$ in physical units
(no scaling), since the test is what amperes per square meter the
device produces. The natural scale is

$$
J_0 = q\, D_0\, C_0 / L_0,
$$

which for Si at $C_0 = 10^{17}\,\mathrm{cm}^{-3}$, $L_0 = 2\,\mu m$,
and $\mu_0 = 1400\,\mathrm{cm^2/(V\,s)}$ comes out around
$2.9 \times 10^4\,\mathrm{A/cm^2}$. The M2 `pn_1d_bias` benchmark
reports $J(V = 0.6\,\mathrm{V}) \approx 1.6 \times 10^3\,\mathrm{A/m^2}$,
which is several decades below $J_0$ as expected for moderate forward
bias.

#### The lambda2 / L_D^2 trap

`Scaling.lambda2` returns the dimensionless ratio $\lambda^2 = L_D^2 /
L_0^2$ (Section 2.3). Early M1 versions of `semi/physics/poisson.py`
used `sc.lambda2` directly as the Laplacian coefficient. Because the
mesh is in physical meters (Invariant 3, not $L_0$ units), this
suppressed diffusion by $L_0^2$ ($\sim 10^{-12}$ on a micron device)
and Newton converged on the Laplace solution. The fix was to use
$L_D^2 \cdot \varepsilon_r = \lambda^2 \cdot L_0^2 \cdot \varepsilon_r$.
This caveat is preserved as a historical note because the same trap
applies to any new scaled form added to the DD block: while the mesh
sits in meters, every $\nabla$ contributes a $1 / L_0$ in scaled units,
so each `div(... grad ...)` operator picks up an explicit $L_0^2$ (for
the continuity rows) or $L_D^2 = \lambda^2 L_0^2$ (for the Poisson row).

See ADR 0002 (`docs/adr/0002-nondimensionalization-mandatory.md`) for
the rationale of the scaling choice and ADR 0004
(`docs/adr/0004-slotboom-variables-for-dd.md`) for why the DD rows are
written in Slotboom form rather than mixed (n, p, psi) primitive
variables.

For the MMS-specific construction of the manufactured forcing terms
on the same scaled equations (variants A / B / C, exact solutions
$\hat\psi_e, \hat\Phi_{n,e}, \hat\Phi_{p,e}$, and the SNES tolerance
note), see `docs/mms_dd_derivation.md`.

## 3. Boundary conditions

### 3.1 Ohmic contacts

At an ohmic contact the carrier densities equal their majority-side
equilibrium values and the semiconductor is locally in quasi-neutral
charge balance. In terms of $\psi$, this gives the Dirichlet value

$$
\psi = V_t\, \operatorname{asinh}\!\left( \frac{N_D - N_A}{2 n_i} \right) \;+\; V_\mathrm{applied}.
$$

For the coupled system, the quasi-Fermi potentials at the contact equal
the applied bias:

$$
\Phi_n = \Phi_p = V_\mathrm{applied}.
$$

See `build_psi_dirichlet_bcs` and `build_dd_dirichlet_bcs` in `semi/bcs.py`.

### 3.2 Gate contacts

At a gate over an oxide, the applied potential is specified as a
Dirichlet condition on $\psi$, offset by the metal-semiconductor work
function difference $\phi_{ms}$ (ignored for ideal gates):

$$
\psi_\mathrm{gate} = V_g - \phi_{ms}.
$$

No carriers are defined in the oxide; continuity equations are solved
only in semiconductor regions. The normal component of
$\varepsilon \nabla \psi$ is continuous across the oxide/semiconductor
interface (natural in the Galerkin form). Gate BCs are implemented for
the 2D MOS capacitor benchmark; see Section 6.

### 3.3 Insulating boundaries

Homogeneous Neumann on every field:

$$
\nabla \psi \cdot \hat{\mathbf{n}} = 0, \qquad
\nabla \Phi_n \cdot \hat{\mathbf{n}} = 0, \qquad
\nabla \Phi_p \cdot \hat{\mathbf{n}} = 0.
$$

This is the natural condition and requires no explicit code.

## 4. Reference values

All values are room temperature ($T = 300\,\mathrm{K}$) unless noted.
Source of truth for material parameters: `semi/materials.py`.

### 4.1 Universal

| Symbol     | Value                                          | Notes                                  |
|------------|------------------------------------------------|----------------------------------------|
| $q$        | $1.602176634 \times 10^{-19}$ C                | 2019 SI redefinition exact             |
| $k_B$      | $1.380649 \times 10^{-23}$ J/K                 | 2019 SI redefinition exact             |
| $\varepsilon_0$ | $8.8541878128 \times 10^{-12}$ F/m        |                                        |
| $V_t(300\,\mathrm{K})$ | $25.852\,\mathrm{mV}$              | $= k_B T / q$                          |

### 4.2 Silicon

| Symbol        | Value                               | Notes                                               |
|---------------|-------------------------------------|-----------------------------------------------------|
| $\varepsilon_r$ | $11.7$                             |                                                     |
| $E_g$         | $1.12$ eV                           |                                                     |
| $\chi$        | $4.05$ eV                           | Electron affinity                                   |
| $N_c$         | $2.86 \times 10^{19}\,\mathrm{cm}^{-3}$ | Sze, 3rd ed. Table 7                            |
| $N_v$         | $3.10 \times 10^{19}\,\mathrm{cm}^{-3}$ |                                                 |
| $n_i$         | $1.0 \times 10^{10}\,\mathrm{cm}^{-3}$ | Altermatt 2003 consensus                         |
| $\mu_n$       | $1400\,\mathrm{cm^2\,V^{-1}\,s^{-1}}$ | Undoped, 300 K                                   |
| $\mu_p$       | $450\,\mathrm{cm^2\,V^{-1}\,s^{-1}}$  |                                                  |

> **Note on $n_i$.** The older textbook value $1.45 \times 10^{10}\,\mathrm{cm}^{-3}$
> is superseded; modern TCAD tools (Sentaurus, Atlas, etc.) use values close
> to $1.0 \times 10^{10}$ following Altermatt's 2003 reassessment. We match
> the modern convention. If you compare against a hand calculation that
> uses 1.45e10, expect built-in voltage offsets of order $V_t \ln(1.45^2)
> \approx 19\,\mathrm{mV}$.

### 4.3 Other semiconductors

See `semi/materials.py` for Ge, GaAs. Insulators (SiO2 $\varepsilon_r = 3.9$,
HfO2 $\varepsilon_r = 25.0$, Si3N4 $\varepsilon_r = 7.5$) carry only a relative
permittivity; semiconductor-specific fields are zero on them.

## 5. Verification & Validation

This section documents the M4 V&V suite and its results. The code
lives under `semi/verification/`, the runner is
`scripts/run_verification.py`, artifacts land under
`results/mms_poisson/`, `results/mesh_convergence/`,
`results/conservation/`, and `results/mms_dd/`. CI runs
`python scripts/run_verification.py all` inside the `docker-fem` job on
every push and uploads the full `results/` tree.

The distinction followed here is the Roache / Oberkampf-Roy one.
**Verification** asks "did we solve the equations we wrote down
correctly?"; it is a property of the code and its discretization.
**Validation** asks "are those equations the right model for the
device?"; it compares the code to experiment or to an accepted
reference model. The four activities below are all verification
activities. The earlier single-point benchmark verifiers (`pn_1d`,
`pn_1d_bias`, `pn_1d_bias_reverse`) stand as validation against
depletion-approximation and Shockley / SNS analytics; they remain in
CI unchanged.

### 5.1 MMS for equilibrium Poisson (Phase 1)

Module: `semi/verification/mms_poisson.py`. Subcommand:
`run_verification.py mms_poisson`.

Manufactured smooth exact solutions (1D sin, 2D sin-sin with distinct
wavenumbers) are injected into the equilibrium Poisson form with UFL
weak-form forcing. The mesh is refined over several levels; the
observed L^2 and H^1 rates are measured on each refinement pair. The
finest-pair rate is gated at L^2 >= 1.85 and H^1 >= 0.85, with
monotone error reduction required across the sweep.

Finest-pair rates on the latest run (N = 320 for 1D, N = 128 for 2D):

| Study          | L^2 rate | H^1 rate |
|----------------|----------|----------|
| 1D linear      | 2.000    | 1.000    |
| 1D nonlinear   | 2.000    | 1.000    |
| 2D triangles   | 1.998    | 0.999    |

The `2d_quad_smoke` sanity check confirms the quad-element L^2 error
is within an order of magnitude of the triangle error at N = 64
(ratio ~0.39, well inside the [0.1, 10] acceptance band).

### 5.2 Mesh convergence on `pn_1d` (Phase 2)

Module: `semi/verification/mesh_convergence.py`. Subcommand:
`run_verification.py mesh_convergence`.

A mesh sweep over N in [50, 100, 200, 400, 800, 1600] on the M1
equilibrium pn junction reports V_bi, peak |E|, depletion width W,
Newton iterations, and solve time. Errors are recorded in two ways:
against the depletion-approximation reference (which is itself an
approximation), and as Cauchy self-convergence ratios (consecutive
mesh differences, which isolate pure FE discretization error).

**Honest flag.** The depletion approximation is a model, not the
truth. The FEM solver converges to the full Poisson-Boltzmann
solution, so relative errors against the depletion references
plateau on fine meshes at the physics-model gap, not at zero.
Monotone reduction is therefore gated only over the first four
refinements (N up to 400), before the plateau; the Cauchy ratios
are gated across the same range at >= 1.8x per doubling, which is
the mathematically meaningful convergence indicator. V_bi is set
exactly by the Ohmic BCs and is mesh-independent to machine
epsilon (reported only, not gated).

Latest run, consecutive Cauchy ratios (N = 100/50, 200/100, 400/200):

| Quantity  | Ratio 100/50 | Ratio 200/100 | Ratio 400/200 |
|-----------|--------------|---------------|---------------|
| E_peak    | 1.99         | 2.10          | 2.31          |
| W         | 2.01         | 4.10          | 4.21          |

All >= 1.8x as required. Residual error vs. depletion at the finest
level (N = 1600): E_peak 4.3% gap, W 3.2% gap, consistent with the
expected physics-model offset.

### 5.3 Conservation checks (Phase 3)

Module: `semi/verification/conservation.py`. Subcommand:
`run_verification.py conservation`.

- **Charge conservation** on `pn_1d` equilibrium. Integrates
  q * (p - n + N_D - N_A) over the device and asserts
  |Q_net| < 1e-10 * Q_ref with
  Q_ref = q * max|N_net| * L_device. Latest run:
  Q_net = 4.8e-19 C/m^2 vs. threshold 3.2e-9 C/m^2 (rel 1.5e-17).

- **Current continuity** on `pn_1d_bias` forward and
  `pn_1d_bias_reverse`. At each target bias the total current
  J_total = J_n + J_p is sampled on ten interior facets; we assert
  max |J - mean| / |mean| < 5% (forward) or 15% (reverse). Latest
  worst-case max_rel across the target set:

  | Benchmark            | V         | worst max_rel | tol    |
  |----------------------|-----------|---------------|--------|
  | pn_1d_bias           | 0.6 V     | 1.91%         | 5%     |
  | pn_1d_bias_reverse   | -0.55 V   | 0.020%        | 15%    |

  All inside tolerance.

### 5.4 MMS for coupled drift-diffusion (Phase 4)

Module: `semi/verification/mms_dd.py`. Subcommand:
`run_verification.py mms_dd`. Derivation and rate-threshold
rationale: `docs/mms_dd_derivation.md`.

Nine studies (three variants x {1D default, 1D nonlinear, 2D}) on
the full three-block Slotboom residual with forcing terms injected
in weak form:

- **Variant A (psi-only).** `phi_n_e = phi_p_e = 0`. The
  continuity diffusion integrals drop out and R_e collapses to
  zero. Exercises the Poisson row of the block residual
  (`build_dd_block_residual`) including sign, scaling, and
  coupling-to-zero-quasi-Fermi plumbing. Variant A's continuity-block
  errors are at machine noise (~1e-35 to 1e-41); these are
  reported but not rate-gated because dividing one noise level
  by another produces meaningless rates. We gate only the psi
  block on Variant A (L^2 >= 1.75, H^1 >= 0.80).
- **Variant B (full coupling, no R).** All three fields nontrivial,
  lifetimes set to 1e+20 so R_e is negligible but the SRH kernel
  is still evaluated. Verifies the drift-diffusion operators in
  their coupled form, independent of recombination numerics.
- **Variant C (full coupling with SRH).** Realistic lifetimes
  tau_n = tau_p = 1e-7 s. Verifies the full residual including
  the cross-block coupling through R_e.

Finest-pair rates (N = 320 for 1D, N = 64 for 2D), default
amplitudes (A_psi, A_n, A_p) = (0.5, 0.3, -0.3):

| Study                | Block  | L^2 rate | H^1 rate |
|----------------------|--------|----------|----------|
| 1D A (linear)        | psi    | 2.000    | 1.000    |
| 1D B (linear)        | psi    | 2.000    | 1.000    |
| 1D B (linear)        | phi_n  | 2.000    | 1.000    |
| 1D B (linear)        | phi_p  | 2.000    | 1.000    |
| 1D C (linear)        | psi    | 2.000    | 1.000    |
| 1D C (linear)        | phi_n  | 2.000    | 1.000    |
| 1D C (linear)        | phi_p  | 2.000    | 1.000    |
| 2D A                 | psi    | 1.990    | 0.996    |
| 2D B                 | psi    | 1.990    | 0.996    |
| 2D B                 | phi_n  | 1.995    | 1.000    |
| 2D B                 | phi_p  | 1.995    | 1.000    |
| 2D C                 | psi    | 1.990    | 0.996    |
| 2D C                 | phi_n  | 1.995    | 1.000    |
| 2D C                 | phi_p  | 1.995    | 1.000    |

Every gated rate clears the L^2 >= 1.75 / H^1 >= 0.80 floor with
>= 0.19 of headroom, matching theoretical P1 Lagrange rates to within
roundoff. The derivation (`mms_dd_derivation.md`) documents the
block-residual scale-disparity issue that forced the SNES
tolerance tweak to `atol = 0.0` with `stol = 1e-12`.

### 5.5 MMS for multi-region Poisson (M6: 2D MOS capacitor)

Guards the Si/SiO2 coefficient-jump assembly used by the MOS
capacitor. A manufactured psi_exact(x, y) is C^0 across the interface
y = y_int and satisfies eps-weighted flux continuity there (no
surface delta source); per-region forcing is supplied via a cellwise
DG0 eps_r Function. The derivation is in `mos_derivation.md` section
7. Convergence sweep on the Si/SiO2 rectangle (t_Si = 70 nm, t_ox =
30 nm, oxide inflated from the device 5 nm to keep a reasonable
aspect ratio on a square uniform mesh):

    N in [20, 40, 80, 160]     finest-pair rate_L^2 = 2.000, rate_H^1 = 1.000

The CLI (`scripts/run_verification.py mms_poisson`) enforces
`rate_L^2 >= 1.99` and `rate_H^1 >= 0.95` on the finest pair; the
pytest gate uses the looser `>= 1.85 / >= 0.85` thresholds common to
the single-region MMS tests.

## 6. 2D MOS capacitor (M6: 2D MOS capacitor)

Condensed reference; the full derivation is in
`docs/mos_derivation.md`.

### 6.1 Device structure

- Silicon substrate, uniform N_A = 1e17 cm^-3, 500 nm thick
- SiO2 gate oxide, 5 nm thick
- Body contact on the bottom face (ohmic)
- Gate contact on the top face (ideal gate, phi_ms = 0)

The benchmark mesh is uniform at 1 nm per cell vertically (505 cells
over 505 nm), which lands the Si/SiO2 interface on a grid line and
gives 5 cells in the oxide. Lateral extent is 500 nm with 4 cells
(the device is translation invariant in x away from the lateral
edges). Either a dense uniform or a graded vertical mesh is
acceptable; the uniform path is used here because CI time is not
the bottleneck at this resolution (the full sweep runs in ~2 s).

### 6.2 Physics model

Equilibrium Poisson with multi-region coefficient, solved once per
gate bias:

- Stiffness L_D^2 eps_r(x) grad psi . grad v over the full mesh with
  cellwise DG0 eps_r (Si: 11.7, SiO2: 3.9)
- Space-charge rho_hat = n_i_hat (exp(-psi) - exp(psi)) + N_net_hat
  restricted to silicon cells via `dx(subdomain_id=semi_tag)`
- Oxide: Laplacian only, no space charge, no Slotboom variables

The body ohmic BC sets psi_body = -phi_F (equilibrium p-type value
under our psi=0-at-intrinsic convention); the gate BC is
psi_gate = V_gate - phi_ms.

Gate and body DOFs cannot overlap: gate contacts live on the oxide
side of the interface where Slotboom variables do not exist, so
`build_dd_dirichlet_bcs` explicitly skips them.

### 6.3 Flatband, threshold, and the BC-convention shift

With our BC convention a p-type substrate at equilibrium has
psi_body = -phi_F = -V_t ln(N_A / n_i). Flatband (psi flat through
silicon) therefore requires psi_gate = psi_body, giving

    V_FB = phi_ms - phi_F

not the textbook V_FB = phi_ms. For the ideal-gate M6 device
(phi_ms = 0, N_A = 1e17 cm^-3, T = 300 K), V_FB = -0.417 V and
V_T = -0.417 + 2 * 0.417 + sqrt(4 eps_s q N_A phi_F) / C_ox = +0.658 V.
The benchmark sweep covers V_gate in [-0.9, +1.2] V (roughly V_FB -
0.5 through V_T + 0.5) and the verifier window sits strictly inside
the depletion regime at [V_FB + 0.2, V_T - 0.1] = [-0.217, +0.558] V.

The 0.2 V low-edge margin is wider than the nominal 0.1 V because at
psi_s < ~2 V_t the carrier tail reaches across a significant fraction
of W_dep and the depletion approximation naturally drifts toward 10%.
Shrinking the window (not loosening the tolerance) is the intended
knob per `mos_derivation.md` section 6.9.

### 6.4 Capacitance extraction and theory

At each V_gate the integrated silicon space charge gives

    Q_gate(V_gate) = -(q / W_lat) * integral_{Omega_Si} rho(x, y) dA

(2D: dividing the charge-per-unit-depth by the lateral extent W_lat
converts to per-area). The simulated capacitance is the centered
finite difference dQ_gate/dV_gate on the sweep grid. Theory
(depletion approximation) is the series combination

    1/C = 1/C_ox + 1/C_dep(psi_s)
    C_ox  = eps_ox / t_ox
    C_dep = sqrt(eps_s q N_A / (2 psi_s))

with psi_s(V_gate) obtained by closed-form inversion of

    V_gate - V_FB = psi_s + a sqrt(psi_s),   a = sqrt(2 eps_s q N_A) / C_ox

(substitute u = sqrt(psi_s) and solve the quadratic
u^2 + a u - V_ov = 0).

### 6.5 Verifier result

Inside the verifier window [V_FB + 0.2, V_T - 0.1] V the worst
relative error |C_sim - C_theory| / C_theory is ~9.3% (at V_gate =
V_FB + 0.217 V, the window edge), inside the 10% tolerance. The
benchmark artifacts are `results/mos_2d/{psi_2d, potentials_1d, cv,
qv}.png`. Monotone non-increasing C in the window is also checked
(C starts near C_ox in accumulation, drops through depletion, rises
back toward C_ox in inversion; inside the window it is strictly
decreasing up to ~1% noise).

## 7. 3D doped resistor (M7: 3D doped resistor)

M7: 3D doped resistor is a dimension extension, not a physics extension. Both the
equilibrium Poisson form (Section 1.1, 2.2) and the coupled Slotboom
drift-diffusion block (Section 1.3, 2.5) are written against the
abstract `nabla` operator and the cell measure `dx`, and the existing
scaled coefficients `L_D^2 eps_r` (Poisson) and `L_0^2 mu_hat`
(continuity) pick up no new dimensional factors in 3D: every
`nabla` still contributes one `1 / L_0` and every `dx` still
contributes `L_0^d`, with the `d`-dependence cancelling once the
weak form is divided through by the reference density rate (see
Section 2.5). In code, `semi/physics/poisson.py` and
`semi/physics/drift_diffusion.py` are untouched on this branch; the
only dimension-sensitive site is the mesh builder, which already
dispatches to `create_interval / create_rectangle / create_box` by
`cfg["dimension"]`, and the builtin box-tagger, which already
indexes centroids over `range(tdim)`.

The full derivation (device geometry, analytical ohmic resistance,
V-I linearity metric and 1% tolerance rationale, 3D slice plot
strategy, gmsh loader test strategy) lives in
`docs/resistor_derivation.md`.

### 7.1 Device

Rectangular silicon bar aligned along the x axis: `L = 1 um`,
cross-section `W x W = 200 nm x 200 nm`, uniform n-type
`N_D = 1e18 cm^-3`. Two ohmic contacts on the `x = 0` and `x = L`
faces; the other four faces are insulating (natural BC). Minority
holes at equilibrium are `n_i^2 / N_D = 1e2 cm^-3`, sixteen orders
of magnitude below N_D, so the device is safely in the ohmic
majority-carrier regime.

### 7.2 Analytical resistance

For a uniform bar with constant mobility in the low-field
low-injection limit,

```
R = L / (q N_D mu_n A),   I = V / R
```

Using `mu_n = 1400 cm^2/(V s)` from `semi/materials.py` and
`A = 4e-14 m^2`:

```
R_theory = 1.0e-6 / (1.602e-19 * 1.0e24 * 0.14 * 4e-14) = 1115 Ohm
```

This is exact in the `V -> 0, uniform-everything` limit; the 1%
verifier tolerance reflects that fact, in contrast to the 10-20%
tolerances on pn-junction and MOS verifiers where the reference
itself carries depletion-approximation modeling error.

### 7.3 Bipolar sweep support

The resistor sweep is symmetric about V = 0
(`V in {-0.010, -0.005, 0.000, +0.005, +0.010} V`). The
pn-junction / MOS benchmarks were always unipolar, so the bias
runner (`semi/runners/bias_sweep.py`) grew a two-leg walk for this
benchmark: when the requested sweep has both negative and positive
entries, it walks `V = 0 -> min(V) -> max(V)` with a fresh
`AdaptiveStepController` on each leg. This is a driver-level change
only; the physics residual, current extraction, and BC code paths
are untouched. A unit test in `tests/test_bipolar_sweep.py` pins
the leg computation by constructing a sweep list crossing zero and
asserting `bipolar_legs` is populated.

### 7.4 V-I linearity verifier

Implemented as `verify_resistor_3d` in `scripts/run_benchmark.py`:

1. Sweep the right contact across 5 equally spaced points in
   `[-0.01, +0.01] V`.
2. Extract `I(V)` as the UFL facet integral of `J_n . n` over the
   `x = L` contact (same postprocessing path as `pn_1d_bias`,
   exercised here on a 3D facet set for the first time).
3. Compute `R_sim(V) = V / I(V)` at each nonzero V.
4. Compare against `R_theory` computed from `mu_n` read from
   `semi/materials.py` and `N_D` read from the benchmark JSON.
5. Assert `max |R_sim - R_theory| / R_theory < 1%`.

Sanity checks: `|I(V=0)|` must be within numerical noise
(`< R_theory * V_t * 1e-6`) and `sign(I(V)) == sign(V)` for every
nonzero V. Either violation indicates a facet-normal orientation
bug or a swapped ohmic BC, not a physics correction.

### 7.5 Gmsh loader

Unstructured tetrahedral meshes are loaded through
`dolfinx.io.gmsh.read_from_msh` in
`semi/mesh.py::_build_from_file`. Physical groups stored in the
`.msh` file are returned verbatim as `cell_tags` and `facet_tags`;
the builtin JSON box-tagger is bypassed for file-source meshes
because the `.msh` tagging is authoritative. A committed fixture
`benchmarks/resistor_3d/fixtures/box.msh` (with its reproducible
`box.geo` source) exercises the same resistor geometry on an
unstructured mesh; the builtin and gmsh variants must each clear
the 1% V-I tolerance and must agree with each other by
transitivity.

No changes to the equilibrium or drift-diffusion residuals were
required for the 3D path. MMS convergence in 3D is deferred: the
1D and 2D MMS suites already cover the dimension-agnostic assembly
(Section 5.1-5.5), and the 3D ohmic-resistor theory match against a
closed-form reference on two distinct mesh paths (builtin Cartesian
tetrahedralization + gmsh unstructured tetrahedra) is sufficient
evidence that the 3D code path and the gmsh loader are correct.
