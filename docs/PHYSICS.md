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

For Day 1 scope (equilibrium, complete ionization) we take
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

> **Cautionary note (Day 1 bug).** An earlier version of
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

**Why bias ramping is necessary.** At V = 0.6 V on the Day 2 pn
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

### 2.5 Scaled continuity (for Day 2)

Under the same scaling, with the time scale $t_0 = L_0^2 / D_0$ and
$D_0 = V_t \mu_0$ (Einstein for the reference mobility), the steady-state
electron continuity equation becomes

$$
-\nabla \cdot \bigl( \mu_n V_t\, n_i \exp(\hat\psi - \hat\Phi_n) \nabla \hat\Phi_n \bigr)
    = q\, R(n, p).
$$

The current-scale becomes $J_0 = q D_0 C_0 / L_0$. Details of the Day 2
scaled forms will be fixed when the code lands.

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

See `_build_ohmic_bcs` in `semi/run.py`.

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
interface (natural in the Galerkin form). Gate BCs are deferred to the
Day 5 MOS benchmark.

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
