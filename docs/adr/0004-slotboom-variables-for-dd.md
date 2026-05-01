# 0004. Slotboom variables for drift-diffusion

- Status: Accepted (carrier-density mapping subsection partially superseded by ADR 0015 for Fermi-Dirac regime)
- Date: 2026-04-20

## Context

The drift-diffusion currents

$$
\mathbf{J}_n = q \mu_n n \mathbf{E} + q D_n \nabla n,
\qquad
\mathbf{J}_p = q \mu_p p \mathbf{E} - q D_p \nabla p
$$

have a mixed advective-diffusive structure. The local Peclet number is
$\mathrm{Pe} = |\mathbf{E}| L_\mathrm{element} / (2 V_t)$. In any
realistic semiconductor device under bias, $|\mathbf{E}|$ at the
junction or in the channel is in the range $10^4$ to $10^6$ V/cm, which
gives $\mathrm{Pe} \gg 1$ for any mesh resolution short of atomistic.

Standard Galerkin FEM is known to be unstable (oscillatory) in the
drift-dominant regime. The textbook fixes are:

- **SUPG / streamline diffusion:** add stabilization terms tuned to
  $\mathrm{Pe}$. Requires careful parameter selection; mass
  conservation becomes approximate; coupling to Poisson is more
  complex.
- **Scharfetter-Gummel box integration:** exact exponential fitting on
  a dual mesh. Excellent for finite volumes, but not a natural fit for
  high-order FEM spaces and does not generalize cleanly to 3D
  unstructured meshes.
- **Slotboom variables (quasi-Fermi potentials):** rewrite $n$ and $p$
  in terms of $\Phi_n$ and $\Phi_p$ using Boltzmann relations. The
  resulting current is a pure gradient,
  $\mathbf{J}_n = -q \mu_n n \nabla \Phi_n$, with a positive
  coefficient. Standard Galerkin on $\Phi_n$ is coercive and stable
  without any stabilization terms.

## Decision

Use Slotboom variables $\Phi_n, \Phi_p$ as the primary continuity
unknowns. The coupled block unknown is $(\hat\psi, \hat\Phi_n,
\hat\Phi_p)$. Electron and hole densities are recovered pointwise from
the Boltzmann relations in `docs/PHYSICS.md` Section 1.2.

## Consequences

Easier:

- No SUPG parameters to tune; no stabilization terms to maintain across
  elements or dimensions.
- Conservation is exact (the continuity equations are in divergence
  form of a pure gradient).
- Ohmic and gate BCs on $\Phi_n$ and $\Phi_p$ are natural (they equal
  the applied bias at contacts).

Harder:

- The coefficient $n = n_i \exp((\hat\psi - \hat\Phi_n)/V_t)$ can
  underflow in deep depletion (large negative argument), or overflow in
  heavy accumulation. The coupled system is stiff. Bias ramping
  continuation (M3: Adaptive continuation) is mandatory; jumping from zero to high forward
  bias will not converge.
- Reverse bias near breakdown pushes the exponent large negative; we
  rely on PETSc SNES line search and bias ramping to stay in the
  feasible region.
- Users familiar with the raw $(n, p)$ formulation have to translate
  when reading the code.

## Related

- `docs/PHYSICS.md` Sections 1.2 and 1.3.
- Invariant 6 in `PLAN.md`.
- ADR 0015 -- Slotboom variables with Fermi-Dirac statistics. ADR 0015
  supersedes the carrier-density mapping subsection of this ADR for
  operation under `physics.statistics.model: "fermi-dirac"`. The choice
  of primary unknowns and the coercivity argument in this ADR remain in
  force unchanged.
- Future ADR if we ever need SUPG for a non-Slotboom subsystem (for
  example, hot-carrier transport).
