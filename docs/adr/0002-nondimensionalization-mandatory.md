# 0002. Nondimensionalization is mandatory

- Status: Accepted
- Date: 2026-04-20

## Context

The raw semiconductor device equations combine quantities that span more
than 40 orders of magnitude:

| Quantity            | Typical value               |
|---------------------|-----------------------------|
| $q$                 | $1.6 \times 10^{-19}$ C     |
| $\varepsilon_0$     | $8.85 \times 10^{-12}$ F/m  |
| $n_i$ (Si, 300 K)   | $1 \times 10^{16}$ m$^{-3}$ |
| $N$ (typical)       | $1 \times 10^{23}$ m$^{-3}$ |
| $\psi$              | $\sim 1$ V                  |
| $L$ (device)        | $\sim 10^{-6}$ m            |

Feeding these directly into Newton gives a Jacobian with condition
number in excess of $10^{30}$. The linear solve at each Newton step is
numerically meaningless at that conditioning, and SNES fails to
converge (residuals stagnate or explode).

This applies equally to Poisson alone and to the coupled
Poisson-continuity system.

## Decision

All FEM forms operate on nondimensional unknowns:

- $\hat\psi = \psi / V_t$
- $\hat n, \hat p, \hat N = \{n, p, N\} / C_0$
- $\hat\Phi_n, \hat\Phi_p = \{\Phi_n, \Phi_p\} / V_t$

The **spatial coordinate is deliberately kept in meters** so that mesh
extents and facet locations in the JSON remain in SI. This is an
asymmetric scaling: fields scaled, geometry not.

The consequence for the scaled Poisson LHS coefficient is

$$
L_D^2 \varepsilon_r = \frac{\varepsilon_0 V_t}{q C_0}\, \varepsilon_r
     = \lambda^2 L_0^2 \varepsilon_r,
$$

not $\lambda^2 \varepsilon_r$. This is spelled out in `docs/PHYSICS.md`
Section 2.

## Consequences

Easier:

- Newton converges quadratically on well-posed problems.
- Boundary conditions simplify ($\hat\psi_\mathrm{BC} =
  \operatorname{asinh}(\hat N / 2 \hat n_i) + V_\mathrm{applied} / V_t$).
- The small parameter $\lambda^2 L_0^2 / L_\mathrm{feature}^2$ makes the
  depletion regime explicit as a singular perturbation.

Harder:

- Every new physics form must be written in scaled variables from the
  start. Mixing scaled and unscaled quantities is a source of bugs.
- The asymmetric scaling choice (fields scaled, x not) is unusual and
  caused the M1 coefficient bug (see `docs/PHYSICS.md` Section 2.3).
  All new physics modules must reference that section and use $L_D^2$,
  not $\lambda^2$, when the coordinate is in meters.

## Related

- `docs/PHYSICS.md` Section 2 for the derivation.
- Invariant 2 and 3 in `PLAN.md`.
