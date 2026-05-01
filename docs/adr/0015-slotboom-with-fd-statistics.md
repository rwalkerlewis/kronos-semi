# ADR 0015: Slotboom variables with Fermi-Dirac statistics

**Date:** 2026-05-01
**Status:** Accepted
**Supersedes (partial):** ADR 0004 (Slotboom variables for DD) -- the
  carrier-density mapping subsection only. ADR 0004 remains in force for
  the choice of primary unknowns and for all Boltzmann-regime operation.
**Milestone:** M16.4 -- Fermi-Dirac statistics

---

## Context

ADR 0004 adopted Slotboom quasi-Fermi potentials `(psi, phi_n, phi_p)` as
the primary DD unknowns and derives electron and hole densities via the
Boltzmann relations:

```
n = n_i * exp((psi - phi_n) / V_t)
p = n_i * exp((phi_p - psi) / V_t)
```

These expressions break down when the doping density approaches or exceeds
the effective density of states (N_C or N_V). For Si this occurs above
approximately 10^19 cm^-3; for GaN and SiC it can occur at much lower
concentrations. Fermi-Dirac statistics correct the carrier densities by
replacing the exponential with the Fermi-Dirac integral of order 1/2.

The question this ADR addresses is: **how does the Slotboom formulation
survive the transition to Fermi-Dirac statistics?** Three options exist:

1. **Abandon Slotboom, switch primary unknowns to (n, p).**
   Rejected: ADR 0004 and ADR 0014 are locked; this path reintroduces
   the Peclet-number instability that Slotboom cures.

2. **Keep Slotboom primary unknowns; redefine the n(psi, phi_n) mapping.**
   Slotboom's coercivity argument depends only on the current being a
   pure gradient `J_n = -q mu_n n grad(phi_n)`; it does not require
   that `n` is an exponential of `(psi - phi_n)`. Any monotone mapping
   from the quasi-Fermi potential to the carrier density preserves the
   coercivity and therefore the stability. Chosen.

3. **Use a different set of quasi-Fermi variables.**
   Possible in principle but introduces a variable transformation that
   is inconsistent with every existing test and benchmark. Rejected.

## Decision

Keep `(psi, phi_n, phi_p)` as the primary DD unknowns (ADR 0004 / ADR
0014 hold unchanged). Replace the Boltzmann carrier-density mapping with
the Fermi-Dirac mapping when `physics.statistics.model == "fermi-dirac"`.

### Fermi-Dirac mapping

Define the reduced quasi-Fermi energy for electrons as:

```
eta_n = (E_Fn - E_C) / kT  =  (psi - phi_n - chi_tilde) / V_t
```

where `chi_tilde = chi / V_t` is the scaled electron affinity (zero for
a homogeneous material; position-dependent in M17 heterojunctions).

Then:

```
n = N_C * F_half(eta_n)
p = N_V * F_half(eta_p)
```

where `F_half` is the Fermi-Dirac integral of order 1/2:

```
F_half(eta) = (2/sqrt(pi)) * integral_0^inf  sqrt(E) / (1 + exp(E - eta)) dE
```

The current retains its pure-gradient form:

```
J_n = -q mu_n n grad(phi_n)
```

because `phi_n` remains the quasi-Fermi potential. The coefficient `n`
is now evaluated via `F_half` instead of the exponential, but the weak
form structure is unchanged.

### Approved approximations for F_half

Two approximations are supported, selected via
`physics.statistics.approximation`:

1. **Blakemore** (default, fast):
   ```
   F_half(eta) ~= 1 / (exp(-eta) + 0.27)
   ```
   Accurate to ~5% for all eta; suitable for production use.
   Pure Python; callable from UFL.

2. **FDI-1/2 via Bednarczyk-Bednarczyk rational approximation**
   (validation):
   ```
   F_half(eta) -- rational approximation, degree (p=4, q=4)
   ```
   Accurate to 4e-4 for all eta in [-20, 20]. Used for MMS validation
   and benchmark comparison; slightly slower but within 10% of Blakemore
   for assembly cost.

The full Fermi integral is not evaluated via numerical quadrature in the
FEM kernel; that would dominate assembly time. The B-B rational
approximation is the production-quality alternative.

### Inverse mapping (for boundary data)

At ohmic contacts, the equilibrium quasi-Fermi level equals the applied
bias: `phi_n = phi_p = V / V_t`. The Dirichlet BC on `psi` is set by
inverting the charge-neutrality condition `p - n + N_D - N_A = 0` with
FD densities. The numerical root-find uses `scipy.optimize.brentq` on
the pure-Python `F_half` expression; it runs once at setup time and is
not in the Newton iteration hot path.

### Generalized Einstein relation

Under Fermi-Dirac statistics the diffusivity is no longer simply
`D = (kT/q) mu`. The generalized Einstein relation is:

```
D_n / mu_n = (kT/q) * (n / (dF_half/d eta_n * N_C / V_t))
            = (kT/q) * F_half(eta_n) / F_minus_half(eta_n)
```

where `F_minus_half` is the Fermi-Dirac integral of order -1/2.

In the Slotboom formulation this ratio does **not** appear explicitly
in the weak form because the current is expressed as `J_n = -q mu_n n
grad(phi_n)` and the diffusion is absorbed into the gradient of `phi_n`.
The generalized Einstein relation therefore only appears in the
**interpretation** of the current (not in the assembled residual or
Jacobian). This is recorded here to prevent future implementors from
adding a spurious D_n correction term to the Slotboom form.

## Consequences

Easier:

- Stability of the Slotboom weak form is preserved; no new stabilization
  is needed.
- Boltzmann remains the default; all existing tests and benchmarks are
  unchanged.
- FD is an opt-in model selection; the Boltzmann path is not touched.

Harder:

- The `n(psi, phi_n)` mapping is no longer a simple UFL expression; it
  requires a user-defined nonlinear function compiled via
  `dolfinx.fem.Expression` or a UFL `conditional`.
- Boundary data (psi_eq at ohmic contacts) requires a numerical root-find
  at setup time.
- The MMS manufactured solution must use the FD inverse to construct
  consistent boundary values.

## Related

- ADR 0004 -- Slotboom variables for steady-state DD (primary-unknown
  choice remains; mapping superseded here for FD regime).
- ADR 0014 -- Slotboom transient (applies unchanged; mapping swap is
  transparent because chain-rule mass term still uses `n_ufl(psi, phi_n)`).
- ADR 0016 -- Heterojunction band alignment (chi_tilde extension of this
  mapping).
- `docs/PHYSICS.md` Section 10 "Statistics and ionization" (to be written
  in Phase J).
- M16.4 implementation: `semi/physics/statistics.py`.
