# Why Slotboom variables for drift-diffusion

Galerkin FEM on raw drift-diffusion is unstable when drift dominates
diffusion, which is almost everywhere in a real device. The Slotboom
transform rewrites the continuity equations in quasi-Fermi potentials
$\Phi_n$, $\Phi_p$:

$$n = n_i\,\exp\!\big((\psi - \Phi_n)/V_t\big), \quad
  p = n_i\,\exp\!\big((\Phi_p - \psi)/V_t\big).$$

Currents become pure gradients of $\Phi$:

$$J_n = -q\,\mu_n\,n\,\nabla\Phi_n, \qquad
  J_p = -q\,\mu_p\,p\,\nabla\Phi_p,$$

so the weak form is well-posed without upwind stabilisation. The
finite-element coupling is symmetric in psi and the two
quasi-Fermi rows.

## Where this lives in the code

- [`semi/physics/drift_diffusion.py`](../../semi/physics/drift_diffusion.py)
  — Slotboom continuity residuals; signs follow the ADR convention
  (electron residual carries `-R*v`, hole residual `+R*v`).
- [`semi/physics/slotboom.py`](../../semi/physics/slotboom.py) —
  conversion helpers between primary-density `(psi, n, p)` and
  quasi-Fermi `(psi, phi_n, phi_p)` representations.
- [`semi/runners/bias_sweep.py`](../../semi/runners/bias_sweep.py) —
  primary unknowns are `(psi, phi_n, phi_p)` (ADR 0014).
- [`semi/runners/transient.py`](../../semi/runners/transient.py) —
  also Slotboom after ADR 0014 superseded ADR 0009.

## Verification

Slotboom DD is verified by the MMS-DD suite (three variants ×
three grids) under
[`semi/verification/mms_dd.py`](../../semi/verification/mms_dd.py).
Finest-pair $L^2$ rates are at theoretical 2.000.

## Caveat: Galerkin Slotboom on coarse grids at low bias

Galerkin Slotboom on P1 Lagrange is OK for coarse Shockley match at
high forward bias ($V \ge 0.5$ V, within 10% of ideal
$J_s(\exp(V/V_t)-1)$). At lower bias the ideality factor is > 1
because SRH depletion-region recombination dominates. For realistic
Si doping (10¹⁷ cm⁻³), $\hat{n}$ varies by ~15 orders of magnitude
across the device; the continuity Jacobian is therefore
ill-conditioned. Scharfetter–Gummel edge-flux assembly
([`semi/fem/scharfetter_gummel.py`](../../semi/fem/scharfetter_gummel.py))
is the proper fix and was wired in M13.1.

## ADRs

- [adr/0004-slotboom-variables-for-dd.md](../adr/0004-slotboom-variables-for-dd.md)
- [adr/0012-scharfetter-gummel.md](../adr/) — SG primitives.
- [adr/0014-slotboom-primary-unknowns.md](../adr/) — supersedes ADR 0009.
