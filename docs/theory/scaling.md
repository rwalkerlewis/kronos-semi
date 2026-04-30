# Why nondimensional scaling

A raw 1 µm device at 10¹⁷ cm⁻³ doping has permittivity ~10⁻¹¹,
elementary charge ~10⁻¹⁹, density ~10²³. Newton on that directly
diverges: the Jacobian condition number is ~10³⁰.

Every field in the engine is therefore scaled so the dimensionless
potential ψ̂ is O(1) and the carrier ratios n/C₀ are O(1). The scaled
Poisson equation has a small parameter

$$\lambda^{2} = \frac{\varepsilon_{r}\,V_{t}}{q\,C_{0}\,L_{0}^{2}} \sim 10^{-4},$$

which is the squared Debye-length-to-device ratio. This is a singular
perturbation, and that is *exactly* the depletion-region physics: the
small parameter is the physical reason space-charge regions are
narrow.

## Where this lives in the code

- [`semi/scaling.py`](../../semi/scaling.py) — `Scaling` dataclass with
  reference units `L_0`, `V_0 = V_t`, `C_0`, and helpers
  `make_scaling_from_config`. Pure Python, no dolfinx.
- [`semi/physics/poisson.py`](../../semi/physics/poisson.py) — the
  Poisson coefficient `L_D^2 * eps_r = lambda2 * L0^2 * eps_r` is
  applied to the stiffness term.
- [`semi/physics/drift_diffusion.py`](../../semi/physics/drift_diffusion.py)
  — the continuity residual coefficient is `L_0^2 * mu_hat` where
  `mu_hat = mu / mu_0`.

## Invariants

- Mesh coordinates stay in physical meters; only field values are
  nondimensionalized.
- Never skip `make_scaling_from_config` in a runner. The raw
  equations have a Jacobian condition number > 10³⁰ and Newton
  diverges immediately.
