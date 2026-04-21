# 1D pn junction benchmark

Symmetric abrupt pn junction, silicon, $N_A = N_D = 10^{17}$ cm⁻³, junction at 1 µm.

## Purpose

Equilibrium-Poisson verification. The depletion approximation gives analytical expressions for:

- Built-in potential $V_{bi} = V_t \ln(N_A N_D / n_i^2) \approx 0.833$ V
- Depletion width $W = \sqrt{2 \varepsilon V_{bi} / (q N_{\text{eff}})} \approx 147$ nm
- Peak field $|E_{\max}| = q N_A x_p / \varepsilon \approx 113.5$ kV/cm

These are verified against the FEM solution in `notebooks/01_pn_junction_1d.ipynb`.

## Physics covered

- Nonlinear Poisson with Boltzmann equilibrium statistics
- Ohmic contact BCs from charge-neutrality / asinh formula
- Nondimensional scaling (see `semi.scaling.Scaling`)

## Physics NOT covered (future benchmarks)

- Applied bias (this is equilibrium only)
- Drift-diffusion continuity equations
- Recombination (generation-recombination is trivially zero at equilibrium)

## Running

```bash
python -c "from semi import schema, run; result = run.run(schema.load('pn_junction.json'))"
```

Or open the Colab notebook via the README badge.
