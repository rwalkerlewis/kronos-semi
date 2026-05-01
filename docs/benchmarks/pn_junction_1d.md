# Benchmark: 1D pn junction (equilibrium)

The simplest benchmark in the repo. A 1D step junction at
$N_a = N_d = 10^{17}$ cm$^{-3}$ in silicon, solved at zero bias, and
compared against the depletion approximation.

## Files

- Input: [`benchmarks/pn_1d/pn_junction.json`](../../benchmarks/pn_1d/pn_junction.json)
- README: [`benchmarks/pn_1d/README.md`](../../benchmarks/pn_1d/README.md)
- Notebook: [`notebooks/01_pn_junction_1d.ipynb`](../../notebooks/01_pn_junction_1d.ipynb)
- Plots: [`results/pn_1d/`](../../results/pn_1d/)

## What it verifies

- $V_{\text{bi}} \approx 0.834$ V for symmetric $10^{17}/10^{17}$ Si junction.
- Peak $|E| \approx 113.5$ kV/cm (within a few percent of depletion approximation).
- Mass-action $np = n_i^2$ in bulk to numerical precision.
- Charge neutrality in quasi-neutral regions to $10^{-15}$.

## Theory

- [`docs/PHYSICS.md`](../PHYSICS.md) — equilibrium Poisson, scaling.
- [`docs/theory/scaling.md`](../theory/scaling.md) — why the scaled form.
- [`docs/theory/slotboom.md`](../theory/slotboom.md) — Slotboom variables (used in the bias sweep follow-up benchmark `pn_1d_bias`).

## Implementation pointers

- Runner: [`semi/runners/equilibrium.py`](../../semi/runners/equilibrium.py).
- Weak form: [`semi/physics/poisson.py`](../../semi/physics/poisson.py).
- Verifier: `verify_pn_1d` in [`scripts/run_benchmark.py`](../../scripts/run_benchmark.py).

## Running

```bash
docker compose run --rm benchmark pn_1d
```

Or from a Colab notebook (links from
[`notebooks/01_pn_junction_1d.ipynb`](../../notebooks/01_pn_junction_1d.ipynb)).
