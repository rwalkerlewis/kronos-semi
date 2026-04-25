# pn_1d_turnon — Transient turn-on benchmark

## Purpose

Verify that the BDF1/BDF2 transient solver captures minority-carrier
injection dynamics correctly. Specifically, the benchmark asserts that
the extracted effective hole lifetime τ_eff from the anode current
transient matches the SRH input lifetime τ_p within 5 %.

## Device

1D pn junction, 20 µm total length:
- p-type left half: N_A = 10¹⁷ cm⁻³
- n-type right half: N_D = 10¹⁷ cm⁻³

Contacts:
- Anode (x = 0): ohmic, V = 0.6 V (forward bias stepped on at t = 0⁺)
- Cathode (x = 20 µm): ohmic, grounded

Physics: Boltzmann statistics, constant mobility (μ_n = 1400, μ_p = 450
cm²/(V·s)), SRH recombination with τ_n = τ_p = 10⁻⁷ s.

## Simulation parameters

| Parameter | Value |
|-----------|-------|
| Mesh resolution | 200 elements |
| dt | 1 ns |
| t_end | 750 ns = 7.5 × 10⁻⁷ s (≈ 7.5 × τ_p) |
| BDF order | 2 |

## Acceptance criterion

Fit J_anode(t) = J_ss × (1 − exp(−t/τ_eff)) in the window
[τ_p, 4τ_p]. Assert |τ_eff − τ_p| / τ_p < 0.05 (5 %).

## Running

```bash
# Full docker run
docker compose run --rm benchmark pn_1d_turnon

# Or from the repo root
python scripts/run_benchmark.py benchmarks/pn_1d_turnon/config.json
python benchmarks/pn_1d_turnon/verify.py results/pn_1d_turnon/iv.json
```

## Physical explanation

At t = 0⁺ the junction is at equilibrium. The applied forward bias of
0.6 V injects minority carriers (holes into the n side, electrons into
the p side). The current rises as the injected charge builds up toward
the steady-state profile. The time constant of this transient is the
minority-carrier lifetime τ_p in the long-diode limit, giving the
characteristic J(t) = J_ss (1 − e^{−t/τ}) shape. Extracting τ from
the transient verifies that the transient solver captures charge storage
correctly.
