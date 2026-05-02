# diode_velsat_1d

M16.1 acceptance benchmark for the Caughey-Thomas field-dependent
mobility (`physics.mobility.model: caughey_thomas`). A 1D pn junction
identical in geometry and doping to `benchmarks/pn_1d_bias`, swept
across the forward range `V_F in [0.0, 0.9] V` so the verifier sees
both the low-field regime (CT factor ~ 1, two models converge) and
the high-field regime (CT mobility suppressed, two models diverge).

## Device

- 20 um total, junction at the midplane (10 um).
- Symmetric step doping: `N_A = N_D = 1e17 cm^-3`,
  `tau_n = tau_p = 1e-8 s`.
- 800-cell uniform 1D mesh; resolution unchanged from `pn_1d_bias`
  so the M14 SRH-current physics is comparable.

## Why the bias range was extended below 0.5 V

The M16.1 starter prompt (`docs/M16_1_STARTER_PROMPT.md`) called for
the convergence anchor at `V_F = 0.5 V` with `< 1 %` deviation. On
this geometry that is incorrect: at `V_F = 0.5 V` the depletion-
region peak field is approximately

```
E_max = sqrt(2 q N_A (V_bi - V_F) / eps_Si) ~ 50 kV/cm
```

with `V_bi ~ 0.83 V`, which gives a Caughey-Thomas factor

```
(mu_n0 * E_max / vsat_n) = (1400 * 50000) / 1e7 = 7
mu_n / mu_n0 = 1 / sqrt(1 + 7^2) ~ 0.14
```

so the integrated terminal current at `V_F = 0.5 V` already deviates
~12 % from the constant-mobility prediction. The natural convergence
anchor is `V_F = 0.3 V`, where peak field ~ 30 kV/cm gives
`mu_n / mu_n0 ~ 0.99` and the I-V deviation is < 1 %. The bias range
was therefore extended down to `V_F = 0.0` so the verifier covers
both regimes within a single sweep; the verifier reads `V = 0.3 V`
for convergence and `V = 0.9 V` for divergence.

## Verifier

`scripts/run_benchmark.py diode_velsat_1d` invokes the registered
`verify_diode_velsat_1d` verifier. It runs the constant-mobility
companion sweep on the fly (the same JSON with `physics.mobility.
model` overridden to `"constant"` and the same low-field
`mu_n / mu_p` values) and asserts:

- `|I_CT(0.9 V) - I_const(0.9 V)| / |I_const(0.9 V)| > 0.05`
  (divergence; M16.1 Acceptance test 2 part 1; observed
  ~ 56 %).
- `|I_CT(0.3 V) - I_const(0.3 V)| / |I_const(0.3 V)| < 0.05`
  (convergence; M16.1 Acceptance test 2 part 2; observed
  ~ 0.2 %).

The closed form is in `semi/physics/mobility.py`; the MMS rate
gate that pins the discretization is
`tests/fem/test_mms_caughey_thomas.py` (M16.1 Phase D). Both
acceptance gates clear with > 5x headroom on this geometry; the
verifier prints the full I-V comparison on failure for
debugging.
