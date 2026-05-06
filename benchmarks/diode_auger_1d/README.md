# diode_auger_1d

M16.3 acceptance benchmark for the Auger recombination kernel
(`physics.recombination.auger: true`). A 1D pn junction in the
moderate-to-high-injection regime where the cubic-in-density Auger
rate adds materially to the Shockley-Read-Hall (SRH) recombination
current.

## Device

- 20 um total, junction at the midplane (10 um).
- Symmetric step doping: `N_A = N_D = 1e17 cm^-3`.
- `tau_n = tau_p = 1e-7 s` (Si default lifetimes).
- 800-cell uniform 1D mesh; resolution unchanged from `pn_1d_bias`
  so the SRH-only path is comparable.
- Forward sweep `V_F` in `[0, 0.9] V` (0.05 V step).

The earlier M16.3 draft of this benchmark used `N_A = N_D = 1e15 cm^-3`
on the assumption the device would reach an asymptotic
`delta(0) ~ n_i exp(V/2V_t) ~ 3.7e17 cm^-3` at `V_F = 0.9 V`.
Empirically the FEM solve with 1e15 doping is bulk-series-resistance
limited and only injects to `n ~ p ~ 3.6e15 cm^-3` at the junction
at `V_F = 0.9 V` (the bulks drop most of the applied bias).
At that injection level the Auger contribution is below 0.01 % of
the SRH current and the divergence acceptance gate cannot be cleared.
Bumping the doping to `1e17 cm^-3` raises the injected `n p`
product enough that Auger is materially exercised at the bias the
device actually reaches.

## Engineered Auger coefficients

- `C_n = C_p = 1.0e-27 cm^6/s`. This is approximately 3000x the Si
  Dziewior-Schmid defaults (`2.8e-31 cm^6/s` and `9.9e-32 cm^6/s`).
  The values are engineered (not physical Si Auger): the kernel's
  correctness is the test, not the coefficient choice. With the
  device actually reaching `n ~ p ~ 1.5e17 cm^-3` at the junction
  the engineered C produces a ~28 % SRH-vs-(SRH+Auger) divergence
  at `V_F = 0.9 V`, comfortably above the 20 % acceptance gate.
  Si-default coefficients on this device produce a sub-percent
  divergence and do not exercise the kernel meaningfully.

## Verifier

`scripts/run_benchmark.py diode_auger_1d` invokes the registered
`verify_diode_auger_1d` verifier. It runs the SRH-only companion
sweep on the fly (the same JSON with `physics.recombination.auger`
overridden to `false`) and asserts:

- `|J_Auger(0.9 V) - J_SRH(0.9 V)| / |J_SRH(0.9 V)| > 0.20`
  (divergence; M16.3 Acceptance test 2 part 1).

The closed-form Hall-Auger long-diode reference
(`semi/diode_analytical.py::shockley_iv_with_auger`) is computed
and printed for diagnostics. It is **not** a hard gate: the closed
form assumes pure high-injection long-diode operation without bulk
series resistance, and the FEM device cannot realize that
idealization at `V_F = 0.9 V`. The simulated junction voltage is
materially below the applied bias because the resistive bulks drop
the difference, so the simulated `|J|` is roughly an order of
magnitude smaller than the asymptotic closed form predicts. The
kernel's correctness is pinned by the divergence gate above and by
the MMS Variant F rate gate (`semi/verification/mms_dd.py`,
`tests/fem/test_mms_auger.py`); the analytical comparison only
fixes the order of magnitude.

The verifier prints the full I-V comparison on failure for
debugging.
