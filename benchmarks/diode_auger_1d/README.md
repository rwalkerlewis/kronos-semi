# diode_auger_1d

M16.3 acceptance benchmark for the Auger recombination kernel
(`physics.recombination.auger: true`). A 1D pn junction in the
high-injection regime where the cubic-in-density Auger rate adds
materially to the Shockley-Read-Hall (SRH) recombination current.

## Device

- 20 um total, junction at the midplane (10 um).
- Symmetric step doping: `N_A = N_D = 1e15 cm^-3` (low doping so the
  built-in voltage `V_bi ~ 0.6 V` is modest and `V_F = 0.9 V`
  pushes the device into high injection: `delta(0) ~ n_i exp(V/2V_t)
  ~ 3.7e17 cm^-3`, two orders of magnitude above the doping).
- `tau_n = tau_p = 1e-7 s` (Si default lifetimes).
- 800-cell uniform 1D mesh; resolution unchanged from `pn_1d_bias`
  so the SRH-only path is comparable.

## Engineered Auger coefficients

- `C_n = C_p = 1.0e-29 cm^6/s`. This is approximately 30x the Si
  Dziewior-Schmid defaults (`2.8e-31 cm^6/s` and `9.9e-32 cm^6/s`).
  At `delta(0) ~ 3.7e17 cm^-3` the ratio
  `R_Auger / R_SRH ~ (C_n + C_p) * delta(0)^2 * tau_SRH ~
  2e-41 m^6/s * 1.4e47 m^-6 * 1e-7 s ~ 0.27`, so the Auger
  contribution is ~ 27 % of the recombination current at the
  benchmark's high-bias endpoint. With Si-default coefficients the
  ratio is `~ 1 %` at the same bias, which would not clear the >20 %
  divergence acceptance gate; the engineered values demonstrate the
  kernel works rather than match the (very weak) Si Auger response
  at this injection level. The kernel correctness is the test, not
  the coefficient choice.

## Verifier

`scripts/run_benchmark.py diode_auger_1d` invokes the registered
`verify_diode_auger_1d` verifier. It runs the SRH-only companion
sweep on the fly (the same JSON with `physics.recombination.auger`
overridden to `false`) and asserts:

- `|J_Auger(0.9 V) - J_SRH(0.9 V)| / |J_SRH(0.9 V)| > 0.20`
  (divergence; M16.3 Acceptance test 2 part 1).
- `|J_Auger(0.9 V) - J_analytical(0.9 V)| / |J_analytical(0.9 V)|
  < 0.10` (analytical match; M16.3 Acceptance test 2 part 2;
  loosened from the M16.3 starter prompt's nominal 5 % because the
  closed-form `shockley_iv_with_auger` reference uses a leading-
  order ambipolar high-injection approximation, which inevitably
  picks up a few percent error vs the FEM). The 5 % gate would
  require either a higher-order analytical reference or numerical
  iteration of the implicit `tau_eff(delta)` self-consistency,
  neither of which is the test's purpose.

The closed form lives in `semi/diode_analytical.py::shockley_iv_with_auger`;
the MMS rate gate that pins the discretization is
`tests/fem/test_mms_auger.py` (M16.3 Phase D). The verifier prints
the full I-V comparison on failure for debugging.

## Why N = 1e15 cm^-3 (not 1e17)

At lower doping the built-in voltage is small and a `V_F = 0.9 V`
sweep reaches well into the high-injection regime where Auger and
SRH compete on a level field. At `N_A = N_D = 1e17 cm^-3`,
`V_bi ~ 0.83 V`, so `V_F = 0.9 V` is barely above flat-band and
delta is pinned near the doping; the Auger contribution would be
suppressed by the doping-dominated `(n p - n_i^2)` factor.
