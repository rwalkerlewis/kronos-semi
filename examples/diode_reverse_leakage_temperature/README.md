# diode_reverse_leakage_temperature

## What this is

A practical 1D pn diode reverse-leakage characterization across three
temperatures (250 K, 300 K, 350 K). 20 um total length with a step
doping at the midpoint (N_A = 1e16 cm^-3 on the p-side, N_D = 1e16
cm^-3 on the n-side), uniform 800 cells, SRH lifetime tau = 100 ns,
constant mobility. The anchor JSON pins T = 300 K and the smoke
verifier reruns the configuration at T = 250 K and T = 350 K via
in-process re-runs that override `physics.temperature`. The intended
use case is "how does the reverse leakage current of my pn diode
change as the package warms up by 50 K," and the config is meant to
be cloned and re-parametrized for your own diode.

This example is the educational complement to
[`examples/schottky_iv_temperature`](../schottky_iv_temperature). The
two examples share the same workflow shape (three-temperature sweep,
shipped anchor at the midpoint, verifier loops via override) but
exercise two distinct temperature-dependent physics mechanisms:

- **Schottky thermionic emission.** Saturation current scales as
  `J_sat = A* T^2 exp(-q phi_B / kT)`. Dominant temperature
  dependence is `exp(-q phi_B / kT)` with the metal-to-semiconductor
  barrier height in the exponent.
- **pn diode reverse leakage (this example).** Generation-limited
  leakage scales as `J_gen propto n_i / tau_eff propto exp(-E_g /
  (2 kT))`. Dominant temperature dependence is `exp(-E_g / (2 kT))`
  with **half** of the band gap in the exponent (the factor of two
  comes from the n_i scaling: n_i^2 = N_C N_V exp(-E_g / kT), so
  n_i propto exp(-E_g / 2 kT)).

That `exp(-E_g/(2 kT))` vs `exp(-q phi_B / kT)` distinction is the
educational pivot: same workflow, different temperature signature,
two different physics mechanisms. Side-by-side these two examples
let users see how the activation energy from the Arrhenius slope
differs by a factor of two between the two mechanisms.

This is an example, not a V&V gate. For the Shockley-saturation
analytical reference at zero or moderate forward bias, see
[`benchmarks/pn_1d_bias_reverse`](../../benchmarks/pn_1d_bias_reverse)
(structural template for this example's geometry).

## Physics features exercised

- **SRH recombination / generation.** `physics.recombination.srh =
  true` with `tau_n = tau_p = 100 ns`. In reverse bias the depletion
  region widens and the local quasi-Fermi-level split drives the SRH
  rate negative (generation rather than recombination); the
  reverse-leakage current density is approximately
  `J_gen = q n_i W_dep / tau_eff` where `tau_eff` aggregates `tau_n`
  and `tau_p`. This is the dominant generation mechanism in the
  depletion region at the doping levels used here; Auger does not
  contribute meaningfully at low injection.
- **Temperature dependence.** `physics.temperature` is the parametric
  knob; the verifier sweeps {250 K, 300 K, 350 K}. Temperature
  enters through the n_i scaling (the dominant factor in the
  generation rate), the thermal voltage `kT/q`, and the
  band-gap narrowing (a small effect at moderate doping).
- **Constant mobility.** `physics.mobility.model = "constant"`. The
  reverse-leakage current is generation-rate-limited rather than
  drift-limited at this geometry, so the mobility model choice is
  largely irrelevant for the reverse-bias regime; the constant branch
  is kept minimal for clarity.
- **Auger off.** `physics.recombination.auger = false`. Auger
  scales as n^2 p (or n p^2) and is negligible at the carrier
  densities reached under reverse bias (deeply depleted; n p << n_i^2
  at strong reverse bias).

## How to run

From the repository root:

```
docker compose run --rm benchmark diode_reverse_leakage_temperature
```

Output lands in `results/diode_reverse_leakage_temperature/`:

- `iv_reverse_leakage_temperature.png`: three semilog-y |I_R| vs V_F
  curves overlaid (color-coded by temperature) plus a companion
  Arrhenius subplot of log(|I_R|) vs 1/T at V_F = -3 V.
- The result JSON written by the engine for the anchor T = 300 K
  run is in the standard `results/diode_reverse_leakage_temperature/`
  location; the companion runs land in `_T_<temperature>K/`
  subdirectories.

## Expected output

The reverse-leakage current density at V_F = -3 V scales roughly as

    |J_R(T)| ~ q * n_i(T) * W_dep(V_F) / tau_eff(T)

with the dominant temperature dependence in n_i:

    n_i(T) propto T^(3/2) exp(-E_g / (2 kT))

For Si with E_g = 1.12 eV, the exponential factor between T = 250 K
and T = 350 K is

    exp(-E_g / (2 k * 250 K)) / exp(-E_g / (2 k * 350 K)) ~ 5e+1

So at fixed reverse bias, |I_R| at 350 K is roughly 30 to 50 times
larger than |I_R| at 250 K. The exact ratio depends on:

- The depletion-region width W_dep(V_F) (which has its own weak
  temperature dependence through V_bi(T)).
- The effective lifetime tau_eff(T) (mostly T-independent for SRH).
- The pre-factor T^(3/2) in n_i (a factor of (350/250)^(3/2) ~ 1.7 on
  top of the exponential).
- Band-gap narrowing at high doping (negligible here at N_A = N_D =
  1e16 cm^-3).

The qualitative gate the smoke verifier checks is the strict ordering
`|I_R(350 K)| > |I_R(300 K)| > |I_R(250 K)|` at V_F = -3 V, plus
finiteness and monotonic |I| vs |V_F| in the reverse-bias range.
**The verifier does not assert a specific ratio**; the
exp(-E_g/(2 kT)) Arrhenius extraction is the user's exercise, not the
verifier's gate. Tightening the gate would create churn whenever
parameters change.

If the reported ordering at V_F = -3 V is wrong (for example
|I_R(250 K)| > |I_R(350 K)|), check:

- The temperature setting (`physics.temperature`). The companion runs
  override this in deepcopies; the anchor run keeps it at 300 K.
- The SRH lifetime (`physics.recombination.tau_n` and `tau_p`). The
  reverse-leakage current scales inversely with tau; if your
  override changes tau, the absolute leakage changes but the
  ordering should be preserved.
- The doping levels. At very heavy doping the Shockley reverse-
  saturation contribution may compete with generation-region
  leakage; for the moderate doping used here (1e16) generation
  dominates at moderate reverse bias.
- Band-gap-narrowing flags. The framework's BGN model adds a small
  T-dependent shift; it should not flip the ordering but may move
  the absolute value.

To extract an activation energy from this example, fit
log(|I_R(V_F = -3 V, T)|) vs 1/T to a straight line; the slope is
`-E_a / k` where the activation energy E_a is approximately
`E_g / 2 ~ 0.56 eV` for the SRH-generation-dominated regime.

## How to adapt

Most users will want to start by changing one of:

- **Temperature set.** Edit the verifier (or override
  `physics.temperature` in your own copy of the JSON) to whatever
  temperature points are relevant to your operating envelope.
- **Doping.** Change `doping[0].profile`. Heavier doping shrinks
  W_dep at fixed V_F and reduces leakage; very light doping
  (sub-1e15 cm^-3) widens W_dep enough that surface states and
  contact effects start to matter.
- **Lifetime.** `physics.recombination.tau_n` / `tau_p`. Shorter
  lifetime increases leakage. The shipped value (100 ns) is moderate
  for clean Si; defect-rich or implant-damaged Si has tau in the
  ns to ps range.
- **Bias range.** `contacts[0].voltage_sweep`. The default
  [0, -5] V at 0.1 V step covers moderate reverse bias; widen to
  characterize approach to breakdown.
- **Geometry.** `mesh.extents` and `doping[0].profile.location`. The
  shipped 20 um geometry with the junction at 10 um gives plenty of
  bulk on each side for the depletion region to grow at deep
  reverse bias.

## Choice of parametric mechanism

Same as `schottky_iv_temperature`: the schema does not have a
built-in "sweep over `physics.temperature`" mechanism, so the
verifier loops via in-process semi.run.run() calls on deepcopies.
The shared `_temperature_overrides` helper in
`scripts/run_benchmark.py` is the same one used by
`schottky_iv_temperature`; future temperature-sweep examples can
register through it.

## Notes on CI runtime

Three small 1D 800-cell bias_sweep runs are cheap; estimate
~30 seconds total wall-clock on a typical CI runner. No runtime
tuning needed.
