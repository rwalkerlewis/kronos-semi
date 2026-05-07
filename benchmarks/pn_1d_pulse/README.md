# pn_1d_pulse

A 1D pn diode under a switching transient. The anode bias is held at
0.0 V for `t < 5 ns` and steps to 0.6 V at `t = 5 ns`, materialized
through the M16.7 `voltage_t.step` schema variant on the anode
contact.

This benchmark demonstrates the `voltage_t` step variant. The formal
verification gate for time-varying contact voltage lives in audit
case 06 (`tests/audit/test_06_transient_fft_vs_ac_sweep.py`); the
verifier here is a lightweight smoke test (finiteness, qualitative
turn-on) so the demonstration stays fast in CI.

## Geometry and physics

Identical to `pn_1d_turnon` so the post-step current is directly
comparable:

- 1D symmetric step junction, length 20 um, junction at 10 um.
- N_A = N_D = 1.0e17 cm^-3 on each side.
- mu_n = 1400, mu_p = 450 cm^2/V/s.
- SRH lifetimes tau_n = tau_p = 1.0e-7 s.
- Boltzmann statistics, T = 300 K.

The mesh is the same 800-cell uniform grid that `pn_1d_turnon` uses.

## Scenario

```text
V_anode(t) = 0.0  V    for t <  5 ns
             0.6  V    for t >= 5 ns
```

The pre-step regime is reverse / zero bias: the diode is at
equilibrium and the terminal current is dominated by SRH leakage
(several orders of magnitude smaller than the post-step injection
current). The post-step regime is the same forward-bias turn-on as
`pn_1d_turnon`, but reached via the M16.7 `voltage_t.step` schema
field instead of a fixed bias applied at `t = 0`.

The 50 ns simulation gives 25 timesteps before the step at
`t = 5 ns` and 225 timesteps after, both at `dt = 200 ps`.

## Difference from `pn_1d_turnon`

`pn_1d_turnon` ships its target bias as `contacts[].voltage = 0.6` and
runs with `bc_ramp_steps = 0`, so the device sees the full forward
bias at `t = 0`. There is no time before the step; only the turn-on
transient is captured.

`pn_1d_pulse` ships `voltage_t = {type: "step", t0: 5e-9, v0: 0.0,
v1: 0.6}` and `bc_ramp_steps = 0`, so the device sits at zero bias
for the first 5 ns and then sees the step. Both runs reach the same
forward-bias steady state by the end of the 50 ns window; the
verifier asserts agreement within 20 %.

## Verifier

`scripts/run_benchmark.py verify_pn_1d_pulse` checks:

1. I(t) is finite and non-NaN at every recorded timestep.
2. The pre-step current at `t = 4.8 ns` (just before the step) is
   at least an order of magnitude smaller than the post-step current
   at `t = 50 ns`.
3. The post-step current at `t = 50 ns` is within 20 % of the
   `pn_1d_turnon` final-time current at the same physical conditions
   (recomputed analytically from the diode forward I-V at V_F = 0.6 V
   to keep this benchmark free of cross-benchmark file dependencies).

The 20 % tolerance reflects the slightly different solver warmup
between the two paths (`pn_1d_turnon` starts at the V_F equilibrium
IC; `pn_1d_pulse` starts at the V = 0 equilibrium IC and the time
loop walks past the step).
