# rc_ac_sweep

M14 small-signal AC validation benchmark.

The device is a **1D pn junction** identical to `pn_1d_bias_reverse`:
20 µm symmetric step junction (N_A = N_D = 1e17 cm^-3), τ = 1e-8 s,
800 cells. The runner is configured with `solver.type = "ac_sweep"` at
DC bias V_DC = -1.0 V and a 41-point logspace frequency sweep from
1 Hz to 1 GHz.

At reverse bias the diode is depletion-dominated, so the small-signal
admittance Y(ω) is approximately Y(ω) ≈ -j ω C_dep(V_DC) over the
swept range (no measurable conduction at this bias, well below the
1/(R_bulk·C_dep) cross-over which is in the 10 GHz range for this
device).

## Verifier

`benchmarks/rc_ac_sweep/verify.py` checks:

1. **Capacitance plateau matches analytical depletion C within 5%** at
   every frequency in [1 Hz, 1 MHz]. The analytical reference is

       C_dep(V_DC) = sqrt(q * eps_Si * N_eff / (2 * (V_bi - V_DC)))
       N_eff = N_A * N_D / (N_A + N_D)
       V_bi  = V_t * ln(N_A * N_D / n_i^2)

2. **Frequency flatness**: max(C) / min(C) over [1 Hz, 1 MHz] must be
   within 2%. Captures the depletion regime: no carrier-transport
   dispersion is expected below ~10 MHz for this device.

3. **High-frequency cap roll-off** is allowed: C above 100 MHz may
   drift modestly as carrier transit times start to enter the picture.
   This is informational, not asserted.

## Run

    docker compose run --rm benchmark rc_ac_sweep
