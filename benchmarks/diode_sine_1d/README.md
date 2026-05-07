# diode_sine_1d

A 1D pn diode under sinusoidal large-signal forward-bias drive. The
anode bias is `V(t) = V_DC + dV * sin(2 pi F t)` with V_DC = 0.4 V,
dV = 50 mV, and F = 1 MHz. The waveform is materialized through the
M16.7 `voltage_t.table` schema variant: the JSON ships 800
`(time, value)` pairs covering 4 cycles at 200 samples per cycle.

This benchmark demonstrates the `voltage_t` table variant. The formal
small-signal AC-vs-FFT verification gate lives in audit case 06
(`tests/audit/test_06_transient_fft_vs_ac_sweep.py`); the verifier
here is a lightweight smoke test (finiteness, fundamental peak,
harmonic content) so the demonstration stays fast in CI.

## Geometry and physics

Identical to `pn_1d_turnon`:

- 1D symmetric step junction, length 20 um, junction at 10 um.
- N_A = N_D = 1.0e17 cm^-3 on each side.
- mu_n = 1400, mu_p = 450 cm^2/V/s.
- SRH lifetimes tau_n = tau_p = 1.0e-7 s.
- Boltzmann statistics, T = 300 K.

The mesh is the same 800-cell uniform grid that `pn_1d_turnon` uses.

## Scenario

```text
V_anode(t) = 0.4 + 0.05 * sin(2 pi * 1e6 * t)   for 0 <= t <= 4 us
```

50 mV around a 0.4 V DC operating point on a Si diode is large-signal
(the small-signal limit would be 1 mV around the same DC point, which
is what audit case 06 drives). The current response therefore has
visible second-harmonic content from the diode-equation
`exp(qV / kT)` nonlinearity; the verifier asserts that the second
harmonic is at most 50 % of the fundamental, which is a generous
gate intended to catch table-construction errors and timestep
issues rather than to constrain the physics.

The simulation runs for 4 cycles at `dt = 5 ns` (200 samples per
cycle). `bc_ramp_steps = 20` ramps the device from V = 0 to V = 0.4 V
before the time loop so the time loop starts near the steady-state
operating point and the transient settles to its sinusoidal limit
within the first cycle.

## Relationship to audit case 06

Audit case 06 drives the same waveform shape on the same
`rc_ac_sweep` device (which has the same geometry as `pn_1d_turnon`)
but with `dV = 1 mV` (small-signal). The audit case computes the
admittance `Y = J_FFT[F] / V_FFT[F]` and asserts agreement with the
M14 AC sweep at the same frequency within 5 %.

`diode_sine_1d` ships at `dV = 50 mV` so the harmonic content is
directly visible to the eye in the FFT plot; it is not a small-
signal Y verifier. The audit case carries the formal V&V gate; this
benchmark is illustrative.

## Verifier

`scripts/run_benchmark.py verify_diode_sine_1d` checks:

1. I(t) is finite and non-NaN at every recorded timestep.
2. The FFT of I(t) (after demeaning and windowing) has its dominant
   peak at the fundamental F = 1 MHz (the input drive frequency).
3. The second-harmonic amplitude `|I_FFT[2F]|` is at most 50 % of
   the fundamental amplitude `|I_FFT[F]|`. A second harmonic that
   dominates indicates either (a) a table-construction bug, (b) a
   timestep too coarse to resolve the fundamental, or (c) a runaway
   nonlinearity at this DC bias.
