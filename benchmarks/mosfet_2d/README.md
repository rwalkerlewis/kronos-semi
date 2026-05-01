# mosfet_2d benchmark

A 2D n-channel MOSFET on a uniform 50 x 21 builtin mesh. P-type Si body
(N_A = 1e16 cm^-3) with a 5 nm SiO2 gate stack; n+ source and drain
implants are Gaussian (peak 5e19 cm^-3, sigma = (0.4 um, 0.15 um)).

The bias_sweep runner sweeps V_GS from 0 to 1.5 V in 0.05 V steps with
V_DS held at 0.05 V (drain ohmic, static). The smaller step size is
needed because the surface-potential transition through onset of
strong depletion (V_GS ~ 0.2 to 0.3 V) is sharp enough to defeat
Newton with a 0.1 V step; the JSON's `solver.continuation` block
also widens `max_halvings` from the default 6 to 10 and tightens
`min_step` to 1e-5 V for the same reason. The verifier compares the
2D-simulated drain current against the Pao-Sah long-channel linear
formula in the depletion-to-strong-inversion transition region.

## Pao-Sah linear-regime reference

In the long-channel limit at V_DS << V_GS - V_T (linear / triode
regime), the drain current per unit channel width is

    I_D / W = (mu_n / L_ch) * C_ox * (V_GS - V_T) * V_DS

where:

- `mu_n` is the channel electron mobility (1400 cm^2/V/s for Si at 300 K
  by default).
- `L_ch` is the channel length, taken here as the lateral extent of the
  gate-oxide region (3.0 um). The implant-tail correction is absorbed
  into the 20 % verifier tolerance.
- `C_ox = eps_ox / t_ox` is the oxide capacitance per unit area
  (eps_ox = 3.9 * eps_0, t_ox = 5 nm; C_ox ~ 6.91 mF/m^2).
- `V_T` is the long-channel threshold voltage from
  `semi.cv.analytical_moscap_params` (Hu Eq. 5.5.1). For this device
  V_T ~ 0.66 V on top of V_FB = -phi_F.

The 2D simulation reports current density per unit z width through the
drain facet line (length 2.005 um). The verifier multiplies the recorded
J_drain by the drain line length to obtain I_D / W and compares.

## Verifier window

V_GS in [V_T + 0.2, V_T + 0.6] V at V_DS = 0.05 V; tolerance 20 %.

- The lower bound (V_T + 0.2) sits clear of subthreshold so the
  Boltzmann-tail correction does not dominate. The strict-Pao-Sah
  square-law form needs V_GS - V_T well above kT/q (~0.026 V at 300 K).
- The upper bound (V_T + 0.6) keeps V_DS / (V_GS - V_T) small so
  saturation and velocity-saturation corrections are still negligible.
  M16.1 (Caughey-Thomas field-dependent mobility) widens the window
  upward; in the meantime the 20 % tolerance absorbs the residual
  high-field deviation.
- The 20 % tolerance further covers the long-channel approximation
  (L_ch as the gate-oxide lateral extent rather than the
  implant-corrected effective channel length), ohmic source-and-drain
  series resistance (small at V_DS = 0.05 V on a 5 um body), and
  the 50 x 21 mesh discretisation error.

## Run

    python scripts/run_benchmark.py mosfet_2d

Exits 0 with the verifier reporting the worst-case relative error in
the window. Plots land under `results/mosfet_2d/`:

- `id_vgs.png`: I_D vs V_GS overlaid with the Pao-Sah analytical
  curve and the verifier window.
