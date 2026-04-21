# pn_1d_bias_reverse

Reverse-bias 1D pn junction. Same device as `pn_1d_bias` (20 um
symmetric 1e17 junction, tau = 1e-8 s) but the anode is swept from
0 V down to -2 V in 0.1 V steps with adaptive growth up to 0.5 V.

Verifier: in the saturation region [-2, -0.5] V the simulated reverse
current |J| must match the Shockley long-diode saturation current
J_s = q n_i^2 (D_n/(L_n N_A) + D_p/(L_p N_D)) within 20%. The
tolerance is deliberately loose because J_s is very small
(~1.5e-7 A/m^2) and the relative error on signed integrals is
numerically noisy. The interval near zero (-0.5, 0) is excluded
because the transition region is not quantitatively Shockley.

Run:

    docker compose run --rm benchmark pn_1d_bias_reverse
