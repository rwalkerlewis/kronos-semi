# pn_1d_bias

Forward-bias 1D pn junction with coupled drift-diffusion (Slotboom form)
and Shockley-Read-Hall recombination. 20 um device, symmetric
N_A = N_D = 1e17 cm^-3, tau_n = tau_p = 1e-8 s.

The benchmark sweeps the anode from 0 V to +0.6 V in 50 mV steps and
records (V, J) at each step. Verification compares the high-bias
current (V >= 0.5 V) against Shockley long-diode theory

    J_s = q * n_i^2 * (D_n/(L_n N_A) + D_p/(L_p N_D)),  J = J_s (exp(V/V_t) - 1)

The match is tightest (within 10%) at V = 0.6 V where diffusion
dominates. At lower bias, depletion-region SRH recombination raises the
ideality factor (Sah-Noyce-Shockley), so Shockley is not a valid
reference below about 0.5 V and the verifier only checks qualitative
properties (monotone J(V), super-linear growth).

Run:

    docker compose run --rm benchmark pn_1d_bias
