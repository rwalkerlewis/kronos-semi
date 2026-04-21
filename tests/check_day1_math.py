"""
Sanity-check the physics math in the Day 1 notebook WITHOUT dolfinx.
Verifies: scaling constants, built-in potential math, depletion widths,
peak field. If any of these are wrong, the notebook will be wrong too.
"""
import numpy as np

# ---- constants ----
Q    = 1.602176634e-19
KB   = 1.380649e-23
EPS0 = 8.8541878128e-12
EPS_R_SI = 11.7
NI_SI    = 1.0e10 * 1.0e6
MU_N_SI  = 1400.0 * 1.0e-4
T  = 300.0
Vt = KB * T / Q
print(f"Vt = {Vt*1000:.3f} mV  (expected ~25.85 mV)  {'OK' if abs(Vt - 0.02585) < 0.001 else 'FAIL'}")

# ---- device ----
L_device = 2.0e-6
x_junction = 1.0e-6
N_A_left  = 1.0e17 * 1.0e6
N_D_right = 1.0e17 * 1.0e6
C0 = max(N_A_left, N_D_right)
L0 = L_device

lambda2 = EPS0 * Vt / (Q * C0 * L0**2)
print(f"lambda^2 (bare) = {lambda2:.3e}")
print(f"lambda^2 * eps_r = {lambda2 * EPS_R_SI:.3e}")
# Debye length sqrt(eps * Vt / (q * N))
L_D = np.sqrt(EPS_R_SI * EPS0 * Vt / (Q * C0))
print(f"Debye length = {L_D*1e9:.2f} nm  (expected ~13 nm for 1e17 Si) "
      f"{'OK' if 10 < L_D*1e9 < 16 else 'FAIL'}")
print(f"(L_D / L0)^2 = {(L_D/L0)**2:.3e}  (should equal lambda^2 * eps_r) "
      f"{'OK' if abs((L_D/L0)**2 - lambda2*EPS_R_SI) / (lambda2*EPS_R_SI) < 0.01 else 'FAIL'}")

# ---- built-in potential ----
# Using psi_eq = Vt * asinh(N_net / (2 n_i))
psi_L = np.arcsinh(-N_A_left / (2.0 * NI_SI))
psi_R = np.arcsinh( N_D_right / (2.0 * NI_SI))
Vbi_from_bc = (psi_R - psi_L) * Vt
Vbi_formula = Vt * np.log(N_A_left * N_D_right / NI_SI**2)
print(f"Vbi via asinh BCs  = {Vbi_from_bc:.4f} V")
print(f"Vbi via Vt*ln(..) = {Vbi_formula:.4f} V")
# These should agree for N >> n_i, since asinh(x) ~ ln(2x) for large x
print(f"Relative diff: {abs(Vbi_from_bc - Vbi_formula) / Vbi_formula * 100:.4f}%"
      f"  {'OK' if abs(Vbi_from_bc - Vbi_formula) / Vbi_formula < 0.01 else 'FAIL'}")

# ---- depletion approximation ----
N_eff = N_A_left * N_D_right / (N_A_left + N_D_right)
eps_Si = EPS_R_SI * EPS0
W = np.sqrt(2.0 * eps_Si * Vbi_formula / (Q * N_eff))
xp = W * N_D_right / (N_A_left + N_D_right)
xn = W * N_A_left  / (N_A_left + N_D_right)
print(f"\nDepletion widths: W = {W*1e9:.1f} nm, xp = {xp*1e9:.1f} nm, xn = {xn*1e9:.1f} nm")

# Peak field: |E_max| = q*N_A*xp/eps = q*N_D*xn/eps
E_peak = Q * N_A_left * xp / eps_Si
print(f"Peak |E| = {E_peak*1e-5:.2f} kV/cm "
      f"(1e17/1e17 Si: ~113 kV/cm expected) "
      f"{'OK' if 90e5 < E_peak < 130e5 else 'FAIL'}")

# Sanity: charge balance, xp * N_A = xn * N_D (symmetric so xp = xn)
print(f"Charge balance check: xp*N_A = {xp*N_A_left:.3e}, xn*N_D = {xn*N_D_right:.3e} "
      f"{'OK' if abs(xp*N_A_left - xn*N_D_right)/(xp*N_A_left) < 1e-6 else 'FAIL'}")

# ---- Bulk carrier check ----
# In p-bulk, psi ~ psi_L, so n = ni * exp(psi_L) and p = ni * exp(-psi_L)
# Check that p_bulk ~ N_A
n_pbulk = NI_SI * np.exp(psi_L)
p_pbulk = NI_SI * np.exp(-psi_L)
n_nbulk = NI_SI * np.exp(psi_R)
p_nbulk = NI_SI * np.exp(-psi_R)
print(f"\np-bulk: n = {n_pbulk*1e-6:.2e} cm^-3, p = {p_pbulk*1e-6:.2e} cm^-3 "
      f"(expected p ~ N_A = 1e17)")
print(f"n-bulk: n = {n_nbulk*1e-6:.2e} cm^-3, p = {p_nbulk*1e-6:.2e} cm^-3 "
      f"(expected n ~ N_D = 1e17)")
# Check mass action: n*p = n_i^2
print(f"n*p (p-bulk) = {n_pbulk*p_pbulk:.3e} m^-6; n_i^2 = {NI_SI**2:.3e} m^-6 "
      f"{'OK' if abs(n_pbulk*p_pbulk - NI_SI**2)/NI_SI**2 < 0.01 else 'FAIL'}")

# Check neutrality: p - n + (N_D - N_A) = 0 in bulk
# In p-bulk: p - n - N_A should be 0
# Using psi_L = asinh(-N_A/(2ni)), exp(-psi_L) - exp(psi_L) = -2 sinh(psi_L) = -(-N_A/ni) = N_A/ni
# So p - n = ni * (exp(-psi_L) - exp(psi_L)) = N_A. Then p - n - N_A = 0. Check.
charge_p = p_pbulk - n_pbulk - N_A_left
print(f"Charge neutrality in p-bulk: {charge_p:.3e} m^-3 "
      f"(should be ~0; normalized: {abs(charge_p)/N_A_left:.3e}) "
      f"{'OK' if abs(charge_p)/N_A_left < 1e-10 else 'FAIL'}")

charge_n = p_nbulk - n_nbulk + N_D_right
print(f"Charge neutrality in n-bulk: {charge_n:.3e} m^-3 "
      f"{'OK' if abs(charge_n)/N_D_right < 1e-10 else 'FAIL'}")

print("\nAll sanity checks complete.")
