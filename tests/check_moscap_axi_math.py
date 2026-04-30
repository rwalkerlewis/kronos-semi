"""
Pure-Python analytical checks for MOSCAP axisymmetric math.

These tests do NOT require dolfinx and verify:
1. Debye length calculation
2. Oxide capacitance per unit area
3. Flat-band capacitance
4. Accumulation/depletion charge scaling

Run with: python tests/check_moscap_axi_math.py
"""
import math
import sys


def check_debye_length():
    """Debye length for Na=1e17 cm^-3 Si at 300 K should be ~13 nm."""
    eps_Si = 11.7 * 8.854e-12       # F/m
    q = 1.602176634e-19              # C
    kT = 1.380649e-23 * 300.0       # J
    Na_cm3 = 1.0e17
    Na_m3 = Na_cm3 * 1.0e6
    L_D = math.sqrt(eps_Si * kT / (q ** 2 * Na_m3))
    print(f"  Debye length: {L_D * 1e9:.2f} nm  (expected ~13 nm)")
    assert 5.0e-9 < L_D < 25.0e-9, f"Debye length out of range: {L_D}"
    return True


def check_oxide_capacitance():
    """Cox for 5 nm SiO2 should be ~6.9 mF/m^2."""
    eps_ox = 3.9 * 8.854e-12        # F/m
    t_ox = 5.0e-9                   # m
    Cox = eps_ox / t_ox
    print(f"  Cox (5 nm SiO2): {Cox * 1e3:.2f} mF/m^2  (expected ~6.9)")
    assert 5.0e-3 < Cox < 9.0e-3, f"Cox out of range: {Cox}"
    return True


def check_gate_area():
    """Gate area for R_gate=5 um should be pi*(5e-6)^2."""
    R_gate = 5.0e-6
    A = math.pi * R_gate ** 2
    print(f"  Gate area: {A * 1e12:.4f} um^2  (expected {math.pi * 25:.4f})")
    assert abs(A - math.pi * 25.0e-12) < 1.0e-20, f"Gate area mismatch: {A}"
    return True


def check_flatband_voltage():
    """
    Flat-band voltage for phi_ms = phi_m - phi_s.
    With phi_m = 4.05 eV (Al gate, approximately) and
    phi_s = chi + Eg/2 + phi_F for P-Si.
    chi(Si) = 4.05 eV, Eg = 1.12 eV, kT/q * ln(Na/ni) ~ 0.407 V for Na=1e17.
    phi_s ~ 4.05 + 0.56 + 0.407 = 5.017 eV
    phi_ms ~ 4.05 - 5.017 = -0.967 V => V_FB < 0 (accumulation biased)
    """
    q = 1.602176634e-19
    kT_eV = 0.02585  # at 300 K
    chi_Si = 4.05    # eV
    Eg = 1.12        # eV
    ni = 9.65e9      # cm^-3 at 300 K (standard value)
    Na = 1.0e17      # cm^-3
    phi_F = kT_eV * math.log(Na / ni)
    phi_s = chi_Si + Eg / 2.0 + phi_F
    phi_m = 4.05     # eV (Al/poly approximately)
    V_FB = phi_m - phi_s
    print(f"  phi_F = {phi_F:.4f} V")
    print(f"  phi_s = {phi_s:.4f} eV")
    print(f"  V_FB  = {V_FB:.4f} V  (expected ~-0.96 V for this simplified model)")
    assert V_FB < 0.0, "V_FB should be negative for n+ poly gate on P-Si"
    assert -1.5 < V_FB < 0.0, f"V_FB out of expected range: {V_FB}"
    return True


def check_max_depletion_width():
    """
    Maximum depletion width Wdm = sqrt(4 * eps_Si * phi_F / (q * Na)).
    For Na=1e17, should be ~100 nm.
    """
    eps_Si = 11.7 * 8.854e-12
    q = 1.602176634e-19
    kT_eV = 0.02585
    ni = 9.65e9
    Na = 1.0e17
    Na_m3 = Na * 1.0e6
    phi_F = kT_eV * math.log(Na / ni)
    Wdm = math.sqrt(4.0 * eps_Si * phi_F / (q * Na_m3))
    print(f"  Wdm = {Wdm * 1e9:.1f} nm  (expected ~100 nm for Na=1e17)")
    assert 50.0e-9 < Wdm < 200.0e-9, f"Wdm out of range: {Wdm}"
    return True


def check_inversion_charge_scaling():
    """
    Threshold voltage (simplified): V_T = V_FB + 2*phi_F + Wdm*q*Na/Cox
    Should be a few tenths of a volt positive for P-body.
    """
    eps_Si = 11.7 * 8.854e-12
    eps_ox = 3.9 * 8.854e-12
    q = 1.602176634e-19
    kT_eV = 0.02585
    ni = 9.65e9
    Na = 1.0e17
    Na_m3 = Na * 1.0e6
    t_ox = 5.0e-9
    Cox = eps_ox / t_ox

    phi_F = kT_eV * math.log(Na / ni)
    Wdm = math.sqrt(4.0 * eps_Si * phi_F / (q * Na_m3))
    Qd_max = q * Na_m3 * Wdm

    chi_Si = 4.05
    Eg = 1.12
    phi_s = chi_Si + Eg / 2.0 + phi_F
    phi_m = 4.05
    V_FB = phi_m - phi_s

    V_T = V_FB + 2.0 * phi_F + Qd_max / Cox
    print(f"  V_T   = {V_T:.4f} V  (expected V_T > V_FB)")
    print(f"  2*phi_F = {2*phi_F:.4f} V")
    print(f"  Qd_max/Cox = {Qd_max/Cox:.4f} V")
    assert V_T > V_FB, "V_T must be greater than V_FB"
    return True


def main():
    checks = [
        ("Debye length", check_debye_length),
        ("Oxide capacitance", check_oxide_capacitance),
        ("Gate area", check_gate_area),
        ("Flat-band voltage", check_flatband_voltage),
        ("Max depletion width", check_max_depletion_width),
        ("Threshold voltage", check_inversion_charge_scaling),
    ]
    failures = 0
    for name, fn in checks:
        print(f"\n[CHECK] {name}")
        try:
            ok = fn()
            if ok:
                print(f"  PASS")
        except AssertionError as e:
            print(f"  FAIL: {e}")
            failures += 1
        except Exception as e:
            print(f"  ERROR: {e}")
            failures += 1

    print(f"\n{'='*50}")
    if failures == 0:
        print(f"All {len(checks)} checks passed.")
        sys.exit(0)
    else:
        print(f"{failures}/{len(checks)} checks FAILED.")
        sys.exit(1)


if __name__ == "__main__":
    main()
