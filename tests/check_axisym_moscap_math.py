"""
Analytical sanity checks for the axisymmetric MOSCAP benchmark.

This script mirrors the style of `tests/check_day1_math.py`: pure
Python (numpy + scipy), no dolfinx. It verifies that the closed-form
parameters in ``semi.cv`` match Hu chapter 5 worked examples, that the
HF and LF C-V curves obey the textbook structural properties, and
that the chosen domain radius R is large enough relative to the
maximum depletion width.

Run as either ``pytest tests/check_axisym_moscap_math.py`` or
``python tests/check_axisym_moscap_math.py``; failed checks print
``FAIL`` and exit non-zero.
"""
import math
import sys

import numpy as np

from semi.constants import EPS0, Q
from semi.cv import (
    EPS_R_SI,
    EPS_R_SIO2,
    analytical_moscap_params,
    compute_hf_cv_depletion_clamp,
    compute_lf_cv_fem,
    hf_cv_depletion_approximation,
    lf_cv_quasistatic,
)

# ---- reference parameters (match benchmarks/moscap_axisym_2d) ----
N_A = 5.0e16          # cm^-3 p-Si body
T_OX = 10.0e-9        # m
PHI_MS = -0.95        # V (N+ poly on p-Si, Hu Table 5-1 ballpark)
T = 300.0
R_DOMAIN = 200.0e-6   # m
R_GATE = 50.0e-6      # m
T_SI = 5.0e-6         # m


def _ok(label: str, cond: bool, note: str = "") -> bool:
    tag = "OK  " if cond else "FAIL"
    extra = f" ({note})" if note else ""
    print(f"[{tag}] {label}{extra}")
    return cond


def main() -> int:
    failures = 0

    p = analytical_moscap_params(
        body_dopant="p",
        N_body_cm3=N_A,
        T_ox_m=T_OX,
        phi_ms=PHI_MS,
        T=T,
    )

    # ---- thermal voltage ----
    failures += not _ok(
        "thermal voltage 25.85 mV at 300 K",
        abs(p.Vt - 0.02585) < 5.0e-4,
        f"got {p.Vt*1000:.3f} mV",
    )

    # ---- |phi_B| for Na = 5e16, n_i = 1e10 cm^-3 ----
    # phi_B = Vt ln(Na/n_i) = 0.02585 * ln(5e6) ~ 0.400 V
    expected_phi_B = p.Vt * math.log(5.0e16 / 1.0e10)
    failures += not _ok(
        "phi_B ~ 0.40 V for Na=5e16, n_i=1e10",
        abs(p.phi_B - expected_phi_B) < 1.0e-3,
        f"phi_B = {p.phi_B*1000:.1f} mV",
    )

    # ---- C_ox for Tox = 10 nm SiO2 ----
    # eps_ox = 3.9 * 8.854e-12 = 3.453e-11; / 1e-8 m = 3.453e-3 F/m^2
    expected_C_ox = EPS_R_SIO2 * EPS0 / T_OX
    failures += not _ok(
        "C_ox = eps_ox / Tox",
        abs(p.C_ox_per_area - expected_C_ox) / expected_C_ox < 1.0e-9,
        f"{p.C_ox_per_area*1e3:.4f} mF/m^2",
    )

    # ---- W_dmax: 2 phi_B ~ 0.80 V => sqrt(2 eps_s phi_s / (q N)) ----
    # eps_s = 11.7 * 8.854e-12 = 1.036e-10
    # W = sqrt(2 * 1.036e-10 * 0.80 / (1.602e-19 * 5e22)) m
    #   = sqrt(1.658e-10 / 8.01e3) = sqrt(2.07e-14) = 1.44e-7 m = 144 nm
    eps_s = EPS_R_SI * EPS0
    N_body = N_A * 1.0e6  # m^-3
    W_expected = math.sqrt(2.0 * eps_s * 2.0 * p.phi_B / (Q * N_body))
    failures += not _ok(
        "W_dmax ~ 145 nm for Na=5e16",
        abs(p.W_dmax - W_expected) < 1.0e-9 and 1.0e-7 < p.W_dmax < 2.0e-7,
        f"W_dmax = {p.W_dmax*1e9:.1f} nm",
    )

    # ---- Threshold voltage signs ----
    # For p-body, V_t > V_fb and V_t > 0 typically (nMOS).
    failures += not _ok(
        "V_t > V_fb for p-body",
        p.V_t > p.V_fb,
        f"V_fb = {p.V_fb:.3f} V, V_t = {p.V_t:.3f} V",
    )

    # ---- C_min < C_ox; expect ~ 0.15 * C_ox for Na=5e16, Tox=10 nm ----
    Cmin_norm = p.C_min_per_area / p.C_ox_per_area
    failures += not _ok(
        "C_min/C_ox in [0.10, 0.30] for these parameters",
        0.10 < Cmin_norm < 0.30,
        f"C_min/C_ox = {Cmin_norm:.3f}",
    )

    # ---- Domain radius R >= 5 W_dmax ----
    failures += not _ok(
        "R_domain >= 5 W_dmax (radial truncation safe)",
        R_DOMAIN >= 5.0 * p.W_dmax,
        f"R = {R_DOMAIN*1e6:.0f} um, 5 W_dmax = {5*p.W_dmax*1e6:.2f} um",
    )

    # ---- LF and HF curve structural properties ----
    Vg = np.linspace(-2.0, 2.0, 81)
    C_HF_norm = hf_cv_depletion_approximation(Vg, p, N_body_cm3=N_A)
    C_LF_norm = lf_cv_quasistatic(Vg, p, N_body_cm3=N_A)

    # 1) Accumulation: both -> 1.0 within tolerance. LF uses FD on Q_s,
    #    so trim the first/last 2 indices to skip one-sided gradient
    #    endpoint effects.
    acc_idx = Vg < (p.V_fb - 0.5)
    failures += not _ok(
        "accumulation: HF ~ 1.0",
        np.max(np.abs(C_HF_norm[acc_idx] - 1.0)) < 0.01,
        f"max dev = {np.max(np.abs(C_HF_norm[acc_idx] - 1.0)):.4f}",
    )
    acc_indices = np.where(acc_idx)[0]
    acc_interior = acc_indices[(acc_indices >= 2) & (acc_indices <= len(Vg) - 3)]
    if acc_interior.size:
        # 10 % rather than 2 % because the analytic LF curve approaches
        # C_ox asymptotically; the residual phi_s contribution at
        # Vg = -2 V is still ~ 5 %. The FEM benchmark sweep extends
        # further into accumulation, where 2 % is achievable.
        failures += not _ok(
            "accumulation (interior FD): LF -> C_ox within 10%",
            np.max(np.abs(C_LF_norm[acc_interior] - 1.0)) < 0.10,
            f"max dev = {np.max(np.abs(C_LF_norm[acc_interior] - 1.0)):.4f}",
        )

    # 2) Strong inversion: HF -> Cmin/Cox within 2%
    inv_idx = Vg > (p.V_t + 0.5)
    if np.any(inv_idx):
        failures += not _ok(
            "strong inversion: HF -> C_min/C_ox within 2%",
            np.max(np.abs(C_HF_norm[inv_idx] - Cmin_norm)) < 0.02,
            f"max dev = {np.max(np.abs(C_HF_norm[inv_idx] - Cmin_norm)):.4f}",
        )
        # 3) Strong inversion: LF -> 1.0; trim FD edges so endpoint
        #    boundary effects of np.gradient do not skew the test.
        inv_indices = np.where(inv_idx)[0]
        interior = inv_indices[(inv_indices > 2) & (inv_indices < len(Vg) - 3)]
        if interior.size:
            # Same asymptotic-saturation rationale as in accumulation.
            failures += not _ok(
                "strong inversion: LF -> C_ox within 10%",
                np.max(np.abs(C_LF_norm[interior] - 1.0)) < 0.10,
                f"max dev = {np.max(np.abs(C_LF_norm[interior] - 1.0)):.4f}",
            )

    # 4) Depletion (V_fb < Vg < V_t): LF and HF coincide within 5%
    dep_idx = (Vg > p.V_fb + 0.05) & (Vg < p.V_t - 0.05)
    if np.any(dep_idx):
        failures += not _ok(
            "depletion: |LF - HF| < 0.05",
            np.max(np.abs(C_LF_norm[dep_idx] - C_HF_norm[dep_idx])) < 0.05,
            f"max diff = {np.max(np.abs(C_LF_norm[dep_idx] - C_HF_norm[dep_idx])):.4f}",
        )

    # ---- compute_hf_cv_depletion_clamp postprocessor sanity ----
    phi_s_synth = np.linspace(-0.3, 1.2, 31)  # synthetic surface potentials
    Vg_synth = phi_s_synth + p.V_fb  # not physical, just shape sanity
    C_hf_post = compute_hf_cv_depletion_clamp(Vg_synth, phi_s_synth, p, N_body_cm3=N_A)
    failures += not _ok(
        "FEM-postprocess HF: monotone non-increasing in depletion regime",
        np.all(np.diff(C_hf_post[(phi_s_synth > 0.0) & (phi_s_synth < 2 * p.phi_B)]) <= 1e-9),
        "",
    )
    failures += not _ok(
        "FEM-postprocess HF: clamps at C_min once phi_s >= 2 phi_B",
        abs(C_hf_post[-1] - p.C_min_per_area) < 0.01 * p.C_ox_per_area,
        f"C[-1] = {C_hf_post[-1]*1e3:.4f}, C_min = {p.C_min_per_area*1e3:.4f}",
    )

    # ---- compute_lf_cv_fem postprocessor: linear Q_s -> constant C ----
    Vg_lin = np.linspace(-1.0, 1.0, 11)
    Q_s_lin = -1.0e-3 * Vg_lin  # Q_s = -slope * Vg => C = slope
    C_lf_lin = compute_lf_cv_fem(Vg_lin, Q_s_lin)
    failures += not _ok(
        "FEM-postprocess LF: -dQ_s/dVg recovers slope",
        np.max(np.abs(C_lf_lin - 1.0e-3)) < 1.0e-9,
        "",
    )

    print()
    if failures:
        print(f"FAILED: {failures} check(s) failed.")
        return 1
    print("All analytical checks passed.")
    return 0


def test_axisym_moscap_math():
    """pytest entry point."""
    assert main() == 0


if __name__ == "__main__":
    sys.exit(main())
