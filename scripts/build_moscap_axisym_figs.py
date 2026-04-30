#!/usr/bin/env python3
"""
Generate reference figures and the canonical LF/HF C-V CSV for the
M14.2 axisymmetric MOSCAP benchmark.

Run inside the kronos-semi-test container::

    python scripts/build_moscap_axisym_figs.py

Outputs land under ``benchmarks/moscap_axisym/figs/`` and
``benchmarks/moscap_axisym/expected/`` in the repository tree.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from semi import run as semi_run
from semi import schema
from semi.physics.cv import analytical_moscap_metrics
from semi.constants import cm3_to_m3

REPO = Path(__file__).resolve().parents[1]
BENCH = REPO / "benchmarks" / "moscap_axisym"
FIGS = BENCH / "figs"
EXPECTED = BENCH / "expected"


def main() -> int:
    FIGS.mkdir(parents=True, exist_ok=True)
    EXPECTED.mkdir(parents=True, exist_ok=True)

    cfg = schema.validate(json.loads((BENCH / "moscap_axisym.json").read_text()))
    result = semi_run.run(cfg)

    V = np.array([r["V"] for r in result.iv])
    Q = np.array([r["Q_gate"] for r in result.iv])
    C_LF = np.array([r["C_LF"] for r in result.iv])
    C_HF = np.array([r["C_HF"] for r in result.iv])

    targets = analytical_moscap_metrics(
        N_body=cm3_to_m3(1.0e17), body_type="p", t_ox=5.0e-9,
        eps_r_si=11.7, eps_r_ox=3.9,
        n_i=cm3_to_m3(1.0e10), T=300.0, phi_ms=-0.977,
    )
    Cox, Cmin, V_FB, V_T = (
        targets["Cox"], targets["Cmin"], targets["V_FB"], targets["V_T"],
    )

    # 1. C-V curve normalised by Cox
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(V, C_LF / Cox, "C0-",  label="LF (quasi-static)", lw=2)
    ax.plot(V, C_HF / Cox, "C3--", label="HF",                lw=2)
    ax.axhline(Cmin / Cox, color="grey", ls=":", lw=1, label=f"$C_{{min}}/C_{{ox}}={Cmin/Cox:.3f}$")
    ax.axvline(V_FB, color="k", ls="--", lw=0.8, alpha=0.6)
    ax.axvline(V_T,  color="k", ls="--", lw=0.8, alpha=0.6)
    ax.text(V_FB, 0.05, "$V_{FB}$", ha="center")
    ax.text(V_T,  0.05, "$V_T$",    ha="center")
    ax.set_xlabel("$V_{gate}$ [V]")
    ax.set_ylabel("$C / C_{ox}$")
    ax.set_title("Axisymmetric MOSCAP, NA=1e17 cm$^{-3}$, t$_{ox}$=5 nm (Hu Fig. 5-18)")
    ax.set_ylim(0, 1.1)
    ax.grid(alpha=0.3)
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(FIGS / "cv_lf_hf.png", dpi=140)
    plt.close(fig)

    # 2. Q(V)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(V, Q, "C2-")
    ax.set_xlabel("$V_{gate}$ [V]")
    ax.set_ylabel("$Q_{gate}$ [C/m$^2$]")
    ax.axvline(V_FB, color="k", ls="--", lw=0.8, alpha=0.6)
    ax.axvline(V_T,  color="k", ls="--", lw=0.8, alpha=0.6)
    ax.set_title("Gate charge per unit area")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIGS / "qv.png", dpi=140)
    plt.close(fig)

    # 3. Reference CSV
    with (EXPECTED / "cv_curve_LF_HF.csv").open("w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow([
            "V_gate", "Q_gate_per_area", "C_LF", "C_HF",
            "C_LF_over_Cox", "C_HF_over_Cox",
        ])
        for v, q, cl, ch in zip(V, Q, C_LF, C_HF):
            wr.writerow([f"{v:+.6f}", f"{q:+.6e}", f"{cl:.6e}",
                         f"{ch:.6e}", f"{cl/Cox:.6f}", f"{ch/Cox:.6f}"])

    # 4. Analytical reference values
    with (EXPECTED / "analytical_targets.json").open("w") as fh:
        json.dump(
            {
                "Cox": Cox, "Cmin": Cmin, "Cdep_min": targets["Cdep_min"],
                "W_dmax": targets["W_dmax"], "phi_F": targets["phi_F"],
                "V_FB": V_FB, "V_T": V_T, "Q_dep_max": targets["Q_dep_max"],
            },
            fh, indent=2, sort_keys=True,
        )

    print("wrote:", FIGS / "cv_lf_hf.png")
    print("wrote:", FIGS / "qv.png")
    print("wrote:", EXPECTED / "cv_curve_LF_HF.csv")
    print("wrote:", EXPECTED / "analytical_targets.json")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
