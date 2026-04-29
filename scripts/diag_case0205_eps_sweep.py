"""
Diagnostic: bias_sweep dI/dV finite-difference convergence study at
both Case 02 (V_DC = -1.0 V, reverse bias) and Case 05 (V_DC = +0.4 V,
forward bias) operating points.

After the AC ac_sweep was rewritten in Slotboom variables (audit
Cases 02/05 closure), the AC Re(Y) at low omega is the answer the
discrete drift-diffusion problem would converge to in the limit of
infinitesimal bias_sweep step. This script runs ``bias_sweep`` at
several EPS_V values to find the converged dI/dV and compares it
against ``ac_sweep`` Re(Y) at f = 1 Hz.

Usage:

    python scripts/diag_case0205_eps_sweep.py
"""
from __future__ import annotations

import copy
import json
import os

from semi.runners.ac_sweep import run_ac_sweep
from semi.runners.bias_sweep import run_bias_sweep


def _ac_at(cfg, V_DC):
    cfg_ac = copy.deepcopy(cfg)
    cfg_ac["solver"]["dc_bias"] = {"contact": "anode", "voltage": V_DC}
    cfg_ac["solver"]["ac"]["frequencies"] = {"type": "list", "values": [1.0]}
    for c in cfg_ac["contacts"]:
        if c["name"] == "anode":
            c["voltage"] = V_DC
    return float(complex(run_ac_sweep(cfg_ac).Y[0]).real)


def _bs_J(cfg, V, step):
    cfg_bs = copy.deepcopy(cfg)
    snes = cfg_bs["solver"].get("snes", {})
    cont = cfg_bs["solver"].get("continuation", {})
    cfg_bs["solver"] = {"type": "bias_sweep", "snes": snes, "continuation": cont}
    for c in cfg_bs["contacts"]:
        if c["name"] == "anode":
            c["voltage_sweep"] = {"start": 0.0, "stop": V, "step": step}
            c["voltage"] = 0.0
    return run_bias_sweep(cfg_bs).iv[-1]["J"]


def main() -> None:
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    cfg_path = os.path.join(repo_root, "benchmarks", "rc_ac_sweep", "rc_ac_sweep.json")
    with open(cfg_path) as fh:
        cfg = json.load(fh)

    for label, V_DC, eps_list in [
        ("Case 02 (reverse, V_DC = -1.0 V)", -1.0,
         [0.05, 0.025, 0.01, 0.005, 0.0025, 0.001]),
        ("Case 05 (forward, V_DC = +0.4 V)",  0.4,
         [0.05, 0.025, 0.01, 0.005, 0.0025, 0.001]),
    ]:
        print(f"=== {label} ===")
        G_ac = _ac_at(cfg, V_DC)
        print(f"AC Re(Y) at f=1 Hz: {G_ac:.6e} S")
        print(f"{'EPS_V (V)':>12s}  {'step (V)':>10s}  {'dI/dV (S)':>14s}  {'rel_err':>10s}")
        for eps_v in eps_list:
            step = eps_v
            Jp = _bs_J(cfg, V_DC + eps_v, step)
            Jm = _bs_J(cfg, V_DC - eps_v, step)
            dIdV = (Jp - Jm) / (2.0 * eps_v)
            rel = abs(G_ac - dIdV) / abs(dIdV) if dIdV != 0 else float("nan")
            print(f"  {eps_v:>10.4f}  {step:>10.4f}  {dIdV:>14.6e}  {rel:>10.4%}")
        print()


if __name__ == "__main__":
    main()
