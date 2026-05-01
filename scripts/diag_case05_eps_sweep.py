"""
Diagnostic: bias_sweep dI/dV finite-difference convergence study at
Case 05's operating point (V_DC = 0.4 V, forward bias).

Hypothesis: the audit test's EPS_V = 0.05 V is too large for the
exponential I-V at forward bias (h/V_t = 0.05/0.026 = 1.92), so the
centered-difference dI/dV is biased by sinh(h/V_t)/(h/V_t) ~ 1.77
relative to the true derivative. If we shrink h, the FD dI/dV should
converge to a value close to the AC Re(Y(omega->0)) at V_DC = 0.4 V.

Run from inside the docker container or any environment with dolfinx:

    python scripts/diag_case05_eps_sweep.py
"""
from __future__ import annotations

import copy
import json
import os

from semi.runners.ac_sweep import run_ac_sweep
from semi.runners.bias_sweep import run_bias_sweep


def main() -> None:
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    cfg_path = os.path.join(repo_root, "benchmarks", "rc_ac_sweep", "rc_ac_sweep.json")
    with open(cfg_path) as fh:
        cfg = json.load(fh)

    V_DC = 0.4

    # AC at low omega.
    cfg_ac = copy.deepcopy(cfg)
    cfg_ac["solver"]["dc_bias"] = {"contact": "anode", "voltage": V_DC}
    cfg_ac["solver"]["ac"]["frequencies"] = {"type": "list", "values": [1.0]}
    for c in cfg_ac["contacts"]:
        if c["name"] == "anode":
            c["voltage"] = V_DC
    res_ac = run_ac_sweep(cfg_ac)
    Y0 = complex(res_ac.Y[0])
    G_ac = float(Y0.real)
    print(f"AC Re(Y) at V_DC={V_DC} V, f=1 Hz: {G_ac:.6e} S")
    print()

    def bs_at(V, step):
        cfg_bs = copy.deepcopy(cfg)
        snes = cfg_bs["solver"].get("snes", {})
        cont = cfg_bs["solver"].get("continuation", {})
        cfg_bs["solver"] = {"type": "bias_sweep", "snes": snes, "continuation": cont}
        for c in cfg_bs["contacts"]:
            if c["name"] == "anode":
                c["voltage_sweep"] = {"start": 0.0, "stop": V, "step": step}
                c["voltage"] = 0.0
        return run_bias_sweep(cfg_bs).iv[-1]["J"]

    print(f"{'EPS_V (V)':>12s}  {'step (V)':>10s}  {'dI/dV (S)':>14s}  {'rel_err':>10s}")
    for eps_v in (0.05, 0.025, 0.01, 0.005, 0.0025, 0.001):
        # Use a step that divides eps_v evenly so endpoints land exactly.
        if eps_v >= 0.025:
            step = 0.025
        elif eps_v >= 0.005:
            step = 0.005
        else:
            step = 0.001
        # Adjust step so endpoints fall on grid.
        # Round step to a divisor of eps_v.
        if eps_v % step != 0:
            step = eps_v
        Jp = bs_at(V_DC + eps_v, step)
        Jm = bs_at(V_DC - eps_v, step)
        dIdV = (Jp - Jm) / (2.0 * eps_v)
        rel = abs(G_ac - dIdV) / abs(dIdV) if dIdV != 0 else float("nan")
        print(f"  {eps_v:>10.4f}  {step:>10.4f}  {dIdV:>14.6e}  {rel:>10.4%}")


if __name__ == "__main__":
    main()
