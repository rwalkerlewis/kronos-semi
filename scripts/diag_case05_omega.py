"""
Diagnostic: AC sweep at Case 05's operating point (V_DC = 0.4 V) over a
range of low frequencies, to see whether Re(Y) is omega-independent
(displacement current is irrelevant to the audit) or omega-dependent
(displacement current contributes a real component at omega = 1 Hz).

Hypothesis: if Re(Y) is essentially constant from 1 Hz down to 1e-3 Hz,
the ~12% Case 05 disagreement with bias_sweep dI/dV cannot be blamed on
displacement current at omega = 1 Hz.  The next hypothesis -- the
recombination linearisation or the operating-point/Jacobian consistency --
is then the prime suspect.

Run from inside the docker container or any environment with dolfinx:

    python scripts/diag_case05_omega.py
"""
from __future__ import annotations

import copy
import json
import os

from semi.runners.ac_sweep import run_ac_sweep


def main() -> None:
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    cfg_path = os.path.join(repo_root, "benchmarks", "rc_ac_sweep", "rc_ac_sweep.json")
    with open(cfg_path) as fh:
        cfg = json.load(fh)

    V_DC = 0.4
    omegas_hz = [1.0, 0.1, 0.01, 0.001]
    rows = []
    for f in omegas_hz:
        cfg_ac = copy.deepcopy(cfg)
        cfg_ac["solver"]["dc_bias"] = {"contact": "anode", "voltage": V_DC}
        cfg_ac["solver"]["ac"]["frequencies"] = {"type": "list", "values": [f]}
        for c in cfg_ac["contacts"]:
            if c["name"] == "anode":
                c["voltage"] = V_DC
        result = run_ac_sweep(cfg_ac)
        Y = complex(result.Y[0])
        rows.append((f, Y.real, Y.imag))
        print(f"f = {f:8.3e} Hz   Re(Y) = {Y.real: .8e}   Im(Y) = {Y.imag: .8e}")

    if rows:
        re0 = rows[0][1]
        if re0 != 0.0:
            print()
            print("Relative drift of Re(Y) from f=1Hz baseline:")
            for f, re, _ in rows:
                drift = abs(re - re0) / abs(re0)
                print(f"  f = {f:8.3e} Hz   drift = {drift:.3e}")


if __name__ == "__main__":
    main()
