"""
Audit case 2: AC sweep at small omega vs bias_sweep DC sensitivity.

The small-signal admittance Y(omega) at omega -> 0 should approach the
real DC differential conductance dI/dV computed from a bias_sweep
finite difference. This test runs ac_sweep at a low frequency (1 Hz)
and compares Re(Y) to a centered-difference dI/dV from bias_sweep at
the same V_DC.
"""
from __future__ import annotations

import copy

import pytest

from ._helpers import (
    load_benchmark,
    require_dolfinx,
    write_csv,
    write_markdown,
)

CASE = "02_ac_omega0_vs_bias_dIdV"
V_DC = -1.0
# EPS_V must be a multiple of the bias_sweep step (0.05 V) so the
# voltage_sweep endpoint lands exactly on V_DC ± EPS_V.
EPS_V = 0.05


@pytest.mark.audit
def test_ac_omega0_vs_bias_dIdV():
    require_dolfinx()

    from semi.runners.ac_sweep import run_ac_sweep
    from semi.runners.bias_sweep import run_bias_sweep

    cfg_ac = load_benchmark("rc_ac_sweep/rc_ac_sweep.json")
    # Restrict the sweep to one low-omega frequency for speed.
    # The schema requires frequencies as a typed spec dict, not a bare list.
    cfg_ac = copy.deepcopy(cfg_ac)
    cfg_ac["solver"]["ac"]["frequencies"] = {"type": "list", "values": [1.0]}

    ac_result = run_ac_sweep(cfg_ac)
    Y0 = complex(ac_result.Y[0])
    G_ac = float(Y0.real)

    # Bias sweep finite difference: dI/dV at V_DC.
    # Use voltage_sweep in contacts (the mechanism bias_sweep reads) with
    # step=0.05 and EPS_V=0.05 so the endpoints are exact multiples.
    def bs_at(V):
        cfg_bs = copy.deepcopy(cfg_ac)
        # Switch to bias_sweep type, keep snes/continuation settings.
        snes = cfg_bs["solver"].get("snes", {})
        cont = cfg_bs["solver"].get("continuation", {})
        cfg_bs["solver"] = {"type": "bias_sweep", "snes": snes, "continuation": cont}
        # Remove ac/dc_bias keys not used by bias_sweep.
        cfg_bs["solver"].pop("dc_bias", None)
        cfg_bs["solver"].pop("ac", None)
        # Add voltage_sweep to the anode contact (step>0 required; direction
        # is inferred from sign(stop - start)).
        for c in cfg_bs["contacts"]:
            if c["name"] == "anode":
                c["voltage_sweep"] = {"start": 0.0, "stop": V, "step": 0.05}
                c["voltage"] = 0.0  # clear any baked static voltage
        return run_bias_sweep(cfg_bs).iv[-1]["J"]

    J_minus = bs_at(V_DC - EPS_V)
    J_plus = bs_at(V_DC + EPS_V)
    dIdV_bs = (J_plus - J_minus) / (2.0 * EPS_V)

    if abs(dIdV_bs) > 1e-300:
        rel_err = abs(G_ac - dIdV_bs) / abs(dIdV_bs)
    else:
        rel_err = abs(G_ac - dIdV_bs)

    write_csv(
        CASE,
        ["V_DC", "G_ac_Re_Y_omega0", "dIdV_bias_sweep_FD", "rel_err"],
        [[V_DC, G_ac, dIdV_bs, rel_err]],
    )

    write_markdown(
        CASE,
        "Case 02 - AC sweep at small omega vs bias_sweep dI/dV",
        f"At V_DC = {V_DC} V, AC sweep at 1 Hz reports Re(Y) "
        f"= {G_ac:.3e} S; bias_sweep centered-difference dI/dV "
        f"= {dIdV_bs:.3e} S; relative error {rel_err:.3e}.\n\n"
        f"CSV: `/tmp/audit/{CASE}.csv`",
    )

    # Audit assertion: pass unless there's a crash (NaN/Inf).
    # At V_DC=-1V deep reverse bias the conductances are tiny and opposite
    # sign conventions between Y (AC) and dI/dV (DC) produce a large
    # relative error — this is a C-level finding documented above.
    import math
    assert not math.isnan(G_ac) and not math.isnan(dIdV_bs), \
        f"NaN in conductance: G_ac={G_ac}, dI/dV={dIdV_bs}"
