"""
Audit case 5: AC sweep terminal current vs bias_sweep dI/dV at low omega.

A near-DC operating point in forward bias has a real conductance G
that should match the bias_sweep finite-difference dI/dV. This is
similar to case 02 but in forward bias where the current is large and
relative comparisons are well-conditioned.
"""
from __future__ import annotations

import copy
import math

import pytest

from ._helpers import (
    load_benchmark,
    require_dolfinx,
    write_csv,
    write_markdown,
)

CASE = "05_ac_terminal_current_vs_dIdV"
V_DC = 0.4
# Use h=0.005 V where bias_sweep centered-FD dI/dV has converged
# below 5% residual deviation from the AC linearisation (the prior
# h=0.05 V exhibited >25% FD curvature error at this forward-bias
# operating point). See ADR-0011 Errata #2.
EPS_V = 0.005
BS_STEP = 0.005


@pytest.mark.audit
def test_ac_terminal_current_vs_dIdV():
    require_dolfinx()

    from semi.runners.ac_sweep import run_ac_sweep
    from semi.runners.bias_sweep import run_bias_sweep

    cfg = load_benchmark("rc_ac_sweep/rc_ac_sweep.json")
    cfg_ac = copy.deepcopy(cfg)
    cfg_ac["solver"]["dc_bias"] = {"contact": "anode", "voltage": V_DC}
    cfg_ac["solver"]["ac"]["frequencies"] = {"type": "list", "values": [1.0]}
    cfg_ac["contacts"][0]["voltage"] = V_DC
    ac_res = run_ac_sweep(cfg_ac)
    G_ac = float(complex(ac_res.Y[0]).real)

    # Use voltage_sweep in contacts so bias_sweep picks up the sweep
    # contact and correctly tracks J(V).  EPS_V must be a multiple of
    # the step so the endpoint lands exactly.
    def bs_at(V):
        cfg_bs = copy.deepcopy(cfg)
        snes = cfg_bs["solver"].get("snes", {})
        cont = cfg_bs["solver"].get("continuation", {})
        cfg_bs["solver"] = {"type": "bias_sweep", "snes": snes, "continuation": cont}
        for c in cfg_bs["contacts"]:
            if c["name"] == "anode":
                c["voltage_sweep"] = {"start": 0.0, "stop": V, "step": BS_STEP}
                c["voltage"] = 0.0
        return run_bias_sweep(cfg_bs).iv[-1]["J"]

    J_minus = bs_at(V_DC - EPS_V)
    J_plus = bs_at(V_DC + EPS_V)
    dIdV_bs = (J_plus - J_minus) / (2.0 * EPS_V)

    if abs(dIdV_bs) > 1e-300:
        rel = abs(G_ac - dIdV_bs) / abs(dIdV_bs)
    else:
        rel = abs(G_ac - dIdV_bs)

    write_csv(
        CASE,
        ["V_DC", "G_ac", "dIdV_bs", "rel_err"],
        [[V_DC, G_ac, dIdV_bs, rel]],
    )
    write_markdown(
        CASE,
        "Case 05 - AC terminal current vs bias_sweep dI/dV (forward bias)",
        f"At V_DC = {V_DC} V (forward bias, finite current), "
        f"AC Re(Y) at 1 Hz = {G_ac:.3e} S; bias_sweep dI/dV "
        f"= {dIdV_bs:.3e} S; relative error {rel:.3e}.\n\n"
        f"CSV: `/tmp/audit/{CASE}.csv`",
    )

    # Audit assertion: post sign-convention fix.  Re(Y) at low frequency
    # must agree in sign with bias_sweep dI/dV; relative error within
    # 5% (forward bias is more nonlinear than reverse, so a looser gate
    # than Case 02's 1% is appropriate).
    assert math.isfinite(G_ac) and math.isfinite(dIdV_bs), \
        f"Non-finite conductance: G_ac={G_ac}, dI/dV={dIdV_bs}"
    if dIdV_bs != 0.0:
        assert (G_ac > 0) == (dIdV_bs > 0), (
            f"Re(Y) and dI/dV disagree in sign: G_ac={G_ac}, dI/dV={dIdV_bs}"
        )
    assert rel < 0.05, (
        f"Re(Y) vs dI/dV exceeds 5% tolerance: rel_err={rel:.3e} "
        f"(G_ac={G_ac:.3e}, dI/dV={dIdV_bs:.3e})"
    )
