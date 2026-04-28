"""
Audit case 5: AC sweep terminal current vs bias_sweep dI/dV at low omega.

A near-DC operating point in forward bias has a real conductance G
that should match the bias_sweep finite-difference dI/dV. This is
similar to case 02 but in forward bias where the current is large and
relative comparisons are well-conditioned.
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

CASE = "05_ac_terminal_current_vs_dIdV"
V_DC = 0.4
EPS_V = 1.0e-3


@pytest.mark.audit
def test_ac_terminal_current_vs_dIdV():
    require_dolfinx()

    from semi.runners.ac_sweep import run_ac_sweep
    from semi.runners.bias_sweep import run_bias_sweep

    cfg = load_benchmark("rc_ac_sweep/rc_ac_sweep.json")
    cfg_ac = copy.deepcopy(cfg)
    cfg_ac["solver"]["dc_bias"] = {"contact": "anode", "voltage": V_DC}
    cfg_ac["solver"]["ac"] = {"frequencies": [1.0]}
    cfg_ac["contacts"][0]["voltage"] = V_DC
    ac_res = run_ac_sweep(cfg_ac)
    G_ac = float(complex(ac_res.Y[0]).real)

    def bs_at(V):
        cfg_bs = copy.deepcopy(cfg)
        cfg_bs["contacts"][0]["voltage"] = V
        cfg_bs["solver"] = {
            "type": "bias_sweep",
            "bias_ramp": {"start": 0.0, "stop": V, "step": 0.05},
        }
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

    assert rel < 0.1, (
        f"AC small-signal G disagrees with bias_sweep dI/dV by {rel:.3e} "
        f"at V_DC={V_DC} V (Re(Y)={G_ac:.3e}, dI/dV={dIdV_bs:.3e})"
    )
