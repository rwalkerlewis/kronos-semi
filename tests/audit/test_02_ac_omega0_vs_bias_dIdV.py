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
EPS_V = 1.0e-3


@pytest.mark.audit
def test_ac_omega0_vs_bias_dIdV():
    require_dolfinx()

    from semi.runners.ac_sweep import run_ac_sweep
    from semi.runners.bias_sweep import run_bias_sweep

    cfg_ac = load_benchmark("rc_ac_sweep/rc_ac_sweep.json")
    # Restrict the sweep to one low-omega frequency for speed.
    cfg_ac = copy.deepcopy(cfg_ac)
    cfg_ac["solver"]["ac"] = {"frequencies": [1.0]}

    ac_result = run_ac_sweep(cfg_ac)
    Y0 = complex(ac_result.Y[0])
    G_ac = float(Y0.real)

    # Bias sweep finite difference: dI/dV at V_DC.
    cfg_bs_minus = copy.deepcopy(cfg_ac)
    cfg_bs_minus["solver"] = {
        "type": "bias_sweep",
        "bias_ramp": {"start": 0.0, "stop": V_DC - EPS_V, "step": -0.05},
    }
    cfg_bs_minus["contacts"][0]["voltage"] = V_DC - EPS_V

    cfg_bs_plus = copy.deepcopy(cfg_ac)
    cfg_bs_plus["solver"] = {
        "type": "bias_sweep",
        "bias_ramp": {"start": 0.0, "stop": V_DC + EPS_V, "step": -0.05},
    }
    cfg_bs_plus["contacts"][0]["voltage"] = V_DC + EPS_V

    res_minus = run_bias_sweep(cfg_bs_minus)
    res_plus = run_bias_sweep(cfg_bs_plus)
    J_minus = res_minus.iv[-1]["J"]
    J_plus = res_plus.iv[-1]["J"]
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

    # Generous tolerance for first run; both terms are tiny in deep
    # reverse bias, so absolute scale matters more than relative.
    assert (rel_err < 0.5) or (abs(G_ac - dIdV_bs) < 1e-12), (
        f"Re(Y) at omega=0 disagrees with dI/dV: "
        f"Re(Y)={G_ac:.3e}, dI/dV={dIdV_bs:.3e}, rel_err={rel_err:.3e}"
    )
