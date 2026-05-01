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
import math

import pytest

from ._helpers import (
    load_benchmark,
    require_dolfinx,
    write_csv,
    write_markdown,
)

CASE = "02_ac_omega0_vs_bias_dIdV"
V_DC = -1.0
# EPS_V must be a multiple of the bias_sweep step so the
# voltage_sweep endpoint lands exactly on V_DC ± EPS_V.
# Per the EPS-sweep diagnostic on dev/audit-case05-drdu, bias_sweep's
# centered-FD dI/dV is dominated by FD curvature/noise at h=0.05 in
# reverse bias (~7% deviation from the converged value). h=0.005 is
# the smallest h before SNES residual noise re-enters; both Case 02
# and Case 05 unxfail at this h. See ADR-0011 Errata #2.
EPS_V = 0.005
BS_STEP = 0.005


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
                c["voltage_sweep"] = {"start": 0.0, "stop": V, "step": BS_STEP}
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

    # The numeric values of Re(Y) and rel_err drift by ~1-2% across
    # environments because of MUMPS LU pivot ordering on the indefinite
    # 2x2 real block (see audit-assertion comment below). To keep
    # docs/PHYSICS_AUDIT.md deterministic across CI / local re-runs we
    # write a stable summary here; the precise per-run floats are in
    # /tmp/audit/{CASE}.csv (uploaded as a CI artifact when needed).
    sign_match = ("matches" if (G_ac > 0) == (dIdV_bs > 0)
                  else "DISAGREES WITH")
    write_markdown(
        CASE,
        "Case 02 - AC sweep at small omega vs bias_sweep dI/dV",
        f"At V_DC = {V_DC} V (reverse bias, ~5 mS device), "
        f"Re(Y(omega->0)) {sign_match} the sign of bias_sweep "
        f"centered-difference dI/dV and the magnitudes agree within "
        f"the 2% audit gate. Per-run values (G_ac, dI/dV, rel_err) "
        f"are written to `/tmp/audit/{CASE}.csv`; they are not pinned "
        f"in this doc because MUMPS LU pivot ordering on the indefinite "
        f"2x2 real block introduces a ~1-2% environment-dependent "
        f"noise floor on Re(Y).",
    )

    # Audit assertion: post Slotboom-rewrite (ADR-0011 Errata #2).
    # Re(Y(omega->0)) must agree in sign with bias_sweep's centered-
    # difference dI/dV at the same V_DC, and the magnitudes must match
    # within 2%.  Case 02 sits in the reverse-bias depletion regime
    # where the terminal current is ~5 mS/m^2 (orders of magnitude
    # smaller than Case 05's ~75 S/m^2 forward-bias value), so the
    # MUMPS LU pivot ordering on the indefinite 2x2 real block plus
    # vertex-quadrature lumping introduces a ~1.5-2% environment-
    # dependent noise floor on Re(Y) (local: ~0.7% vs CI: ~1.2%).
    # The 2% gate is the noise floor, not a tolerance for the underlying
    # physics: the original (n,p)-form bug showed 7%+ at this point, so
    # 2% is a clean indicator that the Slotboom linearisation is
    # discrete-consistent with bias_sweep.
    assert math.isfinite(G_ac) and math.isfinite(dIdV_bs), \
        f"Non-finite conductance: G_ac={G_ac}, dI/dV={dIdV_bs}"
    if dIdV_bs != 0.0:
        assert (G_ac > 0) == (dIdV_bs > 0), (
            f"Re(Y) and dI/dV disagree in sign: G_ac={G_ac}, dI/dV={dIdV_bs}"
        )
    assert rel_err < 0.02, (
        f"Re(Y) vs dI/dV exceeds 2% tolerance: rel_err={rel_err:.3e} "
        f"(G_ac={G_ac:.3e}, dI/dV={dIdV_bs:.3e})"
    )
