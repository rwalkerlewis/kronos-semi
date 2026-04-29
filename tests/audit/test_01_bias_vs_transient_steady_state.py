"""
Audit case 1: bias_sweep vs transient at deep steady state.

Formalizes the M13.1 close-out result (PR #52) by sampling more than
one operating point. The transient runner is allowed to relax to deep
steady state by stepping until residuals stop changing; the result is
compared to `run_bias_sweep` at the same V_F.

The expectation set by M13.1 is "<= 1e-4 relative error at deep
steady state". This test reports the per-bias relative L2 of (psi, n,
p) and the relative IV-current disagreement.
"""
from __future__ import annotations

import copy

import numpy as np
import pytest

from ._helpers import (
    load_benchmark,
    relative_l2,
    require_dolfinx,
    write_csv,
    write_markdown,
)

CASE = "01_bias_vs_transient_steady_state"
BIASES = [0.3, 0.5]  # V_F values; deep SS for forward bias diode


@pytest.mark.audit
def test_bias_sweep_vs_transient_steady_state():
    require_dolfinx()

    from semi.runners.bias_sweep import run_bias_sweep
    from semi.runners.transient import run_transient

    cfg_base = load_benchmark("pn_1d_turnon/config.json")

    rows: list[list[float]] = []
    for V_F in BIASES:
        # Bias sweep solve at V_F.  Use `contacts[*].voltage_sweep` so
        # bias_sweep actually steps the swept contact (otherwise it
        # records J=0 at that contact, defeating the IV comparison).
        cfg_bs = copy.deepcopy(cfg_base)
        for c in cfg_bs["contacts"]:
            if c["name"] == cfg_base["contacts"][0]["name"]:
                c["voltage_sweep"] = {"start": 0.0, "stop": V_F, "step": 0.05}
                c["voltage"] = 0.0
        snes = cfg_bs["solver"].get("snes", {})
        cont = cfg_bs["solver"].get("continuation", {})
        cfg_bs["solver"] = {"type": "bias_sweep", "snes": snes, "continuation": cont}
        bs_result = run_bias_sweep(cfg_bs)

        # Transient with BC ramp + long t_end so the solution settles.
        cfg_tr = copy.deepcopy(cfg_base)
        cfg_tr["contacts"][0]["voltage"] = V_F
        cfg_tr["solver"]["bc_ramp_steps"] = 10
        cfg_tr["solver"]["t_end"] = 5.0e-6
        cfg_tr["solver"]["dt"] = 5.0e-9
        cfg_tr["solver"]["max_steps"] = 2000
        cfg_tr["solver"]["output_every"] = 10000  # only need final
        tr_result = run_transient(cfg_tr)

        # Pull last-snapshot fields.
        #
        # run_bias_sweep returns SimulationResult with .psi_phys / .n_phys /
        # .p_phys as plain numpy arrays (physical units, Volts and m^-3).
        #
        # run_transient returns TransientResult whose .fields dict maps
        # field names to lists of per-snapshot numpy arrays, also in
        # physical units.  We take the final snapshot (index -1).
        tr_fields = getattr(tr_result, "fields", {}) or {}

        try:
            bs_psi = np.asarray(bs_result.psi_phys)
            bs_n = np.asarray(bs_result.n_phys)
            bs_p = np.asarray(bs_result.p_phys)

            if not (tr_fields.get("psi") and len(tr_fields["psi"]) > 0):
                raise KeyError("psi not in transient fields or empty snapshot list")
            tr_psi = np.asarray(tr_fields["psi"][-1])
            tr_n = np.asarray(tr_fields["n"][-1])
            tr_p = np.asarray(tr_fields["p"][-1])
        except (KeyError, AttributeError) as exc:
            pytest.fail(
                f"{CASE}: required field missing from runner output: {exc}. "
                f"bs keys = {[a for a in ('psi_phys','n_phys','p_phys') if getattr(bs_result, a, None) is not None]}; "
                f"tr_fields keys = {list(tr_fields)}"
            )

        psi_err = relative_l2(tr_psi, bs_psi)
        n_err = relative_l2(tr_n, bs_n)
        p_err = relative_l2(tr_p, bs_p)

        # IV: take last current at the swept contact.
        bs_J = bs_result.iv[-1]["J"] if bs_result.iv else float("nan")
        tr_J = tr_result.iv[-1]["J"] if tr_result.iv else float("nan")
        if abs(bs_J) > 1e-300:
            iv_err = abs(tr_J - bs_J) / abs(bs_J)
        else:
            iv_err = abs(tr_J - bs_J)

        rows.append([V_F, psi_err, n_err, p_err, iv_err])

    write_csv(
        CASE,
        ["V_F", "rel_L2_psi", "rel_L2_n", "rel_L2_p", "rel_J"],
        rows,
    )

    summary_lines = ["| V_F (V) | rel_L2 psi | rel_L2 n | rel_L2 p | rel J |", "|---:|---:|---:|---:|---:|"]
    for V, e_psi, e_n, e_p, e_J in rows:
        summary_lines.append(
            f"| {V:.3f} | {e_psi:.3e} | {e_n:.3e} | {e_p:.3e} | {e_J:.3e} |"
        )

    write_markdown(
        CASE,
        "Case 01 - bias_sweep vs transient (deep steady state)",
        "Compares `run_bias_sweep` and `run_transient` (relaxed to deep "
        "steady state) on the pn_1d_turnon device at multiple forward "
        "biases. M13.1's close-out claim: agreement at 1e-4 relative.\n\n"
        + "\n".join(summary_lines)
        + "\n\nCSV: `/tmp/audit/" + CASE + ".csv`",
    )

    # Soft gate: this is a discovery audit, not a regression. Tolerance
    # is intentionally generous; tighten in a follow-up PR if all six
    # cases come in clean.
    for V, e_psi, e_n, e_p, e_J in rows:
        assert e_psi < 1e-2, f"psi disagreement {e_psi:.3e} at V_F={V}"
        assert e_n < 5e-2, f"n disagreement {e_n:.3e} at V_F={V}"
        assert e_p < 5e-2, f"p disagreement {e_p:.3e} at V_F={V}"
        # Now that bias_sweep is configured with a real voltage_sweep,
        # the swept-contact terminal current is non-zero and the IV
        # comparison is meaningful.  Gate at 5% relative.
        assert e_J < 5e-2, f"J disagreement {e_J:.3e} at V_F={V}"
