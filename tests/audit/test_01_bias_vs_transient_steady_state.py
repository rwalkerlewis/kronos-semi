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
    final_field,
    load_benchmark,
    make_bias_sweep_cfg,
    relative_l2,
    require_dolfinx,
    write_csv,
    write_markdown,
)

CASE = "01_bias_vs_transient_steady_state"
BIASES = [0.1, 0.3, 0.5]  # V_F values; deep SS for forward bias diode


@pytest.mark.audit
def test_bias_sweep_vs_transient_steady_state():
    require_dolfinx()

    from semi.runners.bias_sweep import run_bias_sweep
    from semi.runners.transient import run_transient

    cfg_base = load_benchmark("pn_1d_turnon/config.json")

    rows: list[list[float]] = []
    for V_F in BIASES:
        # Bias sweep solve at V_F. The pn_1d_turnon config is authored
        # for the transient runner; reconfigure it for bias_sweep using
        # the helper (rewrites the swept contact's voltage_sweep and
        # injects M12-relaxed SNES tols where needed).
        cfg_bs = make_bias_sweep_cfg(cfg_base, "anode", V_F)
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

        # Pull final fields from each runner. `bias_sweep` exposes
        # `psi_phys`/`n_phys`/`p_phys` attributes; `transient` writes
        # snapshot lists into `.fields`. `final_field` hides the
        # difference.
        try:
            psi_err = relative_l2(final_field(tr_result, "psi"), final_field(bs_result, "psi"))
            n_err = relative_l2(final_field(tr_result, "n"), final_field(bs_result, "n"))
            p_err = relative_l2(final_field(tr_result, "p"), final_field(bs_result, "p"))
        except KeyError as exc:
            pytest.fail(
                f"{CASE}: required field missing from runner output: {exc}"
            )

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
    for V, e_psi, e_n, e_p, _e_J in rows:
        assert e_psi < 1e-2, f"psi disagreement {e_psi:.3e} at V_F={V}"
        assert e_n < 5e-2, f"n disagreement {e_n:.3e} at V_F={V}"
        assert e_p < 5e-2, f"p disagreement {e_p:.3e} at V_F={V}"
