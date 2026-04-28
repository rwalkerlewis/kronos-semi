"""
Audit case 4: equilibrium runner vs bias_sweep at V=0.

The equilibrium solver should produce the same (psi, n, p) field as
the bias_sweep runner halted at V=0 (the V=0 step in the ramp). Any
disagreement points at a sign or BC discrepancy in one of the two
paths.
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

CASE = "04_equilibrium_vs_bias_sweep_V0"


@pytest.mark.audit
def test_equilibrium_vs_bias_sweep_V0():
    require_dolfinx()

    from semi.runners.bias_sweep import run_bias_sweep
    from semi.runners.equilibrium import run_equilibrium

    cfg_base = load_benchmark("pn_1d_bias/pn_junction_bias.json")

    cfg_eq = copy.deepcopy(cfg_base)
    cfg_eq["solver"] = {"type": "equilibrium"}
    eq_res = run_equilibrium(cfg_eq)

    cfg_bs = make_bias_sweep_cfg(cfg_base, "anode", 0.0, relaxed_snes=False)
    bs_res = run_bias_sweep(cfg_bs)

    # `run_equilibrium` references psi to the intrinsic Fermi level
    # (psi = +/- V_t arcsinh(N/2 n_i)); `run_bias_sweep` references psi
    # to the swept-contact built-in potential plus applied bias. The
    # two conventions differ by a constant gauge shift; n, p, and any
    # observable current are gauge-invariant and must agree at V=0.
    # Sort by x-coordinate so we compare matching nodal locations even
    # if dolfinx returns DOFs in different orders.
    psi_eq = final_field(eq_res, "psi")
    psi_bs = final_field(bs_res, "psi")
    n_eq = final_field(eq_res, "n")
    n_bs = final_field(bs_res, "n")
    p_eq = final_field(eq_res, "p")
    p_bs = final_field(bs_res, "p")

    ie = np.argsort(eq_res.x_dof[:, 0])
    ib = np.argsort(bs_res.x_dof[:, 0])
    psi_eq_s = psi_eq[ie]
    psi_bs_s = psi_bs[ib]
    n_eq_s = n_eq[ie]
    n_bs_s = n_bs[ib]
    p_eq_s = p_eq[ie]
    p_bs_s = p_bs[ib]

    # Gauge-aligned psi disagreement: subtract the mean offset before
    # comparing, since psi is only defined up to an additive constant.
    psi_offset = float(np.mean(psi_bs_s - psi_eq_s))
    psi_aligned_err = relative_l2(psi_eq_s + psi_offset, psi_bs_s)
    psi_raw_err = relative_l2(psi_eq_s, psi_bs_s)
    e_n = relative_l2(n_eq_s, n_bs_s)
    e_p = relative_l2(p_eq_s, p_bs_s)

    rows = [
        ["psi_raw", psi_raw_err],
        ["psi_gauge_aligned", psi_aligned_err],
        ["psi_gauge_offset_V", psi_offset],
        ["n", e_n],
        ["p", e_p],
    ]

    write_csv(CASE, ["field", "value"], rows)
    write_markdown(
        CASE,
        "Case 04 - equilibrium vs bias_sweep at V=0",
        "Compares `run_equilibrium` to `run_bias_sweep` halted at "
        "V=0 on the pn_1d_bias device. The two runners use different "
        "psi reference conventions (intrinsic Fermi vs swept-contact "
        "BC), so the raw psi disagreement is a constant gauge shift; "
        "n and p are gauge-invariant.\n\n"
        f"- psi raw rel_L2: {psi_raw_err:.3e}\n"
        f"- psi gauge offset (mean shift): {psi_offset:+.3e} V\n"
        f"- psi gauge-aligned rel_L2: {psi_aligned_err:.3e}\n"
        f"- n   rel_L2: {e_n:.3e}\n"
        f"- p   rel_L2: {e_p:.3e}\n\n"
        f"CSV: `/tmp/audit/{CASE}.csv`",
    )

    # Soft gates: gauge-invariant fields should agree to discretization
    # noise. The raw psi gate is intentionally generous because of the
    # gauge-shift convention difference; the gauge-aligned check is the
    # meaningful one.
    assert psi_aligned_err < 1e-3, (
        f"gauge-aligned psi disagrees: {psi_aligned_err:.3e}"
    )
    assert e_n < 5e-3, f"n disagrees: {e_n:.3e}"
    assert e_p < 5e-3, f"p disagrees: {e_p:.3e}"
