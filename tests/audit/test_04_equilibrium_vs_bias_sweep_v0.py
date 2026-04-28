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
    load_benchmark,
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

    cfg_bs = copy.deepcopy(cfg_base)
    cfg_bs["contacts"][0]["voltage"] = 0.0
    cfg_bs["solver"] = {
        "type": "bias_sweep",
        "bias_ramp": {"start": 0.0, "stop": 0.0, "step": 0.05},
    }
    bs_res = run_bias_sweep(cfg_bs)

    eq_fields = getattr(eq_res, "fields", {}) or {}
    bs_fields = getattr(bs_res, "fields", {}) or {}

    def first(name, fields):
        for k in (name, "potential" if name == "psi" else name):
            if k in fields and len(fields[k]) > 0:
                return np.asarray(fields[k][0])
        raise KeyError(name)

    try:
        e_psi = relative_l2(first("psi", eq_fields), first("psi", bs_fields))
        e_n = relative_l2(first("n", eq_fields), first("n", bs_fields))
        e_p = relative_l2(first("p", eq_fields), first("p", bs_fields))
    except KeyError as exc:
        pytest.fail(
            f"{CASE}: required field missing: {exc}; "
            f"eq keys = {list(eq_fields)}, bs keys = {list(bs_fields)}"
        )

    rows = [["psi", e_psi], ["n", e_n], ["p", e_p]]

    write_csv(CASE, ["field", "rel_L2"], rows)
    write_markdown(
        CASE,
        "Case 04 - equilibrium vs bias_sweep at V=0",
        "Compares `run_equilibrium` to `run_bias_sweep` halted at "
        "V=0 on the pn_1d_bias device.\n\n"
        f"- psi rel_L2: {e_psi:.3e}\n"
        f"- n   rel_L2: {e_n:.3e}\n"
        f"- p   rel_L2: {e_p:.3e}\n\n"
        f"CSV: `/tmp/audit/{CASE}.csv`",
    )

    # The two paths should agree to round-off. 1e-8 is generous; a
    # disagreement bigger than this is a real finding.
    assert e_psi < 1e-8, f"psi disagrees: {e_psi:.3e}"
    assert e_n < 1e-6, f"n disagrees: {e_n:.3e}"
    assert e_p < 1e-6, f"p disagrees: {e_p:.3e}"
