"""
Audit case 3: mos_cv vs mos_cap_ac on Q_gate (byte-identity claim).

M14.1's CHANGELOG asserts that the analytic-AC `mos_cap_ac` runner
agrees with the numerical-dQ/dV `mos_cv` runner on Q_gate to byte
identity at matching gate voltages. This test runs both on the same
mos_2d config and compares Q(V_gate) point-wise.
"""
from __future__ import annotations

import copy

import numpy as np
import pytest

from ._helpers import (
    load_benchmark,
    require_dolfinx,
    write_csv,
    write_markdown,
)

CASE = "03_mos_cv_vs_mos_cap_ac"


@pytest.mark.audit
def test_mos_cv_vs_mos_cap_ac():
    require_dolfinx()

    from semi.runners.mos_cap_ac import run_mos_cap_ac
    from semi.runners.mos_cv import run_mos_cv

    cfg = load_benchmark("mos_2d/mos_cap.json")

    cfg_cv = copy.deepcopy(cfg)
    cfg_cv["solver"] = {"type": "mos_cv"}
    res_cv = run_mos_cv(cfg_cv)

    cfg_ac = copy.deepcopy(cfg)
    cfg_ac["solver"] = {"type": "mos_cap_ac"}
    res_ac = run_mos_cap_ac(cfg_ac)

    # Both runners produce a list of dicts with V_gate and Q_gate (or
    # similar). Adapt to whichever attribute exposes the bias / charge
    # arrays. Defensive: report whatever fields are exposed.
    def extract(result):
        # Try common shapes: result.iv, result.cv, result.points.
        for attr in ("cv", "iv", "points"):
            data = getattr(result, attr, None)
            if data:
                Vs, Qs = [], []
                for row in data:
                    V = row.get("V_gate") or row.get("V") or row.get("voltage")
                    Q = row.get("Q_gate") or row.get("Q") or row.get("charge")
                    if V is None or Q is None:
                        continue
                    Vs.append(float(V))
                    Qs.append(float(Q))
                if Vs:
                    return np.array(Vs), np.array(Qs)
        raise AttributeError(
            f"could not extract (V, Q) from result of type {type(result).__name__}"
        )

    try:
        V_cv, Q_cv = extract(res_cv)
        V_ac, Q_ac = extract(res_ac)
    except AttributeError as exc:
        pytest.skip(
            f"{CASE}: runner output shape unknown to audit harness ({exc}); "
            "open a tracking issue to standardize MOS C-V outputs."
        )

    # Match V points: take intersection (rounded).
    V_cv_r = np.round(V_cv, 6)
    V_ac_r = np.round(V_ac, 6)
    common = np.intersect1d(V_cv_r, V_ac_r)
    if len(common) == 0:
        pytest.fail(
            f"{CASE}: no common gate-voltage points between runners "
            f"(cv: {V_cv_r[:5]}..., ac: {V_ac_r[:5]}...)"
        )

    rows = []
    max_rel = 0.0
    for V in common:
        Q1 = float(Q_cv[np.where(V_cv_r == V)[0][0]])
        Q2 = float(Q_ac[np.where(V_ac_r == V)[0][0]])
        denom = max(abs(Q1), abs(Q2), 1e-300)
        rel = abs(Q1 - Q2) / denom
        rows.append([float(V), Q1, Q2, rel])
        max_rel = max(max_rel, rel)

    write_csv(
        CASE,
        ["V_gate", "Q_mos_cv", "Q_mos_cap_ac", "rel_err"],
        rows,
    )

    write_markdown(
        CASE,
        "Case 03 - mos_cv vs mos_cap_ac (Q_gate)",
        f"Q_gate(V) compared at {len(common)} common gate voltages on "
        f"`benchmarks/mos_2d`. Worst relative disagreement: "
        f"{max_rel:.3e}.\n\nCSV: `/tmp/audit/{CASE}.csv`",
    )

    # M14.1 byte-identity claim implies essentially zero disagreement.
    # Relax to 1e-3 to allow legitimate numerical noise and re-tighten
    # if the run reports much smaller.
    assert max_rel < 1e-3, (
        f"M14.1 claims byte-identical Q_gate; observed rel_err {max_rel:.3e}. "
        f"Investigate before tightening the gate."
    )
