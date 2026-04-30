"""
Runner for axisymmetric MOSCAP LF/HF C-V sweeps.

Dispatched by solver.type == 'cv_lf' or 'cv_hf'.
"""
from __future__ import annotations

from typing import Any

import numpy as np


def run_moscap_axi_cv(cfg: dict[str, Any]):
    """
    Run LF and/or HF C-V for the axisymmetric MOSCAP.

    Returns a SimulationResult with iv rows carrying
    {V, Q_gate, J: 0, C} for LF or {V, Q_gate, J: 0, C_lf, C_hf}
    for the combined mode.
    """
    from ..physics.cv import compute_cv_curve
    from ..run import SimulationResult

    stype = cfg.get("solver", {}).get("type", "cv_lf")

    if stype == "cv_lf":
        result = compute_cv_curve(cfg, mode="LF")
        iv_rows = [
            {"V": float(vg), "Q_gate": float(q), "J": 0.0, "C": float(c)}
            for vg, q, c in zip(result["Vg"], result["Q_sub"], result["C"])
        ]
        return SimulationResult(cfg=cfg, iv=iv_rows, bias_contact="gate")
    elif stype == "cv_hf":
        result_lf = compute_cv_curve(cfg, mode="LF")
        result_hf = compute_cv_curve(cfg, mode="HF")
        iv_rows = [
            {
                "V": float(vg),
                "Q_gate": float(q),
                "J": 0.0,
                "C_lf": float(c_lf),
                "C_hf": float(c_hf),
            }
            for vg, q, c_lf, c_hf in zip(
                result_lf["Vg"], result_lf["Q_sub"],
                result_lf["C"], result_hf["C"]
            )
        ]
        return SimulationResult(cfg=cfg, iv=iv_rows, bias_contact="gate")
    else:
        raise ValueError(f"run_moscap_axi_cv: unexpected solver type {stype!r}")
