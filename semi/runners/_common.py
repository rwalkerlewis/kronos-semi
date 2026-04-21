"""
Helpers shared by `runners/equilibrium.py` and `runners/bias_sweep.py`.

Pure-Python; no dolfinx import.
"""
from __future__ import annotations

from typing import Any


def reference_material(cfg: dict[str, Any]):
    """
    Return the first semiconductor material referenced under
    `cfg["regions"]`. Used by both runners to source `n_i`,
    `epsilon_r`, etc.
    """
    from ..materials import get_material

    for _region_name, region in cfg["regions"].items():
        mat = get_material(region["material"])
        if mat.is_semiconductor():
            return mat
    raise ValueError("No semiconductor region found; nothing to solve.")
