"""
Top-level entry point: run a simulation from a validated config dict.

Supports:
    solver.type == "equilibrium"      equilibrium Poisson (Day 1).
    solver.type == "drift_diffusion"  coupled solve at baked biases.
    solver.type == "bias_sweep"       bias ramp: walk a contact's
                                      voltage_sweep and solve coupled.

Implementation note: this module is a thin dispatcher. The actual
runners live in `semi.runners.equilibrium` and
`semi.runners.bias_sweep`; the post-step current evaluation and IV
recording helpers live in `semi.postprocess`. `SimulationResult` stays
here for backward-compat (it is the public output type and is imported
by external callers, including `scripts/run_benchmark.py`).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np


@dataclass
class SimulationResult:
    """Container for simulation outputs."""
    cfg: dict[str, Any]
    mesh: Any = None
    V: Any = None
    psi: Any = None
    phi_n: Any = None
    phi_p: Any = None
    psi_phys: np.ndarray | None = None
    phi_n_phys: np.ndarray | None = None
    phi_p_phys: np.ndarray | None = None
    n_phys: np.ndarray | None = None
    p_phys: np.ndarray | None = None
    x_dof: np.ndarray | None = None
    N_hat: Any = None
    scaling: Any = None
    solver_info: dict[str, Any] = field(default_factory=dict)
    iv: list[dict[str, float]] = field(default_factory=list)
    bias_contact: str | None = None


def run(cfg: dict[str, Any]) -> SimulationResult:
    """Dispatch on solver.type."""
    from .runners import run_bias_sweep, run_equilibrium, run_mos_cv

    stype = cfg.get("solver", {}).get("type", "equilibrium")
    if stype == "equilibrium":
        return run_equilibrium(cfg)
    if stype in ("drift_diffusion", "bias_sweep"):
        return run_bias_sweep(cfg)
    if stype == "mos_cv":
        return run_mos_cv(cfg)
    raise ValueError(f"Unknown solver.type {stype!r}")


# Backward-compatible re-exports. External callers (V&V suite, tests)
# import these names from `semi.run`; keep them resolvable here.
def __getattr__(name: str):
    if name == "run_equilibrium":
        from .runners.equilibrium import run_equilibrium
        return run_equilibrium
    if name == "run_bias_sweep":
        from .runners.bias_sweep import run_bias_sweep
        return run_bias_sweep
    if name == "_fmt_tag":
        from .postprocess import fmt_tag
        return fmt_tag
    if name == "_resolve_sweep":
        from .runners.bias_sweep import _resolve_sweep
        return _resolve_sweep
    if name == "run_mos_cv":
        from .runners.mos_cv import run_mos_cv
        return run_mos_cv
    raise AttributeError(f"module 'semi.run' has no attribute {name!r}")
