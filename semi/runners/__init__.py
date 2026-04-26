"""
Runner implementations for the public `semi.run.run` dispatcher.

Each runner here owns one solver-type code path. Public modules:

    equilibrium    M1 equilibrium Poisson under Boltzmann statistics.
    bias_sweep     Coupled Slotboom drift-diffusion with adaptive
                   continuation across an applied bias range.
    transient      BDF1/BDF2 time integration in (psi, n, p) form.
    ac_sweep       Small-signal AC analysis around a DC operating point
                   (M14): solves (J + j*omega*M) delta_u = -dF/dV delta_V
                   for a frequency sweep and reports Y, Z, C, G.
    mos_cap_ac     MOS-capacitor differential C(V_gate) via analytic PDE
                   sensitivity (M14.1, replaces `mos_cv` finite-difference
                   dQ/dV).

Both steady-state runners return the `SimulationResult` dataclass defined
in `semi.run`; importing `SimulationResult` from `semi.run` is the
backward-compatible path for callers that already do so. The transient
runner returns `TransientResult` from `semi.results`; the AC sweep
runner returns `AcSweepResult` from `semi.results`.

The `semi.run.run(cfg)` dispatcher is the public entry point and is
unchanged in signature.
"""
from __future__ import annotations

from .ac_sweep import run_ac_sweep
from .bias_sweep import run_bias_sweep
from .equilibrium import run_equilibrium
from .mos_cap_ac import run_mos_cap_ac
from .mos_cv import run_mos_cv
from .transient import run_transient

__all__ = [
    "run_equilibrium", "run_bias_sweep", "run_mos_cv", "run_transient",
    "run_ac_sweep", "run_mos_cap_ac",
]
