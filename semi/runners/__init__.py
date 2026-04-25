"""
Runner implementations for the public `semi.run.run` dispatcher.

Each runner here owns one solver-type code path. Public modules:

    equilibrium    M1 equilibrium Poisson under Boltzmann statistics.
    bias_sweep     Coupled Slotboom drift-diffusion with adaptive
                   continuation across an applied bias range.
    transient      BDF1/BDF2 time integration in (psi, n, p) form.

Both steady-state runners return the `SimulationResult` dataclass defined
in `semi.run`; importing `SimulationResult` from `semi.run` is the
backward-compatible path for callers that already do so. The transient
runner returns `TransientResult` from `semi.results`.

The `semi.run.run(cfg)` dispatcher is the public entry point and is
unchanged in signature.
"""
from __future__ import annotations

from .bias_sweep import run_bias_sweep
from .equilibrium import run_equilibrium
from .mos_cv import run_mos_cv
from .transient import run_transient

__all__ = ["run_equilibrium", "run_bias_sweep", "run_mos_cv", "run_transient"]
