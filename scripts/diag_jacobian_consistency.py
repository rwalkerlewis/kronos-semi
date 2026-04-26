#!/usr/bin/env python3
"""
Jacobian-consistency diagnostic via PETSc -snes_test_jacobian.

Why
---
SNES needs the analytic Jacobian to match a finite-difference reference
to ~1e-6 (relative, in the Frobenius norm) for Newton to converge at
the tightened atol used in M12 (`docs/adr/0008-snes-tolerances.md`).
Issue #17 calls for an audit of three benchmark configs:

    * `benchmarks/pn_1d_bias`        (M2 forward sweep, single-region DD)
    * `benchmarks/mosfet_2d`         (M12 multi-region with Gaussian n+ implants)
    * `benchmarks/mos_2d`            (M14.1 mos_cap_ac equilibrium-Poisson sensitivity)

Usage
-----
    python scripts/diag_jacobian_consistency.py <config.json>

Each Newton iteration's PETSc output gets

    Testing hand-coded Jacobian, ...
    ||J - Jfd||_F/||J||_F = ..., ||J - Jfd||_F = ...
    ||J - Jfd||_inf/||J||_inf = ..., ||J - Jfd||_inf = ...

interleaved with the ordinary `0 SNES Function norm` traces. Filter
post-hoc to focus on a specific operating point.

Implementation note
-------------------
Setting `-snes_test_jacobian` via the global PETSC_OPTIONS env var does
not work because every SNES the runners build is namespaced by a per-
solve `petsc_options_prefix` (so the un-prefixed global option is left
unused). Instead we patch `semi.solver.DEFAULT_PETSC_OPTIONS` so the
key gets merged into the per-prefix options dict each runner already
forwards into NonlinearProblem.

This script is a one-shot diagnostic, not a CI gate. It is not wired
into `.github/workflows/ci.yml` to keep the FD-Jacobian assembly cost
out of every push.
"""
from __future__ import annotations

import argparse
import sys


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument(
        "config", help="path to a benchmark config JSON (semi/schema.py validated)",
    )
    args = parser.parse_args(argv)

    import semi.solver

    semi.solver.DEFAULT_PETSC_OPTIONS["snes_test_jacobian"] = None
    semi.solver.DEFAULT_PETSC_OPTIONS["snes_test_jacobian_threshold"] = 0.0

    from semi import run as semi_run
    from semi import schema

    cfg = schema.load(args.config)
    print(f"[diag_jacobian] config           = {args.config}", flush=True)
    print(f"[diag_jacobian] solver.type      = {cfg['solver']['type']}", flush=True)
    print("[diag_jacobian] DEFAULT_PETSC_OPTIONS keys:", flush=True)
    for k, v in semi.solver.DEFAULT_PETSC_OPTIONS.items():
        if k.startswith("snes_test_jacobian"):
            print(f"    {k} = {v}", flush=True)
    print("=" * 78, flush=True)
    semi_run.run(cfg)
    print("=" * 78, flush=True)
    print(f"[diag_jacobian] done             = {args.config}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
