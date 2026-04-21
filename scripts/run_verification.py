#!/usr/bin/env python3
"""
V&V (Verification & Validation) runner CLI.

Sister script to `scripts/run_benchmark.py`. Each subcommand runs one
verification activity, writes CSV+PNG artifacts under `results/`, and
exits 0 on success.

Subcommands (Phase 1 implements `mms_poisson` and `all`):

    python scripts/run_verification.py mms_poisson
    python scripts/run_verification.py mesh_convergence    # Phase 2 (placeholder)
    python scripts/run_verification.py conservation        # Phase 3 (placeholder)
    python scripts/run_verification.py mms_dd              # Phase 4 (placeholder)
    python scripts/run_verification.py all
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = REPO_ROOT / "results"


def _gate_finest_pair_rate(rows, key: str, threshold: float, label: str) -> tuple[bool, str]:
    """Pass-fail check on the finest-pair rate column."""
    rate = rows[-1].get(key)
    if rate is None or rate != rate:  # NaN check without numpy
        return False, f"{label}: finest-pair {key} is NaN"
    ok = rate >= threshold
    marker = "PASS" if ok else "FAIL"
    return ok, f"[{marker}] {label}: finest-pair {key} = {rate:.3f} (threshold {threshold:.2f})"


def cmd_mms_poisson(args) -> int:
    """Run all MMS-Poisson studies and gate the finest-pair convergence rates."""
    from semi.verification.mms_poisson import report_table, run_cli_study

    out_dir = Path(args.out) if args.out else RESULTS_DIR / "mms_poisson"
    print(f"[run_verification] mms_poisson -> {out_dir.relative_to(REPO_ROOT)}")
    studies = run_cli_study(out_dir)

    all_ok = True
    for label in ("1d_linear", "1d_nonlinear", "2d_triangles"):
        rows = studies[label]
        print()
        print(report_table(rows, header=f"=== mms_poisson {label} ==="))
        ok_l2, msg_l2 = _gate_finest_pair_rate(rows, "rate_L2", 1.85, f"{label} L2")
        ok_h1, msg_h1 = _gate_finest_pair_rate(rows, "rate_H1", 0.85, f"{label} H1")
        print(msg_l2)
        print(msg_h1)
        all_ok = all_ok and ok_l2 and ok_h1

    smoke = studies["2d_quad_smoke"][0]
    print()
    print("=== mms_poisson 2d_quad_smoke (single mesh, N=64) ===")
    print(
        f"  triangle e_L2 = {smoke['e_L2_triangle']:.3e}, "
        f"quad e_L2 = {smoke['e_L2_quad']:.3e}, "
        f"ratio = {smoke['ratio_L2_quad_over_triangle']:.3f}"
    )
    smoke_ok = 0.1 < smoke["ratio_L2_quad_over_triangle"] < 10.0
    print(("[PASS]" if smoke_ok else "[FAIL]") + " quad/triangle ratio in [0.1, 10]")
    all_ok = all_ok and smoke_ok

    return 0 if all_ok else 5


def cmd_not_implemented(name: str):
    def _run(_args) -> int:
        print(f"[run_verification] {name}: not implemented yet (Phase 2/3/4)")
        return 0
    return _run


def cmd_all(args) -> int:
    """Run every implemented V&V activity in sequence."""
    rc = cmd_mms_poisson(args)
    if rc != 0:
        return rc
    # Future: chain mesh_convergence, conservation, mms_dd here.
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run a kronos-semi V&V study and write CSV+PNG artifacts."
    )
    parser.add_argument("--out", help="output root directory (default: results/<study>)")
    sub = parser.add_subparsers(dest="study", required=True)

    sub.add_parser("mms_poisson", help="MMS for equilibrium Poisson (Phase 1)")
    sub.add_parser("mesh_convergence", help="pn_1d mesh convergence (Phase 2)")
    sub.add_parser("conservation", help="current and charge conservation (Phase 3)")
    sub.add_parser("mms_dd", help="MMS for coupled drift-diffusion (Phase 4)")
    sub.add_parser("all", help="run every implemented study")

    args = parser.parse_args(argv)
    if args.study == "mms_poisson":
        return cmd_mms_poisson(args)
    if args.study == "mesh_convergence":
        return cmd_not_implemented("mesh_convergence")(args)
    if args.study == "conservation":
        return cmd_not_implemented("conservation")(args)
    if args.study == "mms_dd":
        return cmd_not_implemented("mms_dd")(args)
    if args.study == "all":
        return cmd_all(args)
    parser.error(f"Unknown study {args.study!r}")
    return 2


if __name__ == "__main__":
    sys.exit(main())
