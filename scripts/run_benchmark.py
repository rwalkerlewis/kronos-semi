#!/usr/bin/env python3
"""
Benchmark runner CLI.

Usage:
    python scripts/run_benchmark.py <benchmark_name>

Where <benchmark_name> matches a directory under benchmarks/ containing
exactly one *.json input file (or a file named <benchmark_name>.json).

For each benchmark, this script:
  - loads and validates the JSON config
  - calls semi.run.run(cfg) to produce a SimulationResult
  - prints SNES diagnostics (iterations, converged reason, solve time, DOFs)
  - invokes a registered verifier that asserts physical correctness
  - writes plots to results/<benchmark_name>/ using matplotlib's Agg backend

Exits 0 if everything passes, nonzero otherwise.
"""
from __future__ import annotations

import argparse
import glob
import os
import sys
import time
import traceback
from pathlib import Path
from typing import Any, Callable

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[1]
BENCHMARKS_DIR = REPO_ROOT / "benchmarks"
RESULTS_DIR = REPO_ROOT / "results"


# --------------------------------------------------------------------------- #
# Verifier registry                                                           #
# --------------------------------------------------------------------------- #

Verifier = Callable[[Any], list[tuple[str, bool, str]]]
_VERIFIERS: dict[str, Verifier] = {}


def register(name: str) -> Callable[[Verifier], Verifier]:
    def _wrap(fn: Verifier) -> Verifier:
        _VERIFIERS[name] = fn
        return fn
    return _wrap


# --------------------------------------------------------------------------- #
# pn_1d verifier                                                              #
# --------------------------------------------------------------------------- #

@register("pn_1d")
def verify_pn_1d(result) -> list[tuple[str, bool, str]]:
    """
    Verify equilibrium 1D pn junction against depletion-approximation theory.

    Checks:
      - V_bi matches V_t * ln(N_A * N_D / n_i^2) within 5%
      - peak |E| matches q * N_A * x_p / eps within 10% (symmetric junction)
      - bulk densities within factor of 2 of majority doping
      - mass-action n*p == n_i^2 in bulk to 1%
    """
    from semi.constants import EPS0, Q, cm3_to_m3
    from semi.materials import get_material

    cfg = result.cfg
    sc = result.scaling
    V_t = sc.V0

    # Doping from JSON (step profile assumed for pn_1d)
    prof = cfg["doping"][0]["profile"]
    if prof["type"] != "step":
        return [("pn_1d: step-doping", False, f"profile type {prof['type']}")]

    N_A = cm3_to_m3(prof["N_A_left"])
    N_D = cm3_to_m3(prof["N_D_right"])
    junction_x = float(prof["location"])

    mat = get_material(cfg["regions"]["silicon"]["material"])
    n_i = mat.n_i
    eps = mat.epsilon

    # Sort dof values by x for 1D analysis
    x = result.x_dof[:, 0]
    order = np.argsort(x)
    x = x[order]
    psi = result.psi_phys[order]
    n = result.n_phys[order]
    p = result.p_phys[order]

    checks: list[tuple[str, bool, str]] = []

    # 1. V_bi
    V_bi_theory = V_t * np.log(N_A * N_D / n_i**2)
    V_bi_sim = float(psi[-1] - psi[0])
    rel_err_Vbi = abs(V_bi_sim - V_bi_theory) / V_bi_theory
    checks.append((
        "V_bi within 5%",
        rel_err_Vbi < 0.05,
        f"sim={V_bi_sim:.4f} V, theory={V_bi_theory:.4f} V, rel_err={rel_err_Vbi:.2%}",
    ))

    # 2. Peak |E| = -dpsi/dx
    E = -np.gradient(psi, x)
    E_peak_sim = float(np.max(np.abs(E)))
    # symmetric junction -> x_p = W/2 where
    #   W = sqrt(2 eps V_bi (N_A + N_D) / (q N_A N_D))
    W_theory = np.sqrt(2.0 * eps * V_bi_theory * (N_A + N_D) / (Q * N_A * N_D))
    x_p_theory = W_theory * N_D / (N_A + N_D)  # asymmetric-safe formula
    E_peak_theory = Q * N_A * x_p_theory / eps
    rel_err_E = abs(E_peak_sim - E_peak_theory) / E_peak_theory
    checks.append((
        "peak |E| within 10%",
        rel_err_E < 0.10,
        f"sim={E_peak_sim/1e5:.2f} kV/cm, theory={E_peak_theory/1e5:.2f} kV/cm, "
        f"rel_err={rel_err_E:.2%}",
    ))

    # 3. Bulk densities within factor of 2 of majority doping
    # Left bulk: x near 0 (p-side), majority = holes
    # Right bulk: x near L (n-side), majority = electrons
    L = x[-1]
    left_mask = x < 0.1 * L
    right_mask = x > 0.9 * L
    p_left_avg = float(np.mean(p[left_mask]))
    n_right_avg = float(np.mean(n[right_mask]))
    ratio_p = p_left_avg / N_A
    ratio_n = n_right_avg / N_D
    ok_p = 0.5 < ratio_p < 2.0
    ok_n = 0.5 < ratio_n < 2.0
    checks.append((
        "p-side bulk hole density within factor 2 of N_A",
        ok_p,
        f"<p>_left={p_left_avg:.2e} m^-3, N_A={N_A:.2e} m^-3, ratio={ratio_p:.2f}",
    ))
    checks.append((
        "n-side bulk electron density within factor 2 of N_D",
        ok_n,
        f"<n>_right={n_right_avg:.2e} m^-3, N_D={N_D:.2e} m^-3, ratio={ratio_n:.2f}",
    ))

    # 4. Mass action in bulk
    np_product = n * p
    np_left = float(np.mean(np_product[left_mask]))
    np_right = float(np.mean(np_product[right_mask]))
    ni2 = n_i**2
    err_left = abs(np_left - ni2) / ni2
    err_right = abs(np_right - ni2) / ni2
    checks.append((
        "mass action n*p = n_i^2 in p-side bulk to 1%",
        err_left < 0.01,
        f"<np>_left={np_left:.3e}, n_i^2={ni2:.3e}, rel_err={err_left:.2%}",
    ))
    checks.append((
        "mass action n*p = n_i^2 in n-side bulk to 1%",
        err_right < 0.01,
        f"<np>_right={np_right:.3e}, n_i^2={ni2:.3e}, rel_err={err_right:.2%}",
    ))

    # Stash analytical references for the plotter
    result._analytical = {
        "V_bi": V_bi_theory,
        "W": W_theory,
        "x_p": x_p_theory,
        "E_peak": E_peak_theory,
        "junction_x": junction_x,
        "N_A": N_A, "N_D": N_D, "n_i": n_i, "eps": eps,
    }
    return checks


# --------------------------------------------------------------------------- #
# Plotting                                                                    #
# --------------------------------------------------------------------------- #

def plot_pn_1d(result, out_dir: Path) -> list[Path]:
    """Write potential/field/carrier density plots for a 1D pn junction."""
    x = result.x_dof[:, 0]
    order = np.argsort(x)
    x = x[order]
    psi = result.psi_phys[order]
    n = result.n_phys[order]
    p = result.p_phys[order]
    E = -np.gradient(psi, x)

    x_um = x * 1e6
    paths = []

    # Potential
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(x_um, psi, label=r"$\psi(x)$")
    ax.set_xlabel("x (um)")
    ax.set_ylabel("potential (V)")
    ax.set_title("Equilibrium potential")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    p1 = out_dir / "potential.png"
    fig.savefig(p1, dpi=130)
    plt.close(fig)
    paths.append(p1)

    # Electric field
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(x_um, E / 1e5, label=r"$E = -d\psi/dx$")
    ax.set_xlabel("x (um)")
    ax.set_ylabel("E (kV/cm)")
    ax.set_title("Electric field")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    p2 = out_dir / "field.png"
    fig.savefig(p2, dpi=130)
    plt.close(fig)
    paths.append(p2)

    # Carrier densities (log scale)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.semilogy(x_um, n * 1e-6, label="n")  # m^-3 -> cm^-3
    ax.semilogy(x_um, p * 1e-6, label="p")
    ax.set_xlabel("x (um)")
    ax.set_ylabel("density (cm^-3)")
    ax.set_title("Carrier densities")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    p3 = out_dir / "carriers.png"
    fig.savefig(p3, dpi=130)
    plt.close(fig)
    paths.append(p3)

    return paths


_PLOTTERS: dict[str, Callable[[Any, Path], list[Path]]] = {
    "pn_1d": plot_pn_1d,
}


# --------------------------------------------------------------------------- #
# Driver                                                                      #
# --------------------------------------------------------------------------- #

def find_json(bench_dir: Path, name: str) -> Path:
    candidate = bench_dir / f"{name}.json"
    if candidate.exists():
        return candidate
    jsons = sorted(glob.glob(str(bench_dir / "*.json")))
    if len(jsons) == 1:
        return Path(jsons[0])
    if not jsons:
        raise FileNotFoundError(f"No JSON input found in {bench_dir}")
    raise RuntimeError(
        f"Multiple JSON files in {bench_dir}; name one {name}.json or pass --input"
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run a kronos-semi benchmark.")
    parser.add_argument("name", help="benchmark directory name under benchmarks/")
    parser.add_argument("--input", help="explicit JSON path (overrides lookup)")
    parser.add_argument(
        "--no-verify", action="store_true",
        help="skip verifier assertions (still runs simulation)",
    )
    args = parser.parse_args(argv)

    name = args.name
    bench_dir = BENCHMARKS_DIR / name
    if not bench_dir.is_dir():
        print(f"ERROR: benchmark directory not found: {bench_dir}", file=sys.stderr)
        return 2

    try:
        json_path = Path(args.input) if args.input else find_json(bench_dir, name)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    print(f"[run_benchmark] benchmark = {name}")
    print(f"[run_benchmark] input     = {json_path}")

    # Imports deferred so --help works without dolfinx
    from semi import run as semi_run
    from semi import schema

    cfg = schema.load(str(json_path))

    t0 = time.perf_counter()
    try:
        result = semi_run.run(cfg)
    except Exception:
        print("[run_benchmark] simulation raised:", file=sys.stderr)
        traceback.print_exc()
        return 3
    dt = time.perf_counter() - t0

    info = result.solver_info or {}
    n_dofs = int(result.V.dofmap.index_map.size_global) if result.V is not None else -1

    print("[run_benchmark] SNES diagnostics:")
    print(f"    converged       = {info.get('converged')}")
    print(f"    iterations      = {info.get('iterations')}")
    print(f"    reason          = {info.get('reason')}")
    print(f"    solve_time      = {dt:.3f} s")
    print(f"    dof_count       = {n_dofs}")
    print(f"    scaling         = {result.scaling}")

    if not info.get("converged", False):
        print("[run_benchmark] FAIL: SNES did not converge", file=sys.stderr)
        return 4

    # Output directory
    out_dir = RESULTS_DIR / name
    out_dir.mkdir(parents=True, exist_ok=True)

    # Plots
    plotter = _PLOTTERS.get(name)
    if plotter is not None:
        paths = plotter(result, out_dir)
        print(f"[run_benchmark] wrote {len(paths)} plot(s):")
        for p in paths:
            print(f"    {p.relative_to(REPO_ROOT)}")
    else:
        print(f"[run_benchmark] no plotter registered for {name!r}; skipping plots")

    # Verification
    if args.no_verify:
        print("[run_benchmark] verification skipped (--no-verify)")
        return 0

    verifier = _VERIFIERS.get(name)
    if verifier is None:
        print(f"[run_benchmark] no verifier registered for {name!r}; skipping checks")
        return 0

    print("[run_benchmark] verification:")
    checks = verifier(result)
    all_ok = True
    for label, ok, detail in checks:
        marker = "PASS" if ok else "FAIL"
        print(f"    [{marker}] {label}: {detail}")
        all_ok = all_ok and ok

    if not all_ok:
        print("[run_benchmark] FAIL: one or more verification checks failed",
              file=sys.stderr)
        return 5

    print("[run_benchmark] OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
