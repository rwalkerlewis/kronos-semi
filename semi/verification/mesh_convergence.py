"""
Mesh convergence study for the 1D pn junction benchmark.

Unlike `mms_poisson`, this study has no exact solution to compare against.
Instead it treats the depletion approximation as an engineering reference
and measures how three canonical observables behave as h -> 0:

    V_bi     = psi(x=L) - psi(x=0)                    built-in potential
    E_peak   = max |-dpsi/dx|                          peak electric field
    W        = extent of |E| > 0.01 * E_peak           depletion width

HONEST FLAG
-----------
The depletion approximation is *approximate*. The FEM solver converges
to the true Poisson-Boltzmann solution (with exp(psi) / exp(-psi)
carrier terms), not to the depletion-approximation step-charge model.
On fine meshes the relative errors against V_bi_theory, E_peak_theory,
and W_theory therefore plateau at a physics-limited floor (typically
a few percent of peak field). That plateau is the difference between
the two *models*, not a discretization bug, and tightening the mesh
cannot remove it.

The pytest subset [100, 200, 400] is chosen so discretization error
still dominates the physics-approximation floor and the error ratio
per mesh doubling is bounded below by a meaningful constant. The full
CLI sweep [50..1600] exhibits the plateau and the PNG shows it.

Boundary conditions, doping, material, and all other physics are taken
from `benchmarks/pn_1d/pn_junction.json`; only `mesh.resolution` is
overridden per level.
"""
from __future__ import annotations

import copy
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_PN_1D_CFG = REPO_ROOT / "benchmarks" / "pn_1d" / "pn_junction.json"


@dataclass(frozen=True)
class MeshConvRow:
    """Per-level observables for the pn_1d mesh sweep."""

    N: int
    h: float
    n_dofs: int
    V_bi_sim: float
    V_bi_theory: float
    err_Vbi_rel: float
    E_peak_sim: float
    E_peak_theory: float
    err_Epeak_rel: float
    W_sim: float
    W_theory: float
    err_W_rel: float
    newton_iters: int
    solve_time_s: float


def _load_base_cfg(cfg_path: str | Path | None) -> dict:
    from semi import schema

    path = Path(cfg_path) if cfg_path else DEFAULT_PN_1D_CFG
    return schema.load(path)


def _depletion_references(cfg: dict) -> dict:
    """Compute V_bi, W, E_peak from the depletion approximation."""
    from semi.constants import Q, cm3_to_m3, thermal_voltage
    from semi.materials import get_material

    prof = cfg["doping"][0]["profile"]
    if prof["type"] != "step":
        raise ValueError(
            f"mesh_convergence expects a step doping profile, got {prof['type']!r}"
        )
    N_A = cm3_to_m3(prof["N_A_left"])
    N_D = cm3_to_m3(prof["N_D_right"])

    region_name = next(iter(cfg["regions"]))
    mat = get_material(cfg["regions"][region_name]["material"])
    T = float(cfg["physics"]["temperature"])
    V_t = thermal_voltage(T)
    n_i = mat.n_i
    eps = mat.epsilon

    V_bi = V_t * np.log(N_A * N_D / n_i**2)
    W = np.sqrt(2.0 * eps * V_bi * (N_A + N_D) / (Q * N_A * N_D))
    x_p = W * N_D / (N_A + N_D)
    E_peak = Q * N_A * x_p / eps
    return {
        "V_bi": float(V_bi),
        "W": float(W),
        "E_peak": float(E_peak),
        "N_A": float(N_A),
        "N_D": float(N_D),
        "n_i": float(n_i),
        "eps": float(eps),
    }


def _extract_observables(result, refs: dict) -> tuple[float, float, float]:
    """
    Return (V_bi_sim, E_peak_sim, W_sim) from the solver output.

    W_sim is extracted via charge-integration rather than threshold on |E|:
    Q = 0.5 * integral |rho/q| dx over the whole mesh equals the total
    depletion charge per unit area, and by charge neutrality across the
    junction W = Q * (N_A + N_D) / (N_A * N_D). This definition is
    threshold-free, matches the depletion approximation exactly on its
    own model, and converges cleanly with mesh refinement because the
    integrand is O(N_A) inside the depletion region and near-zero in
    quasi-neutral bulk.
    """
    sc = result.scaling
    x_raw = result.x_dof[:, 0]
    order = np.argsort(x_raw)
    x = x_raw[order]
    psi = result.psi_phys[order]
    n = result.n_phys[order]
    p = result.p_phys[order]
    N_net = result.N_hat.x.array[order] * sc.C0   # m^-3

    V_bi_sim = float(psi[-1] - psi[0])

    E = -np.gradient(psi, x)
    abs_E = np.abs(E)
    E_peak_sim = float(abs_E.max())

    rho_over_q = p - n + N_net  # m^-3
    Q_abs_total = float(np.trapezoid(np.abs(rho_over_q), x))
    Q = 0.5 * Q_abs_total
    N_A = refs["N_A"]
    N_D = refs["N_D"]
    if N_A > 0.0 and N_D > 0.0:
        W_sim = Q * (N_A + N_D) / (N_A * N_D)
    else:
        W_sim = float("nan")

    return V_bi_sim, E_peak_sim, W_sim


def run_one_level(N: int, *, cfg_path: str | Path | None = None) -> MeshConvRow:
    """
    Run the pn_1d benchmark at resolution N and return one row of
    convergence data. Reuses `semi.run.run_equilibrium` so the solver
    path is the same as the production benchmark.
    """
    from semi.run import run_equilibrium

    base_cfg = _load_base_cfg(cfg_path)
    cfg = copy.deepcopy(base_cfg)
    cfg["mesh"]["resolution"] = [int(N)]

    refs = _depletion_references(cfg)

    t0 = time.perf_counter()
    result = run_equilibrium(cfg)
    dt = time.perf_counter() - t0

    info = result.solver_info or {}
    if not info.get("converged", False):
        raise RuntimeError(
            f"pn_1d equilibrium did not converge at N={N}: "
            f"reason={info.get('reason')}, iters={info.get('iterations')}"
        )

    V_bi_sim, E_peak_sim, W_sim = _extract_observables(result, refs)
    extents = cfg["mesh"]["extents"][0]
    L = float(extents[1] - extents[0])
    h = L / float(N)
    n_dofs = int(result.V.dofmap.index_map.size_global)

    def _rel(sim: float, ref: float) -> float:
        if not np.isfinite(sim) or ref == 0.0:
            return float("nan")
        return abs(sim - ref) / abs(ref)

    return MeshConvRow(
        N=int(N),
        h=h,
        n_dofs=n_dofs,
        V_bi_sim=V_bi_sim,
        V_bi_theory=refs["V_bi"],
        err_Vbi_rel=_rel(V_bi_sim, refs["V_bi"]),
        E_peak_sim=E_peak_sim,
        E_peak_theory=refs["E_peak"],
        err_Epeak_rel=_rel(E_peak_sim, refs["E_peak"]),
        W_sim=W_sim,
        W_theory=refs["W"],
        err_W_rel=_rel(W_sim, refs["W"]),
        newton_iters=int(info.get("iterations", -1)),
        solve_time_s=float(dt),
    )


def run_convergence_study(
    Ns: list[int],
    *,
    cfg_path: str | Path | None = None,
) -> list[MeshConvRow]:
    """Run the benchmark at each N in Ns (in input order)."""
    return [run_one_level(N, cfg_path=cfg_path) for N in Ns]


def _ratio(prev: float | None, cur: float) -> float:
    if prev is None or not np.isfinite(cur) or cur <= 0.0 or not np.isfinite(prev):
        return float("nan")
    return prev / cur


def cauchy_errors(
    results: list[MeshConvRow],
    *,
    reference_index: int = -1,
) -> list[dict[str, float]]:
    """
    Self-convergence errors vs the reference mesh (default: finest).

    Cauchy errors remove the physics-model gap (depletion approx vs
    true Poisson-Boltzmann) from the convergence metric. They are the
    standard V&V choice for problems without an analytical solution
    (Roache, Oberkampf-Roy). Error on the reference row itself is NaN
    because it has no finer reference.
    """
    if not results:
        return []
    ref = results[reference_index]
    out: list[dict[str, float]] = []
    for i, r in enumerate(results):
        if i == reference_index % len(results):
            out.append({
                "err_Epeak_cauchy": float("nan"),
                "err_W_cauchy": float("nan"),
            })
            continue
        dE = abs(r.E_peak_sim - ref.E_peak_sim)
        dW = abs(r.W_sim - ref.W_sim) if np.isfinite(r.W_sim) and np.isfinite(ref.W_sim) else float("nan")
        out.append({
            "err_Epeak_cauchy": dE / abs(ref.E_peak_sim) if ref.E_peak_sim != 0.0 else float("nan"),
            "err_W_cauchy": dW / abs(ref.W_sim) if ref.W_sim and np.isfinite(ref.W_sim) else float("nan"),
        })
    return out


def to_table_rows(results: list[MeshConvRow]) -> list[dict]:
    """
    Convert dataclass rows to dicts with ratio-per-doubling columns and
    Cauchy self-convergence columns against the finest mesh in the sweep.
    """
    rows: list[dict] = []
    prev_Vbi: float | None = None
    prev_E: float | None = None
    prev_W: float | None = None
    prev_E_cauchy: float | None = None
    prev_W_cauchy: float | None = None
    cauchy = cauchy_errors(results, reference_index=-1)
    for r, c in zip(results, cauchy, strict=True):
        rows.append({
            "N": r.N,
            "h": r.h,
            "N_dofs": r.n_dofs,
            "V_bi_sim": r.V_bi_sim,
            "V_bi_theory": r.V_bi_theory,
            "err_Vbi_rel": r.err_Vbi_rel,
            "err_Vbi_ratio": _ratio(prev_Vbi, r.err_Vbi_rel),
            "E_peak_sim": r.E_peak_sim,
            "E_peak_theory": r.E_peak_theory,
            "err_Epeak_rel": r.err_Epeak_rel,
            "err_Epeak_ratio": _ratio(prev_E, r.err_Epeak_rel),
            "err_Epeak_cauchy": c["err_Epeak_cauchy"],
            "err_Epeak_cauchy_ratio": _ratio(prev_E_cauchy, c["err_Epeak_cauchy"]),
            "W_sim": r.W_sim,
            "W_theory": r.W_theory,
            "err_W_rel": r.err_W_rel,
            "err_W_ratio": _ratio(prev_W, r.err_W_rel),
            "err_W_cauchy": c["err_W_cauchy"],
            "err_W_cauchy_ratio": _ratio(prev_W_cauchy, c["err_W_cauchy"]),
            "newton_iters": r.newton_iters,
            "solve_time_s": r.solve_time_s,
        })
        if np.isfinite(r.err_Vbi_rel) and r.err_Vbi_rel > 0.0:
            prev_Vbi = r.err_Vbi_rel
        if np.isfinite(r.err_Epeak_rel) and r.err_Epeak_rel > 0.0:
            prev_E = r.err_Epeak_rel
        if np.isfinite(r.err_W_rel) and r.err_W_rel > 0.0:
            prev_W = r.err_W_rel
        if np.isfinite(c["err_Epeak_cauchy"]) and c["err_Epeak_cauchy"] > 0.0:
            prev_E_cauchy = c["err_Epeak_cauchy"]
        if np.isfinite(c["err_W_cauchy"]) and c["err_W_cauchy"] > 0.0:
            prev_W_cauchy = c["err_W_cauchy"]
    return rows


CSV_COLUMNS = [
    "N", "h", "N_dofs",
    "V_bi_sim", "V_bi_theory", "err_Vbi_rel", "err_Vbi_ratio",
    "E_peak_sim", "E_peak_theory",
    "err_Epeak_rel", "err_Epeak_ratio",
    "err_Epeak_cauchy", "err_Epeak_cauchy_ratio",
    "W_sim", "W_theory",
    "err_W_rel", "err_W_ratio",
    "err_W_cauchy", "err_W_cauchy_ratio",
    "newton_iters", "solve_time_s",
]


def write_artifacts(
    rows: list[dict],
    out_dir: Path,
    *,
    csv_name: str = "pn_1d.csv",
    plot_name: str = "convergence.png",
    title: str = "pn_1d mesh convergence vs depletion approximation",
) -> None:
    """Write the mesh-convergence CSV and log-log convergence plot."""
    from ._convergence import write_convergence_csv, write_loglog_plot

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    write_convergence_csv(out_dir / csv_name, rows, CSV_COLUMNS)

    hs = [r["h"] for r in rows]
    series = {
        "E_peak rel err (vs depletion approx)": [r["err_Epeak_rel"] for r in rows],
        "W rel err (vs depletion approx)": [r["err_W_rel"] for r in rows],
        "E_peak Cauchy err (vs finest FEM)": [r["err_Epeak_cauchy"] for r in rows],
        "W Cauchy err (vs finest FEM)": [r["err_W_cauchy"] for r in rows],
    }
    write_loglog_plot(
        out_dir / plot_name, hs, series,
        title=title,
        theoretical_rates={"E_peak Cauchy err (vs finest FEM)": 2.0},
    )


def report_table(rows: list[dict], header: str = "") -> str:
    """Stdout-friendly table (subset of columns for readability)."""
    from ._convergence import format_table

    cols = [
        "N", "h",
        "err_Epeak_rel", "err_Epeak_ratio",
        "err_Epeak_cauchy", "err_Epeak_cauchy_ratio",
        "err_W_rel", "err_W_ratio",
        "err_W_cauchy", "err_W_cauchy_ratio",
        "newton_iters", "solve_time_s",
    ]
    return format_table(rows, cols, header=header)


# ---------------------------------------------------------------------------
# CLI entry point used by scripts/run_verification.py
# ---------------------------------------------------------------------------
CLI_NS = [50, 100, 200, 400, 800, 1600]


def run_cli_study(out_dir: Path) -> list[dict]:
    """Run the full mesh sweep and write CSV + PNG to `out_dir`."""
    out_dir = Path(out_dir)
    results = run_convergence_study(CLI_NS)
    rows = to_table_rows(results)
    write_artifacts(rows, out_dir)
    return rows
