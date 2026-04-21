"""
Method of Manufactured Solutions (MMS) for the coupled drift-diffusion
block residual.

The production residual (`semi/physics/drift_diffusion.py`) assembles
three coupled blocks in Slotboom form:

    Poisson:   -div( L_D^2 eps_r grad psi_hat )            - (p_hat - n_hat + N_hat) = 0
    Electron:  -div( L_0^2 mu_n_hat n_hat grad phi_n_hat ) - R_hat                   = 0
    Hole:      -div( L_0^2 mu_p_hat p_hat grad phi_p_hat ) + R_hat                   = 0

with Slotboom `n_hat = ni_hat exp(psi - phi_n)`, `p_hat = ni_hat exp(phi_p - psi)`
and SRH `R_hat`. For MMS we set `N_hat = 0` and subtract a weak-form
manufactured source from each block so that a chosen smooth exact triple
`(psi_e, phi_n_e, phi_p_e)` is the solution of the modified coupled
system. Each block's discretization error then decays at the FE rate.

Three variants exercise progressively more of the residual:

    A  `phi_n_e = phi_p_e = 0`, `tau_hat = infinity`:   Poisson block only
    B  full triple, `tau_hat = infinity`:               full coupling, R ~ 0
    C  full triple, Si lifetimes:                       full coupling with SRH

See `docs/mms_dd_derivation.md` for the approved derivation; every weak
source below tracks the sections of that document.
"""
from __future__ import annotations

import math
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .mms_poisson import (
    C_0_REF,
    EPS_R_DEFAULT,
    L_0_REF,
    T_REF,
    _all_boundary_facets,
    _build_mesh,
    build_mms_scaling,
)

# ---------------------------------------------------------------------------
# Module-level constants. See derivation Section 2.
# ---------------------------------------------------------------------------
VARIANTS: tuple[str, ...] = ("A", "B", "C")

#: Default amplitude triple `(A_psi, A_n, A_p)` used by the pytest gate.
#: `A_p = -A_n` ensures Variant C's SRH numerator
#: `ni_hat^2 * (exp(phi_p - phi_n) - 1)` is pointwise nonzero (derivation 1).
DEFAULT_AMPS: tuple[float, float, float] = (0.5, 0.3, -0.3)

#: Larger-amplitude triple reserved for the CLI sweep. Keeps |exp(...)|
#: bounded by exp(1.5) ~ 4.48, comfortably inside PETSc double precision.
NONLINEAR_AMPS: tuple[float, float, float] = (1.0, 0.5, -0.5)

#: Variant C SRH lifetime (Si mid-gap default, seconds).
TAU_SI_S: float = 1.0e-7

#: Variant A/B "off" lifetime in scaled units. Keeps the SRH code path
#: live but drives `R_e` below machine precision.
TAU_OFF_HAT: float = 1.0e+20

#: Mobility ratios pinned across the mesh sweep so convergence measures
#: discretization error only. `mu_n_over_mu0 = 1.0` because mu_0 is the
#: electron mobility itself; `mu_p / mu_n = 0.45 / 1.4` from Si defaults
#: in `semi/materials.py`.
MU_N_OVER_MU0: float = 1.0
MU_P_OVER_MU0: float = 0.45 / 1.4


@dataclass(frozen=True)
class MMSDDCase:
    """Specification for one mesh-level MMS-DD run."""

    dim: int                                 # 1 or 2
    N: int                                   # cells per side
    variant: str                             # "A", "B", or "C"
    L: float = L_0_REF                       # device length, m
    A_psi: float = DEFAULT_AMPS[0]
    A_n: float = DEFAULT_AMPS[1]
    A_p: float = DEFAULT_AMPS[2]
    eps_r: float = EPS_R_DEFAULT
    cell_kind: str = "triangle"              # 2D only


@dataclass(frozen=True)
class MMSDDResult:
    """Per-level numerical result with per-block error norms."""

    dim: int
    N: int
    variant: str
    h: float
    n_dofs: int                              # identical across the three blocks
    e_L2_psi: float
    e_H1_psi: float
    e_L2_phi_n: float
    e_H1_phi_n: float
    e_L2_phi_p: float
    e_H1_phi_p: float
    snes_iters: int
    solve_time_s: float
    cell_kind: str = "triangle"
    A_psi: float = DEFAULT_AMPS[0]
    A_n: float = DEFAULT_AMPS[1]
    A_p: float = DEFAULT_AMPS[2]


# ---------------------------------------------------------------------------
# Exact triple and manufactured weak sources (derivation Sections 1, 3)
# ---------------------------------------------------------------------------


def _exact_triple_ufl(mesh, dim: int, L: float,
                      A_psi: float, A_n: float, A_p: float,
                      variant: str):
    """
    Return `(psi_e, phi_n_e, phi_p_e)` as UFL expressions.

    For Variants B and C we use the full sin-product triple; for Variant
    A the Fermi levels are identically zero. This choice matters for the
    weak-source construction in `_build_weak_sources`: Variant A only
    produces `f_psi_weak` (the other two manufactured sources are zero
    and are not subtracted from the production forms, which avoids
    compiling zero integrands).
    """
    import ufl

    x = ufl.SpatialCoordinate(mesh)
    if dim == 1:
        psi_e = A_psi * ufl.sin(2.0 * ufl.pi * x[0] / L)
        if variant == "A":
            phi_n_e = None
            phi_p_e = None
        else:
            phi_n_e = A_n * ufl.sin(ufl.pi * x[0] / L)
            phi_p_e = A_p * ufl.sin(ufl.pi * x[0] / L)
    elif dim == 2:
        psi_e = (
            A_psi
            * ufl.sin(2.0 * ufl.pi * x[0] / L)
            * ufl.sin(3.0 * ufl.pi * x[1] / L)
        )
        if variant == "A":
            phi_n_e = None
            phi_p_e = None
        else:
            phi_n_e = (
                A_n
                * ufl.sin(ufl.pi * x[0] / L)
                * ufl.sin(ufl.pi * x[1] / L)
            )
            phi_p_e = (
                A_p
                * ufl.sin(ufl.pi * x[0] / L)
                * ufl.sin(ufl.pi * x[1] / L)
            )
    else:
        raise ValueError(f"Unsupported dim {dim}")
    return psi_e, phi_n_e, phi_p_e


def _build_weak_sources(
    mesh,
    spaces,
    sc,
    *,
    variant: str,
    eps_r_value: float,
    mu_n_over_mu0: float,
    mu_p_over_mu0: float,
    tau_n_hat: float,
    tau_p_hat: float,
    psi_e,
    phi_n_e,
    phi_p_e,
):
    """
    Build the three weak-form manufactured sources.

    For Variants B and C this returns `(f_psi_weak, f_n_weak, f_p_weak)`.
    For Variant A only `f_psi_weak` is built and the continuity-block
    entries are `None` -- the caller must leave the production continuity
    forms unmodified (they are already satisfied by `(psi_e, 0, 0)` because
    `grad(phi_n_e) = 0` and `R_e = 0` there).

    The weak-form construction follows Section 3 of
    `docs/mms_dd_derivation.md`. We never form `ufl.div(ufl.grad(...))`
    directly; integrating by parts against the test function uses the
    same gradient discretization as the production residual, so MMS
    verifies the exact discretization in use (compare
    `mms_poisson.manufactured_source_weak`).
    """
    import ufl
    from dolfinx import fem
    from petsc4py import PETSc

    from semi.physics.slotboom import n_from_slotboom, p_from_slotboom

    v_psi = ufl.TestFunction(spaces.V_psi)
    v_n = ufl.TestFunction(spaces.V_phi_n)
    v_p = ufl.TestFunction(spaces.V_phi_p)

    L_D2 = fem.Constant(mesh, PETSc.ScalarType(sc.lambda2 * sc.L0 ** 2))
    L0_sq = fem.Constant(mesh, PETSc.ScalarType(sc.L0 ** 2))
    eps_r = fem.Constant(mesh, PETSc.ScalarType(eps_r_value))
    ni_hat = fem.Constant(mesh, PETSc.ScalarType(sc.n_i / sc.C0))
    mu_n = fem.Constant(mesh, PETSc.ScalarType(mu_n_over_mu0))
    mu_p = fem.Constant(mesh, PETSc.ScalarType(mu_p_over_mu0))
    tau_n = fem.Constant(mesh, PETSc.ScalarType(tau_n_hat))
    tau_p = fem.Constant(mesh, PETSc.ScalarType(tau_p_hat))

    if variant == "A":
        # Variant A: phi_n_e = phi_p_e = 0, so n_e * p_e = ni_hat^2 exactly
        # and R_e vanishes pointwise. Slotboom densities reduce to
        # ni_hat * exp(+/- psi_e), matching the Poisson MMS carrier term.
        n_e = ni_hat * ufl.exp(psi_e)
        p_e = ni_hat * ufl.exp(-psi_e)
        f_psi_weak = (
            L_D2 * eps_r * ufl.inner(ufl.grad(psi_e), ufl.grad(v_psi))
            - (p_e - n_e) * v_psi
        ) * ufl.dx
        return f_psi_weak, None, None

    # Variants B and C: full triple.
    n_e = n_from_slotboom(psi_e, phi_n_e, ni_hat)
    p_e = p_from_slotboom(psi_e, phi_p_e, ni_hat)

    # SRH rate, inlined to share the same ni_hat Constant with the Poisson
    # block (matches production code in `build_dd_block_residual`). E_t = 0
    # is enforced module-wide, so n1 = p1 = ni_hat.
    n1 = ni_hat * math.exp(0.0)
    p1 = ni_hat * math.exp(-0.0)
    R_e = (n_e * p_e - ni_hat * ni_hat) / (
        tau_p * (n_e + n1) + tau_n * (p_e + p1)
    )

    f_psi_weak = (
        L_D2 * eps_r * ufl.inner(ufl.grad(psi_e), ufl.grad(v_psi))
        - (p_e - n_e) * v_psi
    ) * ufl.dx
    f_n_weak = (
        L0_sq * mu_n * n_e * ufl.inner(ufl.grad(phi_n_e), ufl.grad(v_n))
        - R_e * v_n
    ) * ufl.dx
    f_p_weak = (
        L0_sq * mu_p * p_e * ufl.inner(ufl.grad(phi_p_e), ufl.grad(v_p))
        + R_e * v_p
    ) * ufl.dx
    return f_psi_weak, f_n_weak, f_p_weak


# ---------------------------------------------------------------------------
# Per-level solve (derivation Sections 4, 5)
# ---------------------------------------------------------------------------


def run_one_level(case: MMSDDCase, *, sc=None) -> MMSDDResult:
    """
    Solve the MMS-DD problem on a single mesh level and return per-block
    L^2 and H^1 seminorm errors together with SNES diagnostics.

    Reuses `build_dd_block_residual` so future edits to the production
    residual are automatically exercised. BCs are homogeneous Dirichlet
    on all three fields (exact, not interpolated, because every exact
    factor is `sin(k*pi*x/L)` with `k` a positive integer and so vanishes
    on the boundary; see derivation Section 4).
    """
    import ufl
    from dolfinx import fem
    from petsc4py import PETSc

    from semi.physics.drift_diffusion import (
        build_dd_block_residual,
        make_dd_block_spaces,
    )
    from semi.solver import solve_nonlinear_block

    from ._norms import h1_seminorm_error_squared, l2_error_squared

    if case.variant not in VARIANTS:
        raise ValueError(f"variant {case.variant!r} not in {VARIANTS}")
    if sc is None:
        sc = build_mms_scaling(L=case.L)

    mesh = _build_mesh(case.dim, case.N, case.L, case.cell_kind)
    spaces = make_dd_block_spaces(mesh)

    # N_hat = 0 for MMS; see derivation preamble.
    N_hat_fn = fem.Function(spaces.V_psi, name="N_hat_zero")
    N_hat_fn.x.array[:] = 0.0

    # Initial guess: zero everywhere. Trivially satisfies the
    # homogeneous Dirichlet BC and the charge-neutrality identity
    # (n_init = p_init = ni_hat, rho_init = 0 since N_hat = 0).
    for fn in (spaces.psi, spaces.phi_n, spaces.phi_p):
        fn.x.array[:] = 0.0
        fn.x.scatter_forward()

    # Per-variant lifetime choice (derivation Section 2, SRH row).
    if case.variant == "C":
        tau_n_hat = tau_p_hat = TAU_SI_S / sc.t0
    else:
        tau_n_hat = tau_p_hat = TAU_OFF_HAT

    # Production residual, per-block.
    F_prod = build_dd_block_residual(
        spaces, N_hat_fn, sc,
        case.eps_r, MU_N_OVER_MU0, MU_P_OVER_MU0,
        tau_n_hat, tau_p_hat, 0.0,   # E_t / V_t = 0 (mid-gap traps)
    )

    # Manufactured weak sources.
    psi_e, phi_n_e, phi_p_e = _exact_triple_ufl(
        mesh, case.dim, case.L,
        case.A_psi, case.A_n, case.A_p, case.variant,
    )
    f_psi_weak, f_n_weak, f_p_weak = _build_weak_sources(
        mesh, spaces, sc,
        variant=case.variant,
        eps_r_value=case.eps_r,
        mu_n_over_mu0=MU_N_OVER_MU0,
        mu_p_over_mu0=MU_P_OVER_MU0,
        tau_n_hat=tau_n_hat,
        tau_p_hat=tau_p_hat,
        psi_e=psi_e,
        phi_n_e=phi_n_e,
        phi_p_e=phi_p_e,
    )
    F_psi = F_prod[0] - f_psi_weak
    if case.variant == "A":
        # f_n_weak = f_p_weak = 0. Leave continuity forms unmodified; they
        # are already satisfied by (psi_e, 0, 0) (derivation Section 3.4).
        F_phi_n = F_prod[1]
        F_phi_p = F_prod[2]
    else:
        F_phi_n = F_prod[1] - f_n_weak
        F_phi_p = F_prod[2] - f_p_weak

    # Homogeneous Dirichlet BCs on all three fields, exact.
    fdim = mesh.topology.dim - 1
    mesh.topology.create_connectivity(fdim, mesh.topology.dim)
    boundary_facets = _all_boundary_facets(mesh, fdim)
    zero = fem.Constant(mesh, PETSc.ScalarType(0.0))
    bdofs_psi = fem.locate_dofs_topological(spaces.V_psi, fdim, boundary_facets)
    bdofs_phi_n = fem.locate_dofs_topological(spaces.V_phi_n, fdim, boundary_facets)
    bdofs_phi_p = fem.locate_dofs_topological(spaces.V_phi_p, fdim, boundary_facets)
    bcs = [
        fem.dirichletbc(zero, bdofs_psi, spaces.V_psi),
        fem.dirichletbc(zero, bdofs_phi_n, spaces.V_phi_n),
        fem.dirichletbc(zero, bdofs_phi_p, spaces.V_phi_p),
    ]

    # SNES tolerances: rtol and max_it match the derivation (Section 5),
    # but atol is driven below the continuity-block initial residual
    # scale. Section 7.3 identifies the block-residual scale disparity
    # (L_D^2 * ...  vs  L_0^2 * n_i_hat * ...) as the reason atol must
    # reach below every block's native floor; the derivation's nominal
    # 1e-16 is tuned to the psi block alone and prematurely terminates
    # SNES on the continuity blocks in 2D (where the integration-area
    # L^2 ~ 4e-12 shrinks the phi_n / phi_p residuals to ~1e-18 at
    # initial). Driving atol to ~0 with stol controlling termination
    # lets Newton iterate until every block has actually been solved.
    petsc_options = {
        "snes_rtol": 1.0e-14,
        "snes_atol": 0.0,
        "snes_stol": 1.0e-12,
        "snes_max_it": 80,
    }
    prefix = (
        f"mms_dd_dim{case.dim}_N{case.N}_var{case.variant}"
        f"_amp{case.A_psi:g}_"
    )
    t0 = time.perf_counter()
    info = solve_nonlinear_block(
        [F_psi, F_phi_n, F_phi_p],
        [spaces.psi, spaces.phi_n, spaces.phi_p],
        bcs, prefix=prefix, petsc_options=petsc_options,
    )
    dt = time.perf_counter() - t0
    if not info.get("converged", False):
        raise RuntimeError(
            f"MMS-DD did not converge at dim={case.dim}, N={case.N}, "
            f"variant={case.variant}: reason={info.get('reason')}, "
            f"iters={info.get('iterations')}"
        )

    # Per-block error norms. For Variant A the exact Fermi fields are 0,
    # so `phi_n_e_err` below is a UFL Zero and the error reduces to the
    # discrete function's own L^2 and H^1 seminorm.
    psi_e_err = psi_e
    phi_n_e_err = phi_n_e if phi_n_e is not None else 0.0
    phi_p_e_err = phi_p_e if phi_p_e is not None else 0.0

    e_L2_psi = float(np.sqrt(max(l2_error_squared(spaces.psi, psi_e_err), 0.0)))
    e_H1_psi = float(np.sqrt(max(h1_seminorm_error_squared(spaces.psi, psi_e_err), 0.0)))
    e_L2_phi_n = float(np.sqrt(max(l2_error_squared(spaces.phi_n, phi_n_e_err), 0.0)))
    e_H1_phi_n = float(np.sqrt(max(h1_seminorm_error_squared(spaces.phi_n, phi_n_e_err), 0.0)))
    e_L2_phi_p = float(np.sqrt(max(l2_error_squared(spaces.phi_p, phi_p_e_err), 0.0)))
    e_H1_phi_p = float(np.sqrt(max(h1_seminorm_error_squared(spaces.phi_p, phi_p_e_err), 0.0)))

    n_dofs = int(spaces.V_psi.dofmap.index_map.size_global)
    return MMSDDResult(
        dim=case.dim,
        N=case.N,
        variant=case.variant,
        h=case.L / case.N,
        n_dofs=n_dofs,
        e_L2_psi=e_L2_psi, e_H1_psi=e_H1_psi,
        e_L2_phi_n=e_L2_phi_n, e_H1_phi_n=e_H1_phi_n,
        e_L2_phi_p=e_L2_phi_p, e_H1_phi_p=e_H1_phi_p,
        snes_iters=int(info.get("iterations", -1)),
        solve_time_s=dt,
        cell_kind=case.cell_kind,
        A_psi=case.A_psi, A_n=case.A_n, A_p=case.A_p,
    )


# ---------------------------------------------------------------------------
# Convergence sweeps (derivation Section 6)
# ---------------------------------------------------------------------------


def run_convergence_study(
    *,
    dim: int,
    variant: str,
    Ns: list[int],
    L: float = L_0_REF,
    A_psi: float = DEFAULT_AMPS[0],
    A_n: float = DEFAULT_AMPS[1],
    A_p: float = DEFAULT_AMPS[2],
    eps_r: float = EPS_R_DEFAULT,
    cell_kind: str = "triangle",
    sc=None,
) -> list[MMSDDResult]:
    """Run MMS-DD on a sequence of refinements. Returns one result per N."""
    if sc is None:
        sc = build_mms_scaling(L=L)
    results: list[MMSDDResult] = []
    for N in Ns:
        case = MMSDDCase(
            dim=dim, N=N, variant=variant, L=L,
            A_psi=A_psi, A_n=A_n, A_p=A_p,
            eps_r=eps_r, cell_kind=cell_kind,
        )
        results.append(run_one_level(case, sc=sc))
    return results


# Column order used by the tables / CSVs.
_BLOCKS: tuple[str, ...] = ("psi", "phi_n", "phi_p")
_ROW_COLUMNS: tuple[str, ...] = (
    "N", "h", "N_dofs",
    "e_L2_psi", "rate_L2_psi", "e_H1_psi", "rate_H1_psi",
    "e_L2_phi_n", "rate_L2_phi_n", "e_H1_phi_n", "rate_H1_phi_n",
    "e_L2_phi_p", "rate_L2_phi_p", "e_H1_phi_p", "rate_H1_phi_p",
    "snes_iters", "solve_time_s",
)


def to_table_rows(results: list[MMSDDResult]) -> list[dict]:
    """Build per-level rows with observed L^2/H^1 rates for every block."""
    from ._convergence import observed_rates

    hs = [r.h for r in results]
    rates: dict[str, list[float]] = {}
    for b in _BLOCKS:
        rates[f"L2_{b}"] = observed_rates(hs, [getattr(r, f"e_L2_{b}") for r in results])
        rates[f"H1_{b}"] = observed_rates(hs, [getattr(r, f"e_H1_{b}") for r in results])
    rows: list[dict] = []
    for i, r in enumerate(results):
        row = {
            "N": r.N,
            "h": r.h,
            "N_dofs": r.n_dofs,
            "e_L2_psi": r.e_L2_psi,
            "rate_L2_psi": rates["L2_psi"][i],
            "e_H1_psi": r.e_H1_psi,
            "rate_H1_psi": rates["H1_psi"][i],
            "e_L2_phi_n": r.e_L2_phi_n,
            "rate_L2_phi_n": rates["L2_phi_n"][i],
            "e_H1_phi_n": r.e_H1_phi_n,
            "rate_H1_phi_n": rates["H1_phi_n"][i],
            "e_L2_phi_p": r.e_L2_phi_p,
            "rate_L2_phi_p": rates["L2_phi_p"][i],
            "e_H1_phi_p": r.e_H1_phi_p,
            "rate_H1_phi_p": rates["H1_phi_p"][i],
            "snes_iters": r.snes_iters,
            "solve_time_s": r.solve_time_s,
        }
        rows.append(row)
    return rows


def report_table(rows: list[dict], header: str = "") -> str:
    """Pretty multi-line stdout table showing psi L^2 + rate for each block."""
    from ._convergence import format_table

    cols = [
        "N", "h", "N_dofs",
        "e_L2_psi", "rate_L2_psi",
        "e_L2_phi_n", "rate_L2_phi_n",
        "e_L2_phi_p", "rate_L2_phi_p",
    ]
    return format_table(rows, cols, header=header)


def write_artifacts(
    rows: list[dict],
    out_dir: Path,
    *,
    title: str,
    csv_name: str = "convergence.csv",
    plot_name: str = "convergence.png",
) -> None:
    """Write `convergence.csv` and `convergence.png` to `out_dir`."""
    from ._convergence import write_convergence_csv, write_loglog_plot

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    write_convergence_csv(out_dir / csv_name, rows, list(_ROW_COLUMNS))
    hs = [r["h"] for r in rows]
    series = {
        "L2 psi":    [r["e_L2_psi"] for r in rows],
        "L2 phi_n":  [r["e_L2_phi_n"] for r in rows],
        "L2 phi_p":  [r["e_L2_phi_p"] for r in rows],
    }
    write_loglog_plot(
        out_dir / plot_name, hs, series,
        title=title,
        theoretical_rates={
            "L2 psi": 2.0,
            "L2 phi_n": 2.0,
            "L2 phi_p": 2.0,
        },
    )


# ---------------------------------------------------------------------------
# Default mesh sequences
# ---------------------------------------------------------------------------
PYTEST_NS_1D: list[int] = [40, 80, 160]
PYTEST_NS_2D: list[int] = [16, 32, 64]
CLI_NS_1D: list[int] = [40, 80, 160, 320]
CLI_NS_2D: list[int] = [16, 32, 64]


def run_cli_study(out_dir: Path) -> dict[str, list[dict]]:
    """
    Artifact-production sweep used by `scripts/run_verification.py`.

    Runs nine studies in total:
      - `1d_<variant>_linear`     for variant in A,B,C (default amps, Ns=CLI_NS_1D)
      - `1d_<variant>_nonlinear`  for variant in A,B,C (NONLINEAR_AMPS, Ns=CLI_NS_1D)
      - `2d_<variant>`            for variant in A,B,C (default amps, Ns=CLI_NS_2D)

    Each study writes `convergence.csv` and `convergence.png` into
    `out_dir/<label>/`, mirroring the Poisson MMS artifact layout.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    sc = build_mms_scaling()

    studies: dict[str, list[dict]] = {}

    for variant in VARIANTS:
        # 1D default amplitudes
        res = run_convergence_study(
            dim=1, variant=variant, Ns=CLI_NS_1D,
            A_psi=DEFAULT_AMPS[0], A_n=DEFAULT_AMPS[1], A_p=DEFAULT_AMPS[2],
            sc=sc,
        )
        rows = to_table_rows(res)
        label = f"1d_{variant}_linear"
        write_artifacts(
            rows, out_dir / label,
            title=(
                f"MMS-DD 1D variant {variant} (linear, "
                f"A_psi={DEFAULT_AMPS[0]}, A_n={DEFAULT_AMPS[1]}, "
                f"A_p={DEFAULT_AMPS[2]})"
            ),
        )
        studies[label] = rows

        # 1D nonlinear amplitudes
        res = run_convergence_study(
            dim=1, variant=variant, Ns=CLI_NS_1D,
            A_psi=NONLINEAR_AMPS[0], A_n=NONLINEAR_AMPS[1], A_p=NONLINEAR_AMPS[2],
            sc=sc,
        )
        rows = to_table_rows(res)
        label = f"1d_{variant}_nonlinear"
        write_artifacts(
            rows, out_dir / label,
            title=(
                f"MMS-DD 1D variant {variant} (nonlinear, "
                f"A_psi={NONLINEAR_AMPS[0]}, A_n={NONLINEAR_AMPS[1]}, "
                f"A_p={NONLINEAR_AMPS[2]})"
            ),
        )
        studies[label] = rows

        # 2D default amplitudes, triangles
        res = run_convergence_study(
            dim=2, variant=variant, Ns=CLI_NS_2D,
            A_psi=DEFAULT_AMPS[0], A_n=DEFAULT_AMPS[1], A_p=DEFAULT_AMPS[2],
            cell_kind="triangle", sc=sc,
        )
        rows = to_table_rows(res)
        label = f"2d_{variant}"
        write_artifacts(
            rows, out_dir / label,
            title=(
                f"MMS-DD 2D variant {variant} (right-diagonal triangles, "
                f"A_psi={DEFAULT_AMPS[0]}, A_n={DEFAULT_AMPS[1]}, "
                f"A_p={DEFAULT_AMPS[2]})"
            ),
        )
        studies[label] = rows

    return studies
