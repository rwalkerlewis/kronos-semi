"""
Method of Manufactured Solutions (MMS) for the scaled equilibrium
Poisson equation.

The production residual (semi/physics/poisson.py) has the shape

    -div( L_D^2 eps_r grad psi_hat ) - rho_hat = 0
    rho_hat = ni_hat * (exp(-psi_hat) - exp(psi_hat)) + N_hat

For MMS we choose a smooth analytical exact solution psi_exact(x), set
N_hat = 0 (no doping interference), and add a manufactured source f_hat
constructed in UFL so that the modified residual

    F_mms = F_production(psi)  -  f_hat * v * dx
    f_hat = -div( L_D^2 eps_r grad psi_exact )
            - ni_hat * (exp(-psi_exact) - exp(psi_exact))

drives Newton to psi_h -> psi_exact. We then measure

    e_L2 = || psi_h - psi_exact ||_{L^2}
    e_H1 = || grad(psi_h - psi_exact) ||_{L^2}

over a sequence of mesh refinements and report observed convergence
rates p = log(e_prev / e_cur) / log(h_prev / h_cur). For P1 Lagrange on
a smooth problem the theoretical rates are 2 in L^2 and 1 in H^1; the
test gates assert finest-pair p within 0.15 of theoretical.

Boundary conditions are homogeneous Dirichlet enforced via
`dolfinx.fem.dirichletbc(zero_constant, ...)`. Both the 1D and 2D
manufactured solutions vanish identically on the boundary, so this is
the exact-BC choice and avoids polluting the MMS error with boundary
interpolation noise.
"""
from __future__ import annotations

import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Module-level reference scales for the MMS Scaling.
#
# These are deliberately hard-coded so the MMS test is self-contained and
# does not drift with benchmark JSON edits. They are picked so:
#   - V_T_REF = thermal voltage at 300 K (standard);
#   - L_0_REF = 2 um device length matches the pn_1d benchmark length scale;
#   - C_0_REF = 1e16 cm^-3 is the typical reference doping for that benchmark.
#
# The combination yields lambda2 ~ 1.6e-3 << 1, so the dimensionless Debye
# length is much shorter than L_0 and the exp(+/- psi_hat) carrier term is
# the dominant nonlinear contribution to the residual: exactly the regime we
# want MMS to exercise. A sanity assertion below enforces lambda2 < 1 so
# this property is self-checking against future scale changes.
# ---------------------------------------------------------------------------
T_REF: float = 300.0                      # K
L_0_REF: float = 2.0e-6                    # m, matches pn_1d benchmark
C_0_REF: float = 1.0e22                    # m^-3 = 1e16 cm^-3, typical doping
EPS_R_DEFAULT: float = 11.7                # Si


@dataclass(frozen=True)
class MMSPoissonCase:
    """Specification for one mesh-level MMS Poisson run."""

    dim: int                                # 1 or 2
    N: int                                  # cells per side
    L: float = L_0_REF                      # device length, m
    amplitude: float = 0.5                  # psi_exact amplitude, scaled units
    eps_r: float = EPS_R_DEFAULT
    cell_kind: str = "triangle"             # "triangle" or "quadrilateral" (2D only)


@dataclass(frozen=True)
class MMSPoissonResult:
    """Per-level numerical result."""

    dim: int
    N: int
    h: float
    n_dofs: int
    e_L2: float
    e_H1: float
    snes_iters: int
    solve_time_s: float
    cell_kind: str = "triangle"
    amplitude: float = 0.5


def build_mms_scaling(L: float = L_0_REF, C0: float = C_0_REF, T: float = T_REF):
    """
    Construct a minimal `Scaling` for MMS independent of any JSON config.

    Uses Si reference parameters (mu_n, n_i) so the resulting `lambda2`
    and `n_i / C0` are physically representative. Asserts `lambda2 < 1`
    so the carrier nonlinearity remains non-trivial; an MMS where the
    nonlinear term is negligible would only verify the Laplacian, not
    the full residual.
    """
    from semi.materials import get_material
    from semi.scaling import Scaling

    mat = get_material("Si")
    sc = Scaling(L0=L, C0=C0, T=T, mu0=mat.mu_n, n_i=mat.n_i)
    if sc.lambda2 >= 1.0:
        raise ValueError(
            f"MMS scaling has lambda2={sc.lambda2:.3e} >= 1; the carrier "
            "nonlinearity becomes negligible and MMS would not exercise the "
            "exp(+/- psi) terms. Pick smaller L or larger C0."
        )
    return sc


def _build_mesh(dim: int, N: int, L: float, cell_kind: str = "triangle"):
    """Build a unit-cell-uniform mesh on [0, L]^dim."""
    import numpy as np
    from dolfinx import mesh as dmesh
    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    if dim == 1:
        return dmesh.create_interval(comm, N, [0.0, L])
    if dim == 2:
        if cell_kind == "triangle":
            ctype = dmesh.CellType.triangle
        elif cell_kind == "quadrilateral":
            ctype = dmesh.CellType.quadrilateral
        else:
            raise ValueError(f"Unknown cell_kind {cell_kind!r}")
        # dolfinx 0.10 requires a DiagonalType keyword even for quads (it is
        # ignored by the quad code path but the binding does not allow None).
        return dmesh.create_rectangle(
            comm,
            [np.array([0.0, 0.0]), np.array([L, L])],
            [N, N],
            cell_type=ctype,
            diagonal=dmesh.DiagonalType.right,
        )
    raise ValueError(f"Unsupported dim {dim}")


def psi_exact_ufl(mesh, dim: int, L: float, amplitude: float):
    """
    Smooth manufactured exact solution that vanishes on the entire boundary.

    1D:  amplitude * sin(2*pi*x/L)         zeros at x=0, L
    2D:  amplitude * sin(2*pi*x/L) * sin(3*pi*y/L)   zeros on all 4 edges

    The all-sin 2D form is required: sin*cos is nonzero at y=0 and would
    require non-homogeneous Dirichlet BCs, polluting the MMS error.
    """
    import ufl

    x = ufl.SpatialCoordinate(mesh)
    if dim == 1:
        return amplitude * ufl.sin(2.0 * ufl.pi * x[0] / L)
    if dim == 2:
        return (
            amplitude
            * ufl.sin(2.0 * ufl.pi * x[0] / L)
            * ufl.sin(3.0 * ufl.pi * x[1] / L)
        )
    raise ValueError(f"Unsupported dim {dim}")


def manufactured_source_weak(mesh, psi_e_ufl, v, sc, eps_r_value: float):
    """
    Build the weak form of the manufactured source contribution.

    Conceptually we want to subtract `f_hat * v * dx` from the production
    residual where

        f_hat = -div( L_D^2 eps_r grad(psi_exact) )
                - ni_hat * (exp(-psi_exact) - exp(psi_exact))

    Integrated against the test function and using integration by parts on
    the diffusion term (psi_exact vanishes on the boundary), this becomes

        ( L_D^2 eps_r grad(psi_exact) . grad(v)
          - ni_hat * (exp(-psi_exact) - exp(psi_exact)) * v ) dx

    Working in weak form avoids `ufl.div(ufl.grad(...))` which can compile
    to numerically zero on coarse meshes when constant prefactors get
    pulled outside, and it uses the same gradient discretization as the
    production residual (so MMS exactly verifies the discretization in
    use, not a slightly different one).
    """
    import ufl
    from dolfinx import fem
    from petsc4py import PETSc

    L_D2 = fem.Constant(mesh, PETSc.ScalarType(sc.lambda2 * sc.L0 ** 2))
    eps_r = fem.Constant(mesh, PETSc.ScalarType(eps_r_value))
    ni_hat = fem.Constant(mesh, PETSc.ScalarType(sc.n_i / sc.C0))

    return (
        L_D2 * eps_r * ufl.inner(ufl.grad(psi_e_ufl), ufl.grad(v))
        - ni_hat * (ufl.exp(-psi_e_ufl) - ufl.exp(psi_e_ufl)) * v
    ) * ufl.dx


def run_one_level(case: MMSPoissonCase, *, sc=None) -> MMSPoissonResult:
    """
    Solve the MMS-Poisson problem on a single mesh level and return the
    error norms together with SNES diagnostics.

    Reuses `build_equilibrium_poisson_form` so any future change to the
    production residual is automatically exercised by MMS.
    """
    import ufl
    from dolfinx import fem
    from petsc4py import PETSc

    from semi.physics.poisson import build_equilibrium_poisson_form
    from semi.solver import solve_nonlinear

    from ._norms import h1_seminorm_error_squared, l2_error_squared

    if sc is None:
        sc = build_mms_scaling(L=case.L)

    mesh = _build_mesh(case.dim, case.N, case.L, case.cell_kind)
    V = fem.functionspace(mesh, ("Lagrange", 1))
    v = ufl.TestFunction(V)

    # Manufactured exact solution (UFL expression on this mesh).
    psi_e = psi_exact_ufl(mesh, case.dim, case.L, case.amplitude)

    # N_hat = 0 to keep the focus on discretization vs the carrier term.
    N_hat = fem.Function(V, name="N_hat_zero")
    N_hat.x.array[:] = 0.0

    # Initial guess: zero. Trivially satisfies the homogeneous BC.
    psi = fem.Function(V, name="psi_hat_mms")
    psi.x.array[:] = 0.0

    # Homogeneous Dirichlet on the full boundary, exact (not interpolated)
    # because psi_exact vanishes there.
    fdim = mesh.topology.dim - 1
    mesh.topology.create_connectivity(fdim, mesh.topology.dim)
    boundary_facets = _all_boundary_facets(mesh, fdim)
    bdofs = fem.locate_dofs_topological(V, fdim, boundary_facets)
    zero = fem.Constant(mesh, PETSc.ScalarType(0.0))
    bcs = [fem.dirichletbc(zero, bdofs, V)]

    # Residual: production form minus the weak-form manufactured source.
    F_prod = build_equilibrium_poisson_form(V, psi, N_hat, sc, case.eps_r)
    src_weak = manufactured_source_weak(mesh, psi_e, v, sc, case.eps_r)
    F = F_prod - src_weak

    # Tighten SNES tolerances so the carrier nonlinearity is fully iterated.
    # The production default snes_atol = 1e-12 lets Newton stop after one
    # step on this problem (the residual norm is dominated by the small
    # L_D^2 prefactor), which leaves an O(amplitude^2) iteration error
    # that floors the L^2 convergence rate. Driving atol to ~machine zero
    # forces enough iterations that the L^2 error reflects pure FE
    # discretization, which is what MMS is supposed to measure.
    petsc_options = {
        "snes_rtol": 1.0e-14,
        "snes_atol": 1.0e-16,
        "snes_stol": 0.0,
        "snes_max_it": 50,
    }
    t0 = time.perf_counter()
    info = solve_nonlinear(
        F, psi, bcs,
        prefix=f"mms_poisson_dim{case.dim}_N{case.N}_amp{case.amplitude:g}_",
        petsc_options=petsc_options,
    )
    dt = time.perf_counter() - t0
    if not info.get("converged", False):
        raise RuntimeError(
            f"MMS-Poisson did not converge at dim={case.dim}, N={case.N}: "
            f"reason={info.get('reason')}, iters={info.get('iterations')}"
        )

    e_L2 = float(np.sqrt(max(l2_error_squared(psi, psi_e), 0.0)))
    e_H1 = float(np.sqrt(max(h1_seminorm_error_squared(psi, psi_e), 0.0)))

    n_dofs = int(V.dofmap.index_map.size_global)
    return MMSPoissonResult(
        dim=case.dim,
        N=case.N,
        h=case.L / case.N,
        n_dofs=n_dofs,
        e_L2=e_L2,
        e_H1=e_H1,
        snes_iters=int(info.get("iterations", -1)),
        solve_time_s=dt,
        cell_kind=case.cell_kind,
        amplitude=case.amplitude,
    )


def _all_boundary_facets(mesh, fdim: int):
    """Indices of every facet on the exterior boundary."""
    from dolfinx import mesh as dmesh

    return dmesh.locate_entities_boundary(
        mesh, fdim, lambda x: np.full(x.shape[1], True)
    )


# ---------------------------------------------------------------------------
# Convergence sweeps (used by tests and the CLI)
# ---------------------------------------------------------------------------

def run_convergence_study(
    *,
    dim: int,
    Ns: list[int],
    L: float = L_0_REF,
    amplitude: float = 0.5,
    eps_r: float = EPS_R_DEFAULT,
    cell_kind: str = "triangle",
    sc=None,
) -> list[MMSPoissonResult]:
    """
    Run MMS-Poisson on a sequence of refinements. Returns one result per
    mesh size, in input order.
    """
    if sc is None:
        sc = build_mms_scaling(L=L)
    results: list[MMSPoissonResult] = []
    for N in Ns:
        case = MMSPoissonCase(
            dim=dim, N=N, L=L, amplitude=amplitude,
            eps_r=eps_r, cell_kind=cell_kind,
        )
        results.append(run_one_level(case, sc=sc))
    return results


def to_table_rows(results: list[MMSPoissonResult]) -> list[dict]:  # pragma: no cover
    """Build the (h, N_dofs, e_L2, rate_L2, e_H1, rate_H1) rows."""
    from ._convergence import observed_rates

    hs = [r.h for r in results]
    e_L2s = [r.e_L2 for r in results]
    e_H1s = [r.e_H1 for r in results]
    r_L2 = observed_rates(hs, e_L2s)
    r_H1 = observed_rates(hs, e_H1s)
    rows = []
    for i, r in enumerate(results):
        rows.append({
            "N": r.N,
            "h": r.h,
            "N_dofs": r.n_dofs,
            "e_L2": r.e_L2,
            "rate_L2": r_L2[i],
            "e_H1": r.e_H1,
            "rate_H1": r_H1[i],
            "snes_iters": r.snes_iters,
            "solve_time_s": r.solve_time_s,
        })
    return rows


def write_artifacts(  # pragma: no cover
    rows: list[dict],
    out_dir: Path,
    *,
    title: str,
    csv_name: str = "convergence.csv",
    plot_name: str = "convergence.png",
) -> None:
    """Write convergence.csv and convergence.png to `out_dir`."""
    from ._convergence import write_convergence_csv, write_loglog_plot

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    columns = [
        "N", "h", "N_dofs",
        "e_L2", "rate_L2",
        "e_H1", "rate_H1",
        "snes_iters", "solve_time_s",
    ]
    write_convergence_csv(out_dir / csv_name, rows, columns)
    hs = [r["h"] for r in rows]
    series = {
        "L2 error": [r["e_L2"] for r in rows],
        "H1 seminorm error": [r["e_H1"] for r in rows],
    }
    write_loglog_plot(
        out_dir / plot_name, hs, series,
        title=title,
        theoretical_rates={"L2 error": 2.0, "H1 seminorm error": 1.0},
    )


def report_table(rows: list[dict], header: str = "") -> str:  # pragma: no cover
    """Pretty multi-line stdout table."""
    from ._convergence import format_table

    cols = ["N", "h", "N_dofs", "e_L2", "rate_L2", "e_H1", "rate_H1"]
    return format_table(rows, cols, header=header)


# ---------------------------------------------------------------------------
# Default mesh sequences for the CLI (broader than what tests use).
# ---------------------------------------------------------------------------
CLI_NS_1D = [40, 80, 160, 320, 640]
CLI_NS_2D = [16, 32, 64, 128]


def run_cli_study(out_dir: Path) -> dict[str, list[dict]]:  # pragma: no cover
    """
    Run the artifact-production sweep used by `scripts/run_verification.py`.

    Returns a dict mapping study label to row list. Three studies are run:
    1D linear (amplitude 0.5), 1D nonlinear (amplitude 2.0), and 2D
    triangles (amplitude 0.5). A separate 2D quadrilateral smoke test at
    N=64 is run and its single-level error is appended to a comparison
    file.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    sc = build_mms_scaling()

    studies: dict[str, list[dict]] = {}

    # 1D linear regime (amplitude 0.5)
    res = run_convergence_study(
        dim=1, Ns=CLI_NS_1D, amplitude=0.5, sc=sc,
    )
    rows = to_table_rows(res)
    write_artifacts(
        rows, out_dir / "1d_linear",
        title="MMS Poisson 1D linear regime (amplitude 0.5)",
    )
    studies["1d_linear"] = rows

    # 1D nonlinear regime (amplitude 2.0)
    res = run_convergence_study(
        dim=1, Ns=CLI_NS_1D, amplitude=2.0, sc=sc,
    )
    rows = to_table_rows(res)
    write_artifacts(
        rows, out_dir / "1d_nonlinear",
        title="MMS Poisson 1D nonlinear regime (amplitude 2.0)",
    )
    studies["1d_nonlinear"] = rows

    # 2D triangles (amplitude 0.5)
    res = run_convergence_study(
        dim=2, Ns=CLI_NS_2D, amplitude=0.5, sc=sc, cell_kind="triangle",
    )
    rows = to_table_rows(res)
    write_artifacts(
        rows, out_dir / "2d_triangles",
        title="MMS Poisson 2D right-diagonal triangles (amplitude 0.5)",
    )
    studies["2d_triangles"] = rows

    # 2D quadrilateral smoke test at N=64; compare against triangle row at
    # the same N. The quad mesh on a square produces fewer DOFs per cell
    # but the L^2 error should be of comparable magnitude to the triangle
    # case at the same h.
    quad_case = MMSPoissonCase(
        dim=2, N=64, amplitude=0.5, cell_kind="quadrilateral",
    )
    quad_result = run_one_level(quad_case, sc=sc)
    tri_64 = next(r for r in studies["2d_triangles"] if r["N"] == 64)
    smoke_row = {
        "N": quad_result.N,
        "h": quad_result.h,
        "N_dofs_quad": quad_result.n_dofs,
        "e_L2_quad": quad_result.e_L2,
        "e_H1_quad": quad_result.e_H1,
        "N_dofs_triangle": tri_64["N_dofs"],
        "e_L2_triangle": tri_64["e_L2"],
        "e_H1_triangle": tri_64["e_H1"],
        "ratio_L2_quad_over_triangle": (
            quad_result.e_L2 / tri_64["e_L2"] if tri_64["e_L2"] > 0 else float("nan")
        ),
    }
    from ._convergence import write_convergence_csv

    write_convergence_csv(
        out_dir / "2d_quad_smoke.csv",
        [smoke_row],
        list(smoke_row.keys()),
    )
    studies["2d_quad_smoke"] = [smoke_row]

    return studies
