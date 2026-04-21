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


# ---------------------------------------------------------------------------
# Multi-region Poisson MMS (Day 6)
#
# The single-region studies above fix eps_r to a scalar Constant. This
# variant exercises the coefficient-jump assembly with a piecewise DG0
# eps_r(x) matching the Day-6 MOS stack: eps_r_Si = 11.7 over the bottom
# fraction of the square, eps_r_ox = 3.9 over the top. The exact
# solution
#
#     psi_exact(x, y) = sin(pi x / W) * h(y)
#
# uses a piecewise-polynomial h(y) chosen so that (i) psi is C^0 across
# the Si/SiO2 interface and (ii) eps_r * d psi/dy is continuous across
# it (the natural interface condition for the weak form). Per-region
# forcings f_Si, f_ox are manufactured directly from psi_exact and the
# local eps_r_k.
#
# The scaling factor L_D^2 from the production Poisson LHS cancels out
# when paired with a consistent RHS, so this MMS runs in pure scalar
# Poisson form without the Slotboom / doping / carrier terms. That
# isolates the coefficient-jump assembly, which is the only new thing
# in Day-6 Poisson physics. See docs/mos_derivation.md section 7.
# ---------------------------------------------------------------------------

# Material constants used by the multi-region MMS. Hard-coded so the
# test is self-contained.
MR_EPS_R_SI: float = 11.7
MR_EPS_R_OX: float = 3.9

# Geometry per docs/mos_derivation.md section 7.1.
MR_W: float = 1.0e-7
MR_T_SI: float = 0.7e-7
MR_T_OX: float = 0.3e-7
MR_Y_INT: float = MR_T_SI
MR_Y_TOP: float = MR_T_SI + MR_T_OX

# psi_exact amplitudes.
MR_AMPL: float = 1.0
MR_H_TOP: float = 0.5


def _mr_B_ox() -> float:
    return 2.0 * MR_EPS_R_SI / (MR_EPS_R_OX * MR_T_SI)


def _mr_gamma_ox() -> float:
    return (MR_H_TOP - MR_AMPL * (1.0 + _mr_B_ox() * MR_T_OX)) / (MR_T_OX ** 2)


def _mr_psi_exact_ufl(mesh):
    """Piecewise-polynomial exact solution psi_exact(x, y) on the full square.

    Because UFL does not provide a clean per-subdomain branch for the
    exact solution (it lives in a single expression on the mesh), we
    use a `ufl.conditional` on the spatial y coordinate to select the
    correct polynomial in each region. At the interface both branches
    give the same value (by construction), so the conditional is
    unambiguous.
    """
    import ufl

    x = ufl.SpatialCoordinate(mesh)
    sx = ufl.sin(ufl.pi * x[0] / MR_W)

    A = MR_AMPL
    t_Si = MR_T_SI
    y_int = MR_Y_INT
    B_ox = _mr_B_ox()
    gamma_ox = _mr_gamma_ox()

    h_Si = A * (x[1] / t_Si) ** 2
    dy = x[1] - y_int
    h_ox = A + A * B_ox * dy + gamma_ox * dy ** 2

    # conditional(y <= y_int, h_Si, h_ox)
    h = ufl.conditional(ufl.le(x[1], y_int), h_Si, h_ox)
    return sx * h


def _mr_build_mesh(N_x: int, N_y: int):
    """Rectangle on [0, W] x [0, y_top] with y_int as a grid line.

    Requires N_y such that y_int = (N_y * T_SI / Y_TOP) is an integer,
    else the interface falls inside a cell and the coefficient-jump
    assembly becomes ill-defined in a DG0 eps_r field.
    """
    import numpy as _np
    from dolfinx import mesh as dmesh
    from mpi4py import MPI

    expected_layer = N_y * MR_T_SI / MR_Y_TOP
    if abs(expected_layer - round(expected_layer)) > 1.0e-9:
        raise ValueError(
            f"N_y = {N_y} does not land y_int on a grid line "
            f"(y_int / h_y = {expected_layer!r}). Use N_y divisible by 10."
        )
    return dmesh.create_rectangle(
        MPI.COMM_WORLD,
        [_np.array([0.0, 0.0]), _np.array([MR_W, MR_Y_TOP])],
        [N_x, N_y],
        cell_type=dmesh.CellType.triangle,
        diagonal=dmesh.DiagonalType.right,
    )


def _mr_build_region_tags(msh):
    """Cell tags: 1 = silicon (y_mid <= y_int), 2 = oxide (y_mid > y_int)."""
    import numpy as _np
    from dolfinx import mesh as dmesh

    tdim = msh.topology.dim
    msh.topology.create_entities(tdim)
    n = msh.topology.index_map(tdim).size_local + msh.topology.index_map(tdim).num_ghosts
    all_cells = _np.arange(n, dtype=_np.int32)
    msh.topology.create_connectivity(tdim, 0)
    c2v = msh.topology.connectivity(tdim, 0)
    coords = msh.geometry.x
    values = _np.empty(n, dtype=_np.int32)
    for c in range(n):
        verts = c2v.links(int(c))
        y_mid = coords[verts, 1].mean()
        values[c] = 1 if y_mid <= MR_Y_INT else 2
    return dmesh.meshtags(msh, tdim, all_cells, values)


def _mr_build_eps_r_function(msh, cell_tags):
    """DG0 eps_r with eps_r_Si on tag 1 and eps_r_ox on tag 2."""
    from dolfinx import fem
    from petsc4py import PETSc

    V_DG0 = fem.functionspace(msh, ("DG", 0))
    eps_r_fn = fem.Function(V_DG0, name="eps_r")
    eps_r_fn.x.array[:] = PETSc.ScalarType(1.0)
    for tag, eps in ((1, MR_EPS_R_SI), (2, MR_EPS_R_OX)):
        for c in cell_tags.find(tag):
            dof = V_DG0.dofmap.cell_dofs(int(c))[0]
            eps_r_fn.x.array[dof] = PETSc.ScalarType(eps)
    eps_r_fn.x.scatter_forward()
    return eps_r_fn


def _mr_all_boundary_facets(mesh, fdim: int):
    import numpy as _np
    from dolfinx import mesh as dmesh

    return dmesh.locate_entities_boundary(
        mesh, fdim, lambda x: _np.full(x.shape[1], True)
    )


def run_mms_poisson_2d_multiregion(N_x: int, N_y: int) -> "MMSPoissonResult":
    """Solve the multi-region Poisson MMS problem on one (N_x, N_y) level.

    Piecewise eps_r(x, y) with a Si/SiO2 coefficient jump at y = y_int.
    Dirichlet BCs are taken directly from psi_exact (non-homogeneous
    at the top edge y = y_top). The exact solution is C^0 across the
    interface and has continuous eps-weighted normal flux there, so
    the weak form matches psi_exact up to FE discretization error and
    the finest-pair L^2 rate targets 2.
    """
    import time as _time

    import numpy as _np
    import ufl
    from dolfinx import fem
    from petsc4py import PETSc

    from semi.solver import solve_nonlinear

    from ._norms import h1_seminorm_error_squared, l2_error_squared

    msh = _mr_build_mesh(N_x, N_y)
    cell_tags = _mr_build_region_tags(msh)
    eps_r_fn = _mr_build_eps_r_function(msh, cell_tags)

    V = fem.functionspace(msh, ("Lagrange", 1))
    v = ufl.TestFunction(V)
    psi = fem.Function(V, name="psi_hat_mr_mms")
    psi.x.array[:] = 0.0

    psi_e = _mr_psi_exact_ufl(msh)

    # Dirichlet on entire boundary, interpolated from psi_exact via a
    # Function so dolfinx can apply it to the DOFs that coincide with
    # grid nodes. psi_exact is smooth in each region and continuous at
    # the interface, so interpolation is exact for P1 away from y_int
    # and accurate to O(h^2) at the interface row (which is adequate
    # for an L^2 rate of 2).
    fdim = msh.topology.dim - 1
    msh.topology.create_connectivity(fdim, msh.topology.dim)
    boundary_facets = _mr_all_boundary_facets(msh, fdim)
    bdofs = fem.locate_dofs_topological(V, fdim, boundary_facets)

    psi_bc_fn = fem.Function(V, name="psi_bc_mr_mms")
    psi_bc_fn.interpolate(fem.Expression(psi_e, V.element.interpolation_points))
    bcs = [fem.dirichletbc(psi_bc_fn, bdofs)]

    # Weak form: integral eps_r grad(psi) . grad(v) dx - integral f v dx = 0
    # where f_k = -div(eps_r_k grad(psi_exact)). Using integration by
    # parts on the manufactured source term (with psi_exact taking its
    # boundary values through the Dirichlet BC lifting), this reduces
    # to
    #     F(psi) = integral eps_r [grad(psi) - grad(psi_exact)] . grad(v) dx
    # so at convergence grad(psi) == grad(psi_exact) weakly. No L_D^2
    # prefactor: the scaling factor cancels out of both terms (section
    # 7.3 of docs/mos_derivation.md).
    F = eps_r_fn * ufl.inner(ufl.grad(psi) - ufl.grad(psi_e), ufl.grad(v)) * ufl.dx

    petsc_options = {
        "snes_rtol": 1.0e-14,
        "snes_atol": 1.0e-16,
        "snes_stol": 0.0,
        "snes_max_it": 50,
    }
    t0 = _time.perf_counter()
    info = solve_nonlinear(
        F, psi, bcs,
        prefix=f"mms_poisson_mr_Nx{N_x}_Ny{N_y}_",
        petsc_options=petsc_options,
    )
    dt = _time.perf_counter() - t0
    if not info.get("converged", False):
        raise RuntimeError(
            f"MMS multi-region Poisson did not converge at "
            f"N=({N_x}, {N_y}): reason={info.get('reason')}"
        )

    e_L2 = float(_np.sqrt(max(l2_error_squared(psi, psi_e), 0.0)))
    e_H1 = float(_np.sqrt(max(h1_seminorm_error_squared(psi, psi_e), 0.0)))

    n_dofs = int(V.dofmap.index_map.size_global)
    # Report the y-direction spacing as the characteristic h (the
    # coefficient jump is perpendicular to y, so this is the
    # refinement direction that controls the interface error).
    h = MR_Y_TOP / N_y
    return MMSPoissonResult(
        dim=2,
        N=N_y,
        h=h,
        n_dofs=n_dofs,
        e_L2=e_L2,
        e_H1=e_H1,
        snes_iters=int(info.get("iterations", -1)),
        solve_time_s=dt,
        cell_kind="triangle",
        amplitude=MR_AMPL,
    )


def run_mr_convergence_study(Ns: list[int]) -> list["MMSPoissonResult"]:
    """Run the multi-region MMS on a sequence of refinements.

    Uses a square mesh with N_x = N_y = N and requires N divisible by
    10 so that the Si/SiO2 interface lands on a grid line.
    """
    results: list[MMSPoissonResult] = []
    for N in Ns:
        results.append(run_mms_poisson_2d_multiregion(N, N))
    return results


CLI_NS_MR = [20, 40, 80, 160]


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

    # 2D multi-region (Si/SiO2 coefficient jump) -- Day 6.
    mr_results = run_mr_convergence_study(CLI_NS_MR)
    mr_rows = to_table_rows(mr_results)
    write_artifacts(
        mr_rows, out_dir / "2d_multiregion",
        title="MMS Poisson 2D multi-region (Si/SiO2 coefficient jump)",
    )
    studies["2d_multiregion"] = mr_rows

    return studies
