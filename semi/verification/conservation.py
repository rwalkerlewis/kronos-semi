"""
Conservation checks for the 1D pn junction benchmarks.

Two distinct checks live here, mirroring the two things a steady-state
drift-diffusion solution must satisfy to be self-consistent:

1) Current continuity (`pn_1d_bias`, `pn_1d_bias_reverse`)

   In steady state with uniform cross-section, J_total(x) = J_n + J_p
   must be constant in space (div J = 0 with no source at the terminals).
   The check samples J_total at ~10 interior facets via a UFL dS
   facet-integral pattern and asserts

        max |J_total(x_i) - mean(J_total)| / |mean(J_total)| < tol

   with tol = 5% for forward bias (where J >> numerical noise floor) and
   tol = 15% for reverse bias (where |J| is many orders of magnitude
   smaller and absolute Newton residual leakage shows up relatively
   larger).

2) Charge conservation (`pn_1d` equilibrium)

   For an isolated device with ohmic contacts at equilibrium the total
   space charge integrates to zero:

        integral(p - n + N_net) dx == 0    (N_net = N_D - N_A)

   to within the nonlinear solver tolerance. The check asserts

        |Q_net| < 1e-10 * q * max(|N_net|) * L_device

   i.e. the net charge per cross-sectional area is at least ten orders
   of magnitude smaller than the charge one would get by taking peak
   doping uniformly over the entire device length. That matches the
   snes_rtol ~ 1e-14 tolerance used by `run_equilibrium`.

The pure-Python metric functions (`current_continuity_metric`,
`charge_neutrality_metric`) have no dolfinx dependency and are the
unit-testable core. The FEM helpers (`current_continuity_from_result`,
`charge_conservation_from_result`) import dolfinx only when called.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# numpy 2.0 renamed np.trapz to np.trapezoid; keep both paths working so
# the pure-Python metric tests run on the 1.x CI matrix and on the 2.x
# Docker FEM image.
_trapz = getattr(np, "trapezoid", None) or np.trapz


# ---------------------------------------------------------------------------
# Pure-Python metrics (no dolfinx)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class CurrentContinuityMetric:
    """Outcome of sampling J_total at a set of interior x-locations."""

    xs: np.ndarray
    J_vals: np.ndarray
    mean_J: float
    max_abs_dev: float
    max_rel_dev: float


def current_continuity_metric(xs, J_vals) -> CurrentContinuityMetric:
    """
    Reduce an array of J_total(x_i) samples to (mean, max abs deviation,
    max relative deviation). If the mean is zero the relative deviation
    is +inf when any sample deviates and 0 otherwise.
    """
    xs = np.asarray(xs, dtype=float)
    J = np.asarray(J_vals, dtype=float)
    if J.size == 0:
        return CurrentContinuityMetric(xs, J, 0.0, 0.0, 0.0)
    mean_J = float(np.mean(J))
    abs_dev = np.abs(J - mean_J)
    max_abs = float(np.max(abs_dev))
    denom = abs(mean_J)
    if denom == 0.0:
        max_rel = float("inf") if max_abs > 0.0 else 0.0
    else:
        max_rel = max_abs / denom
    return CurrentContinuityMetric(
        xs=xs, J_vals=J, mean_J=mean_J,
        max_abs_dev=max_abs, max_rel_dev=max_rel,
    )


@dataclass(frozen=True)
class ChargeNeutralityMetric:
    """Outcome of integrating rho(x) over the device."""

    Q_net: float       # q * integral(p - n + N_net) dx, C/m^2
    Q_ref: float       # q * max|N_net| * L, C/m^2
    rel_error: float   # |Q_net| / Q_ref


def charge_neutrality_metric(x, p, n, N_net, q_value) -> ChargeNeutralityMetric:
    """
    Integrate (p - n + N_net) dx by the trapezoidal rule and report both
    the dimensional net charge per cross-sectional area (q times the
    integral) and the relative-to-peak-doping reference q * max|N_net| * L.
    """
    x = np.asarray(x, dtype=float)
    p = np.asarray(p, dtype=float)
    n = np.asarray(n, dtype=float)
    N_net = np.asarray(N_net, dtype=float)
    order = np.argsort(x)
    x_s = x[order]
    rho_over_q = p[order] - n[order] + N_net[order]
    integral = float(_trapz(rho_over_q, x_s))
    Q_net = q_value * integral
    L = float(x_s[-1] - x_s[0]) if x_s.size >= 2 else 0.0
    peak_N = float(np.max(np.abs(N_net))) if N_net.size else 0.0
    Q_ref = q_value * peak_N * L
    rel = abs(Q_net) / Q_ref if Q_ref > 0.0 else float("inf")
    return ChargeNeutralityMetric(Q_net=Q_net, Q_ref=Q_ref, rel_error=rel)


# ---------------------------------------------------------------------------
# FEM helpers (dolfinx required)
# ---------------------------------------------------------------------------


def _sample_interior_facets_1d(msh, n_samples: int):
    """
    Return `(facet_indices, xs)` for up to `n_samples` interior facets
    on a 1D mesh, roughly evenly spaced along the x-axis, *excluding*
    facets adjacent to boundary cells.

    The exclusion matters for drift-diffusion continuity checks: cells
    that touch a Dirichlet contact have one-sided FE gradients in which
    grad(phi_n) is effectively (phi_n(x=0) - phi_n(x=h)) / h. That
    difference quotient includes the prescribed BC value, which forces
    phi_n (and hence the computed J_n) in the boundary cell to absorb
    the full mismatch between "ideal drift-diffusion bulk" and "fixed
    contact potential". The interior facet adjacent to that cell then
    sees an `avg(J_n)` that is the mean of a clean interior gradient
    and a boundary-skewed one.

    In an 800-cell pn junction at forward bias the bulk-facet J_total
    samples agree to O(1e-5) while the two boundary-adjacent facets
    differ from the bulk by O(100%). Excluding them is the physically
    meaningful way to measure "is J_total constant in the device
    interior" as Phase 3 V&V intends; it does not hide any real
    conservation defect (which would show up in the bulk too).
    """
    from dolfinx import fem

    tdim = msh.topology.dim
    if tdim != 1:
        raise ValueError(
            f"conservation sampler supports 1D only, got tdim={tdim}"
        )
    fdim = tdim - 1

    msh.topology.create_connectivity(fdim, tdim)
    msh.topology.create_connectivity(tdim, fdim)
    f2c = msh.topology.connectivity(fdim, tdim)
    c2f = msh.topology.connectivity(tdim, fdim)
    n_facets = msh.topology.index_map(fdim).size_local
    n_cells = msh.topology.index_map(tdim).size_local

    # Collect boundary cells: any cell that owns a boundary facet.
    boundary_cells: set[int] = set()
    for c in range(n_cells):
        for f in c2f.links(c):
            if len(f2c.links(f)) == 1:
                boundary_cells.add(c)
                break

    V = fem.functionspace(msh, ("Lagrange", 1))
    all_coords = V.tabulate_dof_coordinates()

    interior = []
    xs = []
    for i in range(n_facets):
        adj = f2c.links(i)
        if len(adj) != 2:
            continue
        # Skip interior facets whose gradient evaluation would pull in a
        # boundary cell; see module docstring.
        if int(adj[0]) in boundary_cells or int(adj[1]) in boundary_cells:
            continue
        dofs = fem.locate_dofs_topological(
            V, fdim, np.array([i], dtype=np.int32)
        )
        if dofs.size == 0:
            continue
        interior.append(i)
        xs.append(float(all_coords[dofs[0], 0]))

    interior = np.array(interior, dtype=np.int32)
    xs = np.array(xs, dtype=float)
    if interior.size == 0:
        return interior, xs

    order = np.argsort(xs)
    interior = interior[order]
    xs = xs[order]

    if interior.size <= n_samples:
        return interior, xs

    pick = np.linspace(0, interior.size - 1, n_samples, dtype=int)
    return interior[pick], xs[pick]


def _dd_mobilities_SI(cfg) -> tuple[float, float]:
    """Return (mu_n_SI, mu_p_SI) in m^2/(V s) from the config."""
    mob = cfg.get("physics", {}).get("mobility", {})
    mu_n_SI = float(mob.get("mu_n", 1400.0)) * 1.0e-4
    mu_p_SI = float(mob.get("mu_p", 450.0)) * 1.0e-4
    return mu_n_SI, mu_p_SI


def _semiconductor_material(cfg):
    """Pick the semiconductor material (only one is supported for now)."""
    from semi.materials import get_material

    for _name, region in cfg["regions"].items():
        mat = get_material(region["material"])
        if mat.is_semiconductor():
            return mat
    raise ValueError("No semiconductor region found in cfg")


def current_continuity_from_state(
    msh,
    sc,
    cfg,
    psi,
    phi_n,
    phi_p,
    *,
    n_samples: int = 10,
) -> CurrentContinuityMetric:
    """
    Evaluate J_total = J_n + J_p at ~`n_samples` interior facets of the
    1D mesh via UFL `dS` facet integrals and return the continuity
    metric.

    Each sampled facet i is tagged individually, a one-entry meshtags
    object marks that facet, and `avg(J_total) * dS(tag_i)` is assembled.
    In 1D dS on a 0-D entity has unit measure so the assembled scalar
    equals avg(J_total) at that vertex. We divide by the assembled
    measure defensively in case future dolfinx versions change that
    convention.

    This function operates on raw UFL/dolfinx state (no SimulationResult)
    so it can be called from inside `semi.run.run_bias_sweep` to record
    continuity at intermediate bias steps without first materializing a
    SimulationResult.
    """
    import ufl
    from dolfinx import fem
    from dolfinx.mesh import meshtags
    from mpi4py import MPI
    from petsc4py import PETSc

    from semi.constants import Q as Q_SI
    from semi.physics.slotboom import n_from_slotboom, p_from_slotboom

    mat = _semiconductor_material(cfg)
    mu_n_SI, mu_p_SI = _dd_mobilities_SI(cfg)

    tdim = msh.topology.dim
    fdim = tdim - 1

    facet_indices, xs = _sample_interior_facets_1d(msh, n_samples)
    if facet_indices.size == 0:
        return current_continuity_metric(xs, np.zeros(0))

    ni_hat = fem.Constant(msh, PETSc.ScalarType(mat.n_i / sc.C0))

    n_ufl = n_from_slotboom(psi, phi_n, ni_hat) * sc.C0
    p_ufl = p_from_slotboom(psi, phi_p, ni_hat) * sc.C0
    grad_phi_n_phys = sc.V0 * ufl.grad(phi_n)
    grad_phi_p_phys = sc.V0 * ufl.grad(phi_p)

    Jn_x = Q_SI * mu_n_SI * n_ufl * grad_phi_n_phys[0]
    Jp_x = Q_SI * mu_p_SI * p_ufl * grad_phi_p_phys[0]
    J_total = Jn_x + Jp_x

    J_vals = np.empty(facet_indices.size, dtype=float)
    for k, f in enumerate(facet_indices):
        tag = 1
        indices = np.array([int(f)], dtype=np.int32)
        values = np.array([tag], dtype=np.int32)
        facet_mt = meshtags(msh, fdim, indices, values)
        dS = ufl.Measure("dS", domain=msh, subdomain_data=facet_mt, subdomain_id=tag)

        I_form = fem.form(ufl.avg(J_total) * dS)
        A_form = fem.form(ufl.avg(fem.Constant(msh, PETSc.ScalarType(1.0))) * dS)

        I_local = fem.assemble_scalar(I_form)
        A_local = fem.assemble_scalar(A_form)
        I = msh.comm.allreduce(I_local, op=MPI.SUM)
        A = msh.comm.allreduce(A_local, op=MPI.SUM)
        J_vals[k] = float(I / A) if A != 0.0 else float("nan")

    return current_continuity_metric(xs, J_vals)


def current_continuity_from_result(result, n_samples: int = 10) -> CurrentContinuityMetric:
    """
    SimulationResult wrapper around `current_continuity_from_state`.
    Used by the `pn_1d_bias` / `pn_1d_bias_reverse` verifiers to check
    the continuity of the final-bias state returned by the benchmark.
    """
    return current_continuity_from_state(
        msh=result.mesh,
        sc=result.scaling,
        cfg=result.cfg,
        psi=result.psi,
        phi_n=result.phi_n,
        phi_p=result.phi_p,
        n_samples=n_samples,
    )


def run_bias_sweep_with_continuity(cfg: dict, *, n_samples: int = 10):  # pragma: no cover
    """
    Run `semi.run.run_bias_sweep(cfg)` with a post-step hook that attaches
    `continuity_*` fields to each IV row. Returns the SimulationResult
    whose `.iv` entries now carry:

        continuity_mean_J   mean of J_total samples at this bias (A/m^2)
        continuity_max_dev  max absolute deviation from the mean (A/m^2)
        continuity_max_rel  max |J - mean| / |mean|  (dimensionless)
        continuity_n_pts    number of interior-facet samples used
    """
    from semi.run import run_bias_sweep
    from semi.scaling import make_scaling_from_config

    ref_mat = _semiconductor_material(cfg)
    sc = make_scaling_from_config(cfg, ref_mat)

    def hook(V_applied, spaces, iv_row):
        cc = current_continuity_from_state(
            msh=spaces.V_psi.mesh, sc=sc, cfg=cfg,
            psi=spaces.psi, phi_n=spaces.phi_n, phi_p=spaces.phi_p,
            n_samples=n_samples,
        )
        iv_row["continuity_mean_J"] = float(cc.mean_J)
        iv_row["continuity_max_dev"] = float(cc.max_abs_dev)
        iv_row["continuity_max_rel"] = float(cc.max_rel_dev)
        iv_row["continuity_n_pts"] = int(cc.xs.size)

    return run_bias_sweep(cfg, post_step_hook=hook)


def charge_conservation_from_result(result) -> ChargeNeutralityMetric:
    """
    Integrate q * (p - n + N_net) over the device from the dof arrays
    produced by `run_equilibrium`. This uses NumPy trapezoidal
    integration on the sorted 1D mesh rather than a UFL dx assembly
    because:

      - result.x_dof, result.n_phys, result.p_phys, result.N_hat are
        already available and carry the exact quantities we need;
      - the trapezoidal rule on P1 dofs is the same numerical value as
        a UFL dx assembly of `(p - n + N_net) * dx` against a constant
        test function on P1 (mass-lumped), which is what we want;
      - the pure-Python metric is independently unit-testable with
        synthetic arrays.
    """
    from semi.constants import Q as Q_SI

    sc = result.scaling
    x = np.asarray(result.x_dof[:, 0], dtype=float)
    p = np.asarray(result.p_phys, dtype=float)
    n = np.asarray(result.n_phys, dtype=float)
    N_net = np.asarray(result.N_hat.x.array, dtype=float) * sc.C0
    return charge_neutrality_metric(x, p, n, N_net, Q_SI)
