"""
Top-level entry point: run a simulation from a validated config dict.

Usage:
    from semi import schema, run
    cfg = schema.load('benchmarks/pn_1d/pn_junction.json')
    result = run.run(cfg)

Day 1 scope: equilibrium Poisson (Boltzmann statistics, no applied bias).
Nonequilibrium drift-diffusion, bias sweeps, multi-region will layer in
subsequent iterations.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np


@dataclass
class SimulationResult:
    """Container for simulation outputs."""
    cfg: dict[str, Any]
    mesh: Any = None
    V: Any = None
    psi: Any = None                         # scaled potential Function
    psi_phys: np.ndarray | None = None      # volts, dof-ordered
    n_phys: np.ndarray | None = None        # m^-3, dof-ordered
    p_phys: np.ndarray | None = None        # m^-3, dof-ordered
    x_dof: np.ndarray | None = None         # (N, 3) dof coordinates, physical
    N_hat: Any = None                       # scaled doping Function
    scaling: Any = None
    solver_info: dict[str, Any] = field(default_factory=dict)


def run(cfg: dict[str, Any]) -> SimulationResult:
    """
    Run a simulation specified by `cfg` (already validated via schema.load).

    Currently implements equilibrium Poisson only.
    """
    from dolfinx import fem

    from .doping import build_profile
    from .mesh import build_mesh
    from .physics.poisson import build_equilibrium_poisson_form
    from .scaling import make_scaling_from_config
    from .solver import solve_nonlinear

    # Reference material: first semiconductor region
    ref_mat = _reference_material(cfg)
    sc = make_scaling_from_config(cfg, ref_mat)

    # Mesh and tags
    msh, cell_tags, facet_tags = build_mesh(cfg)
    V = fem.functionspace(msh, ("Lagrange", 1))

    # Doping: evaluate in physical coords, divide by C0 for scaled form
    N_raw_fn = build_profile(cfg["doping"])

    def N_hat_expr(x: np.ndarray) -> np.ndarray:
        return N_raw_fn(x) / sc.C0

    N_hat_fn = fem.Function(V, name="N_net_hat")
    N_hat_fn.interpolate(N_hat_expr)

    # Build BCs from ohmic contacts (equilibrium potential from charge neutrality)
    psi = fem.Function(V, name="psi_hat")
    bcs = _build_ohmic_bcs(cfg, V, msh, facet_tags, sc, ref_mat, N_raw_fn)

    # Initial guess: set psi pointwise to the local bulk equilibrium potential
    # psi_hat(x) = asinh(N_net(x) / (2 n_i)). This matches the ohmic BC values at
    # the contacts (under zero bias) and is already the correct answer deep in
    # each bulk region, so Newton only needs to smooth the junction. Starting
    # from psi = 0 everywhere would linearize exp(psi) around zero (losing the
    # steep nonlinearity) and the Newton iteration converges spuriously to a
    # linear-Poisson solution that underestimates the peak field.
    two_ni = 2.0 * ref_mat.n_i

    def psi_init_expr(x: np.ndarray) -> np.ndarray:
        return np.arcsinh(N_raw_fn(x) / two_ni)

    psi.interpolate(psi_init_expr)
    # Re-apply Dirichlet values so the boundary dofs match exactly (the
    # interpolation above uses per-point doping and so should agree to
    # roundoff, but BC values are derived from facet-centroid doping; keep
    # them consistent with what SNES will enforce).
    for bc in bcs:
        bc.set(psi.x.array)
    psi.x.scatter_forward()

    # Build residual and solve
    F = build_equilibrium_poisson_form(V, psi, N_hat_fn, sc, ref_mat.epsilon_r)
    info = solve_nonlinear(F, psi, bcs, prefix=f"{cfg['name']}_")

    # Post-process
    x_dof = V.tabulate_dof_coordinates()
    psi_phys = psi.x.array * sc.V0
    # Boltzmann-derived carriers (psi.x is in units of V_t)
    n_phys = ref_mat.n_i * np.exp(psi.x.array)
    p_phys = ref_mat.n_i * np.exp(-psi.x.array)

    return SimulationResult(
        cfg=cfg, mesh=msh, V=V, psi=psi,
        psi_phys=psi_phys, n_phys=n_phys, p_phys=p_phys,
        x_dof=x_dof, N_hat=N_hat_fn, scaling=sc, solver_info=info,
    )


def _reference_material(cfg: dict[str, Any]):
    """Pick the first semiconductor region's material as the reference."""
    from .materials import get_material
    for _region_name, region in cfg["regions"].items():
        mat = get_material(region["material"])
        if mat.is_semiconductor():
            return mat
    raise ValueError("No semiconductor region found; nothing to solve.")


def _build_ohmic_bcs(cfg, V, msh, facet_tags, sc, ref_mat, N_raw_fn):
    """
    Construct Dirichlet BCs at ohmic contacts.

    For equilibrium Poisson, the scaled BC value is
        psi_hat = asinh(N_net / (2 n_i)) + V_applied / V_t
    which corresponds to the contact Fermi level being pinned to the
    majority-carrier bulk quasi-Fermi level plus any external bias.
    """
    from dolfinx import fem
    from petsc4py import PETSc

    tdim = msh.topology.dim
    fdim = tdim - 1

    # Resolve facet tag references (names or ints)
    tag_by_name = {}
    for plane in cfg["mesh"].get("facets_by_plane", []):
        tag_by_name[plane["name"]] = int(plane["tag"])

    bcs = []
    for contact in cfg["contacts"]:
        if contact["type"] == "insulating":
            continue
        if contact["type"] != "ohmic":
            # gate / other types deferred
            continue

        facet_ref = contact["facet"]
        tag = tag_by_name[facet_ref] if isinstance(facet_ref, str) else int(facet_ref)

        if facet_tags is None:
            raise RuntimeError(
                "Mesh has no facet tags; can't apply contact BCs. "
                "Check that mesh.facets_by_plane is populated."
            )

        facets = facet_tags.find(tag)
        if len(facets) == 0:
            raise RuntimeError(
                f"No facets found with tag {tag} for contact {contact['name']!r}"
            )

        # Evaluate doping at a representative point on the facet to get
        # the equilibrium potential. Use the first vertex of the first
        # facet; in 1D a facet is a vertex so this is exact. In higher D,
        # uniform doping on the contact is a common assumption.
        N_net = _evaluate_doping_at_facet(msh, facets, fdim, N_raw_fn)

        psi_eq_hat = float(np.arcsinh(N_net / (2.0 * ref_mat.n_i)))
        V_applied = contact.get("voltage", 0.0)
        psi_bc = psi_eq_hat + V_applied / sc.V0

        dofs = fem.locate_dofs_topological(V, fdim, facets)
        bc = fem.dirichletbc(PETSc.ScalarType(psi_bc), dofs, V)
        bcs.append(bc)

    return bcs


def _evaluate_doping_at_facet(msh, facets, fdim, N_raw_fn) -> float:
    """
    Evaluate net doping at the centroid of the first given facet.

    Works in 1D (facet = vertex), 2D (facet = edge), 3D (facet = triangle/quad).
    Returns a single float (assumes doping is uniform across the contact).
    """
    tdim = msh.topology.dim
    msh.topology.create_connectivity(fdim, 0)
    f2v = msh.topology.connectivity(fdim, 0)
    coords = msh.geometry.x  # (num_nodes, 3)

    verts = f2v.links(int(facets[0]))
    # Centroid in physical space, only first `tdim` components matter
    centroid = coords[verts, :tdim].mean(axis=0)

    # Profile expects shape (dim, N); build single-point array
    pt = centroid.reshape(tdim, 1)
    return float(N_raw_fn(pt)[0])
