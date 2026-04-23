"""
Equilibrium Poisson runner (M1).

Solves -div(L_D^2 eps_r grad psi_hat) = p_hat - n_hat + N_hat under
Boltzmann statistics on the validated config and returns a populated
`SimulationResult`. No applied bias, no continuity equations.
"""
from __future__ import annotations

from typing import Any

import numpy as np


def run_equilibrium(cfg: dict[str, Any]):
    """M1 equilibrium Poisson solver (Boltzmann statistics)."""
    from dolfinx import fem

    from ..bcs import build_psi_dirichlet_bcs, resolve_contacts
    from ..doping import build_profile
    from ..mesh import build_mesh
    from ..physics.poisson import build_equilibrium_poisson_form
    from ..run import SimulationResult
    from ..scaling import make_scaling_from_config
    from ..solver import solve_nonlinear
    from ._common import reference_material

    ref_mat = reference_material(cfg)
    sc = make_scaling_from_config(cfg, ref_mat)

    msh, _cell_tags, facet_tags = build_mesh(cfg)
    V = fem.functionspace(msh, ("Lagrange", 1))

    N_raw_fn = build_profile(cfg["doping"])

    def N_hat_expr(x):
        return N_raw_fn(x) / sc.C0

    N_hat_fn = fem.Function(V, name="N_net_hat")
    N_hat_fn.interpolate(N_hat_expr)

    psi = fem.Function(V, name="psi_hat")
    contacts = resolve_contacts(cfg, facet_tags=facet_tags)
    bcs = build_psi_dirichlet_bcs(V, msh, facet_tags, contacts, sc, ref_mat, N_raw_fn)

    two_ni = 2.0 * ref_mat.n_i

    def psi_init_expr(x):
        return np.arcsinh(N_raw_fn(x) / two_ni)

    psi.interpolate(psi_init_expr)
    for bc in bcs:
        bc.set(psi.x.array)
    psi.x.scatter_forward()

    F = build_equilibrium_poisson_form(V, psi, N_hat_fn, sc, ref_mat.epsilon_r)
    info = solve_nonlinear(F, psi, bcs, prefix=f"{cfg['name']}_")

    x_dof = V.tabulate_dof_coordinates()
    psi_phys = psi.x.array * sc.V0
    n_phys = ref_mat.n_i * np.exp(psi.x.array)
    p_phys = ref_mat.n_i * np.exp(-psi.x.array)

    return SimulationResult(
        cfg=cfg, mesh=msh, V=V, psi=psi,
        psi_phys=psi_phys, n_phys=n_phys, p_phys=p_phys,
        x_dof=x_dof, N_hat=N_hat_fn, scaling=sc, solver_info=info,
    )
