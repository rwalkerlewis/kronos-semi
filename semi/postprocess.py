"""
Post-processing helpers for the bias-sweep runner.

This module owns the FEM-side "what do I do with a converged solution
at one bias step" code path: facet metadata lookup, terminal-current
evaluation through a UFL weak form, and IV-row recording. The runner
in `semi.runners.bias_sweep` calls into here once per ramp step.
"""
from __future__ import annotations

import numpy as np


def resolve_contact_facets(cfg, msh, facet_tags, contact_name):
    """
    Look up the facet entities for the named contact and decide an
    outward sign based on whether the plane sits at the high or low
    end of its mesh extent. Returns a metadata dict consumed by
    `evaluate_current_at_contact` and IV recording.
    """
    from dolfinx import fem

    mesh_cfg = cfg.get("mesh") or {}
    tag_by_name = {p["name"]: int(p["tag"])
                   for p in mesh_cfg.get("facets_by_plane", [])}
    contact = next(c for c in cfg["contacts"] if c["name"] == contact_name)
    facet_ref = contact["facet"]
    tag = tag_by_name[facet_ref] if isinstance(facet_ref, str) else int(facet_ref)

    fdim = msh.topology.dim - 1
    facets = facet_tags.find(tag)

    outward_sign = 1
    extents = mesh_cfg.get("extents", [])
    for p in mesh_cfg.get("facets_by_plane", []):
        if int(p["tag"]) == tag:
            axis = int(p["axis"])
            if axis < len(extents):
                xmin, xmax = extents[axis]
                outward_sign = 1 if float(p["value"]) >= 0.5 * (xmin + xmax) else -1
            break

    V_ref = fem.functionspace(msh, ("Lagrange", 1))
    dofs = fem.locate_dofs_topological(V_ref, fdim, facets)
    return {
        "dofs": np.asarray(dofs, dtype=np.int64),
        "facets": facets,
        "outward_sign": int(outward_sign),
        "tag": int(tag),
    }


def evaluate_current_at_contact(spaces, sc, ref_mat, facet_info,
                                mu_n_SI, mu_p_SI) -> float:
    """
    Evaluate J = J_n + J_p at a contact via a UFL weak form:

        J_n = q mu_n n (grad phi_n . n_outward)
        J_p = q mu_p p (grad phi_p . n_outward)

    Integrated over the contact facet and divided by its measure.
    Positive J = current flowing outward (out of the device through
    the contact).
    """
    import ufl
    from dolfinx import fem
    from dolfinx.mesh import meshtags
    from mpi4py import MPI
    from petsc4py import PETSc

    from .constants import Q
    from .physics.slotboom import n_from_slotboom, p_from_slotboom

    msh = spaces.V_psi.mesh
    tag = int(facet_info["tag"])
    facets = facet_info["facets"]
    fdim = msh.topology.dim - 1

    if len(facets) == 0:
        return 0.0

    values = np.full(len(facets), tag, dtype=np.int32)
    indices = np.asarray(facets, dtype=np.int32)
    sort = np.argsort(indices)
    facet_mt = meshtags(msh, fdim, indices[sort], values[sort])
    ds = ufl.Measure("ds", domain=msh, subdomain_data=facet_mt, subdomain_id=tag)
    n_vec = ufl.FacetNormal(msh)

    ni_hat = fem.Constant(msh, PETSc.ScalarType(ref_mat.n_i / sc.C0))
    psi = spaces.psi
    phi_n = spaces.phi_n
    phi_p = spaces.phi_p

    n_ufl = n_from_slotboom(psi, phi_n, ni_hat) * sc.C0
    p_ufl = p_from_slotboom(psi, phi_p, ni_hat) * sc.C0
    grad_phi_n_phys = sc.V0 * ufl.grad(phi_n)
    grad_phi_p_phys = sc.V0 * ufl.grad(phi_p)

    Jn = Q * mu_n_SI * n_ufl * ufl.dot(grad_phi_n_phys, n_vec)
    Jp = Q * mu_p_SI * p_ufl * ufl.dot(grad_phi_p_phys, n_vec)
    current_form = fem.form((Jn + Jp) * ds)
    area_form = fem.form(1.0 * ds)

    I_local = fem.assemble_scalar(current_form)
    A_local = fem.assemble_scalar(area_form)
    I = msh.comm.allreduce(I_local, op=MPI.SUM)
    A = msh.comm.allreduce(A_local, op=MPI.SUM)

    if msh.topology.dim == 1:
        if A == 0.0:
            A = 1.0
        J_total = I / A
    else:
        if A == 0.0:
            return 0.0
        J_total = I / A

    return float(J_total)


def record_iv(iv_rows, V_applied, spaces, sc, ref_mat,
              sweep_contact, sweep_facet_info, mu_n_SI, mu_p_SI):
    """Append one (V, J) row to `iv_rows`, evaluating J at the sweep contact."""
    if sweep_contact is None or sweep_facet_info is None:
        iv_rows.append({"V": float(V_applied), "J": 0.0})
        return
    J = evaluate_current_at_contact(
        spaces, sc, ref_mat, sweep_facet_info, mu_n_SI, mu_p_SI,
    )
    iv_rows.append({"V": float(V_applied), "J": float(J)})


def fmt_tag(v: float) -> str:
    """
    Format a voltage as a PETSc options-prefix-safe tag.

        +0.250 V -> "p0d2500"
        -0.125 V -> "m0d1250"

    Used by `solve_at` to namespace the per-step SNES options prefix.
    """
    return f"{v:+.4f}".replace("+", "p").replace("-", "m").replace(".", "d")
