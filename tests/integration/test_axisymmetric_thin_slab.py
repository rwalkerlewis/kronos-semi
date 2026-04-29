"""
Cross-validation of axisymmetric vs. cartesian assembly in the limit
``R >> L`` where curvature is negligible (ADR 0015).

The test runs the equilibrium Poisson runner twice on a thin pn slab:

* a planar Cartesian rectangle ``x in [0, L], y in [0, H]``;
* an axisymmetric "thin slab" at large radius ``r in [R, R+L], z in [0, H]``.

In the limit ``R / L -> infinity`` the cylindrical surface element
``r dr dz`` reduces to ``R * dr dz``, a constant scalar factor that
cancels in the resulting potential field. Both runners must therefore
produce the same psi_phys profile to within ~1% relative error at
``R = 100 * L``.

The test catches dispatch leaks: any runner that uses raw ``ufl.dx``
in one assembly site while the rest of the form is r-weighted will
return a different psi than the cartesian comparator at large R, and
this test will fail.
"""
from __future__ import annotations

import numpy as np


def _slab_cfg(name: str, x_min: float, x_max: float, coordinates: str) -> dict:
    """Build a 2D pn-slab equilibrium config.

    The transverse y direction is short (5 cells) so the solution is
    essentially 1D along x. Doping is a step at x = (x_min + x_max)/2:
    n-type to the left, p-type to the right. Two ohmic contacts pin
    psi at the left/right ends; top/bottom are natural-Neumann.
    """
    x_mid = 0.5 * (x_min + x_max)
    H = 1.0e-7  # thin transverse extent, m
    return {
        "schema_version": "1.3.0",
        "name": name,
        "dimension": 2,
        "mesh": {
            "source": "builtin",
            "extents": [[x_min, x_max], [0.0, H]],
            "resolution": [40, 5],
            "regions_by_box": [
                {"name": "silicon", "tag": 1,
                 "bounds": [[x_min, x_max], [0.0, H]]},
            ],
            "facets_by_plane": [
                {"name": "left",  "tag": 1, "axis": 0, "value": x_min},
                {"name": "right", "tag": 2, "axis": 0, "value": x_max},
            ],
        },
        "regions": {
            "silicon": {"material": "Si", "tag": 1, "role": "semiconductor"},
        },
        "doping": [
            {
                "region": "silicon",
                "profile": {
                    "type": "step",
                    "axis": 0,
                    "location": x_mid,
                    "N_D_left": 1.0e16,
                    "N_A_left": 0.0,
                    "N_D_right": 0.0,
                    "N_A_right": 1.0e16,
                },
            },
        ],
        "contacts": [
            {"name": "anode",   "type": "ohmic", "facet": "left",  "voltage": 0.0},
            {"name": "cathode", "type": "ohmic", "facet": "right", "voltage": 0.0},
        ],
        "physics": {
            "temperature": 300.0,
            "statistics": "boltzmann",
        },
        "solver": {
            "type": "equilibrium",
            "coordinates": coordinates,
        },
        "output": {},
    }


def _solve_midline(cfg: dict) -> tuple[np.ndarray, np.ndarray]:
    """Run the equilibrium runner and return ``(x_sorted, psi_phys_sorted)``
    along the transverse midline ``y = H/2``.
    """
    from semi.runners.equilibrium import run_equilibrium

    res = run_equilibrium(cfg)
    coords = res.x_dof
    psi = res.psi_phys
    H = 1.0e-7
    mid_mask = np.isclose(coords[:, 1], H / 2.0, atol=H * 0.05)
    if not mid_mask.any():
        ys = np.unique(coords[:, 1])
        ymid = ys[len(ys) // 2]
        mid_mask = np.isclose(coords[:, 1], ymid)
    xs = coords[mid_mask, 0]
    ps = psi[mid_mask]
    order = np.argsort(xs)
    return xs[order], ps[order]


def test_axisymmetric_thin_slab_matches_cartesian():
    """At ``R = 100 * L`` the axisymmetric solution agrees with cartesian
    to within 1% relative L_inf error along the transverse midline."""
    L = 2.0e-6  # device length, m
    R = 100.0 * L

    cfg_cart = _slab_cfg("slab_cart", 0.0, L, "cartesian")
    cfg_axi = _slab_cfg("slab_axi", R, R + L, "axisymmetric")

    x_cart, psi_cart = _solve_midline(cfg_cart)
    x_axi, psi_axi = _solve_midline(cfg_axi)

    assert x_cart.shape == x_axi.shape, (
        f"thin-slab dof mismatch: cart {x_cart.shape} vs axi {x_axi.shape}"
    )
    # Same transverse layout shifted by R; compare on x' = x - x0.
    x_cart_n = x_cart - x_cart.min()
    x_axi_n = x_axi - x_axi.min()
    assert np.allclose(x_cart_n, x_axi_n, rtol=1.0e-9, atol=1.0e-12), (
        "axisymmetric and cartesian midline x coordinates do not match "
        "after shifting by R; mesh layouts must agree for the comparison."
    )

    # Compare psi profiles. The pn built-in is order ~0.7 V; a 1%
    # relative tolerance on the L_inf difference absorbs the O(L/R)
    # curvature contamination yet catches a dispatch leak (which would
    # produce O(1) discrepancies between the two formulations).
    abs_err = float(np.max(np.abs(psi_axi - psi_cart)))
    scale = float(max(np.max(np.abs(psi_cart)), 1.0e-12))
    rel_err = abs_err / scale
    assert rel_err < 1.0e-2, (
        f"axisymmetric vs cartesian thin-slab rel_err = {rel_err:.4e} "
        f"(expected < 1e-2 at R = 100 * L)"
    )
