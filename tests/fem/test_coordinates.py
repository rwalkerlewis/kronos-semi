"""
Unit tests for the cylindrical-coordinate dispatch helpers
(``semi/fem/coordinates.py``). See ADR 0015.

The tests exercise:

* ``get_volume_measure`` returns the bare ``ufl.dx`` for ``"cartesian"``
  and an ``r``-weighted measure for ``"axisymmetric"`` whose integral of
  a constant ``f`` on a 2D rectangle agrees with the analytic
  ``f * (r_max^2 - r_min^2)/2 * (z_max - z_min)``.
* ``get_surface_measure`` analogously: ``ufl.ds`` for ``"cartesian"``,
  ``r * ufl.ds`` for ``"axisymmetric"``.
* Invalid coordinate strings raise ``ValueError`` from both helpers.
"""
from __future__ import annotations

import math

import numpy as np
import pytest


def _build_rect_mesh(Nx: int, Ny: int, x0: float, x1: float, y0: float, y1: float):
    """Build a structured 2D rectangle mesh on [x0, x1] x [y0, y1]."""
    from dolfinx.mesh import CellType, create_rectangle
    from mpi4py import MPI

    return create_rectangle(
        MPI.COMM_WORLD,
        [np.array([x0, y0]), np.array([x1, y1])],
        [Nx, Ny],
        cell_type=CellType.triangle,
    )


def _assemble_constant(msh, dx_measure, value: float = 1.0) -> float:
    """Assemble ``value * dx_measure`` to a real scalar (MPI-reduced)."""
    from dolfinx import fem
    from mpi4py import MPI
    from petsc4py import PETSc

    c = fem.Constant(msh, PETSc.ScalarType(value))
    form = fem.form(c * dx_measure)
    local = fem.assemble_scalar(form)
    return float(msh.comm.allreduce(local, op=MPI.SUM))


def test_get_volume_measure_cartesian_returns_ufl_dx():
    """Cartesian dispatch returns the bare ``ufl.dx`` object."""
    import ufl

    from semi.fem.coordinates import get_volume_measure

    msh = _build_rect_mesh(2, 2, 0.0, 1.0, 0.0, 1.0)
    dx = get_volume_measure(msh, "cartesian")
    assert dx is ufl.dx, "cartesian volume measure must be the unmodified ufl.dx"


def test_get_volume_measure_cartesian_integrates_to_area():
    """Integrating ``1`` over [0,2] x [0,3] in cartesian gives area 6."""
    from semi.fem.coordinates import get_volume_measure

    msh = _build_rect_mesh(8, 6, 0.0, 2.0, 0.0, 3.0)
    dx = get_volume_measure(msh, "cartesian")
    area = _assemble_constant(msh, dx, 1.0)
    assert math.isclose(area, 6.0, rel_tol=1.0e-12)


def test_get_volume_measure_axisymmetric_picks_up_r_factor():
    """In axisymmetric mode the integrand is multiplied by r, so

        integral_{r=r0..r1, z=z0..z1} 1 * r dr dz
            = (r1^2 - r0^2)/2 * (z1 - z0)

    The test asserts this exactly (both sides are linear in the FE
    polynomial space, so quadrature is exact for triangle CG1).
    """
    from semi.fem.coordinates import get_volume_measure

    r0, r1 = 0.5, 2.5
    z0, z1 = 0.0, 1.5
    msh = _build_rect_mesh(16, 12, r0, r1, z0, z1)
    dx = get_volume_measure(msh, "axisymmetric")
    weighted = _assemble_constant(msh, dx, 1.0)
    expected = 0.5 * (r1 ** 2 - r0 ** 2) * (z1 - z0)
    assert math.isclose(weighted, expected, rel_tol=1.0e-10), (
        f"axisymmetric integral {weighted} != analytic {expected}"
    )


def test_get_volume_measure_axisymmetric_at_axis():
    """When ``r0 = 0`` the r-weighted integral gives r1^2/2 * (z1-z0).

    Catches off-by-one bugs in the SpatialCoordinate handling at the
    symmetry axis. The integrand vanishes on the r=0 boundary, which is
    the natural Neumann zero-flux condition described in ADR 0015.
    """
    from semi.fem.coordinates import get_volume_measure

    r1 = 1.0
    z1 = 1.0
    msh = _build_rect_mesh(16, 16, 0.0, r1, 0.0, z1)
    dx = get_volume_measure(msh, "axisymmetric")
    weighted = _assemble_constant(msh, dx, 1.0)
    expected = 0.5 * r1 ** 2 * z1
    assert math.isclose(weighted, expected, rel_tol=1.0e-10)


def test_get_surface_measure_cartesian_returns_ufl_ds():
    """Cartesian dispatch returns the bare ``ufl.ds`` object."""
    import ufl

    from semi.fem.coordinates import get_surface_measure

    msh = _build_rect_mesh(2, 2, 0.0, 1.0, 0.0, 1.0)
    ds = get_surface_measure(msh, "cartesian")
    assert ds is ufl.ds


def test_get_surface_measure_axisymmetric_picks_up_r_factor():
    """Boundary integral of 1 over the four sides of [r0,r1]x[z0,z1].

    Cartesian gives perimeter = 2*(r1-r0) + 2*(z1-z0).
    Axisymmetric gives:
      bottom (z=z0)  : integral_r0^r1 r dr = (r1^2-r0^2)/2
      top    (z=z1)  : (r1^2-r0^2)/2
      left   (r=r0)  : r0 * (z1-z0)
      right  (r=r1)  : r1 * (z1-z0)
    Total: (r1^2 - r0^2) + (r0+r1)*(z1-z0).
    """
    from semi.fem.coordinates import get_surface_measure

    r0, r1 = 0.5, 2.5
    z0, z1 = 0.0, 1.5
    msh = _build_rect_mesh(16, 12, r0, r1, z0, z1)
    ds = get_surface_measure(msh, "axisymmetric")
    weighted = _assemble_constant(msh, ds, 1.0)
    expected = (r1 ** 2 - r0 ** 2) + (r0 + r1) * (z1 - z0)
    assert math.isclose(weighted, expected, rel_tol=1.0e-10), (
        f"axisymmetric ds integral {weighted} != analytic {expected}"
    )


@pytest.mark.parametrize("bad", ["spherical", "cyl", "", "Cartesian", "AXI"])
def test_invalid_coordinates_raises(bad):
    """Both helpers reject unknown coordinate strings with ValueError."""
    from semi.fem.coordinates import get_surface_measure, get_volume_measure

    msh = _build_rect_mesh(2, 2, 0.0, 1.0, 0.0, 1.0)
    with pytest.raises(ValueError, match="Unknown coordinates"):
        get_volume_measure(msh, bad)
    with pytest.raises(ValueError, match="Unknown coordinates"):
        get_surface_measure(msh, bad)


def test_resolve_coordinates_default_is_cartesian():
    """resolve_coordinates returns 'cartesian' when solver.coordinates absent."""
    from semi.fem.coordinates import resolve_coordinates

    assert resolve_coordinates({}) == "cartesian"
    assert resolve_coordinates({"solver": {}}) == "cartesian"
    assert resolve_coordinates(
        {"solver": {"coordinates": "axisymmetric"}}
    ) == "axisymmetric"


def test_resolve_coordinates_rejects_invalid():
    """resolve_coordinates surfaces the same ValueError for invalid inputs."""
    from semi.fem.coordinates import resolve_coordinates

    with pytest.raises(ValueError, match="Unknown coordinates"):
        resolve_coordinates({"solver": {"coordinates": "spherical"}})
