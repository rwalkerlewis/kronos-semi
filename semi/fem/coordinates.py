"""
Coordinate-system dispatch helpers for FEM volume / surface measures.

Cartesian assembly uses the unmodified ``ufl.dx`` and ``ufl.ds`` measures.
Axisymmetric assembly on a 2D ``(r, z)`` mesh multiplies every volume and
surface integrand by the radial coordinate ``r``, which converts a 2D
mesh integral into the physical cylindrical-coordinate integral
``2*pi*int f r dr dz`` (with the ``2*pi`` factor cancelling in the
conductance / capacitance ratios that the engine reports).

See ``docs/adr/0015-axisymmetric-formulation.md`` for the full derivation
and the rationale for handling the dispatch through measures rather than
through manual residual rewrites.

The helpers here are deliberately thin wrappers around UFL: they take
the mesh and return a measure that can be used anywhere the caller
would otherwise have written ``ufl.dx`` or ``ufl.ds``.
"""
from __future__ import annotations

from typing import Any

CARTESIAN = "cartesian"
AXISYMMETRIC = "axisymmetric"

_VALID_COORDINATES = (CARTESIAN, AXISYMMETRIC)


def _check_coordinates(coordinates: str) -> None:
    if coordinates not in _VALID_COORDINATES:
        raise ValueError(
            f"Unknown coordinates: {coordinates!r}. "
            f"Expected one of {_VALID_COORDINATES}."
        )


def get_volume_measure(
    mesh,
    coordinates: str = CARTESIAN,
    *,
    subdomain_data: Any = None,
    subdomain_id: Any = None,
):
    """Return the appropriate UFL volume measure.

    Parameters
    ----------
    mesh
        The dolfinx mesh on which the form is being assembled. Used only
        to source the ``ufl.SpatialCoordinate`` for the axisymmetric
        case; the cartesian branch ignores it.
    coordinates
        ``"cartesian"`` (default) returns the bare ``ufl.dx`` (or its
        subdomain-restricted equivalent). ``"axisymmetric"`` returns
        ``r * ufl.dx`` where ``r`` is the first spatial coordinate of
        ``mesh``.
    subdomain_data, subdomain_id
        Forwarded to ``ufl.Measure`` so callers can restrict the
        integration to specific tagged regions (e.g. silicon-only
        cells in the multi-region MOS Poisson assembly). When neither
        is provided the bare ``ufl.dx`` is returned.

    Raises
    ------
    ValueError
        If ``coordinates`` is not ``"cartesian"`` or ``"axisymmetric"``.
    """
    import ufl

    _check_coordinates(coordinates)

    if subdomain_data is not None or subdomain_id is not None:
        kwargs: dict[str, Any] = {"domain": mesh}
        if subdomain_data is not None:
            kwargs["subdomain_data"] = subdomain_data
        if subdomain_id is not None:
            kwargs["subdomain_id"] = subdomain_id
        dx = ufl.Measure("dx", **kwargs)
    else:
        dx = ufl.dx

    if coordinates == AXISYMMETRIC:
        r = ufl.SpatialCoordinate(mesh)[0]
        return r * dx
    return dx


def get_surface_measure(
    mesh,
    coordinates: str = CARTESIAN,
    *,
    subdomain_data: Any = None,
    subdomain_id: Any = None,
):
    """Return the appropriate UFL surface measure.

    Mirrors :func:`get_volume_measure` for the boundary measure
    ``ufl.ds``. In axisymmetric mode the measure is multiplied by the
    radial coordinate ``r``, so a contact at ``z=z0`` spanning
    ``r in [r1, r2]`` integrates against ``r dr`` (the cylindrical
    surface measure) automatically. The natural Neumann zero-flux BC
    on the symmetry axis ``r=0`` is enforced by the same ``r``
    weighting (the integrand vanishes on that boundary).

    Raises
    ------
    ValueError
        If ``coordinates`` is not ``"cartesian"`` or ``"axisymmetric"``.
    """
    import ufl

    _check_coordinates(coordinates)

    if subdomain_data is not None or subdomain_id is not None:
        kwargs: dict[str, Any] = {"domain": mesh}
        if subdomain_data is not None:
            kwargs["subdomain_data"] = subdomain_data
        if subdomain_id is not None:
            kwargs["subdomain_id"] = subdomain_id
        ds = ufl.Measure("ds", **kwargs)
    else:
        ds = ufl.ds

    if coordinates == AXISYMMETRIC:
        r = ufl.SpatialCoordinate(mesh)[0]
        return r * ds
    return ds


def resolve_coordinates(cfg: dict) -> str:
    """Read ``solver.coordinates`` from a config dict, defaulting to cartesian.

    Centralised here so every runner uses the same fallback behaviour and
    the schema default is honoured even for configs that omit
    ``solver.coordinates`` entirely.
    """
    coords = cfg.get("solver", {}).get("coordinates", CARTESIAN)
    _check_coordinates(coords)
    return coords
