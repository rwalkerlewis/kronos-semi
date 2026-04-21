"""
Shared L^2 / H^1 error-norm helpers for verification studies.

Both helpers return assembled scalar squared norms so callers compose
into the convergence pipeline without thinking about UFL syntax.
"""
from __future__ import annotations


def l2_error_squared(u_h, u_exact_ufl) -> float:
    """
    Squared L^2 error: integral of (u_h - u_exact)^2 over the mesh.

    Parameters
    ----------
    u_h
        dolfinx.fem.Function holding the discrete solution.
    u_exact_ufl
        Any UFL expression, typically a function of SpatialCoordinate.

    Returns
    -------
    float
        The assembled scalar (already MPI-reduced over the comm).
    """
    import ufl
    from dolfinx import fem

    diff = u_h - u_exact_ufl
    form = fem.form(diff * diff * ufl.dx)
    return float(fem.assemble_scalar(form))


def h1_seminorm_error_squared(u_h, u_exact_ufl) -> float:
    """
    Squared H^1 seminorm error: integral of |grad(u_h - u_exact)|^2.
    """
    import ufl
    from dolfinx import fem

    diff = u_h - u_exact_ufl
    form = fem.form(ufl.inner(ufl.grad(diff), ufl.grad(diff)) * ufl.dx)
    return float(fem.assemble_scalar(form))
