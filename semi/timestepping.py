"""
BDF (backward differentiation formula) time-stepping coefficients.

Pure-Python, no dolfinx at module scope. Provides:

    BDFCoefficients(order)  -- coefficient class for BDF1 and BDF2

Design:
    order=1  (backward Euler / BDF1): coeffs = (1, -1)
    order=2  (BDF2):                  coeffs = (3/2, -2, 1/2)

The general BDF-k formula for the discrete time derivative at step n+1:

    du/dt ~ (1/dt) * sum_{k=0}^{order} alpha_k * u^{n+1-k}

where u^{n+1} is the current (unknown) value and u^{n}, u^{n-1}, ...
are the previous (known) values stored in `u_history`.

The :meth:`apply` method evaluates this sum given a list of PAST values
(not including the current unknown). The caller is responsible for adding
the alpha_0/dt contribution of the current unknown to the UFL residual.
"""
from __future__ import annotations

import numpy as np


class BDFCoefficients:
    """
    Coefficients for backward differentiation formula time integration.

    Parameters
    ----------
    order : int
        BDF order. 1 = backward Euler, 2 = BDF2. Only 1 and 2 are
        supported for M13.

    Attributes
    ----------
    order : int
    coeffs : tuple of float
        (alpha_0, alpha_1, ..., alpha_order). alpha_0 multiplies u^{n+1}
        (the current unknown); alpha_k for k >= 1 multiply the k-th
        previous solution u^{n+1-k}.

    Methods
    -------
    apply(u_history, dt) -> numpy.ndarray
        Compute the full discrete time-derivative approximation
        sum_k(alpha_k * u_history[-k-1]) / dt, where u_history[-1] is
        the most recent value and u_history[-2] is one step back.
    """

    _SUPPORTED_COEFFS = {
        1: (1.0, -1.0),
        2: (1.5, -2.0, 0.5),
    }

    def __init__(self, order: int) -> None:
        if order not in self._SUPPORTED_COEFFS:
            raise ValueError(
                f"BDF order {order!r} is not supported; use 1 or 2."
            )
        self.order: int = order
        self.coeffs: tuple[float, ...] = self._SUPPORTED_COEFFS[order]

    def apply(self, u_history: list, dt: float) -> np.ndarray:
        """
        Evaluate the BDF discrete time derivative approximation.

        Parameters
        ----------
        u_history : list of numpy.ndarray
            Time history of the unknown, ordered from oldest to newest.
            Must contain at least `order + 1` entries.
            u_history[-1] is u^{n+1} (the newest value).
            u_history[-2] is u^n, etc.
        dt : float
            Time step size.

        Returns
        -------
        numpy.ndarray
            The discrete time derivative sum_k(alpha_k * u_history[-k-1]) / dt.
            For BDF1: (u^{n+1} - u^n) / dt
            For BDF2: (3/2 u^{n+1} - 2 u^n + 1/2 u^{n-1}) / dt
        """
        if len(u_history) < self.order + 1:
            raise ValueError(
                f"BDF{self.order} requires at least {self.order + 1} "
                f"history values; got {len(u_history)}."
            )
        result = np.zeros_like(np.asarray(u_history[-1], dtype=float))
        for k, alpha_k in enumerate(self.coeffs):
            result += alpha_k * np.asarray(u_history[-k - 1], dtype=float)
        return result / dt
