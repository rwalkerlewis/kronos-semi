"""
Doping profile evaluation.

Builds callables N_net(x) from the JSON doping specification. The output
convention is SIGNED net doping: N_D - N_A, in m^-3. Positive means donor-
dominated (n-type), negative means acceptor-dominated (p-type).

Supports:
    - uniform: constant in a region
    - step: abrupt change along an axis
    - gaussian: analytic Gaussian (common for source/drain implants)
    - (later) table: interpolated from an external file
"""
from __future__ import annotations

from collections.abc import Callable

import numpy as np

from .constants import cm3_to_m3


def build_profile(doping_list: list[dict]) -> Callable[[np.ndarray], np.ndarray]:
    """
    Combine a list of doping entries into a single callable N_net(x).

    Parameters
    ----------
    doping_list : list of dict
        Each entry has keys 'region' and 'profile'. 'profile' follows the
        JSON schema: uniform / step / gaussian.

    Returns
    -------
    Callable
        Function that takes an array x of shape (dim, N) or (N,) and returns
        an array of shape (N,) with the net doping in m^-3.

    Notes
    -----
    Multiple entries are summed, so you can overlay a Gaussian implant on
    top of a uniform background by using two entries. The 'region' field
    is informational until we support region-restricted evaluation.
    """
    evaluators = [_build_one(d["profile"]) for d in doping_list]

    def net_doping(x: np.ndarray) -> np.ndarray:
        # Accept (dim, N) from dolfinx interpolate or (N,) for 1D scalar x
        if x.ndim == 1:
            x = x[np.newaxis, :]
        total = np.zeros(x.shape[1])
        for ev in evaluators:
            total = total + ev(x)
        return total

    return net_doping


def _build_one(profile: dict) -> Callable[[np.ndarray], np.ndarray]:
    t = profile["type"]
    if t == "uniform":
        return _uniform(profile)
    if t == "step":
        return _step(profile)
    if t == "gaussian":
        return _gaussian(profile)
    raise ValueError(f"Unknown doping profile type {t!r}")


def _uniform(p: dict) -> Callable[[np.ndarray], np.ndarray]:
    N_net = cm3_to_m3(p["N_D"] - p["N_A"])

    def f(x: np.ndarray) -> np.ndarray:
        return np.full(x.shape[1], N_net)

    return f


def _step(p: dict) -> Callable[[np.ndarray], np.ndarray]:
    axis = p["axis"]
    loc = p["location"]
    N_net_left  = cm3_to_m3(p["N_D_left"]  - p["N_A_left"])
    N_net_right = cm3_to_m3(p["N_D_right"] - p["N_A_right"])

    def f(x: np.ndarray) -> np.ndarray:
        xa = x[axis]
        return np.where(xa < loc, N_net_left, N_net_right)

    return f


def _gaussian(p: dict) -> Callable[[np.ndarray], np.ndarray]:
    center = np.asarray(p["center"], dtype=float)
    sigma = np.asarray(p["sigma"], dtype=float)
    peak = cm3_to_m3(p["peak"])
    sign = +1.0 if p["dopant"] == "donor" else -1.0
    bg_net = cm3_to_m3(p.get("background_N_D", 0.0) - p.get("background_N_A", 0.0))

    def f(x: np.ndarray) -> np.ndarray:
        # x is (dim, N); center/sigma are (dim,)
        dim = len(center)
        r2 = np.zeros(x.shape[1])
        for d in range(dim):
            r2 = r2 + ((x[d] - center[d]) / sigma[d]) ** 2
        return bg_net + sign * peak * np.exp(-0.5 * r2)

    return f


def evaluate_at_points(doping_list: list[dict], points: np.ndarray) -> np.ndarray:
    """Convenience: build profile and evaluate at given points."""
    f = build_profile(doping_list)
    return f(points)
