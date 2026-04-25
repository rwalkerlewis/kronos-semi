"""
Result types for transient simulation runs.

TransientResult is the public output type for :func:`run_transient`. It
follows the JSON-as-contract invariant: no PETSc or UFL types appear in
the public API fields; only plain Python scalars, lists, dicts, and numpy
arrays are exposed.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np


@dataclass
class TransientResult:
    """
    Output of :func:`semi.runners.transient.run_transient`.

    Attributes
    ----------
    t : list of float
        Timestamps at which IV data was recorded (one per timestep).
    iv : list of dict
        Per-timestep IV records. Each entry is a dict with keys:
        ``t`` (float), ``contact`` (str), ``V`` (float), ``J`` (float),
        ``J_n`` (float), ``J_p`` (float). One dict per timestep per
        contact.
    fields : dict mapping str to list of numpy.ndarray
        Snapshots of scalar fields (``psi``, ``n``, ``p``) at requested
        output times. Keys are field names; values are lists of numpy
        arrays, one per snapshot time. Only populated at snapshot steps
        (every ``output_every`` timesteps) to keep memory bounded.
    meta : dict
        Run metadata with keys:
        - ``order`` (int): BDF order used.
        - ``dt`` (float): fixed timestep size, seconds.
        - ``n_steps_taken`` (int): number of successful timesteps.
        - ``n_failed_steps`` (int): always 0 for fixed dt (M13).
    x_dof : numpy.ndarray or None
        DOF coordinate array from the potential space; useful for 1D
        plotting without re-loading the mesh.
    """
    t: list[float] = field(default_factory=list)
    iv: list[dict[str, Any]] = field(default_factory=list)
    fields: dict[str, list[np.ndarray]] = field(default_factory=dict)
    meta: dict[str, Any] = field(default_factory=dict)
    x_dof: np.ndarray | None = None
