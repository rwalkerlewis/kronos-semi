"""
Result types for transient and AC small-signal simulation runs.

TransientResult is the public output type for :func:`run_transient`.
AcSweepResult is the public output type for :func:`run_ac_sweep` (M14).
Both follow the JSON-as-contract invariant: no PETSc or UFL types appear
in the public API fields; only plain Python scalars, lists, dicts, and
numpy arrays are exposed. Complex numbers in AcSweepResult are
serialized to JSON as ``{"re": <float>, "im": <float>}`` by the
artifact writer (see :mod:`semi.io.artifact`).
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


@dataclass
class AcSweepResult:
    """
    Output of :func:`semi.runners.ac_sweep.run_ac_sweep` (M14).

    Small-signal AC analysis. The system is linearised around a DC
    operating point u_0 = (psi_0, n_0, p_0); the perturbation
    delta_u(omega) satisfies (J + j*omega*M) * delta_u = -dF/dV * delta_V.
    Terminal current and admittance are evaluated per unit perturbation
    voltage at the swept contact.

    Attributes
    ----------
    frequencies : list of float
        Frequency points (Hz).
    Y : list of complex
        Small-signal admittance per unit voltage perturbation at the
        swept contact, in S/m^2 (1D), S/m (2D), or S (3D), matching the
        per-unit-area convention of the conduction current evaluator.
    Z : list of complex
        Small-signal impedance, 1/Y, with the same per-unit-area
        convention. Entries with |Y| < 1e-300 are reported as 0.0 (no
        terminal response).
    C : list of float
        Effective capacitance C(omega) = +Im(Y) / (2*pi*f). The runner
        reports Y with the "positive = current INTO the device"
        convention used by ``postprocess.evaluate_current_at_contact``
        and by ``bias_sweep`` (i.e. ``Re(Y(omega->0))`` matches
        ``dI/dV`` from a centered finite difference at the same
        ``V_DC``). In that convention an ideal capacitor terminal has
        ``Y = +j*omega*C``, so a positive depletion capacitance maps
        to a positive ``Im(Y)`` and to a positive ``C``. Same
        per-unit convention as Y. Entries at f=0 are reported as 0.0.
    G : list of float
        Effective conductance G(omega) = Re(Y). Same per-unit convention.
    dc_bias : dict
        DC operating point used. Keys: ``contact`` (str), ``voltage``
        (float). Mirrors the relevant slice of ``cfg["solver"]``.
    meta : dict
        Run metadata. Keys include:
        - ``n_freqs`` (int): number of frequency points.
        - ``perturbation_voltage`` (float): finite-difference epsilon
          used for the DC sensitivity (V).
        - ``ac_amplitude`` (float): AC amplitude (V) reported back from
          the cfg; admittance is per unit volt regardless.
        - ``dc_solver_iterations`` (int): SNES iterations of the DC
          solve (max over the V_DC and V_DC+eps solves).
        - ``backend`` (str): ``"real-2x2-block"`` (always, until M16
          enables a complex PETSc build).
        - ``displacement_current_included`` (bool): always True for M14.
    """
    frequencies: list[float] = field(default_factory=list)
    Y: list[complex] = field(default_factory=list)
    Z: list[complex] = field(default_factory=list)
    C: list[float] = field(default_factory=list)
    G: list[float] = field(default_factory=list)
    dc_bias: dict[str, Any] = field(default_factory=dict)
    meta: dict[str, Any] = field(default_factory=dict)
