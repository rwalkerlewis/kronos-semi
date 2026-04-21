"""
Adaptive step-size controller for bias-ramp continuation.

Pure-Python (no dolfinx) so the logic can be unit-tested in isolation
and the same controller can drive any outer continuation parameter
(voltage, doping, temperature, geometry).

Design:

The controller tracks a signed step `step` plus a counter `easy_count`
of how many consecutive SNES solves converged in strictly fewer than
`easy_iter_threshold` iterations. Once `easy_count` reaches the
threshold we multiply `step` by `grow_factor` (capped at `max_step_abs`
in magnitude) and reset the counter. A failed solve halves `step`,
resets the counter, and raises `StepTooSmall` once `|step|` would fall
below `min_step_abs`.

Growth is bounded above by `max_step_abs` so a user who sets
`max_step_abs` equal to the sweep's voltage_sweep.step recovers the
Day 2 halving-only behaviour.
"""
from __future__ import annotations

from dataclasses import dataclass, field


class StepTooSmall(RuntimeError):
    """Raised when the controller would halve below `min_step_abs`."""


@dataclass
class AdaptiveStepController:
    """Signed step-size controller with grow-on-easy, halve-on-failure.

    Parameters
    ----------
    initial_step : float
        Signed initial step. Sign carries the sweep direction (positive
        for forward bias, negative for reverse).
    max_step_abs : float
        Upper bound on |step|. Growth is clamped here.
    min_step_abs : float
        Lower bound on |step|. Halving below this raises StepTooSmall.
    easy_iter_threshold : int
        A solve is "easy" when SNES converges in strictly fewer than
        this many iterations. Default 4 matches `docs/PHYSICS.md`
        bias-continuation guidance.
    grow_factor : float
        Multiplier applied to |step| after easy_iter_threshold easy
        solves in a row. Must be > 1 or growth is a no-op.
    """

    initial_step: float
    max_step_abs: float
    min_step_abs: float
    easy_iter_threshold: int = 4
    grow_factor: float = 1.5
    _step: float = field(init=False)
    _easy_count: int = field(default=0, init=False)

    def __post_init__(self) -> None:
        if self.max_step_abs <= 0.0:
            raise ValueError("max_step_abs must be positive")
        if self.min_step_abs <= 0.0:
            raise ValueError("min_step_abs must be positive")
        if self.min_step_abs > self.max_step_abs:
            raise ValueError("min_step_abs must not exceed max_step_abs")
        if self.easy_iter_threshold < 1:
            raise ValueError("easy_iter_threshold must be at least 1")
        if self.grow_factor <= 1.0:
            raise ValueError("grow_factor must be strictly greater than 1")
        self._step = self._clamp(self.initial_step)

    @property
    def step(self) -> float:
        return self._step

    @property
    def easy_count(self) -> int:
        return self._easy_count

    def on_success(self, n_iter: int) -> None:
        """Record a converged solve; may grow |step|."""
        if n_iter < self.easy_iter_threshold:
            self._easy_count += 1
            if self._easy_count >= self.easy_iter_threshold:
                sign = 1.0 if self._step >= 0.0 else -1.0
                grown = abs(self._step) * self.grow_factor
                if grown > self.max_step_abs:
                    grown = self.max_step_abs
                self._step = sign * grown
                self._easy_count = 0
        else:
            self._easy_count = 0

    def on_failure(self) -> None:
        """Record a SNES failure; halves |step| or raises StepTooSmall."""
        self._easy_count = 0
        sign = 1.0 if self._step >= 0.0 else -1.0
        halved = abs(self._step) * 0.5
        if halved < self.min_step_abs:
            raise StepTooSmall(
                f"step would fall to {halved:g} below min_step_abs="
                f"{self.min_step_abs:g}"
            )
        self._step = sign * halved

    def clamp_to_endpoint(self, remaining: float) -> float:
        """Return the signed step for the next try, clamped so the new
        V does not overshoot the sweep endpoint.

        `remaining = V_endpoint - V_prev` has the sign of the sweep
        direction; the returned value is passed to `V_try = V_prev + s`.
        """
        if remaining == 0.0:
            return 0.0
        if (remaining > 0.0) != (self._step > 0.0):
            raise ValueError(
                "remaining and step have opposite signs; the sweep "
                "appears to have reversed direction"
            )
        if abs(remaining) < abs(self._step):
            return remaining
        return self._step

    def _clamp(self, s: float) -> float:
        sign = 1.0 if s >= 0.0 else -1.0
        a = abs(s)
        if a > self.max_step_abs:
            a = self.max_step_abs
        if a < self.min_step_abs:
            a = self.min_step_abs
        return sign * a
