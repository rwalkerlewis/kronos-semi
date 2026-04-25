# 0008. SNES solver tolerances for the coupled drift-diffusion block

- Status: Accepted
- Date: 2026-04-25

## Context

Through M11 the three-block SNES solver in `semi/runners/bias_sweep.py`
used tolerances copied from the MMS V&V suite:

```
snes_rtol  = 1.0e-14
snes_atol  = 1.0e-14
snes_stol  = 1.0e-14
snes_max_it = 60
```

These values were appropriate for the 1D MMS problems (small meshes,
single-material, uniform doping) and for the existing five benchmarks
(pn_1d, pn_1d_bias, pn_1d_bias_reverse, mos_2d, resistor_3d). They
proved too tight for the M12 MOSFET benchmark because:

1. **Multi-order-of-magnitude residual scale disparity.** The Poisson
   block and the two continuity blocks carry residuals that differ by
   10–12 orders of magnitude in the scaled system at high injection
   levels (near the n+ source/drain implants at 5 × 10¹⁹ cm⁻³). A
   relative tolerance of 1 × 10⁻¹⁴ on the composite residual norm
   demanded that every block drop to machine-precision simultaneously,
   which the continuity blocks cannot achieve while the Poisson block
   is still driving large corrections.

2. **Iteration budget too low for fine MOSFET meshes.** 60 iterations
   was sufficient for the 1D benchmarks with smooth doping but not
   for the 2D MOSFET with abrupt Gaussian implants, where the Newton
   path from the initial guess can require more line-search sub-steps.

SNES termination with "diverged" was reported even when the physical
solution was within engineering accuracy, causing the bias continuation
to trigger unnecessary halvings and ultimately fail.

## Decision

Change the tolerances in the `solve_at` closure of `run_bias_sweep`:

| Parameter     | Before    | After     | Rationale                                      |
|---------------|-----------|-----------|------------------------------------------------|
| `snes_rtol`   | 1.0e-14   | 1.0e-10   | Achievable on all block pairs at high injection |
| `snes_atol`   | 1.0e-14   | 1.0e-7    | Absolute floor matching the continuity block scale |
| `snes_stol`   | 1.0e-14   | 1.0e-14   | Unchanged; keeps displacement convergence tight |
| `snes_max_it` | 60        | 100       | Headroom for fine meshes with abrupt implants  |

The MMS verification suite (`semi/verification/mms_dd.py`) retains its
own tolerance configuration and is **not** affected by this change.

## Validation

- All existing tests pass with the relaxed tolerances because the 1D
  benchmarks converge in 3–15 iterations, well within the new budget,
  and the relaxed `rtol`/`atol` are still far below numerical noise on
  those meshes.
- The M12 MOSFET benchmark (`benchmarks/mosfet_2d/mosfet_2d.json`)
  converges without halvings under the new tolerances.
- The new test `tests/fem/test_bias_sweep_multiregion.py::test_bias_sweep_multiregion_multistep`
  confirms that a 3-step forward ramp on a multi-region config produces
  at least 3 IV points and non-decreasing electron current.

## Consequences

**Positive:**
- The MOSFET benchmark and similar multi-region, high-injection problems
  converge reliably.
- Bias continuation no longer triggers spurious halvings near the n+
  implants.

**Negative / risk:**
- Slightly looser physical accuracy per Newton solve. The 1e-10 / 1e-7
  thresholds are still well below the engineering tolerances used in the
  verifiers (±10–20%).
- If a future kernel change introduces a problem that only converges to
  1e-8, these tolerances would accept a less accurate solution without
  warning.

## Deferred

- Per-block residual monitoring: PETSc SNES does not natively expose
  per-field norms through the `NonlinearProblem` wrapper. A monitoring
  hook could be added in a future milestone to log block-level residuals
  for diagnostic purposes.
- Automatic tolerance selection based on the doping range detected at
  mesh-build time. Deferred to M16 (physics completeness).

## Alternatives

1. **Keep 1e-14 and apply manual scaling.** Rescaling the residual vector
   per block so all blocks start at O(1) would make 1e-14 achievable, but
   requires modifying the UFL assembly and PETSc options, which touches
   `semi/physics/drift_diffusion.py` (off-limits per the M12 scope).

2. **Lower the threshold only for MOSFET configs.** Per-config tolerance
   injection would couple the runner to the schema in a fragile way and
   would not help third-party configs with similar implant profiles.

3. **Switch to a field-split preconditioner.** A proper block-Jacobi or
   Schur preconditioner would make convergence more uniform across blocks,
   but that is a larger infrastructure change deferred to M15.

## References

- PETSc SNES manual: Convergence tests `SNES_CONVERGED_FNORM_RELATIVE`,
  `SNES_CONVERGED_FNORM_ABS`, `SNES_CONVERGED_SNORM_RELATIVE`.
- Selberherr, S. (1984). *Analysis and Simulation of Semiconductor Devices*.
  Springer. Chapter 4: Discretization of the semiconductor equations.
- M4 amendment: `docs/PHYSICS.md` §"Verification & Validation" — note
  that `atol=0` was set for MMS to avoid premature termination on the
  2D continuity block initial residual.
