# ADR 0015: Cylindrical-coordinate axisymmetric formulation in the FEM core

**Date:** 2026-04-29
**Status:** Accepted
**Milestone:** M-axisym (Stage 1 of the MOSCAP roadmap)
**Related:** ADR 0011 (AC small-signal, primary unknowns), ADR 0014 (Slotboom transient primary unknowns)

---

## Context

The 2D Cartesian formulation currently shipped by `bias_sweep`,
`equilibrium`, `ac_sweep`, and `mos_cap_ac` misrepresents circular-gate
MOSCAP test structures. A real MOSCAP is a circular gate dot with a
rectangular cross-section in the `(r, z)` half-plane, and the
perimeter-to-area scaling computed on a planar mesh is wrong by a
device-radius-dependent factor. A full 3D mesh would resolve the
geometry but is wasteful for a structure with exact axial symmetry.

A stakeholder requires axisymmetric support to reproduce
Hu Fig 5-18 in a follow-up PR; this ADR closes Stage 1 by adding the
FEM-core formulation, leaving the MOSCAP-specific configuration and
benchmark to subsequent work.

## Decision

Add a top-level `solver.coordinates` configuration flag with values
`"cartesian"` (default) and `"axisymmetric"`. When `"axisymmetric"`
is selected, every volume and surface integrand assembled by the
four affected runners is multiplied by the radial coordinate `r`,
turning the 2D mesh integral into the cylindrical
`2*pi * int f r dr dz` integral (the `2*pi` cancels in the
conductance / capacitance ratios the engine reports).

The default `"cartesian"` keeps every existing benchmark JSON and
production code path bit-identical: the dispatch returns the bare
`ufl.dx` / `ufl.ds` and the residual is unchanged.

## Implementation

- The mesh's first coordinate is interpreted as `r`; the second as
  `z`. The axis of symmetry is the `r=0` line.
- The volume measure becomes `dx_axisym = r * ufl.dx`. The surface
  measure becomes `ds_axisym = r * ufl.ds`. Both wrapped helpers
  also accept `subdomain_data` / `subdomain_id` so multi-region MOS
  assembly (silicon-only space-charge, oxide-included stiffness) is
  supported.
- The gradient operator is unchanged: `ufl.grad` on a 2D mesh is
  already `(d/dr, d/dz)`, which is the correct cylindrical gradient
  under axial symmetry. Only the measure picks up the `r` factor.
- The Neumann zero-flux BC on the symmetry axis `r=0` is enforced
  automatically by the `r`-weighted measure (the surface contribution
  vanishes because the integrand carries an `r` factor that is zero
  on the axis). No explicit BC needs to be applied on `r=0`; mesh
  authors just place the leftmost facet at `r=0` (or `r=epsilon`
  with `epsilon ~ 1e-6 * device_radius` for numerical robustness).
- The `1/r` factor in the Laplacian is implicit: multiplying the
  Poisson equation through by `r` to derive the weak form removes
  it, leaving an integrand that is well-behaved at `r=0` for
  axially-symmetric smooth solutions. Standard quadrature suffices
  on elements adjacent to the axis.
- Doping, mobility, lifetime, and SRH recombination are pointwise
  fields and need no change.
- The dispatch lives in `semi/fem/coordinates.py`. The four affected
  runners (`bias_sweep`, `equilibrium`, `ac_sweep`, `mos_cap_ac`)
  read `solver.coordinates`, build a measure once via the helper,
  and thread it through `semi.physics.poisson`,
  `semi.physics.drift_diffusion`, and `semi.postprocess`. Threading
  the measure through these helpers keeps the residual builders the
  single source of truth for each PDE.

## Validation

- A cylindrical Poisson MMS on the unit `(r, z)` cylinder using a
  smooth manufactured solution `phi(r, z) = (1 - r^2) cos(pi z / 2)`
  with constant `eps`. The test asserts an L2 convergence rate
  `>= 1.9` on a CG1 mesh sequence
  `h in {1/4, 1/8, 1/16, 1/32}` (theoretical rate 2; some
  degradation near `r=0` is acceptable from the `r`-weighted
  measure).
- A thin-slab cross-validation: in the limit `R >> L` the
  axisymmetric solution must match a planar Cartesian slab at the
  same Dirichlet BCs to within 1% relative error at `R = 100 * L`.
  This catches dispatch bugs in any runner that silently uses raw
  `ufl.dx` for an integral that should be `r`-weighted.
- All M11/M12/M14 benchmarks are run with `solver.coordinates`
  defaulting to `"cartesian"` and must remain numerically
  indistinguishable from before this PR.

## Alternatives considered

- **Full 3D mesh:** rejected. Wasteful for axially-symmetric
  devices and substantially more expensive to assemble and solve.
- **1D radial only:** rejected. Loses the z-direction depletion
  profile that a MOSCAP needs.
- **Manual residual rewrites in each runner:** rejected. Touching
  every assembly site by hand is error-prone and would entangle
  the cylindrical change with the residual definitions in
  `semi.physics.*`. Threading a UFL measure through the existing
  helpers leaves the residuals untouched and lets UFL handle the
  Jacobian correctly via `ufl.derivative`.

## References

- Selberherr, *Analysis and Simulation of Semiconductor Devices*,
  §1.2 (cylindrical-coordinate Poisson and continuity).
- Sze, *Physics of Semiconductor Devices*, §2.2 (MOS in
  axisymmetric form).
- ADR 0011 — AC small-signal, primary unknowns.
- ADR 0014 — Slotboom transient primary unknowns.
- This PR (M-axisym) — Stage 1 of the MOSCAP roadmap.
