# 3D resistor benchmark (planned)

Not yet implemented. Planned for Day 6 of the build.

## Physics scope

Uniformly n-doped silicon bar, rectangular 3D geometry, two ohmic contacts on opposite faces. Apply a small bias, solve full drift-diffusion, extract total current.

## Verification target

- Linear IV in the ohmic (low-field) regime
- Resistance matches $R = \rho L / A$ where $\rho = 1 / (q \mu_n n)$
- Sanity: $J = \sigma E$ held pointwise to FEM accuracy

## Purpose

Demonstrate that the framework scales to 3D. Physics is trivial compared to the MOSFET or FinFET we're deferring, but extends the code to 3D tetrahedra, and validates the current-density post-processing that will be needed for IV curve extraction in later benchmarks.

## Status

Same framework dependencies as `mos_2d/` plus:
- 3D mesh generation (builtin box works; file-based gmsh .msh also supported via schema)
- Full drift-diffusion form (Day 2-3 of the build)
