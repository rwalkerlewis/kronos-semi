# Axisymmetric (cylindrical) weak forms

Schema 1.3.0 added the top-level `coordinate_system` field with
values `"cartesian"` (default, unchanged) or `"axisymmetric"`. The
axisymmetric path solves on the meridian half-plane $(r, z)$ with
$r = x_0 \ge 0$ and $z = x_1$.

## Geometry and the $r$-weighted measure

The physical 3D domain is recovered by revolving the half-plane
about the $z$-axis. The volume measure becomes

$$\mathrm{d}V_{3D} = 2\pi\,r\,\mathrm{d}r\,\mathrm{d}z,$$

so every volume integrand of the cartesian weak forms picks up a
factor $r$. The constant $2\pi$ cancels on both sides of the
residual.

For a scalar PDE $-\nabla\!\cdot\!(A\,\nabla u) = f$, the Galerkin
weak form becomes

$$\int_\Omega A\,\nabla u \cdot \nabla v\;r\,\mathrm{d}r\,\mathrm{d}z
\;=\; \int_\Omega f\,v\,r\,\mathrm{d}r\,\mathrm{d}z,$$

with $A$ possibly piecewise (e.g. $\varepsilon_r$ jumping at the
Si/SiO$_2$ interface). The same weighting applies to the
drift-diffusion continuity rows in Slotboom form.

## Boundary conditions

- **Symmetry axis $r = 0$** is a *natural* boundary. The $r$-weighted
  measure makes test-function weights vanish there automatically; no
  Dirichlet row is needed, and the schema rejects Dirichlet contacts
  on $r = 0$ (`semi.schema._validate_coordinate_system`).
- **Outer radial wall $r = R$** is treated as homogeneous Neumann
  (no flux). The user must choose $R$ large enough that the solution
  near the axis is insensitive to the cutoff. For a MOSCAP this is
  approximately $R \gtrsim 5\,W_{\text{dmax}}$.
- **Top and bottom planar boundaries** follow the usual
  Dirichlet/Neumann contact specification.

## Schema enforcement

Cross-field validation in
[`semi/schema.py`](../../semi/schema.py) enforces:

- `coordinate_system == "axisymmetric"` requires `dimension == 2`.
- Mesh radial extent is non-negative.
- No Dirichlet contact may live on the symmetry axis $r = 0$.

## Where this lives in the code

- [`semi/physics/axisymmetric.py`](../../semi/physics/axisymmetric.py) —
  $r$-weighted Poisson (`build_equilibrium_poisson_form_axisym`,
  `build_equilibrium_poisson_form_axisym_mr`) and the Slotboom
  continuity weak forms.
- [`semi/schema.py`](../../semi/schema.py) — `_validate_coordinate_system`
  cross-field check.
- [`benchmarks/moscap_axisym_2d/`](../../benchmarks/moscap_axisym_2d/) —
  reference benchmark, gmsh `.geo` template, analytical CSV.

## Verification

- [`tests/test_coordinate_system.py`](../../tests/test_coordinate_system.py)
  — schema validation (7 tests).
- [`tests/test_axisym_moscap_math.py`](../../tests/test_axisym_moscap_math.py)
  — pure-Python analytical anchors (15 anchors, all green).
- [`tests/test_moscap_axisym_cv.py`](../../tests/test_moscap_axisym_cv.py)
  — dolfinx-gated FEM smoke test.

## References

- [`docs/PHYSICS.md`](../PHYSICS.md) section 2 (weak forms in scaled units).
- Hu, *Modern Semiconductor Devices for IC*, chapter 5 (MOSCAP).
