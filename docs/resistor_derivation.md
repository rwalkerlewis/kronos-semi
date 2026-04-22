# Resistor Derivation: 3D Doped Resistor (Day 7 Gate)

This document is the mathematical specification for the Day-7 3D doped
resistor benchmark, the V-I linearity verifier, and the gmsh `.msh`
loader test strategy. It is a derivation-lite gate in the same spirit
as `docs/mos_derivation.md` for Day 6, but much shorter because the
physics is unchanged from Days 1-6: equilibrium Poisson plus coupled
Slotboom drift-diffusion with SRH recombination, already MMS-verified
in 1D and 2D.

What is new in Day 7 is structural, not physical:

1. A 3D rectangular-bar device that exercises the solver in 3D for the
   first time.
2. A new verifier, V-I linearity, that checks the simulated current
   scales linearly with applied voltage to within 1 percent.
3. A gmsh `.msh` loader wired into `semi/mesh.py::_build_from_file`
   via `dolfinx.io.gmsh`, with a committed fixture mesh for
   regression testing.

No implementation code (other than PLAN.md and this document) is
written until this document is reviewed and approved.

Scaled symbols follow the established Day-2 / Day-4 / Day-6
conventions (mesh in meters, psi and quasi-Fermi potentials in units
of V_t = k_B T / q, densities in units of C_0).

## Section 1: Device geometry

A 3D rectangular bar aligned with the x axis:

- Length along x: `L = 1 um = 1e-6 m` (the current-carrying direction).
- Cross-section in the (y, z) plane: square, `W = 200 nm = 2e-7 m`.
- Cross-sectional area `A = W^2 = 4e-14 m^2`.
- Bounding box `[0, L] x [0, W] x [0, W]` in meters.
- Single material region: silicon, role `semiconductor`, uniform
  n-type doping `N_D = 1e18 cm^-3 = 1e24 m^-3`. Minority carriers
  are negligible: `p_eq = n_i^2 / N_D = (1e10 cm^-3)^2 / 1e18 cm^-3
  = 1e2 cm^-3`, sixteen orders of magnitude below N_D.
- Two ohmic contacts, one on each x-face:
  - `contact_left`: the `x = 0` face. Swept reference, typically
    `V = 0` or the negative end of the sweep.
  - `contact_right`: the `x = L` face. Swept from `-0.01 V` to
    `+0.01 V` in 5 points (including 0 V).
- All other faces (y = 0, y = W, z = 0, z = W) have no Dirichlet BC;
  they are naturally insulating for both Poisson (zero normal flux of
  eps grad psi) and continuity (zero normal flux of J_n and J_p).

Mesh: the baseline benchmark uses the builtin `create_box` with
resolution `[64, 16, 16]` cells along (x, y, z). The second benchmark
loads a gmsh `.msh` fixture with the same nominal geometry but
unstructured tetrahedra.

Region and facet tagging reuses the existing pure-Python box-tagger
for the builtin path. The gmsh path reads physical groups directly
from the `.msh` file and returns them as `cell_tags` and `facet_tags`
from `_build_from_file`.

## Section 2: Analytical resistance

For a uniform n-type bar with constant mobility in a low-field,
low-injection regime, the steady-state DC current collapses to the
standard textbook expression for an ohmic resistor:

```
R = L / (q * N_D * mu_n * A)
I = V / R
```

with:

- `q = 1.602176634e-19 C` (from `semi/constants.py`).
- `N_D = 1e18 cm^-3 = 1.0e24 m^-3`.
- `mu_n` from `semi/materials.py`: the Si entry uses
  `mu_n = cm2_to_m2(1400.0) = 0.14 m^2/(V s)`. The roadmap prompt
  quoted 1350 cm^2/(V s) as a representative value, but the value
  the solver actually uses is 1400 cm^2/(V s). Theory and simulation
  share the same constant, so the 1% tolerance applies regardless of
  which numerical value is used, as long as both sides are consistent.
- `L = 1e-6 m`, `A = (2e-7 m)^2 = 4e-14 m^2`.

Substituting:

```
q * N_D * mu_n * A = 1.602176634e-19 * 1.0e24 * 0.14 * 4e-14
                   = 8.972e-10  (units: S * m)
R = L / (q * N_D * mu_n * A)
  = 1.0e-6 / 8.972e-10
  = 1114.6 Ohm
```

Predicted currents at the three reference biases:

| V (V)  | I = V/R (A)  |
|-------:|-------------:|
|  0.001 | 8.972e-7     |
|  0.010 | 8.972e-6     |
|  0.100 | 8.972e-5     |

All three are comfortably in the ohmic regime: `|V| << V_bi` (no
built-in junction barrier), `|V| <= 0.1 V = ~4 V_t` (no high-field
nonlinearity), and the bar is uniformly doped n-type (no majority-
to-minority flipping). The Day-7 benchmark caps its sweep at 0.01 V
to stay deep inside this regime, well below `V_t = 0.02585 V`. The
higher-bias values in the table are recorded here only to make the
scaling check plain.

## Section 3: V-I linearity verifier

`verify_resistor_3d` (to be added alongside the existing
`verify_pn_1d_bias`, `verify_pn_1d_bias_reverse`, `verify_mos_2d`
entries in `scripts/run_benchmark.py`) consumes the IV rows produced
by the `bias_sweep` runner and performs the following steps:

1. Sweep grid: bias the `contact_right` contact across 5 equally
   spaced points in `[-0.01, +0.01] V`, i.e.
   `V in {-0.010, -0.005, 0.000, +0.005, +0.010} V`.
2. At each sweep point the runner extracts the total current as the
   UFL facet integral of `J_n . n` over the `x = L` contact. The
   result table rows look like `(contact_name, V, I)` and are stored
   in `result.iv`.
3. For each nonzero bias point, compute `R_sim(V) = V / I(V)`.
4. Compute `R_theory` from the closed form in Section 2, reading
   `mu_n` from `semi/materials.py` and `N_D` from the benchmark JSON
   so the verifier is not a string-match.
5. Linearity metric:
   `err = max_over_nonzero_V( |R_sim(V) - R_theory| / R_theory )`.
6. Tolerance: `err < 0.01` (1 percent). A single PASS / FAIL line is
   emitted with the worst point cited.

Sanity checks emitted as additional PASS lines:

- `I(V = 0) < tol_zero` where `tol_zero = R_theory * 1e-6` (the zero-
  bias current must be numerical noise, not a real offset). An offset
  larger than this indicates a sign error in the facet-integral
  current extraction or a nonzero equilibrium current path through
  the structure (which would be a physics bug).
- Sign check: `sign(I(V)) == sign(V)` for every nonzero V. A
  reversed sign means the contact normals are flipped or the
  Dirichlet BC application swapped the contacts.

### Rationale for the 1 percent tolerance

This tolerance is deliberately tighter than the 5-15 percent used
for the pn-junction benchmarks. In those cases the theory reference
is itself an approximation: depletion-region width, Shockley long-
diode current, SNS composite all have order-15 percent modeling
error baked into the closed form. For the ohmic-resistor regime,
the theory expression `I = V / R` is exact in the limit
`V -> 0, eps_r, mu_n, N_D` uniform, and the benchmark operates
deep inside that limit:

- Maximum applied voltage `|V_max| = 0.01 V << V_t = 0.0259 V`, so
  carrier densities deviate from equilibrium by at most
  `exp(V/V_t) - 1 ~= 0.48` in the worst case. Majority carrier
  density along the bar stays within 1 percent of `N_D` because
  `N_D >> n_i exp(V/V_t)` and the neutrality condition is
  essentially unperturbed.
- No junction, no oxide, no gate: there is no barrier that could
  move with bias, so no second-order `dW/dV` correction.
- No SRH recombination effect in bulk on a 1 um bar at this doping:
  the transit time is `L^2 / (mu_n * V) ~ 1e-12 / (0.14 * 0.01)
  = 7e-10 s`, four orders of magnitude shorter than
  `tau_n = tau_p = 1e-8 s`, so generation-recombination is a
  higher-order correction.

Within this regime, any deviation greater than 1 percent indicates
a bug, not a legitimate physics correction. Common culprits, in
order of likelihood:

1. Current-extraction sign or normal-orientation bug in the facet
   integral (the 2D MOS code path used oriented normals; 3D adds
   a third axis that is easy to mis-orient).
2. Boundary-condition sign error on the ohmic contact psi (the
   equilibrium psi at N_D = 1e18 cm^-3 is not zero, it is
   `psi_eq_hat = asinh(N_D / (2 n_i)) ~= ln(N_D / n_i) = ln(1e8)
   = +18.42` in V_t units, as applied in
   `semi/bcs.py::build_psi_dirichlet_bcs`).
3. Mesh resolution too coarse along x for the linear psi gradient.
   `[64, 16, 16]` is empirically enough, but a coarser mesh can
   inflate the quadrature error on the extracted facet integral.
4. Nondimensionalization bug that mis-scales J by `V_t / L_0` or
   `C_0 mu_hat`; these would show up as a constant prefactor on
   every R_sim point.

The project prompt explicitly flags widening the tolerance as an
anti-pattern. Any non-PASS is a stop-and-diagnose trigger, not a
tolerance-loosening trigger.

## Section 4: 3D slice plotting

Two sanity-check plots are emitted for every run:

- `psi_slice_y_midplane.png`: scalar slice of psi at the
  `y = W / 2` plane, showing psi(x, z) at the midplane. Expected
  shape: approximately linear gradient from `psi_L` at `x = 0` to
  `psi_L - V_applied` at `x = L`, constant along z in the bulk.
- `jn_slice_y_midplane.png`: scalar slice of the electron
  current-density magnitude `|J_n|` on the same plane. Expected
  shape: approximately uniform across the bulk, equal to `I / A`,
  with minor boundary-layer adjustment within a few mesh cells of
  the ohmic contacts.

These are visual regressions only, not gated verifiers. They must
render without error and must show the expected monotone / uniform
shape. If the C-V code path's matplotlib + tricontour is awkward
for a 3D submesh extraction, the fallback is to export psi and J_n
to VTX or XDMF with `dolfinx.io.VTXWriter` and post-process with
pyvista in a separate script that the plotter invokes. The chosen
approach is documented in a short comment in `scripts/run_benchmark.py`.

The 1D and 2D plotters (`plot_pn_1d`, `plot_mos_2d`) are not
modified. The 3D path is added as a `dim == 3` branch in the
plot dispatcher.

## Section 5: Gmsh loader test strategy

The gmsh `.msh` loader is implemented in `_build_from_file` in
`semi/mesh.py`. It reads a `.msh` file via
`dolfinx.io.gmsh.read_from_msh` and returns `(mesh, cell_tags,
facet_tags)`. When the caller is `build_mesh`, the returned tags are
used directly; the box-tagger is skipped for file-sourced meshes
because the gmsh physical groups are authoritative.

A fixture is committed under `benchmarks/resistor_3d/fixtures/`:

- `box.geo`: the reproducible gmsh source. Rectangular box with
  the same dimensions as the builtin (1 um x 200 nm x 200 nm),
  with physical groups for the silicon volume (cell tag) and the
  two ohmic-contact faces (facet tags). Natural-BC faces are left
  untagged.
- `box.msh`: the output of `gmsh -3 box.geo -o box.msh`. Must be
  under 100 KB.

Three tests go in `tests/fem/test_mesh_gmsh.py`:

1. Round-trip: load `box.msh`, assert that the cell count matches
   the builtin `create_box` mesh within 5 percent (gmsh produces a
   slightly different tetrahedralization from the builtin for the
   same box) and that the bounding box matches to 1e-12 m.
2. Physical-group preservation: assert that the silicon volume tag
   and the two contact facet tags are present with the expected
   tag IDs, and that facet-tag entity counts match the contact-face
   counts expected from the geometry (6 triangular facets on each
   200 nm x 200 nm face for a modest mesh resolution; exact count
   depends on gmsh, so we assert `> 0` and consistency between the
   two contact faces rather than a specific integer).
3. End-to-end comparison: run the `resistor_3d` benchmark once with
   the builtin mesh and once with the `.msh` fixture. Both must
   produce V-I linearity within the same 1 percent tolerance, and
   R_sim must agree between the two variants within 1 percent.

The project prompt notes that full MMS verification in 3D is not
required for Day 7: the Poisson and drift-diffusion assembly is
dimension agnostic and has been MMS-verified in 1D and 2D, and the
ohmic-resistor theory match on an unstructured tetrahedral mesh is
sufficient evidence that the 3D code path and the gmsh loader are
correct. A 3D MMS convergence study is deferred to Day 8 or a
follow-on PR if time allows.

## Summary

Day 7 is a structural extension, not a physics extension. The
gate-keeping document is short because the only genuinely new
acceptance criterion is the 1 percent V-I linearity tolerance, and
the only new infrastructure is the gmsh loader. Both are specified
in enough detail here that the implementation in Phases 3 and 4 is
mechanical: define the JSON, add the `dim == 3` branches, register
the verifier, and ship the fixture.
