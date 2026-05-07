# moscap_cv_oxide_thickness

## What this is

A practical MOS capacitor C-V characterization across three gate-oxide
thicknesses (2 nm, 5 nm, 10 nm). Geometry is 2D Cartesian, mirroring
[`benchmarks/mos_2d`](../../benchmarks/mos_2d) structurally: 500 nm Si
depth (uniform p-type body, N_A = 5e16 cm^-3), uniform 1 nm vertical
cells, 4 lateral cells of 125 nm each, SiO2 gate oxide on top. The new
axis is oxide thickness: the anchor JSON pins t_ox = 5 nm and the
smoke verifier reruns the configuration at t_ox = 2 nm and t_ox =
10 nm via in-process re-runs that override the mesh extents,
resolution, oxide region bounds, and gate facet position. The intended
use case is "how does my MOSCAP C-V shift when the oxide thickness
changes," and the config is meant to be cloned and re-parametrized for
your own gate stack.

This is an example, not a V&V gate. For the analytical depletion-
approximation C-V correctness reference, see
[`benchmarks/mos_2d`](../../benchmarks/mos_2d) (which uses the
`mos_cap_ac` runner with analytic PDE sensitivity for noise-free
differential capacitance) and the axisymmetric Hu Fig 5-18 reference
in [`benchmarks/moscap_axisym_2d`](../../benchmarks/moscap_axisym_2d).

This is the first example in the catalogue to exercise the `mos_cv`
runner. The first three examples (PR #85) are all `bias_sweep` or
`transient`. The fourth (`pmos_idvgs`) is also `bias_sweep`. C-V vs
I-V is the new output mode here.

## Physics features exercised

- **`mos_cv` runner.** `solver.type = "mos_cv"` solves multi-region
  equilibrium Poisson at each gate bias and integrates the
  semiconductor space charge to recover Q_gate(V_gate). The verifier
  post-processes the iv table to C(V_gate) = dQ/dV via
  `numpy.gradient`. The runner is documented in
  `semi/runners/mos_cv.py`. See also `mos_cap_ac` (the M14.1 successor
  with analytic PDE sensitivity) which is what
  [`benchmarks/mos_2d`](../../benchmarks/mos_2d) and
  [`benchmarks/moscap_axisym_2d`](../../benchmarks/moscap_axisym_2d)
  use for V&V; this example deliberately uses the legacy `mos_cv`
  path so the catalogue covers both runners.
- **Multi-region mesh.** `regions_by_box` partitions the mesh into a
  silicon region (semiconductor) and an oxide region (insulator); the
  multi-region equilibrium Poisson form is built with a cell-tag
  measure that restricts the space-charge term to the silicon region
  only.
- **Gate contact with workfunction.** `contacts[1].type = "gate"` with
  `workfunction = 0` (mid-gap reference). The full Dirichlet BC at
  the gate facet is psi_metal = V_gate - workfunction in scaled units.
- **No mobility, no recombination.** The C-V output is a pure
  electrostatic measurement at equilibrium; the mobility and
  recombination blocks are present in the JSON for schema
  completeness but the `mos_cv` runner does not consume them.

## How to run

From the repository root:

```
docker compose run --rm benchmark moscap_cv_oxide_thickness
```

(The CLI in `scripts/run_benchmark.py` falls back to `examples/` when
the requested name is not found under `benchmarks/`, so the same
invocation works for both directories.)

Output lands in `results/moscap_cv_oxide_thickness/`:

- `cv_oxide_thickness.png`: three C-V curves overlaid, color-coded by
  oxide thickness, on a linear-y axis.
- The result JSON written by the engine for the anchor t_ox = 5 nm
  run is in the standard `results/moscap_cv_oxide_thickness/`
  location; the companion runs land in `_t_ox_<thickness>nm/`
  subdirectories.

## Expected output

The accumulation capacitance per unit area (the strong-accumulation
plateau, V_gate < V_FB) is set by the oxide alone:

    C_ox = eps_0 * eps_SiO2 / t_ox

For SiO2 (eps_SiO2 = 3.9):

    t_ox =  2 nm:  C_ox ~ 1.73e-2 F/m^2
    t_ox =  5 nm:  C_ox ~ 6.91e-3 F/m^2
    t_ox = 10 nm:  C_ox ~ 3.45e-3 F/m^2

Doubling t_ox from 2 nm to 4 nm halves C_ox; the user can verify this
1/t_ox scaling directly from the strong-accumulation plateaus of the
three curves. In strong accumulation (deep negative V_gate on a p-type
body), the silicon surface is essentially metallic (hole accumulation
layer) and the total capacitance equals C_ox.

In weak inversion / depletion, the silicon depletion-region
capacitance C_dep = eps_Si / W_dep adds in series with C_ox so the
total capacitance drops below C_ox; the depth of the C-V dip is set
by the body doping (heavier doping shrinks W_dep,max and gives a
shallower dip).

If the reported strong-accumulation capacitance ordering is wrong (for
example C(2 nm) is not the largest), check:

- The oxide region bounds (`mesh.regions_by_box[1].bounds`) actually
  cover the right vertical interval. Check by-eye against the
  `extents` and the gate facet position.
- The gate facet position
  (`mesh.facets_by_plane[1].value`). It must equal the top of the
  oxide region.
- The oxide material (`regions.oxide.material = "SiO2"`). The
  framework's SiO2 dielectric constant is 3.9; if your override
  changes the material, the C_ox scaling changes correspondingly.
- The `physics.statistics` setting. For body doping <= 5e16 cm^-3,
  Boltzmann statistics is appropriate; FD makes essentially no
  difference at this doping level but may shift inversion-side
  capacitance by a fraction of a percent at very deep inversion.

The mid-gap (V_T extraction) inflection point of the C-V curve sits
at V_gate ~= V_FB + 2 phi_F where phi_F = kT/q * ln(N_A / n_i). For
N_A = 5e16 cm^-3 at 300 K, 2 phi_F ~ 0.78 V; the inflection is
visible as the steepest slope of the depletion-to-inversion
transition.

For the formal V&V version of the C-V analysis with the depletion-
approximation analytical reference, see
[`benchmarks/mos_2d`](../../benchmarks/mos_2d) (Cartesian) or
[`benchmarks/moscap_axisym_2d`](../../benchmarks/moscap_axisym_2d)
(axisymmetric, with the Hu Fig. 5-18 reference).

## How to adapt

Most users will want to start by changing one of:

- **Oxide thickness.** Edit `mesh.extents[1][1]`,
  `mesh.resolution[1]` (so the cell size stays 1 nm), the oxide
  region's `bounds[1]`, and the gate facet's `value`. Heavier oxide
  scales C_ox down by 1/t_ox.
- **Body doping.** Change `doping[0].profile.N_A`. Heavier body
  doping shifts V_T up, deepens the C-V minimum, and shrinks the
  depletion-region width.
- **Gate workfunction.** `contacts[1].workfunction`. A
  polysilicon-on-Si gate adds a workfunction-difference offset that
  shifts the C-V curve laterally without changing its shape; metal
  gates and high-k stacks add their own offsets.
- **Sweep range.** `contacts[1].voltage_sweep`. The default
  [-2, 2] V at 0.1 V step covers accumulation through inversion for
  N_A = 5e16 cm^-3; widen the range or refine the step for sharper
  inflection-point extraction.
- **Body type.** Swap `N_D` and `N_A` to flip the device to an
  n-body MOSCAP (same C_ox scaling, opposite-sign V_T).

## Choice of parametric mechanism

Oxide thickness is geometry, not a single physics knob, so the
verifier-overrides approach is more invasive than the
`schottky_iv_temperature` `physics.temperature` override. The
`_oxide_thickness_overrides` helper in `scripts/run_benchmark.py`
mutates `mesh.extents[1][1]`, `mesh.resolution[1]`,
`mesh.regions_by_box[1].bounds[1]`, and the gate
`facets_by_plane.value` in deepcopies of the anchor cfg. Cell size
stays at 1 nm so the Si/oxide interface always falls on a grid line.
A purely-physics knob would have been cleaner, but oxide thickness
is the canonical MOSCAP characterization sweep and is worth showing.

## Notes on CI runtime

Three equilibrium-Poisson sweeps (41 V_gate points each, 2008-2040
cells per mesh) on a small Cartesian grid are cheap; estimate
~1 minute total wall-clock on a typical CI runner.
