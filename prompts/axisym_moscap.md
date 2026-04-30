# Axisymmetric 2D MOSCAP with LF vs. HF C–V

> Worker prompt for kronos-semi. Implements a 2D axisymmetric (r, z) MOS
> capacitor benchmark and reproduces Fig. 5-18 of Chenming Hu,
> *Modern Semiconductor Devices for Integrated Circuits*, Ch. 5
> (https://www.chu.berkeley.edu/wp-content/uploads/2020/01/Chenming-Hu_ch5-1.pdf).

---

## 1. Goal

Extend the existing **2D MOSCAP benchmark** in this repo to a **2D axisymmetric (r, z)**
model of a circular MOS capacitor and reproduce the qualitative difference between the
**low-frequency (quasi-static) C–V** and the **high-frequency C–V** as shown in
Hu Fig. 5-18.

In Fig. 5-18:

- **Upper curve (LF / QS C–V):** in inversion, C returns to `Cox` because the inversion
  layer can follow the AC signal.
- **Lower curve (HF C–V):** in inversion, C saturates at `Cox · Cdep,min / (Cox + Cdep,min)`
  because the minority carriers cannot follow the AC signal; only the depletion-edge
  majority charge responds.

The deliverable must reproduce both curves on a single plot, on the same device, with
the LF curve clearly returning to `Cox` and the HF curve saturating at the
Cdep,min-limited plateau.

---

## 2. Scope and non-goals

**In scope**

- New benchmark/example: axisymmetric 2D MOSCAP (oxide + silicon, circular gate over
  a cylindrical body).
- A `dimension: "2D-axisymmetric"` (or equivalent) option threaded through
  `semi.schema`, `semi.mesh`, `semi.physics.poisson`, and `semi.solver`.
- A C–V driver that produces both LF and HF curves from a single bias sweep.
- Geometry plot, mesh plot, region-tag plot, and field plots (ψ, n, p) at
  representative biases (accumulation, flat-band, depletion, threshold, inversion).
  Use `pyvista` for fields/mesh and `matplotlib` for the C–V curve.
- A Jupyter notebook that runs end-to-end on Colab (matching the existing
  `notebooks/01_pn_junction_1d.ipynb` style and the FEM-on-Colab dolfinx install pattern).

**Out of scope**

- Drift-diffusion under bias (Day 2+ on the roadmap). Use the equilibrium
  Poisson–Boltzmann form already in `semi.physics.poisson`.
- Quantum corrections to inversion-layer thickness (Hu §5.9). Stay classical.
- Poly-Si gate depletion (Hu §5.8). Use an ideal metal gate (Dirichlet on the gate facet).
- Interface trap / oxide-charge effects. Ideal interface only.
- 3D.

---

## 3. Where to start in the existing code

Before writing anything new, read and follow the conventions in:

- `semi/schema.py` — JSON schema and `load()`.
- `semi/mesh.py` — builtin mesh generator with axis-aligned region/facet tagging.
  Extend, don't rewrite.
- `semi/physics/poisson.py` — equilibrium Poisson under Boltzmann. The integration
  measure is what changes for axisymmetry; the residual structure stays the same.
- `semi/scaling.py` — nondimensional scaling. Reuse as-is. The axisymmetric weight is
  a geometric factor and does not change the scales.
- `semi/run.py` — top-level entry. Add an axisymmetric path; do not branch the planar path.
- `benchmarks/mos_2d/` — the existing planar 2D MOSCAP benchmark. Mirror its layout
  for the new axisymmetric benchmark.
- `notebooks/01_pn_junction_1d.ipynb` — copy its Colab-install header verbatim for the
  new notebook.

If the existing planar 2D MOSCAP is incomplete, fix what is needed for it to serve as a
sanity check that the new axisymmetric weight reduces correctly, but do not rewrite it.

---

## 4. Physics specification

### 4.1 Geometry (axisymmetric)

The (r, z) computational domain represents a half-cross-section of a cylindrical
MOSCAP rotated about the z-axis (r = 0):

```
        r = 0                              r = R_dev
        |                                  |
  z = 0 +----------- gate contact ---------+   <- top (Dirichlet: V_g)
        |          oxide (region 1)        |
 z = tox+----------------------------------+   <- Si/SiO2 interface (internal)
        |                                  |
        |          silicon (region 2)      |
        |                                  |
z = -Hsi+----------- body contact ---------+   <- bottom (Dirichlet: 0 V)
```

- `r ∈ [0, R_dev]`, `z ∈ [-H_si, t_ox]` (sign convention: oxide above z=0, silicon
  below; pick whatever matches the rest of the codebase and document it).
- Symmetry axis at `r = 0` — natural Neumann (no flux). With the r-weighted measure
  this is automatic; no explicit BC needed.
- Outer radial boundary `r = R_dev` — natural Neumann (sufficiently far from the gate
  edge that fringing is negligible). Pick `R_dev ≥ 5 · max(W_dmax, t_ox)` and verify by
  mesh-independence.
- Top facet `z = t_ox`, `r ∈ [0, R_gate]`: gate Dirichlet
  `ψ = V_g - Vfb_correction` (or `V_g` if `Vfb` is folded into the contact spec —
  match the existing schema).
- Top facet `z = t_ox`, `r ∈ [R_gate, R_dev]`: Neumann (free surface above oxide). For
  a first cut, use `R_gate = R_dev` (fully-covered gate, no fringing) so the curve
  matches 1D theory; add `R_gate < R_dev` as a follow-up case to show fringing.
- Bottom facet `z = -H_si`: body ohmic Dirichlet `ψ = 0`.

**Default device** (matches a textbook Hu Fig. 5-8 / Fig. 5-18 case):

- P-type Si body, `N_A = 1e17 cm⁻³`.
- SiO₂ gate dielectric, `t_ox = 5 nm` (thick enough that poly-depletion / quantum
  corrections are negligible — important since we explicitly excluded those).
- `H_si = 1 µm` (>> `W_dmax ≈ 0.1 µm` for `N_A = 1e17`).
- `R_gate = R_dev = 1 µm` (axisymmetric "1D-like" baseline).
- Ideal metal gate, work function chosen so that `V_t > 0` (e.g. emulate N⁺ poly with
  `φ_M ≈ 4.05 V`).
- T = 300 K.

### 4.2 Axisymmetric weak form

Cylindrical Laplacian with azimuthal symmetry:
```
∇·(ε ∇ψ)  =  (1/r) ∂/∂r ( r ε ∂ψ/∂r ) + ∂/∂z ( ε ∂ψ/∂z )  =  -ρ(ψ)
```

Multiplying by a test function `v` and a volume element `dV = 2π r dr dz`, the **2π
factor cancels** in the residual (or is folded into a per-area normalization at the
end). The weak form becomes — relative to the existing planar `poisson.py`:

```
F(ψ; v)  =  ∫_Ω  ε (∇ψ · ∇v) r dr dz   -   ∫_Ω  ρ(ψ) v r dr dz
                  + (Neumann boundary terms, all zero here)
```

**Implementation:** the cleanest dolfinx pattern is to keep the existing UFL form and
just multiply both volume integrands by the radial coordinate. Use
`ufl.SpatialCoordinate(domain)[0]` to get `r` and absorb it into the measure or into
the integrand. Either:

```python
x = ufl.SpatialCoordinate(domain)
r = x[0]
F = eps * ufl.dot(ufl.grad(psi), ufl.grad(v)) * r * dx \
  - rho(psi) * v * r * dx
```

or define a custom measure once and reuse it. Pick whichever is more idiomatic for the
existing `poisson.py` and document it.

The symmetry-axis BC at `r = 0` is automatic: the `r` factor kills any boundary
integral there, so no explicit Dirichlet/Neumann is needed. **Add an assertion that the
mesh has at least one node on `r = 0`** so this invariant is enforced.

### 4.3 Nondimensionalization

Reuse the existing scaling unchanged. The radial weight `r` is dimensional but the
scaling system already handles a length scale `L0`; just confirm by scaling-test that
the nondimensional residual is well-conditioned. Add a one-line test in
`tests/test_axisym_scaling.py`.

### 4.4 Charge model

Use the equilibrium Poisson–Boltzmann form already in `semi.physics.poisson`:
```
ρ(ψ) = q ( p(ψ) - n(ψ) + N_D⁺ - N_A⁻ )
n(ψ) = n_i exp( (ψ - φ_n) / V_t ),    φ_n = 0  (equilibrium)
p(ψ) = n_i exp( (φ_p - ψ) / V_t ),    φ_p = 0  (equilibrium)
```
No carrier continuity equations are solved — this is fine for both LF and HF C–V
because LF is the equilibrium DC response and HF is a small-signal frequency-domain
modification of the same DC state (see §4.6).

### 4.5 Bias sweep and capacitance computation

Sweep `V_g` from `V_g,min` (deep accumulation, e.g. -2 V relative to flat-band) through
`V_g,max` (strong inversion, e.g. +2 V beyond `V_t`) in steps small enough to resolve
the depletion-region shape change (suggested: 0.05 V; 0.02 V near `V_t`).

At each bias compute the **gate charge per unit gate area**:
```
Q_g(V_g) = (1 / A_gate) · ∫_{Ω_Si} ρ(ψ) · 2π r dr dz       (charge balance)
```
(Or compute as the surface charge implied by `D · n` integrated over the gate — both
are equivalent at convergence; do whichever is numerically cleaner. Reuse helpers if
they exist in `semi.solver`.)

### 4.6 LF vs HF capacitance — the key physics

The two curves in Fig. 5-18 differ in **whether the inversion-layer minority carriers
respond to the AC signal**. In an equilibrium Poisson–Boltzmann simulator, both curves
come from the **same** DC solution at each `V_g`; they differ only in how the
small-signal capacitance is extracted.

**LF (quasi-static) capacitance** — minority carriers respond instantaneously:
```
C_LF(V_g) = - dQ_sub / dV_g          (full numerical derivative of the DC sweep)
```
Compute by **central finite differences** on the bias sweep:
```
C_LF[i] ≈ -( Q_sub[i+1] - Q_sub[i-1] ) / ( V_g[i+1] - V_g[i-1] )
```

**HF capacitance** — minority carriers are frozen at their DC value:
At each `V_g`, after the DC solve, do a **single small-signal solve** of the
linearized Poisson equation in which only the majority-carrier charge perturbation
contributes. For a P-body device the HF small-signal residual is:
```
∫ ε ∇(δψ) · ∇v · r dr dz  -  ∫ (∂ρ_maj/∂ψ) · δψ · v · r dr dz  =  0
```
with Dirichlet `δψ = δV_g` on the gate, `δψ = 0` on the body, and where `∂ρ_maj/∂ψ`
includes **only** the majority-carrier (hole) Boltzmann derivative and the dopant term.
Concretely, for a P-body:
```
∂ρ/∂ψ |_HF = -(q/V_t) · p(ψ_DC)              (holes respond)
∂ρ/∂ψ |_LF = -(q/V_t) · ( n(ψ_DC) + p(ψ_DC) )  (both respond)
```
The HF capacitance is then:
```
C_HF(V_g) = -(1/A_gate) · (1/δV_g) · ∫_{Ω_Si} (ρ(ψ_DC + δψ) - ρ(ψ_DC)) · 2π r dr dz
```
or equivalently, if you express the small-signal charge in terms of `δψ` directly via
the linearized form, just integrate the resulting majority-carrier charge perturbation.
Use a small `δV_g` (e.g. 1 mV) for the linearization-based check; the linearized form
should be machine-zero independent of `δV_g`, which is a useful self-test.

**Sanity checks (must pass before declaring success):**

- `C_LF(V_g << V_fb) → C_ox` (accumulation plateau).
- `C_LF(V_g >> V_t) → C_ox` (inversion plateau — this is the "low frequency" feature).
- `C_HF(V_g << V_fb) → C_ox` (same accumulation plateau as LF).
- `C_HF(V_g >> V_t) → C_min = C_ox · C_dep,min / (C_ox + C_dep,min)`, where
  `C_dep,min = ε_Si / W_dmax` and
  `W_dmax = sqrt( 4 · ε_Si · φ_B · ln(N_A/n_i) / (q · N_A) )`. Compute the analytical
  `C_min` and assert agreement to within ~5%.
- `C_LF` and `C_HF` should **coincide** in accumulation, depletion, and at flat-band;
  they only diverge once `V_g > V_t`. The plot should clearly show this.

---

## 5. JSON schema additions

Add to `semi/schema.py` (with full validation and Draft-07 schema):

```jsonc
{
  "name": "moscap_axisym",
  "dimension": "2D-axisymmetric",          // NEW — alongside existing 1, 2
  "mesh": {
    "source": "builtin",
    "extents": [[0.0, 1.0e-6], [-1.0e-6, 5.0e-9]],   // [r-range, z-range]
    "resolution": [80, 200],
    "axisymmetric_axis": 0,                          // NEW — which axis is r
    "regions_by_box": [
      {"name": "oxide",   "tag": 1, "box": [[0, 1e-6], [0, 5e-9]]},
      {"name": "silicon", "tag": 2, "box": [[0, 1e-6], [-1e-6, 0]]}
    ],
    "facets_by_plane": [
      {"name": "gate",   "tag": 1, "axis": 1, "value":  5.0e-9},
      {"name": "body",   "tag": 2, "axis": 1, "value": -1.0e-6},
      {"name": "axis",   "tag": 3, "axis": 0, "value":  0.0  },
      {"name": "outer",  "tag": 4, "axis": 0, "value":  1.0e-6}
    ]
  },
  "regions": {
    "oxide":   {"material": "SiO2", "tag": 1, "role": "insulator"},
    "silicon": {"material": "Si",   "tag": 2, "role": "semiconductor"}
  },
  "doping": [
    {"region": "silicon",
     "profile": {"type": "uniform", "N_A": 1.0e17, "N_D": 0.0}}
  ],
  "contacts": [
    {"name": "gate", "facet": "gate", "type": "gate",  "voltage": 0.0,
     "work_function_eV": 4.05},
    {"name": "body", "facet": "body", "type": "ohmic", "voltage": 0.0}
  ],
  "sweep": {
    "contact": "gate",
    "values": {"start": -2.0, "stop": 2.5, "step": 0.05}
  },
  "cv_analysis": {                          // NEW
    "modes": ["LF", "HF"],
    "delta_V_small_signal": 1.0e-3,
    "majority_carrier": "holes"             // inferred from doping if omitted
  }
}
```

Validate with `jsonschema` on load and add a `tests/test_schema_axisym.py`.

---

## 6. Required file layout

Create:

```
benchmarks/
  moscap_axisym/
    moscap_axisym.json              # input deck described above
    README.md                       # 1-page description with Fig 5-18 reference
    expected/
      cv_curve_LF_HF.png            # reference image for visual diff
      cv_curve_LF_HF.csv            # numerical reference for CI

semi/
  physics/
    poisson.py                      # extended with axisymmetric weak form
    cv.py                           # NEW — LF + HF small-signal capacitance
  mesh.py                           # extended with axisymmetric_axis support
  schema.py                         # extended schema + validation
  run.py                            # axisymmetric branch wired in

notebooks/
  03_moscap_axisym_cv.ipynb         # end-to-end Colab-runnable notebook

tests/
  test_axisym_scaling.py            # condition-number check
  test_schema_axisym.py             # JSON validation
  test_cv_against_analytical.py     # Cox / Cmin / Vt / Vfb tolerances
```

---

## 7. Visualizations (required — do all of these)

Use `pyvista` for 2D field/mesh plots (the env already has it) and `matplotlib` for
line plots. Save every figure as PNG into `benchmarks/moscap_axisym/figs/` and embed
them inline in the notebook.

1. **Geometry sketch (annotated).** A schematic in (r, z) showing: symmetry axis, oxide
   region, silicon region, gate facet, body facet, outer radial boundary, key
   dimensions. Matplotlib is fine — does not need to be from the mesh. Annotate
   `R_gate`, `t_ox`, `H_si`, `R_dev`.
2. **Mesh plot.** The full triangulation in (r, z), with a zoom-in panel near the
   Si/SiO₂ interface (the mesh should be graded there — show that it is).
3. **Region-tag plot.** Same mesh, colored by region tag (oxide vs silicon).
4. **Boundary/facet-tag plot.** Mesh edges colored by facet tag
   (gate / body / axis / outer).
5. **3D revolved-geometry view.** Use `pyvista.RotationalExtrusionFilter` (or rotate
   the 2D mesh manually about r=0) to render the *full 3D cylinder* the axisymmetric
   model represents. Slice it open on a half-plane so the interior is visible. This is
   the single most useful figure for communicating that the model is axisymmetric —
   make it look good.
6. **Field plots (ψ, n, p) at five biases:** strong accumulation, flat-band,
   mid-depletion, threshold, strong inversion. Use a 5×3 grid (rows = biases,
   cols = fields). Use a log-scaled colormap for n and p (they vary over many decades).
7. **C–V curve — the main result.** Single matplotlib figure with both `C_LF / C_ox`
   and `C_HF / C_ox` vs `V_g`, normalized to `C_ox`. Annotate `V_fb`, `V_t`, the
   `C_ox` plateau, and the `C_min` plateau. Style it to look like Hu Fig. 5-18:
   solid line for LF, dashed for HF, both curves on the same axes, regions labeled
   "Accumulation / Depletion / Inversion" along the x-axis.
8. **Convergence figure.** `C_HF` at `V_g = V_t + 1 V` vs mesh refinement (at least 3
   levels). Demonstrates the result is converged.

---

## 8. Acceptance criteria

A reviewer must be able to:

1. `git pull && pip install -e ".[dev]"` then run
   `python -m semi.run benchmarks/moscap_axisym/moscap_axisym.json` and get the C–V
   CSV + PNG in `benchmarks/moscap_axisym/figs/` with no manual steps.
2. Open `notebooks/03_moscap_axisym_cv.ipynb` on Colab via the badge pattern from
   notebook 01 and run all cells top-to-bottom without edits.
3. See, in the final C–V plot, the two curves coinciding from accumulation through
   depletion and **diverging cleanly past `V_t`**, with LF returning to `C_ox` and HF
   saturating at the analytical `C_min`. The plot should be visually comparable to
   Hu Fig. 5-18.
4. `pytest tests/test_cv_against_analytical.py` passes with the tolerances:
   - `|V_fb_simulated - V_fb_analytical| < 20 mV`
   - `|V_t_simulated - V_t_analytical| < 50 mV`
   - `|C_LF(accumulation) / C_ox - 1| < 2%`
   - `|C_LF(strong-inversion) / C_ox - 1| < 5%`
   - `|C_HF(strong-inversion) / C_min_analytical - 1| < 5%`

---

## 9. Implementation hints

- For the small-signal HF solve, the linear system already has the same stiffness
  matrix structure as the Newton step at convergence — reuse the assembled tangent if
  `semi.solver` exposes it. Do **not** re-derive the form from scratch.
- The radial weight `r` is zero on the symmetry axis. This means UFL warnings about
  `1/r` should not appear (good — we never divide by `r`). If they do, you have made
  an algebraic mistake; go back and re-derive the weak form.
- When computing total charge, integrate `ρ * 2π r dr dz` once with a clear constant
  `2π` factor — easier to debug than folding it into the measure.
- For the Vfb calculation, follow Hu Eq. (5.1.1): `V_fb = φ_g - φ_s` with
  `φ_g = work_function_eV`, `φ_s = χ_Si + (E_c - E_F)/q`. Check the sign convention
  against the existing 1D pn junction's contact handling.
- Mesh: graded toward the Si/SiO₂ interface in `z` (smallest element ≈ 0.1 nm at the
  interface, growing geometrically into the bulk). Uniform in `r` for the baseline
  `R_gate = R_dev` case.
- For Colab: the FEM-on-Colab dolfinx install must come before any `import dolfinx`.
  Copy the cell verbatim from `notebooks/01_pn_junction_1d.ipynb`.

---

## 10. References

- Chenming Hu, *Modern Semiconductor Devices for Integrated Circuits*, **Chapter 5** —
  especially §5.6 (MOS C–V), Fig. 5-17 (small-signal charge response in
  accumulation / depletion / inversion), and **Fig. 5-18 (LF vs HF C–V — the target)**.
  PDF: https://www.chu.berkeley.edu/wp-content/uploads/2020/01/Chenming-Hu_ch5-1.pdf
- Sze & Ng, *Physics of Semiconductor Devices*, 3rd ed., §4.3 — alternative derivation
  of the HF C–V minimum.
- Existing kronos-semi conventions: README.md, `semi/schema.py`,
  `semi/physics/poisson.py`, `notebooks/01_pn_junction_1d.ipynb`.

---

## 11. Deliverable summary (what to commit)

A single PR titled `feat: axisymmetric 2D MOSCAP with LF/HF C-V (Hu Fig 5-18)`
containing:

- Schema, mesh, physics, and run extensions for `dimension: "2D-axisymmetric"`.
- New `semi/physics/cv.py` with LF + HF small-signal capacitance.
- New benchmark folder with JSON, README, expected CSV/PNG.
- New Colab-runnable notebook reproducing Fig. 5-18.
- Three new tests (scaling, schema, C–V vs analytical).
- All eight required figures committed and embedded in the notebook.
- README.md updated: roadmap line for "Day 5: 2D MOS capacitor" marked ✓ (with
  axisymmetric variant).
- CHANGELOG.md entry under Unreleased.

Do not regress any existing test. Run `pytest -q` before pushing.
