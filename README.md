# kronos-semi

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/rwalkerlewis/kronos-semi/blob/main/notebooks/01_pn_junction_1d.ipynb)

kronos-semi is a FEniCSx-based finite-element semiconductor device simulator that mimics the capabilities of the COMSOL Semiconductor Module. Simulations are driven by a single JSON file (optionally referencing an external geometry/mesh artifact) and are otherwise plain text. The project ships 1D, 2D, and 3D benchmark problems and a zero-setup Colab notebook (cloud-hosted, no local install) so reviewers can run everything from a browser.

**Audience:** TCAD engineers, device-physics researchers, and evaluation reviewers who want a reproducible, open-source alternative to a commercial Semiconductor Module solver.

## Quick start on Colab (zero setup, under one minute)

Click the badge above, or open the direct link:

```
https://colab.research.google.com/github/rwalkerlewis/kronos-semi/blob/main/notebooks/01_pn_junction_1d.ipynb
```

The first cell installs FEniCSx on Colab via [FEM on Colab](https://fem-on-colab.github.io/) (~30 s). Subsequent cells clone this repo, load a JSON benchmark file, run the solve, and plot results against analytical curves. No local installation is required.

### Notebooks

| Notebook | Benchmark | Colab |
|---|---|---|
| 01 pn junction 1D | Equilibrium Poisson, depletion approximation | [![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/rwalkerlewis/kronos-semi/blob/main/notebooks/01_pn_junction_1d.ipynb) |
| 02 pn junction bias | Forward Shockley + reverse SNS sweep | [![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/rwalkerlewis/kronos-semi/blob/main/notebooks/02_pn_junction_bias.ipynb) |
| 03 MOS C-V | 2D MOS capacitor C-V sweep | [![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/rwalkerlewis/kronos-semi/blob/main/notebooks/03_mos_cv.ipynb) |
| 04 resistor 3D | 3D doped bar, builtin and gmsh meshes | [![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/rwalkerlewis/kronos-semi/blob/main/notebooks/04_resistor_3d.ipynb) |

## Scope vs. COMSOL Semiconductor Module

kronos-semi covers the quasi-static, steady-state subset of the COMSOL Semiconductor Module over 1D, 2D, and 3D meshes.

**In scope (shipped):**
- Poisson equation with multi-region dielectric (Si/SiO2)
- Drift-diffusion in Slotboom (quasi-Fermi potential) form
- SRH recombination with configurable mid-gap trap energy
- Ohmic contacts and ideal gate contacts with work-function offset
- 1D, 2D, and 3D structured meshes (builtin); 3D unstructured meshes via gmsh .msh files
- Bias sweeps with adaptive step-size continuation (unipolar and bipolar)
- Method-of-Manufactured-Solutions and conservation V&V suite

**Explicitly out of scope (post-submission stretch goals, see [docs/ROADMAP.md](docs/ROADMAP.md)):**
- Caughey-Thomas or Lombardi field-dependent mobility
- Auger and radiative recombination
- Fermi-Dirac statistics (Boltzmann throughout; valid below ~10^19 cm^-3)
- Impact ionization and avalanche generation
- AC small-signal analysis
- Transient (time-dependent) solver (steady-state only)
- Band-to-band or trap-assisted tunneling
- Heterojunction band-offset models and position-dependent band structure
- Schottky contacts, tunnel junctions, or contact resistance models
- Thermal coupling (self-heating, lattice temperature)
- Optical generation and photovoltaic carrier sources
- Full MOSFET, FinFET, or 3D transistor geometries

## Status

**End-of-M7 capability matrix**, shipped across PRs #2-#9:

| Capability                         | Dimensions    | Status  | Verifier                                                     |
|------------------------------------|---------------|---------|--------------------------------------------------------------|
| Equilibrium Poisson                | 1D / 2D / 3D  | shipped | MMS finest-pair rates L2 = 2.0, H1 = 1.0                     |
| Coupled Slotboom drift-diffusion   | 1D / 2D       | shipped | MMS finest-pair rates L2 >= 1.99 across variants             |
| SRH recombination                  | 1D / 2D       | shipped | verified against SNS analytical at reverse bias              |
| Ohmic contact BCs                  | 1D / 2D / 3D  | shipped | Shockley diode within 10% at forward bias                    |
| Gate contact BCs with phi_ms       | 2D            | shipped | MOS C-V within 10% in depletion window                       |
| Multi-region Poisson (Si/SiO2)     | 2D            | shipped | multi-region MMS L2 = 2.0                                    |
| File-sourced gmsh .msh meshes      | 3D            | shipped | builtin vs gmsh R-match within 1%                            |
| Adaptive bias continuation         | uni + bipolar | shipped | pn junction forward + reverse, 3D resistor                   |
| 3D ohmic V-I linearity             | 3D            | shipped | V-I linearity within 1%                                      |
| Benchmarks                         | 5             | shipped | pn_1d, pn_1d_bias, pn_1d_bias_reverse, mos_2d, resistor_3d   |
| Conservation / mesh convergence    | 1D            | shipped | charge neutrality, Cauchy rates >= 1.8/doubling              |
| Test suite                         | pure + FEM    | shipped | 206 tests, 95.58% coverage                                   |
| V&V                                | 10 studies    | shipped | 62/62 PASS                                                   |
| CI                                 | lint+test+FEM | shipped | green on dev and main                                        |

For a breakdown of what each milestone shipped, see the capability matrix above and [CHANGELOG.md](CHANGELOG.md).

## Local install

dolfinx is most easily installed via conda:

```bash
conda create -n kronos-semi -c conda-forge python=3.12 fenics-dolfinx mpich pyvista
conda activate kronos-semi
git clone https://github.com/rwalkerlewis/kronos-semi.git
cd kronos-semi
pip install -e ".[dev]"
```

The pure-Python modules (`schema`, `materials`, `scaling`, `doping`, `constants`) don't need dolfinx and can be installed standalone via `pip install -e .` into any environment; the FEM-heavy modules (`mesh`, `physics.poisson`, `solver`, `run`) need dolfinx available at import time.

## Planning documents

Before contributing (human or AI), read these in order:

- [PLAN.md](PLAN.md): authoritative current state, next task, invariants, and non-goals.
- [docs/PHYSICS.md](docs/PHYSICS.md): governing equations, nondimensionalization, and boundary-condition conventions.
- [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md): five-layer component design and dependency rules.
- [docs/adr/](docs/adr/): architecture decision records; open a new ADR before changing any invariant.
- [docs/ROADMAP.md](docs/ROADMAP.md): capability matrix, scope vs. COMSOL Semiconductor Module, and delivery history.

## Docker

A reproducible dev environment is provided on top of the official `ghcr.io/fenics/dolfinx/dolfinx:stable` image (dolfinx 0.10). The source tree is bind-mounted at `/workspaces/kronos-semi` and the package is installed editable, so host edits take effect immediately.

```bash
# Build the image once (matches your host UID/GID so bind-mount writes stay owned by you)
docker compose build

# Run the test suite
docker compose run --rm test

# Run a benchmark (saves plots to results/<name>/)
docker compose run --rm benchmark pn_1d

# Interactive shell in a long-running dev container
docker compose up -d dev
docker compose exec dev bash

# JupyterLab on http://localhost:8888
docker compose up jupyter
```

## Running benchmarks

The primary interface is the benchmark CLI. Pass a benchmark name to run the full solve, verify, and plot pipeline:

```bash
docker compose run --rm benchmark pn_1d
docker compose run --rm benchmark pn_1d_bias
docker compose run --rm benchmark pn_1d_bias_reverse
docker compose run --rm benchmark mos_2d
docker compose run --rm benchmark resistor_3d
```

Each benchmark is defined by a JSON file under `benchmarks/<name>/`. The CLI loads the JSON, runs the solver, checks verifier criteria, and writes plots to `results/<name>/`.

To call the solver from Python (for scripting or notebook use), load the JSON with `semi.schema.load` and pass the resulting config to `semi.run.run`:

```python
from semi import schema, run

cfg = schema.load("benchmarks/pn_1d/pn_junction.json")
result = run.run(cfg)

# result.psi_phys, result.n_phys, result.p_phys are numpy arrays
# result.x_dof gives the mesh coordinates
# result.scaling exposes the nondimensional scales used
```

## JSON input schema

The schema is defined in `semi/schema.py`. Minimal 1D pn junction example:

```json
{
  "name": "pn_junction_1d",
  "dimension": 1,
  "mesh": {
    "source": "builtin",
    "extents": [[0.0, 2.0e-6]],
    "resolution": [400],
    "facets_by_plane": [
      {"name": "anode",   "tag": 1, "axis": 0, "value": 0.0},
      {"name": "cathode", "tag": 2, "axis": 0, "value": 2.0e-6}
    ]
  },
  "regions": {
    "silicon": {"material": "Si", "tag": 1, "role": "semiconductor"}
  },
  "doping": [
    {
      "region": "silicon",
      "profile": {
        "type": "step", "axis": 0, "location": 1.0e-6,
        "N_D_left": 0.0,  "N_A_left": 1.0e17,
        "N_D_right": 1.0e17, "N_A_right": 0.0
      }
    }
  ],
  "contacts": [
    {"name": "anode",   "facet": "anode",   "type": "ohmic", "voltage": 0.0},
    {"name": "cathode", "facet": "cathode", "type": "ohmic", "voltage": 0.0}
  ]
}
```

Doping densities are specified in cm⁻³ (device-physics tradition); everything else is SI.

## Design notes

### Why Slotboom variables

Naive Galerkin FEM on drift-diffusion is unstable when drift dominates diffusion, which is almost everywhere in a real device. The standard FEM remedy is to rewrite the continuity equations in terms of quasi-Fermi potentials (Slotboom variables), which makes the current a pure gradient and the resulting form well-posed without stabilization. The Slotboom form is implemented in `semi/physics/drift_diffusion.py` and verified by the MMS-DD suite (62/62 PASS). See `docs/adr/0004-slotboom-variables-for-dd.md` for the design rationale.

### Why nondimensional scaling

A raw 1 µm device at 10¹⁷ cm⁻³ doping has permittivity ~10⁻¹¹, charge ~10⁻¹⁹, density ~10²³. Newton on that directly diverges. All fields are scaled so ψ̂ is O(1) and carrier ratios n/C₀ are O(1). The scaled Poisson equation has a small parameter λ² = ε V_t / (q C₀ L₀²) ~ 10⁻⁴, which is the squared Debye-length-to-device ratio.

### Why dolfinx 0.10

The `NonlinearProblem` class in 0.10 wraps PETSc SNES directly (the old `NewtonSolver` is deprecated), which gives us line search, sophisticated convergence criteria, and block-system support for free when we add the coupled (ψ, Φₙ, Φₚ) system.

## Roadmap

- **M1** ✓ Equilibrium Poisson, 1D pn junction, Docker env
- **M2** ✓ Slotboom drift-diffusion, coupled (psi, Phi_n, Phi_p) Newton
- **M3** ✓ Adaptive continuation; forward-bias IV curve vs Shockley diode equation
- **M4** ✓ V&V suite; MMS, mesh convergence, conservation, CI
- **M5** ✓ Refactor and test pass (run.py split, bcs.py extracted, coverage 96%)
- **M6** ✓ 2D MOS capacitor (oxide + silicon multi-region, submesh for carriers)
- **M7** ✓ 3D doped resistor (gmsh loader, bipolar sweep, V-I linearity)
- **M8** (in flight) Submission polish: notebooks, catalog, CHANGELOG

Post-submission stretch goals (field-dependent mobility, Auger, Fermi-Dirac, transient, full MOSFET) are documented in [docs/ROADMAP.md](docs/ROADMAP.md).

## Verification

Every analytical result is regenerated at import time by `tests/check_day1_math.py`, which runs without dolfinx. 9/9 checks pass:

- Thermal voltage V_t = 25.852 mV at 300 K
- Debye length for 10¹⁷ Si = 12.93 nm
- V_bi formulas (asinh vs ln) agree to 10⁻⁴
- Mass-action np = n_i² in bulk
- Charge neutrality in bulk satisfied to 10⁻¹⁵
- Peak |E| for symmetric 10¹⁷/10¹⁷ Si junction = 113.5 kV/cm

## License

MIT. See [LICENSE](LICENSE).

## Author

Robert "Bud" Walker · [Rainbow Seismic](https://github.com/rwalkerlewis) · computational geophysics and device simulation consulting.

Built as an evaluation task for KronosAI; extended from prior finite-element work on [PyLith](https://github.com/geodynamics/pylith) (crustal deformation, poroelastic coupling).
