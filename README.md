# kronos-semi

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/rwalkerlewis/kronos-semi/blob/main/notebooks/01_pn_junction_1d.ipynb)

A JSON-driven finite-element semiconductor device simulator built on [FEniCSx](https://fenicsproject.org/) (dolfinx 0.10+), targeting the workflow of the COMSOL Semiconductor Module — Poisson-coupled drift-diffusion with SRH recombination, multi-region support, ohmic/gate/insulating boundary conditions — in a clean, extensible Python package.

## Planning documents

Before contributing (human or AI), read these in order:

- [PLAN.md](PLAN.md): authoritative current state, next task, invariants, and non-goals.
- [docs/PHYSICS.md](docs/PHYSICS.md): governing equations, nondimensionalization, and boundary-condition conventions.
- [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md): five-layer component design and dependency rules.
- [docs/adr/](docs/adr/): architecture decision records; open a new ADR before changing any invariant.
- [docs/ROADMAP.md](docs/ROADMAP.md): full per-day deliverables and verification criteria.

## Status

**End-of-M7 capability matrix**, shipped across PRs #2–#9:

| Capability                         | Dimensions    | Status  | Verifier                                                     |
|------------------------------------|---------------|---------|--------------------------------------------------------------|
| Equilibrium Poisson                | 1D / 2D / 3D  | shipped | MMS finest-pair rates L2 = 2.0, H1 = 1.0                     |
| Coupled Slotboom drift-diffusion   | 1D / 2D       | shipped | MMS finest-pair rates L2 ≥ 1.99 across variants              |
| SRH recombination                  | 1D / 2D       | shipped | verified against SNS analytical at reverse bias              |
| Ohmic contact BCs                  | 1D / 2D / 3D  | shipped | Shockley diode within 10% at forward bias                    |
| Gate contact BCs with φ_ms         | 2D            | shipped | MOS C–V within 10% in depletion window                       |
| Multi-region Poisson (Si/SiO₂)     | 2D            | shipped | multi-region MMS L2 = 2.0                                    |
| File-sourced gmsh .msh meshes      | 3D            | shipped | builtin vs gmsh R-match within 1%                            |
| Adaptive bias continuation         | uni + bipolar | shipped | pn junction forward + reverse, 3D resistor                   |
| 3D ohmic V–I linearity             | 3D            | shipped | V–I linearity within 1%                                      |
| Benchmarks                         | 5             | shipped | pn_1d, pn_1d_bias, pn_1d_bias_reverse, mos_2d, resistor_3d   |
| Conservation / mesh convergence    | 1D            | shipped | charge neutrality, Cauchy rates ≥ 1.8/doubling               |
| Test suite                         | pure + FEM    | shipped | 206 tests, 95.58% coverage                                   |
| V&V                                | 10 studies    | shipped | 62/62 PASS                                                   |
| CI                                 | lint+test+FEM | shipped | green on dev and main                                        |

### Scope

kronos-semi is an **eight-milestone FEM-framework demonstration** targeting the quasi-static, steady-state subset of the COMSOL Semiconductor Module: equilibrium Poisson and coupled drift-diffusion with SRH recombination over 1D/2D/3D meshes, multi-region dielectrics, and ohmic/gate/insulating boundary conditions. The verifier suite proves the numerics against analytical (Shockley, SNS, depletion-approx C–V, ohmic V–I) and manufactured (MMS) solutions.

What it deliberately does **not** implement, relative to a full commercial Semiconductor Module:

- Transient (time-dependent) solver — steady-state only.
- Impact ionization and avalanche generation.
- Heterojunction band-offset models and position-dependent band structure.
- Trap-assisted and band-to-band tunneling.
- Fermi–Dirac statistics (Boltzmann throughout; valid below ~10¹⁹ cm⁻³).
- Mobility models beyond constant — no Caughey–Thomas, no Lombardi surface mobility.
- Non-equilibrium Green's function (NEGF) or quantum-corrected transport.
- Schottky contacts, tunnel junctions, or contact resistance models.
- Thermal coupling (self-heating, lattice temperature).
- Optical generation and photovoltaic carrier sources.

These are deferred as out-of-scope for the eight-milestone framework demonstration; the architecture accommodates them as follow-on work (see [docs/ROADMAP.md](docs/ROADMAP.md) Post-submission).

For a breakdown of what each day shipped, see the capability matrix above and [CHANGELOG.md](CHANGELOG.md).

## Quick start on Colab (zero setup)

Click the badge above. The first cell installs dolfinx on Colab via [FEM on Colab](https://fem-on-colab.github.io/) (~30 s); subsequent cells clone this repo, load a JSON, run the solve, and plot results against analytical curves.

If you prefer a direct link:
```
https://colab.research.google.com/github/rwalkerlewis/kronos-semi/blob/main/notebooks/01_pn_junction_1d.ipynb
```

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

## Running the benchmark

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

### Why Slotboom variables (planned for M2+)

Naive Galerkin FEM on drift-diffusion is unstable when drift dominates diffusion, which is almost everywhere in a real device. The standard FEM remedy is to rewrite the continuity equations in terms of quasi-Fermi potentials (Slotboom variables), which makes the current a pure gradient and the resulting form well-posed without stabilization. The equilibrium case implemented here is a special case where Φₙ = Φₚ = 0 and the continuity equations are trivially satisfied.

### Why nondimensional scaling

A raw 1 µm device at 10¹⁷ cm⁻³ doping has permittivity ~10⁻¹¹, charge ~10⁻¹⁹, density ~10²³. Newton on that directly diverges. All fields are scaled so ψ̂ is O(1) and carrier ratios n/C₀ are O(1). The scaled Poisson equation has a small parameter λ² = ε V_t / (q C₀ L₀²) ~ 10⁻⁴, which is the squared Debye-length-to-device ratio.

### Why dolfinx 0.10

The `NonlinearProblem` class in 0.10 wraps PETSc SNES directly (the old `NewtonSolver` is deprecated), which gives us line search, sophisticated convergence criteria, and block-system support for free when we add the coupled (ψ, Φₙ, Φₚ) system.

## Roadmap

- **M1** ✓ Equilibrium Poisson, 1D pn junction, depletion-approx verification
- **M2** ✓ Slotboom drift-diffusion, coupled (ψ, Φₙ, Φₚ) Newton
- **M3** ✓ Adaptive continuation; forward-bias IV curve vs Shockley diode equation
- **M4** ✓ V&V suite; MMS, mesh convergence, conservation, CI
- **M5** ✓ Refactor and test pass (run.py split, bcs.py extracted, coverage 96%)
- **M6** ✓ 2D MOS capacitor (oxide + silicon multi-region, submesh for carriers)
- **M7** ✓ 3D doped resistor (gmsh loader, bipolar sweep, V-I linearity)
- **Post-submission** Full 2D MOSFET with source/drain/gate/body contacts; 3D FinFET; field-dependent mobility; Auger recombination; Schottky contacts

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
