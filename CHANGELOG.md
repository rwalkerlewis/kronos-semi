# Changelog

## [0.1.0] - Day 1

### Added
- JSON schema with Draft-07 validation (`semi.schema`)
- Material database: Si, Ge, GaAs, SiO₂, HfO₂, Si₃N₄ (`semi.materials`)
- Physical constants module with SI/cm conversions (`semi.constants`)
- Nondimensional scaling class with auto-inference from config (`semi.scaling`)
- Doping profile evaluators: uniform, step, Gaussian (`semi.doping`)
- Builtin mesh generation with region/facet tagging (`semi.mesh`)
- Equilibrium Poisson form under Boltzmann statistics (`semi.physics.poisson`)
- PETSc SNES solver wrapper for dolfinx 0.10+ (`semi.solver`)
- Top-level `run(cfg)` entry point (`semi.run`)
- 1D pn junction benchmark with JSON input (`benchmarks/pn_1d/`)
- Colab notebook: clone repo, import package, run, plot, verify (`notebooks/01_pn_junction_1d.ipynb`)
- Test suite: 36 tests covering materials, schema, scaling, doping, constants
- Offline math sanity check script (runs without dolfinx)
- GitHub Actions CI for pure-Python tests across Python 3.10/3.11/3.12
- Ruff lint configuration

### Verified
- V_bi = 0.833 V for symmetric 10¹⁷/10¹⁷ Si junction (matches analytical)
- Peak |E| = 113.5 kV/cm (within few % of depletion approximation)
- Mass-action np = n_i² in bulk (to numerical precision)
- Charge neutrality in quasi-neutral regions (to 10⁻¹⁵)

### Not yet implemented (see Roadmap in README)
- Drift-diffusion under applied bias
- Coupled (ψ, Φₙ, Φₚ) block Newton
- Bias sweep continuation
- Multi-region (oxide + semiconductor)
- 2D MOS capacitor benchmark
- 3D resistor benchmark
- MOSFET and FinFET
