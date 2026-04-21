# Contributing

## Before you start

Read these documents before writing code. They are the source of truth
for what is done, what is next, and which design decisions are locked.

- [PLAN.md](PLAN.md) for current state, the single in-flight task, and invariants.
- [docs/PHYSICS.md](docs/PHYSICS.md) for equations, scaling conventions, and boundary conditions.
- [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) for the layered design and the imports-allowed rules per layer.
- [docs/adr/](docs/adr/) for accepted architecture decisions. If your change conflicts with any invariant in PLAN.md, open an ADR first.

This is an early-stage evaluation project. If you want to run it or extend it:

## Development setup

```bash
git clone https://github.com/rwalkerlewis/kronos-semi.git
cd kronos-semi
# Option A: conda for dolfinx (full FEM capability)
conda create -n kronos-semi -c conda-forge python=3.12 fenics-dolfinx mpich
conda activate kronos-semi
pip install -e ".[dev]"

# Option B: pure-Python only (no FEM, just schema/materials/scaling/doping)
pip install -e ".[dev]"
```

## Before committing

```bash
pytest tests/ -v       # all 36 should pass
ruff check semi/ tests/
python tests/check_day1_math.py
```

## Structure

- `semi/` — the package. Pure-Python modules (constants, materials, scaling, doping, schema) never import dolfinx. FEM modules (mesh, physics, solver, run) import dolfinx lazily.
- `benchmarks/` — self-contained JSON inputs plus README per benchmark.
- `notebooks/` — thin Colab drivers that clone + import the package.
- `tests/` — pytest suite for pure-Python modules + offline math checks.
- `scripts/` — notebook generators and dev utilities (not installed).

## Conventions

- Internal units are SI. JSON input follows device-physics tradition where doping is in cm⁻³ and mobility in cm²/(V·s), converted at ingest.
- Physics-style variable names (`N_A`, `V`, `F`, `psi_L`, `eps_Si`) are fine; PEP 8 naming is explicitly relaxed via ruff config.
- Any new physics module should include (a) a unit test in `tests/` if it has pure-Python logic, (b) an analytical verification in the corresponding benchmark.
