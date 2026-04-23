# Contributing

## Before you start

Read these documents before writing code. They are the source of truth
for what is done, what is next, and which design decisions are locked.

**Required reading (~45 minutes):**

- [PLAN.md](PLAN.md) — current state, single in-flight task, and invariants.
- [docs/IMPROVEMENT_GUIDE.md](docs/IMPROVEMENT_GUIDE.md) — M9+ roadmap with explicit acceptance tests per milestone. Do not start a milestone without reading its section.
- [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) — five-layer design and per-layer import rules.

**Depending on what you're doing:**

- **New physics or solver code:** [docs/PHYSICS.md](docs/PHYSICS.md) for equations and scaling conventions; then [docs/adr/](docs/adr/) for locked decisions. If your change conflicts with any invariant in `PLAN.md`, open an ADR first.
- **Code navigation and the data flow:** [docs/WALKTHROUGH.md](docs/WALKTHROUGH.md) traces `semi.run.run(cfg)` end-to-end with file:line anchors.
- **No device-physics background:** [docs/PHYSICS_INTRO.md](docs/PHYSICS_INTRO.md) is a tutorial that assumes calculus and programming but not TCAD. Read this before `PHYSICS.md`.
- **Picking up an M9+ milestone:** the matching section of `IMPROVEMENT_GUIDE.md` has your acceptance tests. For M9 specifically, `docs/M9_STARTER_PROMPT.md` is a ready-to-paste prompt for coding agents.

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
pytest tests/ -v               # 206 tests (148 pure-Python + 58 FEM under Docker)
ruff check semi/ tests/
python tests/check_day1_math.py
```

If you added FEM code, also run:

```bash
docker compose run --rm test   # full pytest including FEM tests
python scripts/run_verification.py   # V&V suite (62 studies)
```

## Structure

- `semi/` — the package. Pure-Python modules (`constants`, `materials`, `scaling`, `doping`, `schema`) never import dolfinx. FEM modules (`mesh`, `physics`, `solver`, `run`, `runners`, `bcs`, `postprocess`) import dolfinx lazily inside function bodies, not at module top level.
- `benchmarks/` — self-contained JSON inputs plus README per benchmark.
- `notebooks/` — thin Colab drivers that clone the repo and import the package.
- `tests/` — pytest suite. `tests/*.py` are pure-Python tests; `tests/fem/*.py` need dolfinx.
- `scripts/` — notebook generators, benchmark runner CLI, V&V driver. Not installed.
- `docs/` — physics reference, walkthrough, architecture, ADRs, improvement guide.
- `schemas/` — (M9+) JSON schemas for input and manifest contracts.

## Conventions

- Internal units are SI. JSON input follows device-physics tradition where doping is in cm⁻³ and mobility in cm²/(V·s), converted at ingest.
- Physics-style variable names (`N_A`, `V`, `F`, `psi_L`, `eps_Si`) are fine; PEP 8 naming is explicitly relaxed via ruff config.
- Any new physics module requires (a) a unit test in `tests/` if it has pure-Python logic, (b) an MMS verifier following ADR 0006, (c) an analytical verification in its benchmark.
- Never skip `make_scaling_from_config`. The raw equations have a Jacobian condition number >10³⁰.
- For M9+ work, don't cross milestone boundaries in one PR. The acceptance tests in `IMPROVEMENT_GUIDE.md` are non-negotiable gates.
