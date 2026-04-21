# 0005. Docker for dev, conda and FEM-on-Colab for users

- Status: Accepted
- Date: 2026-04-20

## Context

FEniCSx is not pip-installable from PyPI in any practical sense: it
depends on PETSc, MPI, Basix (C++ core), and FFCx compiled against
specific library versions. Users and contributors need a reliable
install path, and the right path differs by audience:

- **Development:** primary author plus any future contributors need
  a tight edit-test loop with deterministic dolfinx behavior and no
  risk of environment drift.
- **Evaluation reviewers:** want to run the benchmarks with zero setup
  cost.
- **Local users on their own machines:** mix of Linux, macOS, and WSL,
  some with conda preferences, some without.

Candidate install paths considered:

- Docker image (official `ghcr.io/fenics/dolfinx/dolfinx:stable`).
- conda-forge (`fenics-dolfinx` package).
- FEM-on-Colab install script (precompiled wheels + deb packages,
  Colab-specific).
- Building dolfinx from source via pip.
- System packages (apt, brew): available but versions lag.

## Decision

Three supported install paths, in this order:

1. **Docker** is the primary development environment. The repo ships a
   `Dockerfile` layered on `ghcr.io/fenics/dolfinx/dolfinx:stable`, with
   `docker-compose.yml` services for interactive dev, test, benchmark,
   and Jupyter. This is what CI and the primary author use.
2. **conda-forge** is the alternative local install for contributors
   who prefer conda or need native performance (no container
   overhead). Instructions in README.
3. **FEM-on-Colab** is the zero-setup path for evaluation reviewers:
   the Colab notebook includes a first cell that runs the
   FEM-on-Colab install script and then clones the repo. Reviewers
   click the Colab badge, wait roughly 30 seconds, and run cells.

Do not attempt to support a "pip install dolfinx" path.

## Consequences

Easier:

- Three clearly documented install paths, each validated by the
  project: Docker (CI), conda (README), Colab (notebook).
- Docker is bit-exact reproducible for development; conda is fast on
  native hardware; Colab is accessible with no local tooling at all.
- When dolfinx 0.11 lands, we update one base image, bump the Colab
  install URL, and update the conda channel spec.

Harder:

- Three install paths is more documentation surface than one.
- The Docker-on-Apple-Silicon path is slow (dolfinx base image is
  x86_64; on M-series Macs it runs under emulation). Contributors on
  Apple Silicon are encouraged to use conda-forge instead.
- First-time Docker build pulls several GB of base image.

Rejected alternatives:

- **pip-from-source dolfinx:** multi-hour build, brittle dependence on
  the user's PETSc, MPI, and compiler versions, breaks silently on
  PETSc major version bumps.
- **apt / brew / pacman system packages:** version lag makes the API
  drift documented in ADR 0003 unworkable.

## Related

- `Dockerfile`, `docker-compose.yml` at repo root.
- Invariant 9 in `PLAN.md`.
- README's "Quick start on Colab" and "Docker" sections.
