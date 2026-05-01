# 3D Poisson GPU Benchmark (M15)

A uniform-doped 3D Poisson box of ~500k DOFs (80×80×80 mesh = 531441
nodes) used as the **A2 acceptance test** for the GPU linear-solver
path: the GPU solve must beat CPU-MUMPS by ≥5× wall-clock on this
problem.

This benchmark is intended for the nightly GPU CI job. It is not part
of the default `make verify` loop because it requires PETSc-with-CUDA
(or HIP) and an AMGX or hypre install.

## Run

```bash
# CPU reference
python scripts/run_benchmark.py benchmarks/poisson_3d_gpu/input.json --override solver.backend=cpu-mumps

# GPU
python scripts/run_benchmark.py benchmarks/poisson_3d_gpu/input.json
```

The wall-clock for the linear solve is reported in the manifest under
`solver.linear_solve_wall_s`. The acceptance criterion is

    cpu_wall_s / gpu_wall_s >= 5.0

on a single A100 / MI250 with PETSc 3.21+, AMGX 2.4+ or hypre 2.32+.
