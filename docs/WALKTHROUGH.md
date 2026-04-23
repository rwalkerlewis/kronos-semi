# WALKTHROUGH: What Happens When You Call `semi.run.run(cfg)`

**Audience.** A contributor or coding agent who needs to understand the
full pipeline from JSON input to numerical output, with pointers into
the actual code at each step. If you want the physics first, read
[PHYSICS_INTRO.md](PHYSICS_INTRO.md). If you want equations, read
[PHYSICS.md](PHYSICS.md). This document is about the *mechanics*:
which file does what, in what order, with what side effects.

**Input.** A validated config dict (from `semi.schema.load`).

**Output.** A `SimulationResult` dataclass, documented below.

All line numbers below refer to the v0.8.0 state of the codebase and may
drift slightly after M9+ lands. The *structure* is stable; grep for the
symbol name if a line number is off.

---

## 0. Ten-thousand-foot view

```
JSON file
    │
    ▼  semi.schema.load()            (validate, fill defaults, record source path)
validated cfg dict
    │
    ▼  semi.run.run(cfg)             (dispatch on cfg["solver"]["type"])
    │
    ├─► runners.equilibrium          (equilibrium Poisson, no bias)
    ├─► runners.bias_sweep           (coupled DD, optional bias ramp)
    └─► runners.mos_cv               (MOS capacitor C-V)
             │
             ▼  build_mesh → build_profile → make_scaling → build_forms → solve_nonlinear_block
             │                                                                  │
             │                                                                  ▼  PETSc SNES + MUMPS LU
             ▼
        SimulationResult dataclass (held in memory; no file output yet)
```

Everything below expands on one of these arrows.

---

## 1. JSON → validated config

**Entry point.** `semi/schema.py::load(path)` at roughly line 350.

1. Open the file, `json.loads` it into a plain Python dict.
2. Call `validate(cfg)` which runs `jsonschema.Draft7Validator` against
   `schema.SCHEMA` defined earlier in the same file. On failure, raises
   `SchemaError` with every violation listed.
3. Call `_fill_defaults(cfg)`, which mutates the dict in place to set
   physics defaults (temperature 300 K, SRH on, constant mobility,
   Boltzmann statistics), solver defaults (SNES rtol 1e-8, atol 1e-10,
   direct LU via MUMPS), continuation defaults (min step 1e-4, max
   halvings 6, grow factor 1.5), and output defaults (write XDMF on,
   fields = [potential, n, p]).
4. Stash `_source_path` and `_source_dir` into the dict so that mesh
   files referenced by relative path can be resolved later.

The invariant after this step: any code downstream can rely on every
documented config field being present with a sensible value.

**What validation does not catch.** Schema validation only checks
structure and types. It does not check that:

- Facet names referenced in `contacts` actually exist in
  `mesh.facets_by_plane`.
- Region names referenced in `doping` exist in `regions`.
- Doping profiles are physically reasonable.
- The mesh is fine enough to resolve the Debye length.

Those are runtime errors, reported by the mesh and BC builders with
messages that point you back at the offending JSON field.

---

## 2. Dispatch to a runner

**Entry point.** `semi/run.py::run(cfg)` at line 47.

The function is a six-line dispatcher:

```python
stype = cfg.get("solver", {}).get("type", "equilibrium")
if stype == "equilibrium":              return run_equilibrium(cfg)
if stype in ("drift_diffusion",
             "bias_sweep"):             return run_bias_sweep(cfg)
if stype == "mos_cv":                   return run_mos_cv(cfg)
raise ValueError(f"Unknown solver.type {stype!r}")
```

There are three runners. They share a lot of infrastructure but have
distinct responsibilities:

| Runner         | Solver type(s)                     | What it does                                                 |
|----------------|------------------------------------|--------------------------------------------------------------|
| `equilibrium`  | `"equilibrium"`                    | Single Poisson solve, V_applied = 0 everywhere               |
| `bias_sweep`   | `"drift_diffusion"`, `"bias_sweep"` | Coupled (ψ, Φ_n, Φ_p) solve, with optional adaptive bias ramp |
| `mos_cv`       | `"mos_cv"`                         | Same as bias_sweep but on multi-region mesh with submesh     |

The rest of this walkthrough follows `run_bias_sweep` because it is the
most general. The equilibrium runner is a strict subset (single form
instead of three, no continuation loop), and the MOS runner is a
superset (submesh for carriers, multi-region Poisson).

---

## 3. Mesh construction

**Entry point.** `semi/mesh.py::build_mesh(cfg)` at line 27.

Branches on `cfg["mesh"]["source"]`:

### 3.1 `"builtin"` path (most benchmarks)

`semi/mesh.py::_build_builtin` at line 56.

- 1D: `dolfinx.mesh.create_interval` with the given extents and resolution.
- 2D: `create_rectangle`.
- 3D: `create_box`.

Then `_tag_regions` (line 122) walks `mesh.regions_by_box`, computes
centroids of every cell, and tags cells whose centroid lies inside one
of the axis-aligned boxes. Later boxes override earlier ones, which is
how you get nested tagging (e.g., oxide inside silicon). The result is
a `dolfinx.mesh.MeshTags` on cells.

Then `_tag_facets` (line 156) walks `mesh.facets_by_plane`, finds
boundary facets whose coordinate on the named axis matches the given
value to within `tol` (default 1e-12 m). Returns MeshTags on facets.

### 3.2 `"file"` path (resistor_3d gmsh variant)

`semi/mesh.py::_build_from_file` at line 89.

- Resolve the mesh path relative to the JSON's `_source_dir`.
- For format `"gmsh"`, call `dolfinx.io.gmsh.read_from_msh`. Physical
  groups stored in the `.msh` are returned verbatim as `cell_tags` and
  `facet_tags`; the JSON's `regions_by_box` / `facets_by_plane` are
  *not* consulted in this path.
- For format `"xdmf"`, raise `NotImplementedError`. This is on the M12
  roadmap.

**Return value.** `(mesh, cell_tags, facet_tags)`. Any of the tags may
be `None` if the corresponding JSON section is missing.

---

## 4. Doping profile construction

**Entry point.** `semi/doping.py::build_profile(doping_cfg)` at line 1.

Returns a callable `f(x)` where `x` is a `(dim, N)` numpy array of
coordinates and `f(x)` returns net doping `N_D - N_A` at those points
in SI units (m⁻³).

Supported profile types, all as sums:

- `uniform`: constant.
- `step`: piecewise-constant, sharp transition at `location` on the
  given `axis`.
- `gaussian`: `peak * exp(-|r - center|² / (2·sigma²))`, signed positive
  for donor dopant, negative for acceptor.

The returned closure is passed into `fem.Function.interpolate`, which
evaluates it at every DOF of the P1 doping-field function. That array
is then divided by C₀ (peak doping scale) to get the scaled `N_hat_fn`.

**JSON quirk.** Doping densities in JSON are in cm⁻³ (device-physics
convention); `semi.doping.build_profile` converts to m⁻³ at the
boundary. Everywhere else in the codebase, densities are in m⁻³.

---

## 5. Scaling setup

**Entry point.** `semi/scaling.py::make_scaling_from_config(cfg, ref_mat)`
at line 60.

Returns a `Scaling` dataclass with:

- `V0` = thermal voltage V_t = k_B T / q, ≈ 25.85 mV at 300 K.
- `L0` = characteristic length (today, always 1 m — we keep coordinates
  in meters).
- `C0` = doping scale, equals `max(|N|)` across all profiles, floored
  at 10¹⁶ cm⁻³ to avoid div-by-zero on intrinsic devices.
- `mu0` = mobility scale, set from the reference material.
- `t0` = time scale, `L0² / (mu0 · V0)`.
- `lambda2` = ε_vac · V_t / (q · C0 · L0²), the squared Debye ratio.
- `n_i` = intrinsic carrier density, copied from the reference material.

The invariant: everywhere in the FEM code, `psi_hat = psi / V_t`,
`n_hat = n / C0`, and spatial coordinates are in meters. The scaled
Poisson LHS coefficient is `lambda2 * L0²`, which simplifies to a
dimensional value; this is how we keep the mesh in SI while still
benefiting from scaled variables. See ADR 0002 for the rationale.

---

## 6. Function space and initial guess construction

**In bias_sweep runner, lines 59–70 of `semi/runners/bias_sweep.py`.**

```python
spaces = make_dd_block_spaces(msh)           # three P1 Lagrange spaces
V_psi = spaces.V_psi
N_hat_fn = fem.Function(V_psi, name="N_net_hat")
N_hat_fn.interpolate(lambda x: N_raw_fn(x) / sc.C0)

two_ni = 2.0 * ref_mat.n_i
spaces.psi.interpolate(lambda x: np.arcsinh(N_raw_fn(x) / two_ni))
spaces.phi_n.x.array[:] = 0.0
spaces.phi_p.x.array[:] = 0.0
for fn in (spaces.psi, spaces.phi_n, spaces.phi_p):
    fn.x.scatter_forward()
```

Three things to notice:

1. Three separate P1 scalar spaces, not a MixedElement. Keeping them
   separate makes the block residual assembly explicit and lets us
   reuse the equilibrium Poisson solver for initial guesses on simpler
   problems.

2. Initial guess for ψ is `V_t · asinh(N / 2n_i)`, which is the
   local-charge-neutrality potential. Analytically: at equilibrium in
   the bulk, `n - p = N_net`, which combined with `np = n_i²` gives
   `ψ_bulk = V_t · asinh(N_net / 2n_i)`. That's a much better starting
   point than zero, and it's what lets the SNES solver converge in 4–6
   iterations instead of 20+.

3. Initial guess for Φ_n, Φ_p is zero. That's correct at thermal
   equilibrium, and it's still the right seed for biased solves because
   we always start the bias ramp at V = 0.

**`scatter_forward()`** distributes the local array to ghost DOFs on
other MPI ranks. Skipping it causes subtle wrong-answer bugs in
parallel.

---

## 7. Boundary condition resolution

**Entry point.** `semi/bcs.py::resolve_contacts(cfg, facet_tags, voltages)`
at line 1.

Walks `cfg["contacts"]`, resolves each `facet` field (which may be a
string name or an int tag) to an actual integer facet tag. Looks up the
voltage from the `voltages` dict (passed in per bias step) or from the
contact's baked-in `voltage` field. Returns a list of resolved contact
records.

**Then** `semi/bcs.py::build_dd_dirichlet_bcs(spaces, msh, facet_tags, contacts, sc, ref_mat, N_raw_fn)`
at line 100ish.

For each ohmic contact, builds three Dirichlet BCs:

- ψ = V_applied + ψ_built-in(x_contact), where ψ_built-in is `V_t ·
  asinh(N_net(x_contact) / 2n_i)`.
- Φ_n = V_applied.
- Φ_p = V_applied.

All three in *scaled* units, i.e., divided by V_t. For gate contacts
(multi-region only), ψ is set to `V_applied + φ_ms` (work-function
difference) on the oxide-side facet, and Φ_n, Φ_p get no BC because
they don't live on the oxide.

**Return value.** A flat list of `dolfinx.fem.DirichletBC` objects. Each
BC is tied to the subspace it was constructed on, which is how they're
routed to the right block of the coupled solver.

---

## 8. Weak-form construction

**Entry point.** `semi/physics/drift_diffusion.py::build_dd_block_residual`
at line 66.

Builds three UFL residual forms:

```python
F_psi = L_D² · eps_r · inner(grad(psi), grad(v_psi)) · dx  -  rho_hat · v_psi · dx
F_phi_n = L0² · mu_n · n_hat · inner(grad(phi_n), grad(v_n)) · dx  -  R · v_n · dx
F_phi_p = L0² · mu_p · p_hat · inner(grad(phi_p), grad(v_p)) · dx  +  R · v_p · dx
```

where:

- `L_D² = lambda2 * L0²` is the dimensional squared Debye length.
- `L0²` is the squared reference length (1 m² in the current scaling).
- `n_hat = (n_i / C0) * exp(psi - phi_n)` is the Slotboom expression
  for electron density.
- `p_hat = (n_i / C0) * exp(phi_p - psi)` is the Slotboom expression
  for hole density.
- `rho_hat = p_hat - n_hat + N_hat_fn` is the scaled space charge.
- `R` is the scaled SRH rate:
  `(n·p - n_i²) / (tau_p(n + n_1) + tau_n(p + p_1))`.
- `eps_r` is a scalar `fem.Constant` in single-region problems or a
  cellwise DG0 `fem.Function` in multi-region (silicon + oxide)
  problems.

Homogeneous Neumann on the non-contact boundaries is implicit: the
boundary term from integration by parts is just dropped. No `ds`
integrals appear. This is consistent with "no current leaks through
the sides of the device."

**Return value.** A three-element list `[F_psi, F_phi_n, F_phi_p]`.

---

## 9. The solve itself

**Entry point.** `semi/solver.py::solve_nonlinear_block(F_list, u_list, bcs, prefix, petsc_options, kind, entity_maps)`
at line 79.

Wraps `dolfinx.fem.petsc.NonlinearProblem` (dolfinx 0.10+) which itself
wraps PETSc SNES. The PETSc options set in `DEFAULT_PETSC_OPTIONS` are:

```python
{
    "snes_type": "newtonls",                    # Newton with line search
    "snes_linesearch_type": "bt",               # backtracking (robust)
    "snes_rtol": 1.0e-10,
    "snes_atol": 1.0e-12,
    "snes_max_it": 50,
    "snes_monitor": None,                       # print residuals
    "ksp_type": "preonly",                      # no Krylov iteration, direct solve
    "pc_type": "lu",
    "pc_factor_mat_solver_type": "mumps",       # parallel direct LU
}
```

Single call to `problem.solve()`, then pull `getConvergedReason()` and
`getIterationNumber()` off the SNES solver. Return a dict with
`iterations`, `reason`, `converged`.

**Why direct LU?** For problems up to ~50k DOFs (1D, 2D, small 3D),
MUMPS is fast enough and bulletproof. The block Jacobian of the
coupled DD system has terrible conditioning because of the exponential
coupling, so Krylov + AMG tends to stall without heavy tuning. Direct
LU just works. The cost is that it doesn't scale to large 3D. M15
(GPU + iterative Krylov) is the future.

---

## 10. Bias continuation (bias_sweep only)

**Entry point.** `semi/runners/bias_sweep.py` starting around line 87.

For a `voltage_sweep` on one contact, the runner:

1. Resolves the sweep into a list of voltages: e.g., `[0.0, 0.1, 0.2, ..., 0.6]`.
2. If the list spans zero (has both negative and positive values),
   splits it into two legs: `0 → most_negative → most_positive`.
3. For each leg, runs `ramp_leg(V_start, V_target)` which drives an
   `AdaptiveStepController`:
   - Start with `initial_step = nominal_step = abs(sweep_step)`.
   - Try to solve at `V_prev + step`. If SNES converges, record an
     IV row, advance `V_prev`, and potentially grow the step
     (if the last `easy_iter_threshold` solves all took fewer than
     `easy_iter_threshold` Newton iterations).
   - If SNES fails, halve the step and retry. After `max_halvings`
     halvings at the same target, raise.
4. At every successful step, `semi/postprocess.py::record_iv` evaluates
   the terminal current via a UFL weak form on the swept contact's
   facet:
   ```
   J_n = q · mu_n · n · (grad phi_n · n_outward)
   J_p = q · mu_p · p · (grad phi_p · n_outward)
   I = ∫ (J_n + J_p) dS
   ```
   and appends `{V: V_applied, J: I / A}` to the IV list.

**Solution advection.** Between bias steps, the solution from the
previous step *is* the initial guess for the next. This is what makes
continuation work at all — a 0.1 V step of a forward-biased diode at
0.5 V is a small perturbation of an already-well-converged state.

---

## 11. The result object

**Returned by every runner.** `semi.run.SimulationResult` at
`semi/run.py` line 25. A dataclass holding:

| Field           | Type                  | Meaning                                                           |
|-----------------|-----------------------|-------------------------------------------------------------------|
| `cfg`           | `dict`                | The validated input config (copy)                                 |
| `mesh`          | `dolfinx.mesh.Mesh`   | The mesh                                                          |
| `V`             | `FunctionSpace`       | The scalar P1 space used for ψ                                   |
| `psi`           | `fem.Function`        | Scaled electrostatic potential (internal unit = V_t)             |
| `phi_n`, `phi_p`| `fem.Function`        | Scaled quasi-Fermi potentials (`None` for equilibrium runs)      |
| `psi_phys`      | `np.ndarray`          | Physical ψ in volts, one value per DOF                           |
| `phi_n_phys`    | `np.ndarray`          | Physical Φ_n in volts                                             |
| `phi_p_phys`    | `np.ndarray`          | Physical Φ_p in volts                                             |
| `n_phys`        | `np.ndarray`          | Electron density in m⁻³                                           |
| `p_phys`        | `np.ndarray`          | Hole density in m⁻³                                               |
| `x_dof`         | `np.ndarray`          | Mesh DOF coordinates, shape (N, 3)                                |
| `N_hat`         | `fem.Function`        | Scaled doping field                                               |
| `scaling`       | `Scaling`             | The scaling object used                                           |
| `solver_info`   | `dict`                | `{iterations, reason, converged}` from the last SNES solve        |
| `iv`            | `list[dict]`          | Bias-sweep IV pairs `[{"V": v, "J": j}, ...]` (empty for equilibrium) |
| `bias_contact`  | `str \| None`         | Which contact was swept                                           |

**Important caveat.** `psi`, `phi_n`, `phi_p`, `mesh`, `V`, `N_hat` are
live dolfinx objects with PETSc handles. They cannot cross a process
boundary. They cannot be pickled. They are meaningful only in the
Python process that ran the solve.

The numpy arrays (`_phys` versions, `x_dof`) are the portable
representation. Today, if you want to consume results in a separate
process or in a UI, you have to extract those arrays yourself and
serialize them. **M9 (result artifact writer) exists to fix exactly
this.** Read [IMPROVEMENT_GUIDE.md](IMPROVEMENT_GUIDE.md) §3 for the
proposed on-disk contract.

---

## 12. Where the data flows, in one picture

```
JSON file
  │
  │  semi.schema.load
  ▼
cfg dict  ───────────────────────────────────────────────────────┐
  │                                                               │
  │  semi.run.run → dispatches to runner                          │
  ▼                                                               │
runner reads cfg, calls: build_mesh ────────► (mesh, cell_tags, facet_tags)
                         build_profile ─────► N_raw_fn(x)                  │
                         make_scaling ──────► Scaling object              │
                         make_dd_block_spaces (or equivalent)            │
                         interpolate initial guess                       │
                         resolve_contacts + build_dirichlet_bcs (per step)
                         build_dd_block_residual ──► [F_psi, F_phi_n, F_phi_p]
                         AdaptiveStepController (per step, bias ramp)
                           │
                           │  solve_nonlinear_block
                           ▼
                         PETSc SNES + MUMPS LU
                           │
                           ▼
                         converged ψ, Φ_n, Φ_p (in-place in spaces)
                           │
                           │  postprocess.record_iv (if sweeping)
                           ▼
                         SimulationResult populated
                                                                          │
  ┌───────────────────────────────────────────────────────────────────────┘
  │
  ▼  returned to caller

SimulationResult
```

---

## 13. What to read after this

- To understand *why* the equations look the way they do, read
  [PHYSICS_INTRO.md](PHYSICS_INTRO.md) and then [PHYSICS.md](PHYSICS.md).
- To extend the solver (new physics, new solver type, new contact
  type), find the nearest existing component in the list above and
  model your change on it.
- To consume the output from a UI or external process, the right move
  is to wait for or contribute to M9 (artifact writer). See
  [IMPROVEMENT_GUIDE.md](IMPROVEMENT_GUIDE.md) §4.
- For the code-level invariants the five-layer architecture enforces
  (Layer 3 never imports dolfinx, JSON is the contract, etc.), read
  [ARCHITECTURE.md](ARCHITECTURE.md).
