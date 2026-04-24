# 2D long-channel n-MOSFET (M12)

First benchmark that uses the M12 `geometry` input block to drive gmsh
from a `.geo` file. Also the first benchmark with four contacts
(source, drain, body, gate) of mixed type on a single device.

## Device

- Silicon body: 5 um wide (x), 2 um deep (y), uniform p-type, N_A = 5e17 cm^-3
- Gate oxide: 5 nm SiO2 from x = 1.5 um to x = 3.5 um
- Source / drain contacts (ohmic): y = 2 um, 1 um wide each
- Body contact (ohmic): y = 0, full width
- Gate contact (ideal, phi_ms = 0): y = 2 um + 5 nm over the oxide
- Channel length L = 2 um (gate x-extent); drain contact length L_drain = 1 um

## Mesh

Parametric: `geometry.source = "gmsh_geo"` points at `mosfet_2d.geo`.
Per-point characteristic lengths shape the mesh so the channel and
oxide are finer than the bulk body. The `characteristic_length`
JSON field caps cells uniformly; the convergence test walks it through
`[1e-7, 5e-8, 2.5e-8]` m.

At the nominal `characteristic_length = 5e-8`, gmsh produces on the
order of 13k triangles. The engine content-hashes the `.geo` bytes
plus the `clmax` value and caches the resulting `.msh` in
`$KRONOS_MESH_CACHE` or `~/.cache/kronos-semi/mesh/`, so repeat runs
skip gmsh.

## Physics

- Coupled Slotboom drift-diffusion (M2+). The bias_sweep runner
  currently treats multi-region meshes as single-material silicon
  (see "Known caveat" below).
- Ohmic BCs pin psi and both quasi-Fermi potentials. Gate BC pins
  psi only, at the top of the oxide with the workfunction offset.
- Constant mobilities (mu_n = 1400 cm^2 / V s, mu_p = 450 cm^2 / V s).
- SRH is off.

## Verification (what this benchmark gates)

Hard gates (must pass):
1. **Sweep reaches V_DS = 0.1 V.** The geometry + solver end-to-end run.
2. **V_T analytic value is consistent** with the engine's psi=0-at-
   intrinsic BC convention: V_FB = phi_ms - phi_F, V_T = V_FB + 2 phi_F
   + sqrt(2 q eps_Si N_A (2 phi_F)) / C_ox. Reported, not tolerance-
   asserted (the value is documentary rather than compared to sim).
3. **Triode-regime I_D within +/-20% of long-channel theory** at
   V_GS = 1.5 V, V_DS = 0.1 V. Long-channel formula
   `I_D / W = mu_n_eff * C_ox / L * (V_GS - V_T) * V_DS` with
   mu_n_eff = 600 cm^2 / V s (conservative surface mobility, no M16
   field-dependent model). **Currently fails loudly** until the
   multi-region runner is wired; see "Known caveat" below.
4. **Simulated I_D is finite and non-negative.** Catches NaN or sign
   bugs in the facet-current integration.

## Known caveat: simulated I_D vs triode theory fails gate 3

The simulated drain current at V_GS = 1.5 V, V_DS = 0.1 V is ~3e+02
A/m per unit depth, while long-channel triode theory predicts ~9.9
A/m. The ~30x over-prediction (ratio ~2920 percent) is real and
reproducible; it is not a mesh or Newton-convergence artifact (the
Cauchy-sequence test in `tests/fem/test_mosfet_mesh_convergence.py`
is intentionally not gated on the drain-current rate for exactly this
reason: see that module's docstring).

The root cause is that `semi/runners/bias_sweep.py` does not wire
multi-region coupled drift-diffusion. The runner calls
`make_dd_block_spaces(msh)` on the full mesh and assembles
`build_dd_block_residual(..., ref_mat.epsilon_r, ...)` with a scalar
`epsilon_r` from the reference material (silicon, 11.7). The oxide
region therefore uses Si's permittivity rather than SiO2's (3.9), and
continuity equations are assembled on oxide cells that should be pure
insulator.

`run_mos_cv` (M6) handles multi-region for the equilibrium Poisson
case via `build_equilibrium_poisson_form_mr` and
`build_submesh_by_role`. Lifting the same pattern into `bias_sweep`
(i.e., a `build_dd_block_residual_mr` path driven off `cell_tags` plus
`build_submesh_by_role(..., role="semiconductor")`) is the obvious
follow-up milestone. That plumbing exists in
`semi/physics/drift_diffusion.py::DDBlockSpacesMR` already; only the
runner wiring is missing.

Until that lands, the MOSFET benchmark exercises the M12 geometry
pipeline (gmsh subprocess -> content-hashed cache -> dolfinx load ->
tag resolution) end to end, and the verifier's triode-regime gate
fails loudly so CI correctly blocks M12 merge until the multi-region
runner wiring lands.

## Running

```
docker compose run --rm benchmark mosfet_2d
```

Produces two PNGs in `results/mosfet_2d/`:
- `psi_2d.png`: tricontourf of psi(x, y) at the final bias.
- `id_vds.png`: simulated `|I_D|(V_DS)` vs analytical triode theory.

## References

- `docs/PHYSICS.md` (multi-region Poisson and DD derivation; see also
  the Si/SiO2 interface section in §6 for the MOS-capacitor precedent).
- `semi/geometry.py` (subprocess + cache contract).
- `tests/test_geometry.py` (cache-hit and GeometryError tests).
- `tests/fem/test_mosfet_mesh_convergence.py` (geometry-pipeline sweep).
