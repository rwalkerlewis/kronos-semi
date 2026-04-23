# 2D MOS capacitor benchmark (M6: 2D MOS capacitor)

First 2D benchmark and first multi-region device in kronos-semi.

## Device

- Silicon substrate: 500 nm x 500 nm, uniform N_A = 1e17 cm^-3 (p-type)
- SiO2 gate oxide: 500 nm x 5 nm on top
- Body contact: ohmic, y = 0
- Gate contact: ideal-gate Dirichlet (phi_ms = 0), y = 505 nm
- Left / right edges: natural (Neumann) boundaries

## Mesh

Uniform 2D rectangle, triangles: 4 cells laterally, 505 cells
vertically (1 nm per cell). Interface at y = 500 nm is a grid line;
the oxide gets exactly 5 cells.

Total DOFs: 2530 (P1 Lagrange on 4 x 505 rectangle -> 5 x 506 nodes).

## Physics

Multi-region equilibrium Poisson (`solver.type == "mos_cv"`):

- Full-mesh stiffness `L_D^2 * eps_r(x) * grad psi . grad v` with
  cellwise DG0 `eps_r` (Si: 11.7, SiO2: 3.9)
- Boltzmann space-charge `n_i (exp(-psi) - exp(psi)) + N_net`
  restricted to silicon via `dx(subdomain_id=semi_tag)`
- Oxide: Laplacian only (no space charge, no Slotboom variables)

At each V_gate the gate charge per area is

    Q_gate(V_gate) = -(q / W_lat) * integral_{Omega_Si} rho dA

and `C_sim = dQ_gate/dV_gate` via centered finite difference.

## Verification target

Depletion-approximation C-V curve within 10% in the window
[V_FB + 0.2, V_T - 0.1] V. For the ideal-gate M6 device under our
psi=0-at-intrinsic BC convention:

    V_FB = phi_ms - phi_F = 0 - 0.417 V = -0.417 V
    V_T  = V_FB + 2 phi_F + sqrt(4 eps_s q N_A phi_F) / C_ox = +0.658 V
    window = [-0.217, +0.558] V

Worst observed error in the window: 9.25% at V_gate = -0.200 V.

Monotone non-increasing C_sim across the window is also checked (C
drops from near C_ox in accumulation through a minimum at strong
inversion onset).

## Excluded regimes

- **Accumulation** (V_gate < V_FB): holes pile up at the interface;
  the depletion approximation does not model the accumulation
  surface charge. C_sim recovers toward C_ox but is not gated here.
- **Strong inversion** (V_gate > V_T): minority-carrier inversion
  layer forms. Quasi-static C-V recovery toward C_ox is not well
  modelled without a transient solver. Not gated.

The full V_gate sweep [-0.9, +1.2] V is rendered on `cv.png` for
reference (simulation line crosses through both excluded regions);
the verifier window is shaded.

## Running

```
docker compose run --rm benchmark mos_2d
```

Produces four PNGs in `results/mos_2d/`:
- `psi_2d.png`: tricontourf of psi(x, y) at the last V_gate
- `potentials_1d.png`: central-column psi(y) slice
- `cv.png`: C_sim vs depletion-approximation theory
- `qv.png`: Q_gate vs V_gate

## References

- `docs/mos_derivation.md`: full derivation gate (device, interface
  conditions, submesh formulation, gate BC, C-V theory, MMS
  construction).
- `docs/PHYSICS.md` Section 6: condensed physics reference including
  the V_FB shift from the BC convention.
