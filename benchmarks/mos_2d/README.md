# 2D MOS capacitor benchmark (planned)

Not yet implemented. Planned for Day 5 of the build.

## Physics scope

- Two regions: silicon substrate + SiO₂ gate oxide
- Poisson on both regions (no carriers in oxide, full space-charge in semiconductor)
- Gate contact: Dirichlet BC on top of oxide with work-function offset
- Body contact: ohmic on bottom of silicon
- Boltzmann equilibrium carriers (no DD needed for C-V)

## Verification target

MOS C-V sweep from accumulation to inversion, comparing:
- Flat-band voltage
- Threshold voltage (onset of inversion)
- Oxide capacitance at accumulation (analytical: $C_{ox} = \varepsilon_{ox}/t_{ox}$)
- Inversion capacitance

## Status

Framework elements already in place:
- Multi-region JSON schema supports this (`regions` map with distinct materials)
- Material database includes SiO₂ and HfO₂
- Mesh tagging supports axis-aligned box subregions

What's missing:
- Submesh / entity-map handling so that Φₙ, Φₚ live only on semiconductor
- Gate BC implementation (with work-function offset)
- Two-material Poisson form (piecewise εᵣ via a cell-tag-indexed Constant)
