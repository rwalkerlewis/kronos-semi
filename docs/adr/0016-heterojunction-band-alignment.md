# ADR 0016: Heterojunction band alignment via cellwise DG0 fields

**Date:** 2026-05-01
**Status:** Accepted
**Milestone:** M17 -- Heterojunctions

---

## Context

The current engine assumes a single semiconductor material throughout each
semiconductor region. The Poisson and continuity kernels hard-code the
electron affinity `chi` and bandgap `E_g` as scalar constants loaded from
`semi/materials.py` at setup time. There is no mechanism for:

- Multiple semiconductor materials in contact (heterojunction).
- Composition-graded alloys (e.g. Al_x Ga_{1-x} As with varying x).
- Position-dependent band offsets that affect carrier transport.

This ADR establishes how `chi(x)` and `E_g(x)` are represented and how
they enter the existing Poisson and DD weak forms.

## Decision

### Representation

`chi` and `E_g` are represented as **DG0 (piecewise constant) scalar
functions** on the mesh. Each cell is assigned a single value; there is
no sub-cell variation. This is consistent with how material properties
(permittivity, mobility, doping) are already handled throughout the engine.

The DG0 representation has three consequences:

1. Inside each cell, `chi` and `E_g` are constant; the gradient `grad(chi)`
   is zero in the distributional sense on each element interior. The
   heterojunction discontinuity is captured by the jump in `chi` across
   the cell boundary, which the DG0 function represents exactly.

2. No jump condition on the current across the heterojunction needs to be
   added explicitly. The Poisson and DD weak forms in their standard
   Galerkin formulation automatically enforce the natural boundary condition
   that the normal component of `J_n` is continuous across any interface
   at which no explicit BC is applied. This is the physically correct
   thermionic-recombination-free heterojunction condition.

3. For composition-graded regions, the alloy fraction is interpolated
   cellwise at setup time. The grain of the interpolation is the mesh
   element size; a fine mesh reproduces a smooth composition profile.

### Equilibrium psi-reference under band alignment

The equilibrium potential `psi_0` is defined relative to the intrinsic
Fermi level. For a homogeneous material this is simply
`psi_0 = phi_T ln(N_net / n_i)` at ohmic contacts, where `phi_T = kT/q`.

For a heterojunction, the intrinsic Fermi level position depends on the
local bandgap and effective densities of states, which vary with position.
The conduction band minimum is:

```
E_C(x) = E_ref - q chi(x) - q psi(x)
```

where `E_ref` is an arbitrary global energy reference. The valence band
maximum is:

```
E_V(x) = E_C(x) - E_g(x)
```

The Fermi-Dirac equilibrium condition `phi_n = phi_p = 0` at an ohmic
contact implies that the electrostatic potential at the contact satisfies:

```
psi_contact = psi_ref_homo  +  Delta chi / V_t
```

where `Delta chi` is the difference in affinity between the contact's
semiconductor material and the reference semiconductor (the material at
the first region processed in setup). In the M16.4 FD mapping (ADR 0015),
`chi_tilde = chi / V_t` enters the `eta_n` expression directly, so the
ohmic BC root-find already accounts for the band offset.

The Boltzmann path does not support heterojunctions (the non-degenerate
approximation breaks at most heterointerfaces). Schema cross-validation
rejects `M17` composition fields unless `physics.statistics.model ==
"fermi-dirac"`.

### Schema

Regions gain an optional `composition` object:

```jsonc
"regions": {
  "<name>": {
    "material": "<base_material>",
    "composition": {
      "<alloy_param>": {
        "profile": "constant" | "step" | "linear",
        "axis":   <int>,
        "values": [<float>, ...]
      }
    }
  }
}
```

- `"constant"`: the alloy parameter is uniform across the region.
- `"step"`: the region is split at a plane normal to `axis`; `values`
  is a two-element list `[x_low_half, x_high_half]`.
- `"linear"`: the alloy parameter interpolates linearly from the first
  cell to the last cell along `axis`.

Material-property functions for alloy materials (e.g. AlGaAs) are
registered in `semi/materials.py` as callables of the alloy parameter.

### Weak form changes

**Poisson kernel** (`semi/physics/poisson.py`): the intrinsic carrier
density `n_i(x)` becomes position-dependent. It is computed cellwise from
`chi(x)`, `E_g(x)`, `N_C(x)`, `N_V(x)` using the relation:

```
n_i^2 = N_C(x) N_V(x) exp(-E_g(x) / (kT))
```

All four quantities are DG0 fields; `n_i` is therefore also DG0.

**Continuity kernels** (`semi/physics/drift_diffusion.py`): the
`eta_n` and `eta_p` arguments to `F_half` (ADR 0015) include `chi_tilde`:

```
eta_n = (psi - phi_n - chi_tilde(x)) / V_t
eta_p = (phi_p - psi + chi_tilde(x) - E_g_tilde(x)) / V_t
```

where `chi_tilde(x) = chi(x) / V_t` and `E_g_tilde(x) = E_g(x) / V_t`.
These are DG0 fields passed into the UFL form as `dolfinx.fem.Function`
objects.

No jump term is added at heterojunction facets. The standard H1
conforming test functions together with the DG0 material coefficients
give the correct physical behavior: the quasi-Fermi potentials `phi_n`
and `phi_p` are continuous across the heterojunction (because they live
in H1), while the carrier densities are discontinuous (because `eta_n`
jumps due to the jump in `chi_tilde`). This is exactly the thermionic
equilibrium condition at an abrupt heterojunction in the limit of zero
interface recombination.

## Consequences

Easier:

- No new FEM infrastructure is required: DG0 fields already exist for
  permittivity and doping.
- No interface numerical flux terms are needed; the standard Galerkin
  form handles the heterojunction naturally.
- Composition-graded profiles are set up entirely in the Python mesh
  setup layer, with no changes to the UFL form assembly.

Harder:

- `n_i(x)` is no longer a global scalar; any code that uses it as a
  scalar must be updated.
- The ohmic-contact `psi` root-find (ADR 0015) must receive the local
  `chi` and `E_g` values for the material at the contact.
- Benchmarks for M17 (HEMT, HBT) are required; no existing benchmark
  exercises position-dependent `chi`.

## Related

- ADR 0004 -- Slotboom variables (primary-unknown choice unchanged).
- ADR 0015 -- FD statistics with Slotboom mapping (chi_tilde enters here).
- `docs/PHYSICS.md` Section 12 "Heterojunctions" (to be written in
  Phase J).
- M17 implementation: `semi/materials.py` alloy functions, updated
  `semi/physics/poisson.py` and `semi/physics/drift_diffusion.py`.
