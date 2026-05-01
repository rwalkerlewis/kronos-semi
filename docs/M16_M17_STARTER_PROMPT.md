# M16 + M17 starter plan: complete the drift-diffusion physics surface

Handoff plan for GitHub Copilot. Mirrors the convention of
`docs/M9_STARTER_PROMPT.md`. Commit it as `docs/M16_M17_STARTER_PROMPT.md`
after Phase A lands.

This plan covers **M16 (physics completeness)** plus **M17
(heterojunctions)** as enumerated in `docs/IMPROVEMENT_GUIDE.md` section 4.
Together these close the drift-diffusion gap relative to the
COMSOL-Semiconductor surface that kronos-semi targets. Today the engine
ships Slotboom DD with constant mobility, SRH recombination, Boltzmann
statistics, ohmic and ideal-gate contacts, and a single-material
assumption. After this plan: field- and surface-dependent mobility,
Auger and radiative recombination, Fermi-Dirac statistics, incomplete
dopant ionization, Schottky contacts, band-to-band and trap-assisted
tunneling, and heterojunctions with position-dependent chi and E_g.

## Honest scope

Per the M16 entry in `docs/IMPROVEMENT_GUIDE.md`: **2 to 4 days per
sub-item, one PR per item, do not batch**. Plus M17 at ~4 days, plus
documentation closeout. Realistic total: **17 to 30 working days**. If
the calendar is tighter than that, drop scope (skip M16.6 tunneling
and M16.7 incomplete ionization; those are the lowest-yield items per
COMSOL-feature parity), do not compress the per-item discipline.

## Required reading (Copilot must read these before any phase)

Per `CONTRIBUTING.md` "Before you start":

1. `PLAN.md` in full. Confirm the package version and that no in-flight
   PR conflicts with physics changes. If M15 is in flight, **stop and
   wait**; physics changes that touch `semi/solver.py` or
   `semi/runners/*.py` should land after M15's CPU/GPU equivalence
   tests are green.
2. `docs/IMPROVEMENT_GUIDE.md` sections 1, 2, 4 (M16), 4 (M17), 8 (anti-goals).
3. `docs/PHYSICS.md` end-to-end. Especially Sections 2 (scaled DD),
   2.5 (full derivation), 5.5 (multi-region MMS), 6 (MOS reference),
   and 7 (3D extension).
4. `docs/PHYSICS_INTRO.md` if you do not have a TCAD background.
5. `docs/ARCHITECTURE.md` for the five-layer rule. Constitutive models
   (mobility, recombination, statistics, ionization) are Layer 3 (pure
   Python core, no dolfinx); FEM kernels that consume them are Layer 4.
6. `docs/adr/`: 0001, 0002, 0003, 0004 (Slotboom variables, locked),
   0006 (V&V strategy, governs MMS for every new model), 0007, 0011,
   0014. Skim the rest.
7. The five existing physics modules in `semi/physics/`:
   `slotboom.py`, `recombination.py`, `drift_diffusion.py`,
   `poisson.py`, `axisymmetric.py`. The patterns there are what every
   new module must mirror.

## Triage

The M16 specification in `docs/IMPROVEMENT_GUIDE.md` section 4 lists six
sub-items but does not number them. This plan numbers them M16.1
through M16.7 (adding incomplete-ionization as M16.7 since it appears
in `PLAN.md` "Non-goals" but is naturally part of the DD completeness
pass). M17 follows.

Order matters because of dependencies:

1. M16.1 Caughey-Thomas mobility. **No dependencies.** Start here.
2. M16.2 Lombardi surface-mobility. Depends on M16.1's mobility framework.
3. M16.3 Auger and radiative recombination. **No dependencies.** Can run in parallel with M16.1/M16.2.
4. M16.4 Fermi-Dirac statistics. **No dependencies on prior M16 items**, but is itself a dependency for M16.7 and M17.
5. M16.5 Schottky contacts. **No dependencies.** Can run in parallel.
6. M16.6 Tunneling (Kane band-to-band, Hurkx trap-assisted). Depends on M16.3 (recombination framework reuse) and M16.4 (carrier-density expressions).
7. M16.7 Incomplete dopant ionization. Depends on M16.4 (occupation statistics).
8. M17 Heterojunctions. Depends on M16.4 (degeneracy at the barrier).

Practical Copilot execution: run M16.1, M16.3, M16.4, M16.5 in any
order (or in parallel branches); then M16.2, M16.6, M16.7; then M17;
then closeout.

`PLAN.md` "Next task" currently lists M16.1 as the first sub-item.
The owner has explicitly re-prioritized full DD ahead of M15 if M15
is not yet in. Update `PLAN.md` in Phase A.

---

## Phase A: refresh dev guide and stage the milestone numbering

No code. Documentation only.

1. `docs/IMPROVEMENT_GUIDE.md` section 1: refresh the version banner and the
   "What exists / does not exist" lists if M15 has shipped since the
   last refresh. If section 1 still says "as of v0.8.0" or any version older
   than the current `pyproject.toml`, fix it. Source: `PLAN.md`
   "Completed work log".
2. `docs/IMPROVEMENT_GUIDE.md` section 4 M16: split the existing six-bullet
   M16 entry into seven numbered sub-milestones M16.1 through M16.7
   (insert M16.7 incomplete ionization after M16.6 tunneling, since
   it depends on M16.4). Each sub-milestone gets:
   - One-paragraph deliverable
   - Schema extension snippet
   - Acceptance test (analytical or textbook reference)
   - LOC estimate (lift from the existing M16 entry)
   - Dependencies
3. `PLAN.md` "Next task": replace whatever is there with "M16.1
   Caughey-Thomas mobility (first of seven M16 sub-items + M17;
   tracked in `docs/M16_M17_STARTER_PROMPT.md`)".
4. `PLAN.md` "Roadmap" table: split the M16 row into M16.1-M16.7 plus
   keep M17. Status `Planned` for all.
5. `PLAN.md` "Non-goals": remove the items that this plan promotes
   into scope (Caughey-Thomas, Auger and radiative, Fermi-Dirac,
   incomplete ionization, tunneling, Schottky, heterojunctions).
   Leave the GUI item.
6. New ADR `docs/adr/0015-slotboom-with-fd-statistics.md`: drafts the
   convention for mapping Slotboom variables under Fermi-Dirac
   (the eta = (E_F - E_C)/kT mapping replaces the Boltzmann
   exp(psi - phi) form). This is needed before M16.4 starts; if it
   conflicts with ADR 0004, ADR 0015 supersedes the relevant
   subsection. Mark ADR 0004 with the supersession reference.
7. New ADR `docs/adr/0016-heterojunction-band-alignment.md`:
   establishes the convention for chi(x) and E_g(x) as cellwise DG0
   fields and how the affinity discontinuity at a heterojunction
   propagates through the equilibrium psi-reference. Required before
   M17 starts.
8. Run `pytest tests/ -v && ruff check semi/ tests/`. No code paths
   changed; tests should be green.

**Commit message:** `docs: stage M16.1-M16.7 and M17; ADR 0015 (FD statistics) and ADR 0016 (heterojunction band alignment)`

---

## Phase B: M16.1 Caughey-Thomas field-dependent mobility

**Why:** constant mobility breaks above ~10^4 V/cm, which is everywhere
in real device operation. Velocity saturation is the first-order fix.

1. New module `semi/physics/mobility.py` (Layer 3 / Layer 4 split: the
   parameter database and the closed-form `mu(F, N)` expression are
   pure Python; the UFL form that wires it into the continuity
   equation is Layer 4 and lives in this same module since
   `drift_diffusion.py` already does this).
2. Schema additions to `solver.physics.mobility`:
   ```jsonc
   {
     "model": {"enum": ["constant", "caughey-thomas"], "default": "constant"},
     "caughey_thomas": {
       "v_sat_n":  {"type": "number", "default": 1.07e7, "description": "cm/s"},
       "v_sat_p":  {"type": "number", "default": 8.37e6, "description": "cm/s"},
       "beta_n":   {"type": "number", "default": 2.0},
       "beta_p":   {"type": "number", "default": 1.0}
     }
   }
   ```
   Material-default values for Si live in `semi/materials.py`.
3. The CT formula is `mu(F_par) = mu_0 / (1 + (mu_0 |F_par| / v_sat)^beta)^(1/beta)`
   where `F_par` is the field component along the current. Implement it
   as a UFL expression that consumes `grad(psi)`. The DD weak form in
   `drift_diffusion.py` gets a new branch that calls the CT mu instead
   of the constant mu_0. Slotboom variables stay; this is a coefficient
   change, not a stabilization change (ADR 0004 holds).
4. Schema 1.4.0 -> 1.5.0. Add cross-field validation: CT requires a
   semiconductor region.
5. Benchmark `benchmarks/pn_1d_velocity_saturation/`: 1D pn diode at
   high reverse bias (V = -10 V), where the field in the depletion
   region exceeds the CT critical field. Verifier compares the
   simulated saturation current to the textbook closed-form
   `J_sat ~= q n_dep v_sat` within 15%.
6. MMS verifier per ADR 0006: extend `semi/verification/mms_dd.py`
   with a CT variant. The MMS forcing has to account for the
   field-dependent diffusivity; document the manufactured solution in
   `docs/mms_dd_derivation.md` Section 4 (new).
7. Tests: pure-Python `tests/test_mobility_caughey_thomas.py` for the
   closed-form expression; FEM-gated `tests/fem/test_pn_velocity_sat.py`
   for the benchmark.

**Commit message:** `feat(physics): Caughey-Thomas field-dependent mobility (M16.1)`

**Acceptance:** existing benchmarks pass with `mobility.model:
"constant"` byte-identically to today; new benchmark passes; MMS rate
>=1.95 in L2 and >=0.95 in H1 at the finest pair.

---

## Phase C: M16.2 Lombardi surface-mobility

**Why:** required for realistic MOSFET I-V in the inversion regime.
Surface roughness and Coulomb scattering at the Si/SiO2 interface
degrade mu by 2-5x relative to bulk; CT alone misses this.

1. Extend `semi/physics/mobility.py` with the Lombardi surface model.
   The Lombardi formula combines bulk mu_b, acoustic-phonon surface mu_ac,
   and surface-roughness mu_sr via Matthiessen's rule:
   `1/mu_total = 1/mu_b + 1/mu_ac + 1/mu_sr` where mu_ac and mu_sr depend on
   E_perp, the field component perpendicular to the Si/SiO2 interface.
2. Add `model: "lombardi"` to the mobility schema. Lombardi requires a
   declared interface: schema cross-field validation rejects the model
   unless at least one facet has `role: "semi_insulator_interface"`.
3. Computing E_perp: at each cell adjacent to the interface, take
   `E . n_hat` where `n_hat` is the outward normal of the interface facet
   pulled back into the cell. Implement via a DG0 facet-to-cell
   projection in a new helper `semi/physics/interface_projection.py`.
4. Benchmark: extend the existing `benchmarks/mos_2d/` with a strong-
   inversion sweep (V_gate from V_T to V_T + 2 V) and verify that the
   peak channel mu is 2-5x lower than bulk Si mu at the same N_A. New
   benchmark `benchmarks/mosfet_2d_lombardi/` if a separate device is
   needed for the verifier.
5. MMS variant in `mms_dd.py` that builds an interface and verifies the
   E_perp projection converges at first order.

**Commit message:** `feat(physics): Lombardi surface-mobility model (M16.2)`

**Depends on:** M16.1.

---

## Phase D: M16.3 Auger and radiative recombination

**Why:** SRH alone misses the high-injection regime (Auger) and the
direct-gap radiative limit (radiative). Both are needed for any LED,
laser, or high-current diode model.

1. Extend `semi/physics/recombination.py` with two new model kernels:
   - `R_Auger = (C_n n + C_p p)(n p - n_i^2)` with material-dependent
     C_n, C_p (Si: 2.8e-31, 9.9e-32 cm^6/s).
   - `R_rad   = B (n p - n_i^2)` with material-dependent B (GaAs:
     7.2e-10 cm^3/s; Si: ~1e-14, effectively negligible).
2. The recombination schema becomes a list, allowing multiple models
   to sum:
   ```jsonc
   "recombination": [
     {"model": "srh", "tau_n0": 1.0e-7, "tau_p0": 1.0e-7, "E_t": 0.0},
     {"model": "auger", "C_n": 2.8e-31, "C_p": 9.9e-32},
     {"model": "radiative", "B": 1.0e-14}
   ]
   ```
   The existing scalar form auto-migrates to a single-element list at
   schema-load time so existing benchmarks keep validating.
3. Benchmark `benchmarks/pn_1d_high_injection/`: 1D pn diode at very
   high forward bias (V > V_bi) where Auger dominates. Verifier
   compares J(V) slope above the Auger crossover to the analytical
   `J_Auger ~ exp(qV/3kT)` (slope shifts from kT/q to 3kT/q).
4. Benchmark `benchmarks/gaas_pn_1d/`: GaAs pn diode where radiative
   recombination dominates. Reuses the Si pn benchmark structure with
   GaAs material card. Verifier checks the radiative-limited diode
   ideality factor.
5. MMS variants for each new kernel.

**Commit message:** `feat(physics): Auger and radiative recombination (M16.3)`

**Depends on:** none.

---

## Phase E: M16.4 Fermi-Dirac statistics

**Why:** Boltzmann breaks down at degenerate doping (>=10^19 cm^-3 in Si)
and at heterojunction barriers. Required for M16.7 and M17.

This is the largest M16 sub-item (~400 LOC) and it touches the core
Slotboom mapping. Read ADR 0015 (drafted in Phase A) before starting.

1. New module `semi/physics/statistics.py`:
   - `boltzmann_n(psi, phi_n, scaling)` and `boltzmann_p(...)` lifted
     out of `slotboom.py` so the API is symmetric with FD.
   - `fermi_dirac_n(psi, phi_n, scaling, approximation)` where
     `approximation in {"blakemore", "fdi-1half", "fdi-3half"}`.
   - Blakemore: `n = N_C / (exp(-eta) + 0.27)` for fast iteration.
   - Full FDI-1/2: implemented via the Bednarczyk-Bednarczyk rational
     approximation (accurate to 4e-4 for all eta). Pure-Python; can be
     called from UFL via `dolfinx.fem.Expression`.
2. Schema:
   ```jsonc
   "physics": {
     "statistics": {
       "model":         {"enum": ["boltzmann", "fermi-dirac"], "default": "boltzmann"},
       "approximation": {"enum": ["blakemore", "fdi-1half"], "default": "blakemore"}
     }
   }
   ```
3. The DD weak form in `drift_diffusion.py` learns to dispatch on
   `statistics.model`. Slotboom variables remain primary; the change
   is in the n(psi, phi_n) and p(psi, phi_p) mapping. Per ADR 0015, the
   Slotboom-FD form replaces `n = n_i exp(psi - phi_n)` with
   `n = N_C F_half((psi - phi_n - chi_tilde) / V_t)` where chi_tilde is
   the scaled affinity.
4. Benchmark `benchmarks/pn_1d_degenerate/`: heavily doped Si pn
   junction (N_D = N_A = 5e19 cm^-3). Verifier compares V_bi between
   Boltzmann and FD; FD should produce a ~30 mV smaller V_bi at this
   doping (textbook reference: Sze 4th ed. Fig. 1.13).
5. MMS for FD in `mms_dd.py`. Manufactured solution must use the FD
   inverse to set boundary data consistently.
6. Tests: pure-Python `tests/test_statistics.py` for the FD formulas
   (compare Blakemore and FDI-1/2 to a tabulated reference); FEM-gated
   `tests/fem/test_pn_degenerate.py`.

**Commit message:** `feat(physics): Fermi-Dirac statistics with Slotboom mapping (M16.4)`

**Depends on:** ADR 0015 (Phase A).

---

## Phase F: M16.5 Schottky contacts

**Why:** Schottky diodes, metal-semiconductor contacts in MESFETs, and
gate leakage in MOS structures all need thermionic-emission BCs.

1. Extend `semi/bcs.py` with a Schottky BC class:
   - Boundary condition is Robin-type: `J_n . n = q v_n_th (n - n_eq)`
     where `v_n_th = A_star T^2 / (q N_C)` is the thermionic velocity and
     `n_eq = N_C exp(-q phi_B / kT)` with phi_B the Schottky barrier
     height.
   - Same for holes with `phi_Bp = E_g - phi_Bn`.
2. Schema: extend `contacts[].type` enum with `"schottky"` plus
   required `barrier_height_n` (eV) and optional `barrier_height_p`.
   Cross-field validation: Schottky requires a metal `material` card.
3. Add metals to `semi/materials.py`: Al, W, Ti with work functions.
4. Benchmark `benchmarks/schottky_diode_1d/`: Al-Si Schottky diode
   reverse-saturation current. Verifier against
   `J_s = A_dstar T^2 exp(-q phi_B / kT)` within 20%.
5. ADR `docs/adr/0017-schottky-thermionic-emission.md`: documents the
   Robin BC formulation and the sign convention (current INTO
   semiconductor, mirroring the M14 errata convention).

**Commit message:** `feat(bcs): Schottky thermionic-emission contact (M16.5)`

**Depends on:** none.

---

## Phase G: M16.6 Band-to-band and trap-assisted tunneling

**Why:** required for any meaningful Zener diode, floating-gate, or
ESD-protection simulation. Largest single sub-item (~600 LOC).

1. New module `semi/physics/tunneling.py`:
   - Kane band-to-band: `G_BTB = A_K E^2 exp(-B_K / E)` where E is the
     local field magnitude, A_K and B_K are material parameters.
   - Hurkx trap-assisted: a field-enhancement factor multiplying the
     SRH rate, `R_SRH_enhanced = (1 + Gamma(E)) R_SRH` with
     `Gamma(E) = 2 sqrt(3 pi) (E / E_0) exp((E/E_0)^2)`.
2. Both kernels appear as additional generation/recombination terms
   added to the right-hand side of the continuity equations. Wire
   them in via the recombination list from M16.3.
3. Schema:
   ```jsonc
   {"model": "kane",  "A_K": 4.0e14, "B_K": 1.9e7},
   {"model": "hurkx", "E_0": 1.0e7, "tau_n0": 1.0e-7, "tau_p0": 1.0e-7}
   ```
4. Benchmark `benchmarks/zener_diode_1d/`: heavily doped pn junction
   reverse-biased into Zener breakdown. Verifier checks the breakdown
   voltage matches the Sze closed-form to within 25% (BTB tunneling is
   famously sensitive to material parameters; tighter tolerance is
   not warranted).
5. Benchmark `benchmarks/trap_assisted_diode_1d/`: a moderately doped
   pn diode at intermediate reverse bias where Hurkx enhancement
   dominates over pure SRH. Verifier shows J(V) deviating from the
   pure-SRH curve.

**Commit message:** `feat(physics): Kane BTB and Hurkx trap-assisted tunneling (M16.6)`

**Depends on:** M16.3 (recombination list infrastructure), M16.4
(carrier-density expressions in the degenerate regime where BTB
dominates).

---

## Phase H: M16.7 Incomplete dopant ionization

**Why:** at low temperature (T < 100 K) or for deep dopants, fully
ionized doping is wrong by 1-3 orders of magnitude. Required for any
cryogenic device or wide-gap (SiC, GaN) work.

1. New module `semi/physics/ionization.py`:
   - `f_D(N_D, E_D, E_F, T) = 1 - 1/(1 + g_D exp((E_F - E_D)/kT))` for
     donor ionization, similar for acceptors.
   - Material database extension: dopant ionization energies E_D, E_A
     and degeneracy factors g_D=2, g_A=4 in `semi/materials.py`.
2. The Poisson source term changes from `q (p - n + N_D - N_A)` to
   `q (p - n + f_D N_D - f_A N_A)`. The change happens in
   `semi/physics/poisson.py` and `axisymmetric.py`.
3. Schema:
   ```jsonc
   "physics": {
     "ionization": {
       "model": {"enum": ["full", "incomplete"], "default": "full"}
     }
   }
   ```
4. Benchmark `benchmarks/cryogenic_resistor_1d/`: 1D doped Si bar at
   T = 50 K. Verifier compares simulated carrier density to the
   textbook freeze-out formula
   `n ~= sqrt(N_C N_D / 2) exp(-E_D/2kT)` within 10%.
5. Schema cross-validation: `ionization.model: "incomplete"` requires
   `physics.statistics.model: "fermi-dirac"` (Boltzmann + freeze-out
   is inconsistent at the relevant temperatures).

**Commit message:** `feat(physics): incomplete dopant ionization (M16.7)`

**Depends on:** M16.4.

---

## Phase I: M17 Heterojunctions (position-dependent chi and E_g)

**Why:** AlGaAs/GaAs HEMTs, SiGe/Si HBTs, InGaAs/InP photodiodes, none
of which the engine can model today. Largest single phase (~1 week).

Read ADR 0016 (Phase A) before starting.

1. Material model extension (`semi/materials.py`): every material card
   gains `electron_affinity` (eV) and a `bandgap_at_300K` field. Add
   AlGaAs (with Al fraction parameter), InGaAs, GaN, SiGe.
2. Mesh / region extension: regions can declare composition gradients,
   producing cellwise chi(x) and E_g(x) as DG0 functions. Schema:
   ```jsonc
   "regions": {
     "channel": {
       "material": "AlGaAs",
       "composition": {"x_Al": {"profile": "step", "axis": 1, "values": [0.3, 0.0]}}
     }
   }
   ```
3. Poisson and continuity kernels in `poisson.py` and
   `drift_diffusion.py` learn to consume cellwise chi and E_g.
   Continuity of `E_C(x) = -q psi(x) - q chi(x)` and
   `E_V(x) = E_C(x) - E_g(x)` across heterojunctions is enforced by
   the equations naturally; do not add an interface jump condition.
4. The Slotboom mapping with FD statistics from M16.4 already accommodates
   chi(x); confirm by re-running the M16.4 MMS with a constant chi != 0.
5. Benchmark `benchmarks/hemt_2d/`: AlGaAs/GaAs HEMT. Verifier checks
   the 2DEG sheet density at the heterojunction interface against
   the Sze 4th ed. Section 7.5 closed-form within 30%.
6. Benchmark `benchmarks/hbt_1d/`: SiGe/Si HBT. Verifier checks the
   collector current enhancement vs the homojunction reference.

**Commit message:** `feat(physics): heterojunctions with position-dependent chi and E_g (M17)`

**Depends on:** M16.4, ADR 0016 (Phase A).

---

## Phase J: closeout

After all sub-phases land:

1. `PLAN.md`:
   - Move M16.1 through M16.7 and M17 in section Roadmap from `Planned` to
     `Done` with dates.
   - Append entries to section "Completed work log" following the M14.2
     entry's format. Append-only; do not edit prior entries.
   - Refresh section "Current state". Bump `pyproject.toml` from the
     post-M15 version (likely 0.15.x) to 0.16.0 then 0.17.0 if you
     bump per major-physics surface, or do a single bump at closeout
     to reflect the cumulative DD completeness pass.
2. `docs/IMPROVEMENT_GUIDE.md`:
   - Mark M16.1 through M16.7 and M17 Done in section 4 with one-line
     summaries each. Move detail into the section 10 "Shipped milestone
     detail" appendix.
   - Append a section 9 changelog entry: "Drift-diffusion physics surface
     completed: CT and Lombardi mobility, Auger and radiative
     recombination, Fermi-Dirac statistics, Schottky contacts, BTB
     and Hurkx tunneling, incomplete ionization, heterojunctions."
3. `docs/PHYSICS.md`:
   - New Section 8 "Mobility models" covering CT and Lombardi.
   - New Section 9 "Recombination model surface" covering SRH +
     Auger + radiative + tunneling.
   - New Section 10 "Statistics and ionization" covering Boltzmann
     vs FD and full vs incomplete ionization.
   - New Section 11 "Schottky contacts".
   - New Section 12 "Heterojunctions".
   - Each section cross-references the relevant ADR and benchmark.
4. `README.md`:
   - "What works today" list grows the bullets:
     `Field- and surface-dependent mobility (Caughey-Thomas, Lombardi)`,
     `Recombination: SRH, Auger, radiative, BTB, Hurkx`,
     `Boltzmann and Fermi-Dirac statistics; incomplete ionization`,
     `Ohmic, ideal-gate, and Schottky contacts`,
     `Heterojunctions with position-dependent chi and E_g (HEMT, HBT)`.
   - Add a one-line link from "Project status" to the new
     `docs/PHYSICS.md` sections.
5. `docs/schema/reference.md`: regenerate the schema reference for
   the bumped schema version. The new fields are: `physics.mobility`,
   `physics.statistics`, `physics.ionization`, the recombination
   list, the `schottky` contact type, region `composition`.
6. `CHANGELOG.md`: a single section for the cumulative DD-completeness
   release. Sub-bullets for each M16.x and M17.
7. New notebook `notebooks/06_dd_complete_tour.ipynb`: a single Colab
   notebook that exercises every new model on a small representative
   problem. Auto-generated from the new benchmarks where possible.

**Commit message:** `docs: close out M16.1-M16.7 and M17; full DD physics surface complete`

---

## Invariants checklist (re-verify before each commit, every phase)

- [ ] No em dashes anywhere in prose touched (PLAN.md invariant 8).
- [ ] Pure-Python core remains dolfinx-free. New modules `mobility`,
      `recombination`, `statistics`, `ionization`, `tunneling` import
      only stdlib and numpy at module level (PLAN.md invariant 4 /
      ADR 0007 / Layer 3).
- [ ] Slotboom variables remain the only DD primary unknowns (ADR
      0004 / 0014). No SUPG, no streamline diffusion.
- [ ] `make_scaling_from_config` still on every solve path (M1-locked).
- [ ] Existing benchmarks continue to pass byte-identically with
      default model selections (`mobility.model: "constant"`,
      `statistics.model: "boltzmann"`, etc.). New models are opt-in.
- [ ] Schema bumped per minor when fields are added; major-gate
      logic still passes (P4 / M11). Coordinate every phase's bump
      so they don't collide; if multiple phases ship in parallel,
      one wins the minor bump and the others rebase onto it.
- [ ] Every new physics model has an MMS verifier per ADR 0006.
- [ ] No new ADR introduced unless a locked invariant is broken; if
      so, draft an ADR in `docs/adr/` first and stop.
- [ ] No PETSc / UFL types leak into the pure-Python modules.
- [ ] `docs/PHYSICS.md` section count grows monotonically; do not
      renumber existing sections.

## Anti-goals

- Do not batch sub-phases into a single PR. The dev guide is explicit:
  one PR per model.
- Do not reformulate Slotboom as a different primary-variable choice
  to "make FD easier". ADR 0004 and ADR 0014 are locked. ADR 0015
  describes how Slotboom and FD coexist; that is the only allowed
  path.
- Do not introduce SUPG or streamline-diffusion stabilization at any
  point. Slotboom is the stabilization mechanism (ADR 0004).
- Do not add field-dependent diffusivity that contradicts the
  Einstein relation; CT modifies mu, and the corresponding D is
  derived from the Einstein relation (D = (kT/q) mu) in the
  non-degenerate regime. In the degenerate regime under FD, use the
  generalized Einstein relation; document it in ADR 0015.
- Do not hard-code material parameters in the FEM kernels. Every
  numerical constant lives in `semi/materials.py` and is loaded via
  the existing material-card path.
- Do not weaken the V&V acceptance gates. Every benchmark gets a
  closed-form or textbook reference; if no reference is available,
  do not ship the model in this milestone.
- Do not skip the MMS verifier for any new model "because the
  benchmark passed". The benchmark proves you reproduce a known
  device; the MMS proves the kernel converges at the right rate.
  Both are required (ADR 0006).
- Do not claim COMSOL parity in marketing. The README's "What works
  today" list grows; the project tagline does not change.

## Working method

Per the dev guide:

1. Read the required reading list. Do not skim.
2. Pick exactly one M16.x sub-phase. Do not start a second until the
   first is merged.
3. Write the acceptance test first (the benchmark verifier or the
   MMS check). Make it fail. Then make it pass.
4. After each phase: `pytest tests/ -v`,
   `python tests/check_analytical_math.py` (renamed from
   `check_day1_math.py`), `python scripts/run_verification.py all`,
   `ruff check semi/ tests/`. Do not advance with red tests.
5. Re-run all five existing shipped benchmarks (`pn_1d`,
   `pn_1d_bias`, `pn_1d_bias_reverse`, `mos_2d`, `resistor_3d`,
   plus any added since) and confirm byte-identity to the pre-phase
   baseline. Default model selections must preserve historical
   behavior.

## Hand-off

When all phases land, the obvious next picks are:

- M14.2.x Cartesian-2D MOSCAP variant and gate-driven HF C-V (the
  remaining items from the AC backlog).
- 3D MOSFET / FinFET geometry. The DD physics is now complete; the
  remaining gap is geometry and gmsh-driven mesh authoring at scale.
- Validation suite Phase 2, which is now substantially unblocked
  since the engine has the models that Sze and Nicollian-Brews
  reference.
- M18 UI work, in the separate frontend repo.

Pick exactly one per the project rule; surface the choice to the
owner.
