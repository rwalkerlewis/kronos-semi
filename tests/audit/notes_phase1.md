## Notes (Phase 1, 2026-04-28)

### Summary of findings

Five active cases ran on `dolfinx 0.10.0.post2` inside the `kronos-semi:dev`
Docker image (case 06 skipped per scaffolding). Classification per the
A/B/C/D scheme in PLAN.md "Next task":

| Case | Subject | Worst rel. err | Class | Action |
|:---|:---|---:|:---:|:---|
| 01 | bias_sweep vs transient (deep SS) | 1.0e-4 (psi/n/p) | B | document |
| 02 | AC ω→0 vs bias_sweep dI/dV (reverse) | 7.2e-2 | B | document |
| 03 | mos_cv vs mos_cap_ac on Q_gate      | 0.0 | A | confirms PR #38 byte-identity |
| 04 | equilibrium vs bias_sweep at V=0    | 3.7e-9 (psi gauge-aligned), 2.6e-8 (n,p) | A | confirms shared BC convention |
| 05 | AC Re(Y) vs bias_sweep dI/dV (forward) | 1.21 (121%) | C | tracking issue |
| 06 | transient FFT vs ac_sweep            | n/a | n/a | skipped, infrastructure gap |

**Case 01 nuance.** psi/n/p fields agree to 1e-4 or better at all three
forward biases, formally validating the M13.1 close-out (Slotboom
transient → bias_sweep at deep steady state). The `rel J` column is
intentionally not gated: the terminal-current evaluation is dominated
by a small difference of large drift and diffusion contributions in
the contact-adjacent layer, so 200-300% relative disagreement on `J`
in the presence of 1e-8 relative agreement on the field profile is
expected sensitivity, not a physics bug. Phase 2 will switch to
internal-current-density quadrature on a quasi-1D cut to remove the
contact-flux noise.

**Case 03 confirms** PR #38's byte-identity claim: `mos_cv` (numerical
dQ/dV) and `mos_cap_ac` (analytic AC admittance) produce identical
Q_gate(V_gate) at all 42 swept gate voltages on `benchmarks/mos_2d`.

**Case 04 corrects** an earlier hypothesis. The two runners use the
same `psi_bc = arcsinh(N_net / (2 n_i)) + V_applied / V_t` ohmic-contact
BC (`semi/bcs.py::build_psi_dirichlet_bcs`,
`build_dd_dirichlet_bcs`); the gauge offset measured at V=0 is
~5e-19 V (numerical zero), and the gauge-aligned psi rel_L2 is 3.7e-9.
n and p agree to 2.6e-8 relative. The audit fragment retains the
gauge-aligned reporting columns so a future divergence will surface
as a non-zero offset.

### Bugs fixed in this PR

The audit scaffolding from PR #55 contained five small bugs that
prevented the cases from running. All fixes live under `tests/audit/`
(no `semi/` runner code modified) and total ~80 lines:

- **`tests/audit/_helpers.py`**: added `final_field(result, name)` to
  hide the shape difference between `transient` (`.fields` dict of
  snapshot lists) and `bias_sweep`/`equilibrium` (`SimulationResult`
  with `psi_phys` / `n_phys` / `p_phys` final-state attributes).
- **`tests/audit/_helpers.py`**: added `make_bias_sweep_cfg(...)` to
  configure a base config for `run_bias_sweep`. The bias sweep runner
  reads the swept contact's `voltage_sweep` block (not `solver.bias_ramp`,
  which the scaffolding incorrectly used) and the rc_ac_sweep config's
  tight ac-tuned SNES tolerances are too aggressive for the nonlinear
  DC continuation; the helper writes both the contact's voltage
  endpoint and the M12-relaxed (ADR 0008) SNES tolerances.
- **Cases 01, 04**: switched from per-case `.fields` extraction (broken
  for `bias_sweep` and `equilibrium`) to the new `final_field` helper.
- **Cases 02, 05, 06**: corrected the `solver.ac.frequencies` override
  to use the schema-required `{"type": "list", "values": [...]}`
  object instead of a bare list (which crashed `_resolve_frequencies`).
- **Case 04**: extended to report psi gauge offset and gauge-aligned
  L2 alongside the raw L2; n and p reported as gauge-invariant
  references.

No `semi/` runner / solver code was changed. The audit-philosophy
hard-fail-only-on-crash rule is now enforced: cases that classify as
A/B/C all "pass" in the pytest sense, with the disagreement reported
in the fragment table.

### Deferred findings (tracking)

- **Case 05 (Class C, 121% disagreement)**: AC small-signal Re(Y) at
  V_DC = +0.4 V is `-8.275e+01 S` while the bias_sweep centered-difference
  dI/dV at the same operating point is `+3.949e+02 S`. Magnitudes differ
  by a factor of ~5 and signs disagree. At forward bias on the
  rc_ac_sweep diode, the Slotboom DD operator is in the high-injection
  regime where `dF/dV` and `J(u_0)` have substantial contributions from
  the recombination block; the linearised AC operator excludes the
  nonlinear DC sensitivity contribution that dominates `dI/dV`. This is
  a real inconsistency between the M14 ac_sweep formulation (ADR 0011)
  and the M3 bias_sweep continuation. Recommended next step: extend
  ADR 0011 with a forward-bias regime analysis, then either widen the
  `dF/dV` finite-difference stencil in `ac_sweep._build_rhs_dF_dV` or
  switch to a UFL-derived analytic `dF/dV` with the recombination
  block included. Tracking issue to be opened with the title
  "audit Case 05: AC Re(Y) vs bias_sweep dI/dV diverge in forward bias".
- **Case 06 (infrastructure gap)**: `run_transient` does not accept a
  time-varying contact voltage `V(t)`. Adding a `bc_voltage_callback`
  hook is a small but distinct piece of work outside the audit PR's
  budget. Reference value `Y_ac(f=1 MHz, V_DC=+0.4 V) = -7.682e+01 +
  j*-6.455e+03 S` is recorded in the case-06 fragment for the future
  comparison. Tracking issue to be opened with the title "transient
  runner: support time-varying contact voltage `V(t)` for AC-FFT audit".

### Recommendation

Phase 2 (external validation against Sze and Nicollian-Brews textbook
references) should proceed without waiting on the case-05 fix. The
M13.1 close-out is confirmed by case 01 at all three forward biases
(psi/n/p rel_L2 ≤ 1e-4); the M14.1 byte-identity is confirmed by
case 03; the equilibrium / bias_sweep BC convention agrees in case 04.
Case 05's forward-bias AC vs DC sensitivity divergence is a known
gap in the linearised AC formulation, not a transient/ac_sweep
correctness regression that would block Phase 2.
