## Context

You are working in `kronos-semi` at `v0.23.0` (post-examples
catalogue PR #85, commit `3a23383`). Schema is `2.7.0`,
package version `0.23.0`. The repo is at
`https://github.com/rwalkerlewis/kronos-semi`.

Your assignment is **the examples catalogue extension** — three
more examples added to the existing `examples/` directory. This
is **not** an M-numbered milestone; it is a continuation of the
documentation / demonstration layer that PR #85 established.

The first examples PR shipped:

- `examples/nmos_idvgs/` (NMOS Id-Vgs at four V_DS values;
  Caughey-Thomas + Lombardi + FD)
- `examples/schottky_iv_temperature/` (Pt-on-n-Si at three
  temperatures; thermionic emission)
- `examples/power_diode_reverse_recovery/` (rectifier turn-off
  transient; voltage_t + Auger)

This PR adds three more, deliberately chosen to widen the
runner / output-mode coverage without re-treading the first
three:

- `examples/pmos_idvgs/` — the complement to NMOS Id-Vgs.
  Same M16.1 + M16.2 features, opposite-polarity doping, the
  band-edge convention flip. Demonstrates the
  symmetry between n- and p-channel MOSFETs.
- `examples/moscap_cv_oxide_thickness/` — MOSCAP C-V at oxide
  thickness 2 nm, 5 nm, 10 nm. Exercises the `mos_cv` runner
  (untouched by the first three examples) and a different
  output mode (capacitance, not current).
- `examples/diode_reverse_leakage_temperature/` — pn diode
  reverse leakage at three temperatures. Extends the
  temperature-dependence theme from `schottky_iv_temperature`
  to a different physics mechanism: SRH generation in the
  depletion region scales as exp(-E_g/(2 kT)), distinct from
  the thermionic emission's exp(-q phi_B / kT).

Three "future examples" candidates from the first PR's
`examples/README.md` are explicitly **not** in scope here and
remain in the future-examples list:

- **NPN BJT Gummel plot.** Three-terminal infrastructure
  exists but the parameter tuning to get the Gummel plot
  qualitatively right is a multi-day task.
- **SiC Schottky power diode.** Requires SiC material
  parameters (E_g = 3.26 eV, m*, dielectric constant) which
  are not in `semi/materials.py`.
- **Tunnel diode forward I-V.** Forward-bias BBT under M16.6
  has not been independently validated; ship a tunnel-diode
  example only after a future PR confirms forward-bias Kane
  tunneling is captured.

This prompt **is** the starter. Save it as
`docs/EXAMPLES_EXTENSION_STARTER_PROMPT.md` as the first commit
of the PR (Phase 0 below). The descriptive (non-M-numbered)
filename matches the precedent set by
`docs/EXAMPLES_CATALOGUE_STARTER_PROMPT.md` from PR #85.

## Branch and PR rules

- Work on a fresh branch off `main`:
  `git checkout -b dev/examples-extension`. Not an M-numbered
  branch; parallel to the M17/M19 next-task decision.
- One scope, one PR. Do not bundle with M17 / M19 work.
- Push every phase commit immediately. Open the PR after
  Phase 0 so reviewers can watch examples land.
- **Run all five phases consecutively without pausing for
  confirmation, status reports, or "should I continue"
  prompts.** Do not call `ask_user_input` or any equivalent
  between phases. Run Phase 0 through Phase D end to end.
  The only legitimate reasons to stop are: (a) a hard blocker
  that re-reading this prompt and the codebase cannot
  resolve; (b) a `pytest` or `ruff` failure that persists
  after a reasonable repair attempt; (c) the Stop conditions
  section reports green and the PR is ready for review.
  "Should I move on to the next example" is not a legitimate
  reason. Default behavior: pick the most defensible
  interpretation, document inline, continue.
- Title the PR `examples: pmos, moscap C-V, reverse leakage
  temperature`. Use the PR description template at the bottom
  of this prompt.
- **No AI-assistant credits in any shipped artifact.** The
  `includeCoAuthoredBy: false` settings fix is presumed in
  effect; verify with
  `git log --format=%B dev/examples-extension | grep -i
  'co-authored-by: claude'` before opening the PR.

## Required reading (in order; ~25 minutes)

1. `PLAN.md` "Current state" and "Next task". Confirm `main`
   is at v0.23.0 and that no other PR is in flight on
   `dev/examples-*`. The "Next task" is still the M17/M19
   choice; this PR does not change that.
2. `docs/EXAMPLES_CATALOGUE_STARTER_PROMPT.md` end to end.
   That document established the conventions this PR inherits:
   the README structure, smoke-verifier rule, no-tight-tolerances
   rule, no-AI-credit rule, CI matrix entry rule, and the
   five-phase shape. **Do not re-state those conventions in
   this PR's commits.** Refer back to PR #85's prompt; this
   PR is a continuation, not a re-establishment.
3. `examples/README.md` (current state, post PR #85). The
   table-of-three becomes table-of-six in Phase 0.
4. The three shipped examples for stylistic reference:
   - `examples/nmos_idvgs/` (the structural sibling for PMOS)
   - `examples/schottky_iv_temperature/` (the structural
     sibling for the diode reverse-leakage temperature
     example)
   - `examples/power_diode_reverse_recovery/` (the structural
     sibling for any future transient-mode example, not in
     scope here)
5. `benchmarks/mos_2d/` and the `mos_cv` runner
   (`semi/runners/mos_cv.py`). The MOSCAP C-V example uses
   this runner; existing benchmark is the structural template.
6. `benchmarks/pn_1d_bias_reverse/` (the closest existing
   reverse-bias config; structural template for the diode
   reverse-leakage example).
7. `scripts/run_benchmark.py` near the existing `verify_*`
   registrations. Three new entries in this PR.

## Conventions inherited from PR #85

The conventions block from
`docs/EXAMPLES_CATALOGUE_STARTER_PROMPT.md` applies verbatim.
Specifically: no schema change, no new physics modules,
README structure (load-bearing five sections), smoke-test
verifiers (not V&V gates), CI matrix entries (no
allow-failure), naming
(`examples/<device_class>_<test_type>/`), one example per
commit, no em dashes, no AI-credit, coverage gate at 95 in
`docker-fem-tests`, existing-benchmark byte-identity to v0.23.0.

The only convention worth restating here, because it has bitten
the M16.x prompts repeatedly: **smoke-test verifiers, not V&V
gates.** A reverse-leakage example at three temperatures is a
pointed temptation to assert a specific E_g extraction from the
slope of log(I_R) vs 1/T; do not assert that. The example
demonstrates the workflow; quantitative E_g extraction is the
user's exercise, not the verifier's. Same caution applies to
PMOS V_T extraction and MOSCAP C_ox-vs-thickness scaling.

## Phases, one commit per phase

Do not bundle. After each phase, run `ruff check` (Python
touched only — the example JSONs are pure config) and
`pytest tests/`. Do not advance with red tests. **Do not pause
between phases for confirmation.**

---

### Phase 0: ship this starter prompt and update the
examples directory README

1. Save this entire file as
   `docs/EXAMPLES_EXTENSION_STARTER_PROMPT.md`. Strip nothing.
2. Update `examples/README.md`:
   - Expand the table of examples from three rows to six
     (the three from PR #85 plus the three from this PR).
     Keep the existing column layout (name / device class /
     runner / M16 features).
   - Update the "future examples" section: remove the three
     candidates this PR ships (none of the existing future-
     examples list overlaps with this PR; the existing four
     stay as they were).
3. Append a one-line "Author examples extension starter
   prompt" entry to `docs/IMPROVEMENT_GUIDE.md` § 9 changelog
   under `[Unreleased]`.
4. Push the branch and open the PR.

**Commit message:** `docs: ship examples extension starter prompt; expand examples README`

**Acceptance:** the file exists at
`docs/EXAMPLES_EXTENSION_STARTER_PROMPT.md` on
`dev/examples-extension`; the README table has the new
three-row block; CI green on docs-only diff.

---

### Phase A: `examples/pmos_idvgs/` — PMOS Id-Vgs

The complement to NMOS Id-Vgs from PR #85. Same channel length
(100 nm), same gate oxide (2 nm), opposite-polarity doping
(n-type body, p+ source/drain), same Lombardi + Caughey-Thomas
+ FD physics. Voltage sweeps run over negative V_GS (PMOS
inversion is at V_GS < V_T_p, V_T_p < 0).

This example is deliberately simple. The intent is "show the
NMOS-PMOS symmetry"; users adapt the NMOS config to PMOS by
flipping doping signs and band-edge conventions, and the
example makes the precise edits visible.

1. `examples/pmos_idvgs/pmos.json`. Geometry mirrors NMOS:
   2D MOSFET, channel length 100 nm, gate oxide 2 nm.
   N-type body N_D = 5e17 cm^-3 (sign-flipped from the NMOS
   N_A = 5e17). P+ source/drain N_A = 1e20 cm^-3 (degenerate;
   FD applies; `physics.statistics: "fermi_dirac"`).
   `physics.mobility = {model: "lombardi", bulk_model:
   "caughey_thomas", interface_facet_tag: <oxide-silicon
   facet tag>}`. Voltage sweep: parametric over V_DS in
   [-0.05, -0.5, -1.0, -1.8] V with V_GS sweep [-1.8, 0] V
   step 0.05 V at each V_DS.
2. `examples/pmos_idvgs/README.md` per the README structure
   convention. The "What this is" section explicitly notes
   the relationship to `examples/nmos_idvgs/` and walks the
   user through the four edits required to convert NMOS to
   PMOS (doping polarity, source/drain doping polarity,
   V_GS sign, V_DS sign). The "Expected output" section gives
   an order-of-magnitude estimate for |I_D| at V_GS = -1.8 V,
   V_DS = -1.8 V using the textbook square-law formula
   adjusted for hole mobility (mu_p ~ 0.4 mu_n in Si, so
   |I_D_pmos| at the same |V| is roughly 40 % of I_D_nmos).
   Pointer back at `examples/nmos_idvgs/` for direct
   comparison.
3. Verifier `verify_pmos_idvgs` in `scripts/run_benchmark.py`:
   - Run completed (no exception, no NaN in output).
   - For every V_DS, |I_D| is monotonically non-decreasing in
     |V_GS - V_T_p| (the transistor turns on as V_GS goes
     more negative).
   - I_D at V_GS = -1.8 V, V_DS = -1.8 V is **negative** (PMOS
     conducts source-to-drain in the negative direction; the
     verifier asserts the sign as well as finiteness).
4. Plotter `plot_pmos_idvgs`: four |I_D| vs |V_GS| curves
   on a semilog-y plot (subthreshold visible) plus the same
   curves on a linear-y plot (saturation visible). Two-panel
   layout. Mirror `plot_nmos_idvgs`.
5. CI matrix entry near `nmos_idvgs`. Estimate the runtime: a
   PMOS solve is comparable to an NMOS solve at the same
   geometry (the asymmetry is in mobility and the band-edge
   offset, neither of which materially changes solver
   conditioning). If the NMOS runtime in CI was ~2 minutes,
   PMOS will be similar.

**Commit message:** `examples: pmos_idvgs (PMOS complement to nmos_idvgs)`

**Acceptance:** smoke verifier passes; CI green.

---

### Phase B: `examples/moscap_cv_oxide_thickness/` — MOSCAP C-V at three oxide thicknesses

A practical MOSCAP C-V characterization across oxide thickness
(2 nm, 5 nm, 10 nm). Exercises the `mos_cv` runner — which
none of the first three examples touched — and demonstrates
a different output mode (capacitance vs voltage, not current).

1. `examples/moscap_cv_oxide_thickness/moscap.json`. Geometry
   is a 2D axisymmetric MOSCAP (template:
   `benchmarks/moscap_axisym_2d/`) on a p-type silicon body
   N_A = 5e16 cm^-3. SiO2 gate oxide. Body radius 50 um, body
   depth 5 um; oxide thickness is the parametric knob.
   The JSON layout mirrors the schottky_iv_temperature
   pattern: either three sibling JSONs
   (`moscap_2nm.json`, `moscap_5nm.json`,
   `moscap_10nm.json`) or one parametric config that the
   verifier loops over. Match whichever shape
   `schottky_iv_temperature` shipped (consistency).
2. `solver.type: "mos_cv"` with V_GS sweep [-2, 2] V step
   0.05 V (81 points; covers accumulation through inversion
   for a typical Si MOSCAP).
   `physics.statistics: "boltzmann"` is fine (the
   bulk doping 5e16 is below the FD-required threshold).
   `physics.mobility: {model: "constant"}` (mobility does not
   enter the C-V curve under the existing mos_cv extraction;
   keep it minimal).
3. `examples/moscap_cv_oxide_thickness/README.md` per the
   README structure convention. The "Expected output" section:
   the accumulation capacitance per unit area is
   C_ox = eps_0 eps_ox / t_ox, so doubling t_ox from 2 nm to
   4 nm halves C_ox; the user can verify this scaling
   directly from the three curves at strong accumulation.
   Mid-gap V_T extraction guide using the existing C-V curve
   inflection point. Pointer at `benchmarks/mos_2d/` for the
   formal V&V version of the C-V analysis.
4. Verifier `verify_moscap_cv_oxide_thickness` in
   `scripts/run_benchmark.py`:
   - All three runs completed.
   - For each t_ox, C-V curve is well-formed (finite,
     non-NaN, has at least one extremum).
   - In strong accumulation (V_GS = -2 V, the swept
     extremum), C(2 nm) > C(5 nm) > C(10 nm) (thinner oxide
     gives larger accumulation capacitance; smoke check).
5. Plotter `plot_moscap_cv_oxide_thickness`: three C-V
   curves overlaid, color-coded by oxide thickness. Linear-y
   axis (C is positive throughout the sweep).
6. CI matrix entry. MOSCAP C-V is moderately expensive;
   estimate ~1-3 minutes per oxide thickness on a typical
   runner. Three runs total, so ~3-9 minutes wall-clock; if
   that exceeds the runner budget, drop the V_GS step from
   0.05 V to 0.1 V (40 points) or drop one of the oxide
   thicknesses.

**Commit message:** `examples: moscap_cv_oxide_thickness (mos_cv runner; three oxide thicknesses)`

**Acceptance:** smoke verifier passes; CI green.

---

### Phase C: `examples/diode_reverse_leakage_temperature/` — pn diode reverse leakage at three temperatures

A practical pn diode reverse-leakage characterization across
temperature (250 K, 300 K, 350 K). The temperature dependence
of SRH generation in the depletion region scales as
exp(-E_g/(2 kT)) at moderate temperatures (generation-
limited regime); this is distinct from the
exp(-q phi_B / kT) thermionic-emission scaling
demonstrated by `schottky_iv_temperature`. Educational
complement: two examples show two distinct temperature
dependences for two distinct physics mechanisms, side by
side.

1. `examples/diode_reverse_leakage_temperature/diode.json`.
   Geometry: 1D pn junction, 20 um total, N_A = 1e16 cm^-3,
   N_D = 1e16 cm^-3 (moderate doping; the depletion region
   widens enough at reverse bias that SRH generation
   dominates). 800 cells. tau = 1e-7 s (a typical Si
   defect-limited lifetime). V_F sweep [-5, 0] V step 0.1 V
   (51 points; covers the reverse-leakage regime). Three
   sibling JSONs (`diode_250K.json`, `diode_300K.json`,
   `diode_350K.json`) or one parametric config; match the
   `schottky_iv_temperature` shape.
2. `physics.recombination: {srh: true, auger: false}` (SRH is
   the dominant generation mechanism in the depletion region
   at the doping levels used here; Auger does not contribute
   meaningfully at low injection). `physics.mobility:
   {model: "constant"}` (mobility does not affect the
   leakage current substantially in this regime; keep it
   minimal).
3. `examples/diode_reverse_leakage_temperature/README.md` per
   the README structure convention. The "What this is" section
   explicitly relates this to `schottky_iv_temperature`: same
   workflow shape (three-temperature sweep), different physics
   mechanism (SRH generation in the depletion region instead of
   thermionic emission). The "Expected output" section: at
   reverse bias V_F = -3 V, |I_R| at T = 350 K should be
   roughly 30-50 x larger than |I_R| at T = 250 K (the
   exp(-E_g/(2 kT)) scaling with E_g = 1.12 eV gives a factor
   of ~50 between those two temperatures; the scaling is
   approximate because the depletion region also widens with
   reverse bias and there are pre-factor temperature
   dependences). Activation-energy extraction guide pointing
   at the slope of log(|I_R|) vs 1/T at fixed V_F.
4. Verifier `verify_diode_reverse_leakage_temperature` in
   `scripts/run_benchmark.py`:
   - All three runs completed.
   - For each T, the I-V curve is finite, non-NaN, and the
     reverse-bias side has |I| monotonically non-decreasing
     in |V_F| (depletion-region widens with reverse bias;
     leakage grows correspondingly).
   - At V_F = -3 V, |I_R(350 K)| > |I_R(300 K)| > |I_R(250 K)|
     monotonically (the temperature-dependent SRH generation
     ordering must be physically correct).
   - Resist the temptation to assert a specific
     `|I_R(350K)/I_R(250K)|` ratio: the ratio depends on
     pre-factor T-dependences and the depletion-region
     width, both of which are example-knob-dependent.
5. Plotter `plot_diode_reverse_leakage_temperature`: three
   semilog-y |I_R| vs V_F curves overlaid, color-coded by
   temperature. Companion subplot of log(|I_R|) vs 1/T at
   V_F = -3 V (the activation-energy extraction visual).
6. CI matrix entry. Three small 1D runs are cheap; no
   runtime tuning needed.

**Commit message:** `examples: diode_reverse_leakage_temperature (SRH generation T-dependence)`

**Acceptance:** smoke verifier passes; CI green.

---

### Phase D: closeout

1. `examples/README.md`: re-confirm the table-of-six matches
   what shipped (Phases A-C added three rows). The
   "future examples" section is unchanged from PR #85 (this
   PR did not consume any future-examples candidates; the
   three shipped here were not previously listed because they
   are obvious-and-cheap follow-ons rather than candidates
   that needed deferred-scope discussion).
2. `tests/test_examples_register.py`: append three pytests
   asserting the new verifier names are registered in
   `_VERIFIERS`. Mirror the existing three pytests verbatim.
3. `PLAN.md` "Backlog": the existing one-line entry from
   PR #85 (the examples catalogue) gets a sub-bullet noting
   the three additions in this PR. The "Next task" field
   is unchanged (still M17 or M19, maintainer's choice).
4. `docs/IMPROVEMENT_GUIDE.md` § 9 changelog: append an
   "examples extension: pmos_idvgs, moscap_cv_oxide_thickness,
   diode_reverse_leakage_temperature" entry under
   `[Unreleased]`. Move the existing `[Unreleased]` block (the
   Phase 0 entry from this PR plus whatever was already there
   from PR #85's closeout) into `### Released` as part of the
   closeout.
5. `CHANGELOG.md`: this PR does not bump the package version
   (no API or schema change). Add an entry under `[Unreleased]`
   noting the three new examples; the next physics PR will
   roll this into its own bump.
6. `docs/PHYSICS_INTRO.md` § 6: the existing
   "practical examples" bullet from PR #85 already points
   at `examples/`. No edit needed unless the bullet's example
   list is exhaustive (in which case append the three new
   examples). If the bullet says "see `examples/` for
   self-contained device configs (NMOS Id-Vgs, Schottky
   temperature dependence, power diode reverse recovery)
   demonstrating M16 features", expand it to "see `examples/`
   for six self-contained device configs ... demonstrating
   M16 features."

**Commit message:** `docs: close out examples extension (PHYSICS_INTRO, CHANGELOG, IMPROVEMENT_GUIDE)`

---

## Invariants checklist (re-verify before each commit)

- [ ] No em dashes in any new prose or code comment.
- [ ] No mention of Claude, Claude Code, Anthropic, or any AI
      assistant in any shipped artifact.
- [ ] `git log --format=%B dev/examples-extension | grep -i
      'co-authored-by: claude'` returns empty.
- [ ] No schema change (v2.7.0 stays current).
- [ ] No `semi/` source change. Only `scripts/run_benchmark.py`
      (verifier registrations + plot helpers) and
      `tests/test_examples_register.py` are Python-touched.
- [ ] Every existing benchmark and every example from PR #85
      is bit-identical to v0.23.0.
- [ ] Coverage gate holds at 95 in the gated `docker-fem-tests`
      job. The new pytests in `tests/test_examples_register.py`
      add three lines each; the coverage delta is bounded.
- [ ] Each example has a `README.md` matching the load-bearing
      structure convention from PR #85.
- [ ] Each example has a CI matrix entry; none is marked
      `allow-failure: "true"`.
- [ ] Smoke-test verifiers, not V&V gates. Especially: do
      not assert E_g extraction in
      diode_reverse_leakage_temperature; do not assert exact
      C_ox-vs-thickness ratios in moscap_cv_oxide_thickness;
      do not assert exact V_T_p in pmos_idvgs.

## Anti-goals

- Do not bump the schema. Examples use only v2.7.0 features.
- Do not start M17 (heterojunctions) or M19 (3D MOSFET) in
  this PR. The "Next task" pointer is unchanged.
- Do not add a fourth example. Three is enough; the
  future-examples list (NPN BJT, NMOS C-V at multiple body
  biases, SiC Schottky, tunnel diode forward I-V) stays in
  the README for the next examples PR.
- Do not assert tight numerical correctness in the example
  verifiers. Smoke checks only; restated explicitly because
  the temptation is acute on a temperature-dependence
  example (Arrhenius extraction is irresistible) and a
  C-V-vs-oxide example (C_ox = eps/t scaling is irresistible).
  The user's exercise, not the verifier's.
- Do not refactor `scripts/run_benchmark.py` to share
  infrastructure between benchmarks and examples. PR #85
  established the convention that examples register through
  the same dispatcher but live in `examples/` rather than
  `benchmarks/`; do not deepen that coupling.
- Do not add notebooks. PR #85 explicitly anti-goaled them
  for examples; this PR inherits that decision.
- Do not retire the `mosfet_2d` `allow-failure: "true"` flag.
- Do not add `Co-Authored-By` trailers, generated-by footers,
  or any other AI-credit marker.

## Stop conditions

You are done when:

1. The PR is opened on branch `dev/examples-extension`.
2. All three new example smoke verifiers pass in CI:
   - `examples/pmos_idvgs/`
   - `examples/moscap_cv_oxide_thickness/`
   - `examples/diode_reverse_leakage_temperature/`
3. `tests/test_examples_register.py` passes in the gated
   `docker-fem-tests` job; coverage gate holds at 95 without
   a follow-up commit.
4. `examples/README.md` table-of-examples matches what
   shipped (six rows: three from PR #85, three from this PR).
5. PLAN.md, IMPROVEMENT_GUIDE.md (§ 9), PHYSICS_INTRO.md
   (§ 6), CHANGELOG.md reflect the examples extension ship.
6. No commit, PR description, code comment, or doc in this PR
   mentions Claude, Claude Code, Anthropic, or any AI
   assistant. Squash-merge body is also clean.
7. Every existing benchmark and every example from PR #85
   remains bit-identical to v0.23.0 (cite anchors in the PR
   description: pn_1d_bias J(V=0.6 V) = 1.635e+03 A/m^2;
   diode_velsat_1d 56.27 % @ 0.9 V, 0.19 % @ 0.3 V;
   schottky_1d worst-case <10 %; zener_1d worst-case <20 %
   Kane match; nmos_idvgs / schottky_iv_temperature /
   power_diode_reverse_recovery smoke verifiers green).
8. PR reviewed, CI green (modulo the documented `allow-
   failure` on `mosfet_2d`, unchanged in scope from v0.23.0),
   merged.

## PR description template

```
## Summary

Examples catalogue extension: three more practical device
configs added to the `examples/` directory established by
PR #85. Distinct from `benchmarks/` (which gates V&V); these
are illustrative clone-and-modify starting points for users.

Three examples ship in this PR:

- `examples/pmos_idvgs/` — PMOS Id-Vgs at four V_DS values,
  the natural complement to nmos_idvgs from PR #85.
  Demonstrates the NMOS-PMOS symmetry. Caughey-Thomas (M16.1)
  + Lombardi (M16.2) + FD (M16.4) for the p+ source/drain.
- `examples/moscap_cv_oxide_thickness/` — MOSCAP C-V at
  oxide thickness 2 nm, 5 nm, 10 nm. Exercises the mos_cv
  runner (untouched by PR #85 examples) and demonstrates
  the C-V output mode plus the C_ox = eps/t_ox scaling.
- `examples/diode_reverse_leakage_temperature/` — pn diode
  reverse leakage at T = 250, 300, 350 K. Educationally
  complements schottky_iv_temperature: SRH generation in the
  depletion region scales as exp(-E_g/(2 kT)), distinct from
  the thermionic emission's exp(-q phi_B / kT).

No schema bump (v2.7.0 stays current). No `semi/` source
change beyond verifier registrations and plot helpers in
`scripts/run_benchmark.py`. Every existing benchmark and
every example from PR #85 is bit-identical to v0.23.0.

## Acceptance

- [ ] examples/pmos_idvgs/ smoke verifier passes
- [ ] examples/moscap_cv_oxide_thickness/ smoke verifier passes
- [ ] examples/diode_reverse_leakage_temperature/ smoke
      verifier passes
- [ ] Existing benchmarks and PR #85 examples bit-identical
      to v0.23.0
      (anchors: pn_1d_bias J(V=0.6 V) = 1.635e+03 A/m^2;
      diode_velsat_1d 56.27 % @ 0.9 V, 0.19 % @ 0.3 V;
      schottky_1d worst-case <10 %; zener_1d worst-case
      <20 % Kane match)

## Test plan

- [ ] ruff check semi/ tests/
- [ ] pytest tests/
- [ ] pytest --cov=semi --cov-fail-under=95 (in the gated
      docker-fem-tests job)
- [ ] docker compose run --rm benchmark pmos_idvgs
- [ ] docker compose run --rm benchmark moscap_cv_oxide_thickness
- [ ] docker compose run --rm benchmark diode_reverse_leakage_temperature
- [ ] docker compose run --rm benchmark nmos_idvgs
      (PR #85 example, smoke-test bit-identical)
- [ ] docker compose run --rm benchmark pn_1d_bias
      (existing benchmark, bit-identical)

## Notes

- This PR does not advance the M-numbered roadmap. Next-task
  pointer (M17 heterojunctions or M19 3D MOSFET capstone) is
  unchanged.
- Future examples candidates (NPN BJT, NMOS C-V at multiple
  body biases, SiC Schottky, tunnel diode forward I-V) remain
  in the `examples/README.md` future-examples list. None of
  them is consumed by this PR.
- No package version bump. The next physics PR (M17 or M19)
  will roll the docs-only changes from this PR and PR #85
  into its own version.
```

## Hand-off

When this PR lands, the maintainer's M17 / M19 choice is the
remaining open question. The examples catalogue can keep
growing on dev/examples-* branches in parallel; the next
contributor picking up an examples PR adds entries from the
README "future examples" list (or proposes new ones) and
follows the conventions established in
`docs/EXAMPLES_CATALOGUE_STARTER_PROMPT.md` and this prompt.

The two prompts together (`EXAMPLES_CATALOGUE_STARTER_PROMPT.md`
and `EXAMPLES_EXTENSION_STARTER_PROMPT.md`) are the durable
reference for "how do I add an example to this repo." Future
examples-PR starter prompts can shrink to a delta-against-PR-#85
shape rather than restating conventions.
