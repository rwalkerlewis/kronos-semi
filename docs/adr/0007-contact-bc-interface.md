# 0007. Contact BC interface and `semi/bcs.py` layer

- Status: Accepted
- Date: 2026-04-21

## Context

Through Day 4, ohmic contact BC construction lived inline in
`semi/run.py` as the private helpers `_build_ohmic_bcs_psi` and
`_build_dd_ohmic_bcs`. Both helpers read `cfg["contacts"]`, walked
`cfg["mesh"].get("facets_by_plane", [])` to resolve string facet
references to integer tags, evaluated the doping at each contact's
facet centroid, and built `dolfinx.fem.dirichletbc` objects. The bias
sweep called `_build_dd_ohmic_bcs` with a `voltages` dict that
overrode each contact's static `voltage` field at the current ramp
step.

That code worked but it had three issues:

1. **The schema-walk was duplicated.** Both helpers had identical
   `tag_by_name = {p["name"]: int(p["tag"]) for p in ...}` blocks
   and identical `if contact["type"] != "ohmic": continue` filters,
   plus identical raises on missing facet tags. Any future contact
   kind (gate, Schottky) would have to copy the same scaffolding.
2. **The BC builders mixed dolfinx-aware and pure-Python concerns.**
   The "is this contact recognized? does the facet exist? what is the
   applied voltage right now?" decisions are pure-Python (Invariant
   4-eligible) but were buried inside helpers that imported dolfinx
   at module scope.
3. **The Day 5 `run.py` split needed somewhere clean to host the
   BC code.** Pulling these helpers into the per-runner files would
   re-duplicate them between `equilibrium.py` and `bias_sweep.py`.

## Decision

Introduce `semi/bcs.py` as a new module in the pure-Python core tier
(no dolfinx import at module scope), with a two-step API:

1. `resolve_contacts(cfg, facet_tags=None, voltages=None) -> list[ContactBC]`
   walks the JSON config, validates each contact's kind, resolves
   string facet references via `mesh.facets_by_plane`, optionally
   verifies the resolved tag has facets when `facet_tags` is provided,
   and returns a list of pure-data `ContactBC` records. Bias-sweep
   drivers pass `voltages={contact_name: V_step}` to push the next
   ramp value through the same resolver without mutating the config.
2. `build_psi_dirichlet_bcs(...)` and `build_dd_dirichlet_bcs(...)`
   consume the resolved `list[ContactBC]` and a mesh / facet_tags /
   scaling / reference material, and produce the
   `dolfinx.fem.dirichletbc` lists for the equilibrium Poisson row
   and the coupled (psi, phi_n, phi_p) block respectively.

Specific design choices and their reasons:

- **Dataclass over dict for `ContactBC`.** A dict has no schema; a
  field typo (`facet_tag` vs. `facetTag`) silently produces wrong
  Dirichlet values without raising. The dataclass gives static names
  the type checker can see, makes equality and `repr` free, and is
  a cheap upgrade path to `slots=True` if profiling later shows
  resolver allocations matter.
- **`voltages` as an override dict, not a config mutation.** The
  bias sweep solves at one bias per call. Mutating
  `cfg["contacts"][i]["voltage"]` between calls works but couples
  the resolver to the caller's bookkeeping (every caller has to know
  to write the value back, every test has to remember to clean it
  up). A keyword-only `voltages` argument keeps the resolver
  side-effect-free, which is what the V&V tests rely on.
- **Insulating contacts produce no `ContactBC`.** They are a natural
  Neumann BC, so there is nothing for `build_*_dirichlet_bcs` to
  build. Returning a `ContactBC(kind="insulating")` would force
  every consumer to filter; skipping at resolution time is simpler
  and matches the existing `if contact["type"] != "ohmic": continue`
  semantics in the legacy helpers.
- **`gate` and `schottky` kinds resolve but are skipped by the
  builders.** The schema accepts `"gate"` (and we anticipate
  `"schottky"` for Day 6+); validating the kind at `resolve_contacts`
  time means a typo in `cfg["contacts"][i]["type"]` raises during
  config processing rather than during the builder's silent
  skip-loop. The builders themselves only handle ohmic for now;
  gate / Schottky get dedicated builders when their physics lands,
  with no further changes to the resolver.
- **Facet-tag verification is optional (default off).** Pure-Python
  callers (tests, schema-only validators, Phase 1 dry runs) cannot
  produce a real `dolfinx.mesh.MeshTags`. Allowing
  `facet_tags=None` for those callers and asserting non-emptiness
  when `facet_tags` is provided keeps both call paths honest
  without forcing the pure-Python tier to depend on dolfinx.

## Consequences

- `semi/bcs.py` is now the single owner of contact BC construction.
  Adding gate / Schottky support in Day 6+ means one new
  `build_gate_dirichlet_bcs` (or similar) here plus dispatch from
  the runner; no edits to `run.py`, `runners/equilibrium.py`, or
  `runners/bias_sweep.py`.
- `tests/test_bcs.py` covers the resolver entirely on the pure-Python
  CI matrix; `tests/fem/test_bcs.py` covers the dolfinx builders by
  comparison against an embedded copy of the legacy inline code so
  the byte-identical numerics are gated even after the legacy code
  was deleted.
- The pure-Python core tier now has one more module
  (`constants`, `materials`, `scaling`, `doping`, `schema`,
  `continuation`, `diode_analytical`, `bcs`).
  The Invariant 4 statement in PLAN.md is updated to list it.
- The `ContactBC` dataclass is a public type. Future ADRs that change
  the resolver's contract (for example, adding workfunction defaults
  or per-region BC overrides) need to consider backward compat for
  external callers that may instantiate it directly.

## Alternatives considered

- **Inline the helpers in each runner file.** Rejected because it
  re-duplicates the schema walk and prevents Invariant 4 cleanup.
- **Build BCs directly inside `runners/_common.py`.** Rejected
  because `_common.py` is for tiny shared helpers; mixing dolfinx
  imports there would force `_common.py` out of the pure-Python
  tier.
- **Use a dict instead of `ContactBC`.** Rejected: see above.
- **Mutate `cfg["contacts"]` for bias-step voltages.** Rejected: see
  above.
