# §7 — Closing: Roadmap, Honest Gaps, and Takeaways

**Suggested slide title:** "What's Next, What's Missing, and Three Things to Take Away"
**Target time:** 2–3 minutes

---

## Slide 7.1 — Honest Gaps

**[33:00]**

The best way to earn trust is to tell you what the code does not do well
before you find out yourself.

**Thin 3D semiconductor coverage.** The shipped 3D benchmarks are a doped
resistor and a pure-Poisson box used for the GPU speedup test. There is
no real 3D semiconductor device. M19 closes this with a 3D MOSFET on a
gmsh unstructured mesh; it depends on M16.1 (Caughey-Thomas mobility) so
the saturation region has any physical meaning.

**Boltzmann statistics throughout.** The carrier statistics use the
Boltzmann approximation, which is valid below about 10¹⁹ cm⁻³. Every
modern MOSFET source/drain extension is above that threshold. M16.4
replaces Boltzmann with Fermi-Dirac.

**Physics catalog is incomplete.** No Auger or radiative recombination.
No Schottky contacts. No band-to-band or trap-assisted tunneling. No
field-dependent mobility (Caughey-Thomas or Lombardi surface mobility).
These are all on the M16 roadmap, each as a separate PR with its own
analytical acceptance test.

**Key points**
- No 3D semiconductor device yet (M19 is the capstone).
- Boltzmann only, valid < 10¹⁹ cm⁻³ (M16.4 adds Fermi-Dirac).
- M16.1–M16.7 add the missing physics, one ADR and one verifier each.

---

## Slide 7.2 — What's Planned

**[34:00]**

The forward roadmap in priority order:

1. **M16.1–M16.7: Physics completeness pass.** Caughey-Thomas
   field-dependent mobility, Lombardi surface mobility, Auger
   recombination, Fermi-Dirac statistics, Schottky contacts, BBT/TAT
   tunneling, time-varying transient voltages.
2. **M17: Heterojunctions.** Position-dependent band structure,
   enabling HEMTs and heterostructure diodes.
3. **M19: 3D MOSFET capstone.** Full 3D FET on an unstructured mesh with
   CPU and GPU verification.
4. **M19.1: MPI parallel.** Multi-process runners for large 3D problems.
5. **M20: HTTP server hardening.** Authentication, rate limiting, API keys.

Each milestone has a written acceptance test in `docs/IMPROVEMENT_GUIDE.md`
before any code is written. The acceptance test is the definition of done.

**Key points**
- Physics completeness (M16) before 3D (M19) — physical dependencies drive order.
- Every milestone: acceptance test written before implementation begins.
- Roadmap is public; contributions welcome.

---

## Slide 7.3 — Three Takeaways

**[35:00]**

If you remember nothing else from this talk, remember three things:

**1. The three numerical fixes are the foundation.**
Nondimensionalization kills the 10³⁰ condition number. Slotboom
variables make the FEM stable without stabilization. Adaptive
continuation makes biased solves converge. Without all three, the solver
diverges on every real problem. They are not optimizations; they are the
minimum viable set.

**2. JSON as contract, not Python as contract.**
Making JSON the only input format — and enforcing it with a strict schema
with `additionalProperties: false` — means a typo in the input is a
validation error, not a silent wrong answer. Any tool can generate,
inspect, and store a device spec without a Python runtime. This is the
right call for a tool that will be driven by UIs, AI agents, and web APIs.

**3. Verification before validation, and both before shipping.**
A benchmark that matches one analytical formula at one resolution is
not enough. MMS tells you whether the discretization converges correctly.
Conservation checks tell you whether the physics is self-consistent. CI
gates tell you whether a change breaks either. The V&V suite is not a
formality; it has caught actual bugs that the benchmarks missed.

---

## Slide 7.4 — Q&A Seed Questions

**[36:00]** These are questions the audience is likely to ask. Have
answers ready.

**"Why FEniCSx over FEniCS Classic, devsim, or COMSOL?"**
FEniCSx 0.10 provides a stable, actively maintained Python API for UFL
residual forms with PETSc backend. devsim uses a different discretization
(box-integration) that is efficient but not straightforward to extend. We
chose FEniCSx to stay in the finite-element world where MMS verification
is natural. COMSOL is not open source.

**"Why not use the Scharfetter-Gummel discretization?"**
SG is excellent for finite-volume codes. On a standard P1 Lagrange FEM
mesh in 3D it requires special quadrature on edges, which does not
generalize cleanly. Slotboom gives Galerkin stability with no mesh
dependence and second-order convergence confirmed by MMS. ADR 0004
has the full reasoning.

**"What happens at the Debye-length scale?"**
The small parameter λ² means the space-charge layers are thin compared
to the mesh cells in most simulations. The solver resolves them
approximately, not exactly. This is acceptable for computing integrated
quantities (terminal current, total charge); it becomes a problem if you
need field profiles with sub-nanometer resolution. Adaptive mesh
refinement at the junction is on the future roadmap.

**"How hard is it to add a new device?"**
Write a JSON file (existing ones are good templates), optionally write a
new runner if the physics is genuinely different, and write an MMS
verifier. The `CONTRIBUTING.md` walkthrough covers the benchmark layout.
A new 1D device from JSON to passing verifier is typically one day of work
for a contributor who knows the physics.

**"Is this production-ready?"**
No, and the roadmap says so honestly. The HTTP server has no
authentication or rate limiting (M20). The 3D semiconductor coverage is
thin (M19). Boltzmann statistics break above 10¹⁹ cm⁻³ (M16.4). It is
a research simulator and a teaching tool, not a replacement for Sentaurus
or COMSOL.
