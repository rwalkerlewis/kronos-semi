# §5 — Capabilities and Live Demos

**Suggested slide title:** "What the Engine Can Do: A Tour of Shipped Capabilities"
**Target time:** 5–6 minutes

---

## Slide 5.1 — Devices and Solvers

**[22:00]**

Let me walk through what ships in v0.16. The key point is that each
capability is verified against an analytical reference; nothing is
"implemented but untested."

**Equilibrium Poisson (1D / 2D / 3D)**
The foundation. Solve Poisson's equation with fixed doping, recover the
depletion approximation. Verified by MMS: finest-pair L² rate = 2.0,
H¹ rate = 1.0 on 1D, 2D triangles, and 2D multi-region (silicon + oxide).

**pn Junction Bias Sweep (1D)**
Coupled Slotboom drift-diffusion with SRH recombination. Forward-bias IV
matches the Shockley diode equation within 10%. Reverse bias recovers SRH
generation current against the SNS analytical formula.

**2D MOSFET**
Multi-region mesh (silicon + gate oxide) with Gaussian n+ source/drain
implants. The Pao-Sah linear-regime formula provides the analytical
reference in the window [V_T + 0.2, V_T + 0.6] V at V_DS = 0.05 V. The
verifier passes within a 20% tolerance — appropriate for a first-principles
simulation vs. a formula derived from a textbook approximation.

**3D Doped Resistor**
Uniform-doped silicon bar with builtin axis-aligned mesh or gmsh
unstructured mesh. V-I linearity within 1%. Demonstrates 3D and
file-sourced mesh paths.

**Key points**
- Every shipped device has a registered analytical verifier.
- Equilibrium Poisson → pn bias → MOSFET: physics builds cumulatively.
- 3D path works; the 3D semiconductor gap closes with M19.

---

## Slide 5.2 — Transient, AC, and Axisymmetric

**[24:00]**

**Transient solver (M13/M13.1)**
BDF1 (backward Euler) and BDF2 time integration in Slotboom primary
unknowns (ψ, Φ_n, Φ_p). The steady-state limit test matches the bias
sweep result within 1×10⁻⁴ relative error. The `pn_1d_turnon` benchmark
tracks a diode from equilibrium to forward bias and compares the transient
approach to steady state within 5%.

**AC small-signal (M14 / M14.1)**
Two modes. For pn diodes, the small-signal linearised system
(J + jωM)δu = −(dF/dV)δV is solved as a real 2×2 block at each
frequency, recovering depletion capacitance within 0.4% of the
analytical formula from 1 Hz to 1 MHz. For MOS capacitors, the
`mos_cap_ac` runner solves the linearised Poisson sensitivity dQ/dV
analytically at each gate bias, replacing noisy numerical differentiation
and agreeing with it to machine precision.

**Axisymmetric 2D MOSCAP (M14.2)**
Cylindrical symmetry: the meridian half-plane is used with r-weighted
inner products. Reproduces the Hu textbook Figure 5-18 MOSCAP LF/HF
C-V split. Worst error vs. the depletion-approximation reference: 6.79%
in the verifier window.

**Key points**
- Transient: BDF1/BDF2 in Slotboom unknowns; verified by rate test and steady-state limit.
- AC: two-terminal depletion C within 0.4%; MOSCAP analytic dQ/dV to machine precision.
- Axisymmetric: r-weighted FEM; validates against Hu Fig. 5-18.

---

## Slide 5.3 — GPU Path and Mesh Flexibility

**[25:30]**

**GPU linear-solver path (M15)**
Set `solver.backend` to `gpu-amgx`, `gpu-hypre`, or `auto`. The linear
solve at each Newton step runs on the device when PETSc-CUDA or PETSc-HIP
is available. The CPU-MUMPS path is bit-identical to the pre-GPU version.
The acceptance test for M15 requires ≥5× linear-solve speedup at 500k
DOFs on the 3D Poisson benchmark. The `GET /capabilities` endpoint reports
available backends so a UI can gate the option.

**Mesh flexibility**
Builtin axis-aligned 1D/2D/3D boxes are enough for the physics benchmarks.
For real device geometry, the engine ingests gmsh `.msh` files and XDMF
files. Physical groups propagate verbatim as cell and facet tags.

**Result artifacts**
Every solve writes a `runs/<run_id>/` directory: a schema-validated
`manifest.json`, field files (XDMF), IV CSVs, and convergence logs. Any
tool that can read JSON and XDMF can consume the results without importing
the Python package.

**Key points**
- GPU path: no code change required; set `solver.backend` in JSON.
- Mesh: builtin, gmsh .msh, XDMF — physical groups pass through verbatim.
- Artifacts: schema-validated manifest + XDMF + CSV; readable by any tool.

**Transition:** Let me now talk about how we got here — the development
process and the design decisions that shaped the code.
