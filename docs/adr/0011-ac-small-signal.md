# ADR 0011: Small-signal AC sweep — formulation, primary variables, and PETSc-real 2x2 block reformulation

- **Status:** Accepted
- **Date:** 2026-04-26
- **Milestone:** M14 — small-signal AC analysis

---

## Context

`docs/IMPROVEMENT_GUIDE.md` §M14 lists small-signal AC analysis as the
next solver-type addition after the M13 transient runner. The user-
facing motivation:

- **Capacitance-voltage extraction.** True C(V) curves require either a
  transient ramp + numerical differentiation of Q(V) or a small-signal
  AC measurement of -Im(Y)/omega at each V. The latter is ~10x cheaper
  per V-point (one linear solve per frequency vs hundreds of nonlinear
  time-steps) and is what every commercial DD simulator exposes.
- **Admittance / S-parameter prediction.** RF MOSFETs and HEMTs
  characterise the device by Y- (or S-) parameters at GHz frequencies.
  The DD-level AC sweep is the engine that produces these.
- **Noise analysis groundwork.** Small-signal noise ((1/f), shot,
  generation-recombination) reduces to a related linear AC system
  driven by stochastic sources at each frequency. Standing up the AC
  infrastructure now positions us for that analysis later (M16+).

The math reuses the steady-state Jacobian J of the bias_sweep solver
and the lumped mass matrix M shipped with the M13 transient runner.

## Decision

### Formulation

Around a converged DC operating point u_0 = (psi_0, n_0, p_0) at a
chosen contact bias V_DC, the small-signal perturbation u(t) = u_0 +
delta_u * exp(j*omega*t) under V(t) = V_DC + delta_V * exp(j*omega*t)
satisfies the linear frequency-domain system

    (J + j*omega*M) * delta_u = -dF/dV * delta_V

where:

- **J** is the steady-state DD Jacobian dF/du at u_0, in primary-density
  form (psi, n_hat, p_hat).
- **M** is the lumped mass matrix on the carrier-density rows. The
  Poisson row has no time derivative (Maxwell's equations are quasi-
  electrostatic in the DD model used throughout this engine), so
  M_psi = 0 and only M_n and M_p are nonzero. The lumped diagonal is
  the same construct shipped with M13 (`semi/fem/mass.py`).
- **dF/dV** is the right-hand side arising from the small change in the
  contact's Dirichlet BC. The runner forms it implicitly via finite
  difference of the DC operating point: solving the steady-state DD
  problem at V_DC and at V_DC + eps_V, the difference of the two
  solutions delta_u_DC = (u(V_DC + eps_V) - u(V_DC)) / eps_V satisfies
  J * delta_u_DC = -dF/dV at u_0 by Newton's lemma. The AC RHS is then
  b = J * delta_u_DC (one mat-vec).

### Primary unknowns: (psi, n, p) primary-density

The bias_sweep runner solves in Slotboom variables (psi, phi_n, phi_p)
because steady-state DD is coercive in that form (see ADR 0004). The
transient runner switched to (psi, n_hat, p_hat) primary-density form
to make the time-derivative term natural and to share a mass matrix
with the carrier-density continuity equations (see ADR 0009).

The AC runner adopts the same primary-density choice as the transient
runner. Two reasons:

1. **Mass matrix reuse.** `semi.fem.mass.assemble_lumped_mass` is
   defined on the (n, p) spaces. Reformulating M for Slotboom variables
   would require a chain-rule mass matrix that couples the Poisson and
   continuity rows (since dn/dt depends on dpsi/dt and dphi_n/dt) — see
   the discussion in ADR 0009. The (psi, n, p) form avoids this entirely.
2. **Ecosystem coherence.** Future M14.1 (transient C-V extraction
   matching AC C-V), M16 (Auger and field-dependent mobility) are all
   easier to extend in (psi, n, p) form.

The DC operating point is obtained by running `run_bias_sweep` (which
gives Slotboom variables) and then converting to (n, p) via Boltzmann:

    n_hat = (n_i / C_0) * exp(psi_hat - phi_n_hat)
    p_hat = (n_i / C_0) * exp(phi_p_hat - psi_hat)

This is exact under Boltzmann statistics (the only statistics M14
supports; Fermi-Dirac is M16 work).

### PETSc complex-vs-real choice — real 2x2 block reformulation

The `dolfinx-real` PETSc build (the only build in the engine's CI
Docker image) uses real `PetscScalar` (float64). Solving (J + j*omega*M)
delta_u = b directly would require a `dolfinx-complex` build, which
neither CI nor any benchmark environment currently has.

We therefore use the standard **real 2x2 block reformulation**: write
delta_u = x + j*y and the equation as

    [   J     -omega*M ] [ x ]   [ Re(b) ]
    [ omega*M    J     ] [ y ] = [ Im(b) ]

(BC rows of J are replaced by identity by the dolfinx assembler; M is
zeroed at BC rows so the BC perturbation is imposed identically at
every frequency.) The block size is 2 * (3N) = 6N where N is the
per-block DOF count. A direct LU factorisation (MUMPS) is used at each
frequency. For the M14 benchmark mesh sizes (sub-1k DOFs in 1D, sub-
10k in 2D MOS C-V) this is sub-second per frequency.

Switching to a complex PETSc build in the future would shrink the
solved matrix to 3N complex entries and roughly halve the per-frequency
LU cost, but the runner contract (`AcSweepResult` shape, the
`backend` meta tag) does not depend on this choice. Adding a complex
backend would only require:

1. A complex PETSc build in the Docker image (`PETSC_USE_COMPLEX=1`).
2. A branch in `_solve_2x2_real_block` that builds (J + j*omega*M) as
   a single complex matrix when `PETSc.ScalarType` is complex.

Both are deferred (M16 or beyond). Documenting the build-time switch
here keeps the option open without committing to a Docker image rebuild
during M14.

### Sign convention for terminal admittance

The runner reports Y with the **"positive = current INTO the device"**
convention used throughout `semi/postprocess.py` for the steady-state
runners (`evaluate_current_at_contact` integrates
`q μ_n n ∇φ_n · n_FE`, which by Slotboom algebra is the negative of
the (n,p)-primary outward current density at the contact, i.e. it
reports the current INTO the device). Under this convention an
ideal capacitor terminal has

    Y = +j * omega * C

so the runner reports

    C(omega) = +Im(Y) / (2 * pi * f)

and a positive physical capacitance maps to a positive C in
`AcSweepResult.C`. The same sign convention applies to G = Re(Y):
positive G corresponds to a positive `dI/dV` from `bias_sweep` at
the same operating point.

> **Note (2026-04-28):** The original M14 implementation had a sign
> bug — it reported Y in the OUT-of-device convention while this
> ADR's prose described the same OUT convention. Audit Cases 02 and
> 05 surfaced the inconsistency with `bias_sweep`. Both the
> implementation and this section were corrected to the
> INTO-device convention. See the **Errata** section at the bottom
> of this ADR for diagnosis and fix details.

## Validation

### MMS-style consistency tests
`tests/mms/test_ac_consistency.py` checks three internal invariants of
the formulation that do not depend on a closed-form analytical
solution:

1. Y(omega = 0) is purely real (the 2x2 block decouples; imaginary
   solution component is identically zero).
2. Re(Y) at omega = 0 and omega = 2*pi*1e-3 Hz agree to better than 1%
   (linearisation-error scale of the eps_V = 1e-3 V finite-difference
   DC sensitivity).
3. C(f) is omega-independent at low frequency in the depletion regime
   (tested at f = 1 Hz vs f = 2 Hz; matches to better than 1e-6).

### Acceptance benchmark — RC depletion capacitance
`benchmarks/rc_ac_sweep` runs a 1D pn diode at V_DC = -1.0 V across
41 logspace frequencies from 1 Hz to 1 GHz. The verifier checks the
capacitance plateau matches the analytical depletion C within 5% over
[1 Hz, 1 MHz]; in current measurements it matches to **0.41%** worst-
case (well under tolerance). Plateau flatness max(C)/min(C) is **1.000**
across 27 samples in band.

### Acceptance test — pn depletion-C agreement across bias
`tests/fem/test_ac_dc_limit.py` evaluates the M14 acceptance criterion
2 from `docs/IMPROVEMENT_GUIDE.md`: low-frequency C(omega) matches
analytical depletion C within 5% at three biases V_DC in {-2.0, -1.0,
-0.5} V. Current results match to within 1% at V = -1.0 V and within
0.4% at V = -2.0 V.

## Alternatives considered

### Time-domain transient + FFT (rejected)
Drive a sinusoidal small-amplitude perturbation in time, FFT the
terminal current, extract Y(omega). Conceptually simpler — reuses the
existing M13 transient runner verbatim — but requires hundreds to
thousands of timesteps per frequency point and re-solves the full
nonlinear system at every step. Per-frequency cost is two to three
orders of magnitude higher than the M14 linear-solve approach, and
spectral leakage / windowing introduces extra accuracy concerns.
Rejected as the primary path; remains available as a sanity-check
mechanism in M16.

### Harmonic balance (rejected)
Express V(t) as a Fourier series over the steady period, solve the
coupled multi-harmonic nonlinear system. Standard for circuit
simulators (e.g. SPECTRE), genuinely needed for large-signal AC
behaviour. Overkill for small-signal where the linearisation is
already exact. Rejected for M14; potentially M18+ work.

### Complex PETSc build (deferred)
See "PETSc complex-vs-real choice" above. The real-block reformulation
is well-understood, has a known 2x condition-number penalty (the
largest eigenvalue of the 2x2 block is sqrt(2) larger than the original
J's), and works without a Docker rebuild. Deferred.

### Slotboom-variable AC formulation (rejected)
Could reuse the bias_sweep machinery directly without converting to
density form. Requires a chain-rule mass matrix that couples Poisson
with continuity blocks (see ADR 0009 §"Why density form for transient");
mass matrix would no longer be a simple lumped diagonal on the
continuity rows. Discarded for the same reasons that drove ADR 0009.

## Consequences

### Positive
- Small-signal AC analysis becomes a first-class solver type; the
  schema, runner, result type, and CLI driver all support it.
- Mass matrix infrastructure shipped with M13 is exercised by a second
  consumer (transient + AC), validating the M13 contract.
- The runner's `backend` meta field documents the real-block path,
  making the future complex-PETSc switch a discoverable upgrade.

### Negative / risk
- The real-block formulation doubles the linear-system size per
  frequency. Acceptable for M14 mesh sizes; would require the M15 GPU
  path (or a complex PETSc build) for ~100k+ DOF problems.
- The DC sensitivity is computed by finite difference (eps_V = 1e-3 V)
  rather than directly assembling dF/dV. This couples the AC accuracy
  to the SNES tolerance of the bias_sweep runner. The choice of eps_V
  trades linearisation error against finite-difference noise; 1e-3 V
  is well below the V_t ≈ 0.026 V scale and well above the SNES atol
  floor of 1e-7 in scaled units.

### Deferred
- Complex PETSc backend (see above).
- Direct dF/dV assembly via UFL action on a BC-perturbation function.
  Cleaner than finite difference; deferred to M16 along with the
  Auger/field-mobility refresh.
- `mos_2d` C-V upgrade to use AC admittance instead of `mos_cv`'s
  numerical dQ/dV. Listed as M14 deliverable in `IMPROVEMENT_GUIDE.md`
  but is **not** done in this PR — `mos_cv` is left untouched. Tracked
  as M14.1.

## Errata (2026-04-28) — terminal admittance sign convention fix

The original M14 implementation reported `Y` and `C` with the
**opposite** sign to the convention this ADR claimed. The Phase 1
physics validation audit (Cases 02 and 05 in `docs/PHYSICS_AUDIT.md`)
flagged this as a Class C finding: `Re(Y(omega→0))` at the same DC
operating point came out with the same magnitude but opposite sign
to `bias_sweep`'s centered-difference `dI/dV`. Since `bias_sweep` is
the convention authority used by all M11/M12 benchmarks, the fix is
a global sign correction in `run_ac_sweep`, not in `bias_sweep`.

### Convention (post-fix)

`Y` is reported with the **"positive = current INTO the device"**
convention used by `bias_sweep` /
`semi/postprocess.evaluate_current_at_contact` (i.e. the same sign
that `iv_row["J"]` carries on the steady-state IV record). Under
this convention an ideal capacitor terminal has

    Y = +j * omega * C

and the runner reports

    C(omega) = +Im(Y) / (2 * pi * f)

so a positive depletion capacitance maps to a positive
`AcSweepResult.C[k]`, matching the "Validation" section above and
the `benchmarks/rc_ac_sweep` verifier. `G = Re(Y)` is positive when
the linearised terminal current increases with applied bias (i.e.
`dI/dV > 0`).

### Why the original implementation had it backwards

In Slotboom variables (`n = n_i exp(psi - phi_n)`), `bias_sweep`
integrates `+q*mu_n*n*grad(phi_n).n_outward` at the contact. After
substituting Slotboom into the (n,p)-primary outward flux
`q*mu_n*V_t*grad(n) - q*mu_n*n*grad(psi)`, algebra gives

    bias_sweep_J = -[(n,p)-primary outward J]

The M14 ac_sweep runner computes the linearised `(n,p)-primary
outward J` (and an analogous outward displacement current); the
result is the small-signal terminal current flowing OUT of the
device. Bias_sweep's reported value is the negative of that, so to
make `Re(Y(omega→0))` agree with `dI/dV` we must negate the
assembled total current. The fix is one negation at the assembly
site in `run_ac_sweep` plus a matching sign-flip in the `C` formula
so that `result.C` (and hence the RC depletion-C verifier) is
numerically unchanged.

### Verification

After the fix, audit Cases 02 and 05 are tightened from NaN-guards
to `sign(Re(Y)) == sign(dI/dV)` plus a 1% (Case 02) / 5% (Case 05)
relative-error gate. `tests/mms/test_ac_consistency.py` and
`tests/fem/test_ac_dc_limit.py` are unchanged and pass bit-
identically (they compare via absolute values and ratios).
`benchmarks/rc_ac_sweep` is unchanged: the verifier reads
`result.C`, which is bit-identical because both `Im(Y)` and the
read-out formula flipped sign together.


  Devices*. Springer. Chapter 7 §"Small-signal AC analysis".
- Schenk, A. (1998). *Advanced Physical Models for Silicon Device
  Simulation*. Springer. Chapter 5 §"Linearised continuity equations
  for AC analysis".
- Snowden, C. (1989). *Semiconductor Device Modelling*. Peter
  Peregrinus. §"Small-signal extraction from drift-diffusion".
- ADR 0004 — `docs/adr/0004-slotboom-variables-for-dd.md` (steady-state
  variable choice).
- ADR 0009 — `docs/adr/0009-transient-formulation.md` (primary-density
  form for time-dependent solves).
- ADR 0010 — `docs/adr/0010-bdf-time-integration.md` (lumped mass matrix
  rationale).
