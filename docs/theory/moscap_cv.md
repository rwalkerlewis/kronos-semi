# MOSCAP capacitance–voltage: LF and HF derivations

This note maps every analytical formula in
[`semi/cv.py`](../../semi/cv.py) to its reference in Hu,
*Modern Semiconductor Devices for IC*, chapter 5, and explains why
the current high-frequency implementation is a depletion-approximation
clamp rather than a true AC small-signal solve.

## Setup

A MOS capacitor under bias $V_g$ has gate stack
gate / oxide ($T_{\text{ox}}$, $\varepsilon_{\text{ox}}$) / Si body
($N_a$ for nMOS or $N_d$ for pMOS). The flat-band condition zeros the
band bending; the threshold condition pins the surface potential at
strong inversion.

## Closed-form anchors

Implemented in `analytical_moscap_params`
([`semi/cv.py`](../../semi/cv.py#L62)):

| Quantity | Formula | Hu reference |
|---|---|---|
| Thermal voltage | $V_t = kT/q$ | — |
| Bulk Fermi potential | $|\phi_B| = V_t\,\ln(N_a/n_i)$ | Eq. 5.2.4 |
| Oxide capacitance | $C_{\text{ox}} = \varepsilon_{\text{ox}}/T_{\text{ox}}$ | Eq. 5.4.1 |
| Flat-band voltage | $V_{\text{fb}} = \phi_{ms} - Q_f/C_{\text{ox}}$ | Eq. 5.4.3 |
| Strong-inversion surface pot. | $\psi_s = 2|\phi_B|$ | Eq. 5.5.1 |
| Max depletion width | $W_{\text{dmax}} = \sqrt{2\varepsilon_s\,(2\phi_B)/(qN_a)}$ | Eq. 5.6.4 |
| Threshold voltage | $V_t^* = V_{\text{fb}} + 2\phi_B + \sqrt{2\varepsilon_s qN_a (2\phi_B)}/C_{\text{ox}}$ | Eq. 5.5.1 |
| HF minimum capacitance | $C_{\min} = (1/C_{\text{ox}} + W_{\text{dmax}}/\varepsilon_s)^{-1}$ | Eq. 5.6.5 |

## Reference numbers (Hu Fig. 5-18 parameters: $N_a = 5\times10^{16}$ cm$^{-3}$, $T_{\text{ox}} = 10$ nm, $\phi_{ms} = -0.95$ V)

These numbers are produced by `analytical_moscap_params` and asserted
by [`tests/test_axisym_moscap_math.py`](../../tests/test_axisym_moscap_math.py):

- $V_{\text{fb}} = -0.950$ V
- $V_t = +0.181$ V
- $|\phi_B| = 0.399$ V
- $W_{\text{dmax}} = 144$ nm
- $C_{\min}/C_{\text{ox}} = 0.173$

## Low-frequency C–V

Quasi-static: minority carriers track the bias. The total gate
charge $Q_g(V_g) = -Q_s(V_g)$ where $Q_s$ is the surface charge
including the inversion sheet. Capacitance is
$C_{\text{LF}}(V_g) = -\mathrm{d}Q_s/\mathrm{d}V_g$.

Implementation:
- Pure-Python analytical reference: `lf_cv_quasistatic`
  ([`semi/cv.py`](../../semi/cv.py#L210)).
- FEM postprocessor: `compute_lf_cv_fem`
  ([`semi/cv.py`](../../semi/cv.py#L288)) takes
  $(V_g, Q_s)$ from a steady-state MOSCAP sweep and finite-differences
  $-\mathrm{d}Q_s/\mathrm{d}V_g$.

LF should asymptote to $C_{\text{ox}}$ in deep accumulation
($V_g \le V_{\text{fb}} - 0.5$ V) and deep inversion
($V_g \ge V_t + 1.0$ V).

## High-frequency C–V (depletion-approximation clamp)

In the strict definition (Hu §5.6), the inversion layer cannot
follow at high frequency, so the small-signal capacitance saturates
at $C_{\min}$ once the surface potential reaches $2|\phi_B|$. The
correct way to extract this from a FEM solve is an AC small-signal
sweep at the operating point: solve $(J + j\omega M)\,\delta u
= -\mathrm{d}F/\mathrm{d}V\,\delta V$ at high $\omega$ and report
$C(\omega) = -\Im(Y)/(2\pi f)$.

The current implementation **does not** do that. It implements the
classical depletion-approximation clamp:

$$C_{\text{HF}}(V_g) = \begin{cases}
C_{\text{LF}}(V_g) & \text{if } \psi_s(V_g) < 2|\phi_B| \\
C_{\min} & \text{otherwise}
\end{cases}$$

Implementation: `hf_cv_depletion_approximation`
([`semi/cv.py`](../../semi/cv.py#L159)) and the FEM postprocessor
`compute_hf_cv_depletion_clamp`
([`semi/cv.py`](../../semi/cv.py#L313)). This matches the textbook
$C_{\text{HF}}$ shape to within a few percent for the standard
$N_a = 5\times10^{16}$, $T_{\text{ox}} = 10$ nm anchor case but is
not a true AC simulation.

A rigorous AC small-signal HF method is tracked as a follow-up; see
the `[follow-up]` GitHub issue list referenced from the post-merge
cleanup PR description.

## Where this lives in the code

- [`semi/cv.py`](../../semi/cv.py) — analytical helpers and FEM
  post-processors.
- [`semi/physics/axisymmetric.py`](../../semi/physics/axisymmetric.py) —
  $r$-weighted weak forms used by the MOSCAP benchmark.
- [`benchmarks/moscap_axisym_2d/`](../../benchmarks/moscap_axisym_2d/)
  — reference CSV and gmsh meridian mesh.
- [`notebooks/05_moscap_axisym_cv.ipynb`](../../notebooks/05_moscap_axisym_cv.ipynb)
  — end-to-end demonstration.

## See also

- [theory/axisymmetric.md](axisymmetric.md) — geometry and $r$-weighted weak forms.
- [benchmarks/moscap_axisym_2d.md](../benchmarks/moscap_axisym_2d.md) — benchmark landing page.
