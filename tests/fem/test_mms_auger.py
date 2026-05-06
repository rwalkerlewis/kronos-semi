"""
MMS convergence tests for the M16.3 Auger recombination kernel
(Variant F in `semi/verification/mms_dd.py`).

ADR 0006 mandates an MMS verifier for every new physics module. The
Auger kernel is cubic in carrier density (versus SRH's
bilinear-over-linear); the P1 Galerkin discretization should
preserve the standard polynomial-order rates: finest-pair L2 rate
>= 1.99 and H1 rate >= 0.99 on every gated block.

Variant F engineers `MMS_F_C_*_HAT_FOR_FORM` so the Auger
contribution is comparable in magnitude to SRH at the typical
manufactured amplitudes (~30 % reduction target, matching M16.1
Variant D and M16.2 Variant E); the path is materially exercised,
not numerically dormant.

Mesh sequences mirror Variants D and E (one level shy of the full
4-level 1D sweep, full 3-level 2D sweep with [32, 64, 128]) so the
FEM test budget stays inside docker-fem's 15-minute ceiling. The
CLI driver (`scripts/run_verification.py mms_dd`) runs Variant F on
the same sequence and gates the same rates.
"""
from __future__ import annotations

import math

PYTEST_NS_1D = [40, 80, 160]
PYTEST_NS_2D = [32, 64, 128]

# M16.3 acceptance floor (Acceptance test in
# docs/IMPROVEMENT_GUIDE.md § M16.3):
RATE_L2_FLOOR = 1.99
RATE_H1_FLOOR = 0.99


def _finest_pair_rate(hs, errors):
    """Log slope of the last two (h, err) pairs."""
    h_prev, h_cur = hs[-2], hs[-1]
    e_prev, e_cur = errors[-2], errors[-1]
    return math.log(e_prev / e_cur) / math.log(h_prev / h_cur)


def _assert_monotone_and_rates(results, blocks):
    """Strict monotone L^2 reduction plus finest-pair rate gates."""
    hs = [r.h for r in results]
    for b in blocks:
        eL2 = [getattr(r, f"e_L2_{b}") for r in results]
        eH1 = [getattr(r, f"e_H1_{b}") for r in results]
        assert all(eL2[i + 1] < eL2[i] for i in range(len(eL2) - 1)), (
            f"L^2 error on block {b!r} should decrease monotonically: {eL2}"
        )
        rate_L2 = _finest_pair_rate(hs, eL2)
        rate_H1 = _finest_pair_rate(hs, eH1)
        assert rate_L2 >= RATE_L2_FLOOR, (
            f"block {b!r}: finest-pair L^2 rate {rate_L2:.3f} < "
            f"{RATE_L2_FLOOR:.2f} (M16.3 acceptance gate)"
        )
        assert rate_H1 >= RATE_H1_FLOOR, (
            f"block {b!r}: finest-pair H^1 rate {rate_H1:.3f} < "
            f"{RATE_H1_FLOOR:.2f} (M16.3 acceptance gate)"
        )


def test_mms_dd_1d_variant_F_auger():
    """1D Variant F: full coupling with SRH + Auger; all three blocks
    gated at L2 >= 1.99, H1 >= 0.99."""
    from semi.verification.mms_dd import run_convergence_study

    results = run_convergence_study(
        dim=1, variant="F", Ns=PYTEST_NS_1D,
    )
    _assert_monotone_and_rates(results, blocks=("psi", "phi_n", "phi_p"))


def test_mms_dd_2d_variant_F_auger():
    """2D Variant F: same as the 1D test on right-diagonal triangles."""
    from semi.verification.mms_dd import run_convergence_study

    results = run_convergence_study(
        dim=2, variant="F", Ns=PYTEST_NS_2D, cell_kind="triangle",
    )
    _assert_monotone_and_rates(results, blocks=("psi", "phi_n", "phi_p"))
