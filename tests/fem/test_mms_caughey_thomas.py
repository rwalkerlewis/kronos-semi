"""
MMS convergence tests for the M16.1 Caughey-Thomas field-dependent
mobility (Variant D in `semi/verification/mms_dd.py`).

ADR 0006 mandates an MMS verifier for every new physics module. The
Caughey-Thomas closed form is a smooth nonlinearity of the carrier-
specific quasi-Fermi gradient, so the P1 Galerkin discretization
should preserve the standard polynomial-order rates: finest-pair
L2 rate >= 1.99 and H1 rate >= 0.99 on every gated block.

The Variant D MMS engineers `MMS_D_VSAT_*_FOR_FORM` so the
dimensionless ratio `(mu0 * F_par / vsat)^beta` is O(0.3) at the
typical manufactured gradient (mu reduction ~7%); the path is
materially exercised, not numerically dormant.

Mesh sequences here mirror Variants B and C (3 levels in 1D, 3 in
2D) so the FEM test budget stays inside docker-fem's 15-minute
ceiling. The CLI driver (`scripts/run_verification.py mms_dd`) runs
the longer 4-level 1D sweep and gates the same rates.
"""
from __future__ import annotations

import math

PYTEST_NS_1D = [40, 80, 160]
# 2D needs [32, 64, 128]: the [16, 32, 64] sequence used by Variants
# A/B/C bottoms out at finest-pair rate 1.990 from triangle-mesh
# boundary-layer effects, just under the M16.1 acceptance floor of
# 1.99. One extra refinement level pulls the rate cleanly to ~1.997
# without breaking the docker-fem 15-minute budget.
PYTEST_NS_2D = [32, 64, 128]

# M16.1 acceptance floor (Acceptance test 1 in
# docs/IMPROVEMENT_GUIDE.md § M16.1):
RATE_L2_FLOOR = 1.99
RATE_H1_FLOOR = 0.99


def _finest_pair_rate(hs, errors):
    """Log slope of the last two (h, err) pairs."""
    h_prev, h_cur = hs[-2], hs[-1]
    e_prev, e_cur = errors[-2], errors[-1]
    return math.log(e_prev / e_cur) / math.log(h_prev / h_cur)


def _assert_monotone_and_rates(results, blocks):
    """
    Strict monotone L^2 reduction plus finest-pair rate gates on
    every block in `blocks` (subset of psi / phi_n / phi_p).
    """
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
            f"{RATE_L2_FLOOR:.2f} (M16.1 acceptance gate)"
        )
        assert rate_H1 >= RATE_H1_FLOOR, (
            f"block {b!r}: finest-pair H^1 rate {rate_H1:.3f} < "
            f"{RATE_H1_FLOOR:.2f} (M16.1 acceptance gate)"
        )


def test_mms_dd_1d_variant_D_caughey_thomas():
    """1D Variant D: full coupling with SRH and Caughey-Thomas
    mobility; all three blocks gated at the milestone L2 >= 1.99,
    H1 >= 0.99 floor."""
    from semi.verification.mms_dd import run_convergence_study

    results = run_convergence_study(
        dim=1, variant="D", Ns=PYTEST_NS_1D,
    )
    _assert_monotone_and_rates(results, blocks=("psi", "phi_n", "phi_p"))


def test_mms_dd_2d_variant_D_caughey_thomas():
    """2D Variant D: same as the 1D test on right-diagonal triangles."""
    from semi.verification.mms_dd import run_convergence_study

    results = run_convergence_study(
        dim=2, variant="D", Ns=PYTEST_NS_2D, cell_kind="triangle",
    )
    _assert_monotone_and_rates(results, blocks=("psi", "phi_n", "phi_p"))
