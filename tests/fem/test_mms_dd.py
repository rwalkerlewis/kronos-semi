"""
MMS convergence tests for the coupled drift-diffusion block residual.

Mesh sequences here are the pytest-budget subset (3 levels in 1D and 2D)
so the full FEM test suite stays under the docker-fem budget. The CLI in
`scripts/run_verification.py` runs the fuller sweep used for the
published convergence plots. See `docs/mms_dd_derivation.md` for the
approved derivation and rate targets.

Variant A gates the psi block only: `phi_n_e = phi_p_e = 0` identically,
so the Fermi-level errors are at Newton-noise level and do not carry a
meaningful observed rate. Variants B and C gate all three blocks.
"""
from __future__ import annotations

import math


PYTEST_NS_1D = [40, 80, 160]
PYTEST_NS_2D = [16, 32, 64]

# Derivation Section 5: L^2 >= 1.75 and H^1 >= 0.80 per gated block.
RATE_L2_FLOOR = 1.75
RATE_H1_FLOOR = 0.80


def _finest_pair_rate(hs, errors):
    """Log slope of the last two (h, err) pairs."""
    h_prev, h_cur = hs[-2], hs[-1]
    e_prev, e_cur = errors[-2], errors[-1]
    return math.log(e_prev / e_cur) / math.log(h_prev / h_cur)


def _assert_monotone_and_rates(results, blocks):
    """
    Assert strict monotone L^2 error reduction and finest-pair rates on
    every block in `blocks`. `blocks` is a subset of ("psi", "phi_n",
    "phi_p").
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
            f"{RATE_L2_FLOOR:.2f}"
        )
        assert rate_H1 >= RATE_H1_FLOOR, (
            f"block {b!r}: finest-pair H^1 rate {rate_H1:.3f} < "
            f"{RATE_H1_FLOOR:.2f}"
        )


def test_mms_dd_1d_variant_A():
    """
    Variant A: Poisson block only, phi_n_e = phi_p_e = 0. Only the psi
    block has a meaningful convergence rate; the continuity blocks solve
    the trivially-zero exact solution and their errors are at Newton
    noise level.
    """
    from semi.verification.mms_dd import run_convergence_study

    results = run_convergence_study(
        dim=1, variant="A", Ns=PYTEST_NS_1D,
    )
    _assert_monotone_and_rates(results, blocks=("psi",))


def test_mms_dd_1d_variant_B():
    """
    Variant B: full three-block coupling with tau_hat = infinity so the
    SRH term is numerically negligible but the kernel path is still
    exercised. All three blocks must converge.
    """
    from semi.verification.mms_dd import run_convergence_study

    results = run_convergence_study(
        dim=1, variant="B", Ns=PYTEST_NS_1D,
    )
    _assert_monotone_and_rates(results, blocks=("psi", "phi_n", "phi_p"))


def test_mms_dd_1d_variant_C():
    """
    Variant C: full coupling with physical Si SRH lifetimes; R_e is
    pointwise nonzero because A_p = -A_n forces phi_p_e - phi_n_e != 0.
    All three blocks must converge at the gate rate.
    """
    from semi.verification.mms_dd import run_convergence_study

    results = run_convergence_study(
        dim=1, variant="C", Ns=PYTEST_NS_1D,
    )
    _assert_monotone_and_rates(results, blocks=("psi", "phi_n", "phi_p"))


def test_mms_dd_2d_variant_A():
    """2D right-diagonal triangles, Variant A (psi block only)."""
    from semi.verification.mms_dd import run_convergence_study

    results = run_convergence_study(
        dim=2, variant="A", Ns=PYTEST_NS_2D, cell_kind="triangle",
    )
    _assert_monotone_and_rates(results, blocks=("psi",))


def test_mms_dd_2d_variant_B():
    """2D Variant B: full coupling, R ~ 0; all three blocks must converge."""
    from semi.verification.mms_dd import run_convergence_study

    results = run_convergence_study(
        dim=2, variant="B", Ns=PYTEST_NS_2D, cell_kind="triangle",
    )
    _assert_monotone_and_rates(results, blocks=("psi", "phi_n", "phi_p"))


def test_mms_dd_2d_variant_C():
    """2D Variant C: full coupling with SRH; all three blocks gated."""
    from semi.verification.mms_dd import run_convergence_study

    results = run_convergence_study(
        dim=2, variant="C", Ns=PYTEST_NS_2D, cell_kind="triangle",
    )
    _assert_monotone_and_rates(results, blocks=("psi", "phi_n", "phi_p"))
