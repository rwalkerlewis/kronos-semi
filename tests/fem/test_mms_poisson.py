"""
MMS-Poisson convergence tests.

Mesh sequences here are the pytest-budget subset (4 levels in 1D, 3 in
2D) so the full FEM test suite stays under one minute on the docker-fem
job. The CLI in `scripts/run_verification.py` runs the fuller sweep
used for the published convergence plot.
"""
from __future__ import annotations

import math

import pytest

PYTEST_NS_1D = [40, 80, 160, 320]
PYTEST_NS_2D = [16, 32, 64]


def _finest_pair_rate(hs, errors):
    """log slope of the last two (h, err) pairs."""
    h_prev, h_cur = hs[-2], hs[-1]
    e_prev, e_cur = errors[-2], errors[-1]
    return math.log(e_prev / e_cur) / math.log(h_prev / h_cur)


def test_mms_poisson_1d_linear_regime():
    """
    amplitude = 0.5: exp(+/- psi) ~ 1 + psi + psi^2/2, so the residual is
    dominated by the Laplacian perturbation. L^2 rate must be >= 1.85 on
    the finest pair; H^1 rate must be >= 0.85.
    """
    from semi.verification.mms_poisson import run_convergence_study

    results = run_convergence_study(dim=1, Ns=PYTEST_NS_1D, amplitude=0.5)
    hs = [r.h for r in results]
    eL2 = [r.e_L2 for r in results]
    eH1 = [r.e_H1 for r in results]

    # Monotone error reduction across the sweep
    assert all(eL2[i + 1] < eL2[i] for i in range(len(eL2) - 1)), \
        f"L^2 error should decrease monotonically: {eL2}"

    rate_L2 = _finest_pair_rate(hs, eL2)
    rate_H1 = _finest_pair_rate(hs, eH1)
    assert rate_L2 >= 1.85, (
        f"finest-pair L^2 rate {rate_L2:.3f} < 1.85 for amplitude=0.5"
    )
    assert rate_H1 >= 0.85, (
        f"finest-pair H^1 rate {rate_H1:.3f} < 0.85 for amplitude=0.5"
    )


def test_mms_poisson_1d_nonlinear_regime():
    """
    amplitude = 2.0: exp(2) ~ 7.4, so the carrier-nonlinearity dominates
    the residual. This is the regime physical devices operate in
    (real pn junctions see ~20 V_T of band bending). L^2 rate must
    still be >= 1.85 on the finest pair. Rate degradation here is a
    legitimate finding to document; the prompt explicitly says do not
    lower the tolerance to make the test pass.
    """
    from semi.verification.mms_poisson import run_convergence_study

    results = run_convergence_study(dim=1, Ns=PYTEST_NS_1D, amplitude=2.0)
    hs = [r.h for r in results]
    eL2 = [r.e_L2 for r in results]
    eH1 = [r.e_H1 for r in results]

    assert all(eL2[i + 1] < eL2[i] for i in range(len(eL2) - 1)), \
        f"L^2 error should decrease monotonically: {eL2}"

    rate_L2 = _finest_pair_rate(hs, eL2)
    rate_H1 = _finest_pair_rate(hs, eH1)
    assert rate_L2 >= 1.85, (
        f"finest-pair L^2 rate {rate_L2:.3f} < 1.85 for amplitude=2.0"
    )
    assert rate_H1 >= 0.85, (
        f"finest-pair H^1 rate {rate_H1:.3f} < 0.85 for amplitude=2.0"
    )


def test_mms_poisson_2d_convergence():
    """
    2D right-diagonal triangles; L^2 rate >= 1.85 on the finest pair.
    """
    from semi.verification.mms_poisson import run_convergence_study

    results = run_convergence_study(
        dim=2, Ns=PYTEST_NS_2D, amplitude=0.5, cell_kind="triangle",
    )
    hs = [r.h for r in results]
    eL2 = [r.e_L2 for r in results]
    eH1 = [r.e_H1 for r in results]

    assert all(eL2[i + 1] < eL2[i] for i in range(len(eL2) - 1)), \
        f"L^2 error should decrease monotonically: {eL2}"

    rate_L2 = _finest_pair_rate(hs, eL2)
    rate_H1 = _finest_pair_rate(hs, eH1)
    assert rate_L2 >= 1.85, (
        f"finest-pair 2D L^2 rate {rate_L2:.3f} < 1.85"
    )
    assert rate_H1 >= 0.85, (
        f"finest-pair 2D H^1 rate {rate_H1:.3f} < 0.85"
    )


def test_mms_poisson_2d_quad_smoke():
    """
    Single-mesh-size quadrilateral mesh at N=64; error must be of the
    same order as the triangle case at the same N. This is a smoke
    test, not a convergence study: we only verify the operator works
    on a quad mesh, not that convergence rates are identical (the two
    element families differ in a non-trivial way).
    """
    from semi.verification.mms_poisson import (
        MMSPoissonCase,
        run_convergence_study,
        run_one_level,
    )

    quad = run_one_level(MMSPoissonCase(
        dim=2, N=64, amplitude=0.5, cell_kind="quadrilateral",
    ))
    tri = run_convergence_study(
        dim=2, Ns=[64], amplitude=0.5, cell_kind="triangle",
    )[0]
    ratio = quad.e_L2 / tri.e_L2
    # Quad on a square has half the DOFs of right-diagonal triangles at
    # the same N, so some error inflation is expected; a ratio in
    # [0.1, 10] range means the two are comparable, not pathological.
    assert 0.1 < ratio < 10.0, (
        f"quad/triangle L^2 ratio {ratio:.3f} outside [0.1, 10] at N=64"
    )


def test_mms_poisson_2d_multiregion_convergence():
    """
    2D multi-region (Si/SiO2) coefficient-jump assembly.

    The exact solution is C^0 across the interface y = y_int and
    satisfies eps_r * (d psi / d y) continuity by construction, so the
    FE problem is well-posed and the finest-pair L^2 rate must hit
    theoretical 2.0. Pytest uses a 1.85 floor so solver-tolerance
    jitter does not destabilise the gate; the stricter 1.99 floor
    lives in run_verification.py per the M6 acceptance criterion
    in docs/mos_derivation.md section 7.
    """
    from semi.verification.mms_poisson import run_mr_convergence_study

    # N_y must be divisible by 10 so y_int lands on a grid line.
    results = run_mr_convergence_study([10, 20, 40, 80])
    hs = [r.h for r in results]
    eL2 = [r.e_L2 for r in results]
    eH1 = [r.e_H1 for r in results]

    assert all(eL2[i + 1] < eL2[i] for i in range(len(eL2) - 1)), (
        f"L^2 error should decrease monotonically: {eL2}"
    )

    rate_L2 = _finest_pair_rate(hs, eL2)
    rate_H1 = _finest_pair_rate(hs, eH1)
    assert rate_L2 >= 1.85, (
        f"multi-region finest-pair L^2 rate {rate_L2:.3f} < 1.85"
    )
    assert rate_H1 >= 0.85, (
        f"multi-region finest-pair H^1 rate {rate_H1:.3f} < 0.85"
    )


@pytest.mark.parametrize("dim", [1, 2])
def test_mms_poisson_snes_converges(dim):
    """
    Newton converges in a small number of iterations on a representative
    mesh level for each dim. Zero iterations is allowed (initial guess
    can be exact in the linearized regime) but more than 20 indicates
    something has regressed in the residual or initial guess.
    """
    from semi.verification.mms_poisson import MMSPoissonCase, run_one_level

    N = 80 if dim == 1 else 32
    result = run_one_level(MMSPoissonCase(dim=dim, N=N, amplitude=0.5))
    assert result.snes_iters < 20, (
        f"SNES took {result.snes_iters} iterations on dim={dim} N={N}; "
        "expected single digits"
    )
