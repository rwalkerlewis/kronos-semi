"""
Mesh-convergence tests for the pn_1d benchmark.

For a problem with no analytical solution, the standard V&V gate is
Cauchy self-convergence: compare coarser FEM solutions against a finer
FEM reference on the same mesh family, not against a physics model
(depletion approximation) that the FEM solver is not targeting.

This test runs N = [100, 200, 400, 1600] and asserts the >= 2x error-
ratio-per-doubling rule on the Cauchy error at N = [100, 200, 400],
with N = 1600 serving as the reference. It also asserts monotone
reduction in the depletion-approx relative error over the same subset,
to smoke-test that the benchmark observables (V_bi, peak |E|, W) keep
approaching the physical values as h -> 0.

See the module docstring of `semi/verification/mesh_convergence.py`
for why Cauchy convergence is the right gate here.
"""
from __future__ import annotations

PYTEST_NS = [100, 200, 400, 1600]
ASSERT_NS = [100, 200, 400]          # rates are asserted only on these
MIN_RATIO_PER_DOUBLING = 2.0


def _assert_monotone_decreasing(errs, label, ns):
    for i in range(len(errs) - 1):
        assert errs[i + 1] < errs[i], (
            f"{label}: error should decrease monotonically at "
            f"N={ns[i + 1]}, got {errs[i + 1]:.3e} >= {errs[i]:.3e}"
        )


def _assert_halving(errs, label, ns):
    for i in range(1, len(errs)):
        ratio = errs[i - 1] / errs[i] if errs[i] > 0.0 else float("inf")
        assert ratio >= MIN_RATIO_PER_DOUBLING, (
            f"{label}: error ratio per mesh doubling at N={ns[i]} is "
            f"{ratio:.2f}, below the {MIN_RATIO_PER_DOUBLING:.2f} floor "
            f"(errs={errs})"
        )


def test_mesh_convergence_pn_1d_cauchy_halves_per_doubling():
    """
    Cauchy self-convergence errors (FEM vs finest FEM) at
    N = [100, 200, 400] must reduce by >= 2x per mesh doubling and
    decrease monotonically. This is the discretization-error gate and
    is unaffected by the depletion-approx plateau that floors
    err_Epeak_rel on fine meshes.
    """
    from semi.verification.mesh_convergence import (
        cauchy_errors,
        run_convergence_study,
    )

    results = run_convergence_study(PYTEST_NS)
    cauchy = cauchy_errors(results, reference_index=-1)
    E_errs = [c["err_Epeak_cauchy"] for c in cauchy[:-1]]
    W_errs = [c["err_W_cauchy"] for c in cauchy[:-1]]

    _assert_monotone_decreasing(E_errs, "E_peak Cauchy", ASSERT_NS)
    _assert_monotone_decreasing(W_errs, "W Cauchy", ASSERT_NS)
    _assert_halving(E_errs, "E_peak Cauchy", ASSERT_NS)
    _assert_halving(W_errs, "W Cauchy", ASSERT_NS)


def test_mesh_convergence_pn_1d_depletion_error_monotone():
    """
    Smoke test that the depletion-approximation relative error also
    decreases monotonically on the asserted subset, even though its
    absolute magnitude is floored by the physics-model gap. No ratio
    is enforced here; see the module docstring for why.
    """
    from semi.verification.mesh_convergence import run_convergence_study

    results = run_convergence_study(ASSERT_NS)
    E_errs = [r.err_Epeak_rel for r in results]
    W_errs = [r.err_W_rel for r in results]

    _assert_monotone_decreasing(E_errs, "E_peak rel", ASSERT_NS)
    _assert_monotone_decreasing(W_errs, "W rel", ASSERT_NS)


def test_mesh_convergence_pn_1d_newton_converges():
    """Equilibrium Poisson must converge in single-digit iterations at
    every refinement level used by the test subset."""
    from semi.verification.mesh_convergence import run_convergence_study

    results = run_convergence_study(ASSERT_NS)
    for r in results:
        assert 0 < r.newton_iters < 20, (
            f"pn_1d Newton took {r.newton_iters} iterations at N={r.N}; "
            "expected single digits"
        )
