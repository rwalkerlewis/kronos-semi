"""
FEM C-V regression test for the axisymmetric MOSCAP benchmark.

Compares an FEM-extracted C-V curve at
``benchmarks/moscap_axisym_2d/fem_cv.csv`` against the analytical
reference at ``benchmarks/moscap_axisym_2d/reference_cv.csv`` on:

  1. C_HF_min within 2 % of analytical C_min.
  2. C_LF within 2 % of C_ox in deep accumulation
     (V_g <= V_fb - 0.5 V) and deep inversion (V_g >= V_t + 1.0 V).
  3. LF and HF coincide within 1 % for V_g < V_t - 0.1 V (depletion).

The test SKIPS cleanly (does not fail) when ``fem_cv.csv`` is absent
so contributors without dolfinx are not blocked. The FEM CSV is
produced by running ``notebooks/05_moscap_axisym_cv.ipynb`` end-to-end;
this is currently a follow-up (see the post-merge cleanup PR's Open
Questions).

Pure-Python: only requires numpy and the analytical helpers in
``semi.cv``. The CSV format is the same as ``reference_cv.csv``:
columns ``Vg_V, C_LF_F_per_m2, C_HF_F_per_m2, C_LF_norm, C_HF_norm``
with leading ``#``-prefixed comment lines.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
BENCH_DIR = REPO_ROOT / "benchmarks" / "moscap_axisym_2d"
FEM_CSV = BENCH_DIR / "fem_cv.csv"
REF_CSV = BENCH_DIR / "reference_cv.csv"


def _load_cv(path: Path) -> np.ndarray:
    """Load the (Vg, C_LF, C_HF, ...) CSV, skipping ``#``-prefixed lines."""
    with path.open() as fh:
        lines = [line for line in fh if not line.lstrip().startswith(("#", '"#'))]
    # The first non-comment line is the header; skip it.
    if lines and not lines[0][0].isdigit() and not lines[0][0] in "+-.":
        lines = lines[1:]
    data = np.loadtxt(lines, delimiter=",")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data


def _moscap_anchors() -> dict[str, float]:
    """Recompute the analytical anchors that the benchmark targets."""
    from semi.cv import analytical_moscap_params

    m = analytical_moscap_params(
        body_dopant="p",
        N_body_cm3=5.0e16,
        T_ox_m=10.0e-9,
        phi_ms=-0.95,
    )
    return {
        "V_fb": m.V_fb,
        "V_t": m.V_t,
        "C_ox": m.C_ox_per_area,
        "C_min": m.C_min_per_area,
    }


@pytest.fixture(scope="module")
def fem_cv() -> np.ndarray:
    if not FEM_CSV.exists():
        pytest.skip(
            f"FEM C-V CSV not found at {FEM_CSV}; run "
            "notebooks/05_moscap_axisym_cv.ipynb end-to-end to produce it."
        )
    return _load_cv(FEM_CSV)


@pytest.fixture(scope="module")
def anchors() -> dict[str, float]:
    return _moscap_anchors()


def test_chf_min_matches_analytical(fem_cv, anchors):
    """C_HF_min from the FEM curve should be within 2 % of analytical C_min."""
    C_HF = fem_cv[:, 2]
    chf_min = float(np.min(C_HF))
    rel_err = abs(chf_min - anchors["C_min"]) / anchors["C_min"]
    assert rel_err < 0.02, (
        f"C_HF_min={chf_min:.4e} F/m^2 deviates from analytical "
        f"C_min={anchors['C_min']:.4e} by {100 * rel_err:.2f}% (> 2%)."
    )


def test_clf_plateau_in_accumulation_and_inversion(fem_cv, anchors):
    """C_LF should saturate to C_ox in deep accumulation and deep inversion."""
    Vg = fem_cv[:, 0]
    C_LF = fem_cv[:, 1]
    C_ox = anchors["C_ox"]

    acc_mask = Vg <= anchors["V_fb"] - 0.5
    inv_mask = Vg >= anchors["V_t"] + 1.0

    if not acc_mask.any():
        pytest.skip("FEM sweep does not extend deep enough into accumulation.")
    if not inv_mask.any():
        pytest.skip("FEM sweep does not extend deep enough into inversion.")

    acc_err = float(np.max(np.abs(C_LF[acc_mask] - C_ox) / C_ox))
    inv_err = float(np.max(np.abs(C_LF[inv_mask] - C_ox) / C_ox))

    assert acc_err < 0.02, (
        f"C_LF plateau in accumulation deviates from C_ox by "
        f"{100 * acc_err:.2f}% (> 2%)."
    )
    assert inv_err < 0.02, (
        f"C_LF plateau in inversion deviates from C_ox by "
        f"{100 * inv_err:.2f}% (> 2%)."
    )


def test_lf_hf_coincide_in_depletion(fem_cv, anchors):
    """LF and HF should agree to within 1 % for V_g < V_t - 0.1 V."""
    Vg = fem_cv[:, 0]
    C_LF = fem_cv[:, 1]
    C_HF = fem_cv[:, 2]

    mask = Vg < anchors["V_t"] - 0.1
    if not mask.any():
        pytest.skip("FEM sweep does not include any depletion-regime points.")

    diff = np.abs(C_LF[mask] - C_HF[mask]) / np.maximum(C_LF[mask], 1e-30)
    worst = float(np.max(diff))
    assert worst < 0.01, (
        f"LF/HF coincidence in depletion fails: worst relative "
        f"difference {100 * worst:.2f}% (> 1%)."
    )


def test_reference_csv_loads_and_anchors_match():
    """Sanity: the analytical reference CSV is parseable and aligns with semi.cv."""
    assert REF_CSV.exists(), f"reference_cv.csv missing at {REF_CSV}"
    data = _load_cv(REF_CSV)
    assert data.shape[1] >= 5
    anchors = _moscap_anchors()
    # The reference is normalised; final HF point at strong inversion must
    # equal C_min within tight tolerance.
    C_HF_min_ref = float(np.min(data[:, 2]))
    rel = abs(C_HF_min_ref - anchors["C_min"]) / anchors["C_min"]
    assert rel < 1e-3, (
        f"reference_cv.csv C_HF_min ({C_HF_min_ref:.4e}) and "
        f"semi.cv C_min ({anchors['C_min']:.4e}) disagree by "
        f"{100 * rel:.4f}%."
    )
