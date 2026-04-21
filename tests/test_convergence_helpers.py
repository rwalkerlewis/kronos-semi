"""
Pure-Python tests for `semi.verification._convergence` helpers.

Covers `observed_rates`, `format_table`, `write_convergence_csv`, and
`write_loglog_plot`. These helpers run on every V&V CLI invocation but
are not exercised by the convergence-study tests (which inspect the
`.h` and `.e_L2` results directly), so they need their own coverage.
"""
from __future__ import annotations

import math
from pathlib import Path

import pytest

from semi.verification._convergence import (
    format_table,
    observed_rates,
    write_convergence_csv,
    write_loglog_plot,
)


def test_observed_rates_two_pairs_recovers_quadratic():
    # Halve h; quarter the error -> rate 2 by construction.
    rates = observed_rates([1.0, 0.5], [1.0, 0.25])
    assert math.isnan(rates[0])
    assert rates[1] == pytest.approx(2.0)


def test_observed_rates_three_pairs():
    rates = observed_rates([1.0, 0.5, 0.25], [1.0, 0.25, 0.0625])
    assert math.isnan(rates[0])
    assert rates[1] == pytest.approx(2.0)
    assert rates[2] == pytest.approx(2.0)


def test_observed_rates_handles_nonpositive_error_with_nan():
    rates = observed_rates([1.0, 0.5, 0.25], [1.0, 0.0, 0.25])
    # Pairs touching the zero sample produce NaN (cannot take log).
    assert math.isnan(rates[1])
    assert math.isnan(rates[2])


def test_observed_rates_single_point_returns_nan():
    rates = observed_rates([1.0], [1.0])
    assert len(rates) == 1
    assert math.isnan(rates[0])


def test_observed_rates_length_mismatch_raises():
    with pytest.raises(ValueError):
        observed_rates([1.0, 0.5], [1.0])


def test_observed_rates_equal_h_returns_nan():
    rates = observed_rates([1.0, 1.0], [1.0, 0.5])
    # log(h_prev/h_cur) = 0 -> NaN (avoid div-by-zero).
    assert math.isnan(rates[1])


def test_format_table_renders_header_columns_and_numbers():
    rows = [
        {"N": 16, "h": 1.0e-2, "e_L2": 5.0e-3, "rate_L2": float("nan")},
        {"N": 32, "h": 5.0e-3, "e_L2": 1.25e-3, "rate_L2": 2.0},
    ]
    out = format_table(rows, ["N", "h", "e_L2", "rate_L2"], header="MMS sweep")
    lines = out.splitlines()
    assert lines[0] == "MMS sweep"
    assert "N" in lines[1] and "rate_L2" in lines[1]
    # First-row rate is NaN -> rendered as "---".
    assert "---" in lines[2]
    # Second-row rate is 2.0 -> rendered with the rate formatter (3 dp).
    assert "2.000" in lines[3]
    # Scientific notation for floats that are not rates.
    assert "5.000e-03" in lines[2]


def test_format_table_no_header_omits_first_line():
    rows = [{"a": 1.0, "b": "x"}]
    out = format_table(rows, ["a", "b"], header="")
    assert "x" in out
    assert "a" in out.splitlines()[0]


def test_write_convergence_csv_round_trips_columns(tmp_path: Path):
    rows = [
        {"N": 16, "h": 1.0e-2, "e_L2": 5.0e-3},
        {"N": 32, "h": 5.0e-3, "e_L2": 1.25e-3},
    ]
    csv_path = tmp_path / "out" / "sweep.csv"
    write_convergence_csv(csv_path, rows, ["N", "h", "e_L2"])
    text = csv_path.read_text().splitlines()
    assert text[0] == "N,h,e_L2"
    assert text[1] == "16,0.01,0.005"
    assert text[2] == "32,0.005,0.00125"


def test_write_convergence_csv_creates_parent_directory(tmp_path: Path):
    csv_path = tmp_path / "deeply" / "nested" / "out.csv"
    write_convergence_csv(csv_path, [{"a": 1}], ["a"])
    assert csv_path.exists()


def test_write_convergence_csv_omits_unknown_columns(tmp_path: Path):
    rows = [{"a": 1, "extra": "ignored"}]
    csv_path = tmp_path / "out.csv"
    write_convergence_csv(csv_path, rows, ["a", "missing"])
    lines = csv_path.read_text().splitlines()
    assert lines[0] == "a,missing"
    assert lines[1] == "1,"


def test_write_loglog_plot_writes_png_with_reference_overlay(tmp_path: Path):
    png = tmp_path / "rates.png"
    write_loglog_plot(
        png,
        hs=[1.0, 0.5, 0.25, 0.125],
        series={"L2": [1.0, 0.25, 0.0625, 0.015625]},
        title="quadratic",
        theoretical_rates={"L2": 2.0},
    )
    assert png.exists()
    assert png.stat().st_size > 0


def test_write_loglog_plot_handles_no_finite_anchor(tmp_path: Path):
    png = tmp_path / "all_nan.png"
    nan = float("nan")
    write_loglog_plot(
        png,
        hs=[1.0, 0.5],
        series={"L2": [nan, nan]},
        title="no anchor",
        theoretical_rates={"L2": 2.0},
    )
    assert png.exists()


def test_write_loglog_plot_no_theoretical_rates(tmp_path: Path):
    png = tmp_path / "noref.png"
    write_loglog_plot(
        png,
        hs=[1.0, 0.5],
        series={"L2": [1.0, 0.25]},
        title="no overlay",
    )
    assert png.exists()
