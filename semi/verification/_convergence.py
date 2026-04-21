"""
Convergence-rate bookkeeping shared across verification studies.

Given a sequence of (h, error) pairs from a mesh-refinement sweep,
compute observed rates, format a stdout table, and write a CSV plus
log-log PNG plot. None of these are physics-specific.
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np


def observed_rates(hs: list[float], errors: list[float]) -> list[float]:
    """
    Element-wise log slope between consecutive (h, error) pairs.

    Returned list has length len(hs); the first entry is NaN because
    no previous pair exists. This matches the human-specified DataFrame
    convention where rate_L2[0] = NaN and tests assert on .iloc[-1].
    """
    if len(hs) != len(errors):
        raise ValueError("hs and errors must have equal length")
    if len(hs) < 2:
        return [float("nan")] * len(hs)
    rates: list[float] = [float("nan")]
    for i in range(1, len(hs)):
        h_prev, h_cur = hs[i - 1], hs[i]
        e_prev, e_cur = errors[i - 1], errors[i]
        if e_prev <= 0.0 or e_cur <= 0.0 or h_prev == h_cur:
            rates.append(float("nan"))
        else:
            rates.append(math.log(e_prev / e_cur) / math.log(h_prev / h_cur))
    return rates


def write_convergence_csv(
    path: Path,
    rows: list[dict],
    columns: list[str],
) -> None:
    """Write a CSV with an explicit column order. No pandas dependency."""
    import csv

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({c: row.get(c, "") for c in columns})


def write_loglog_plot(
    path: Path,
    hs: list[float],
    series: dict[str, list[float]],
    *,
    title: str,
    theoretical_rates: dict[str, float] | None = None,
) -> None:
    """
    Write a log-log error-vs-h PNG. `series` maps label to error sequence;
    optional `theoretical_rates` overlays guide lines anchored at the
    first finite (h, error) point of the matching series.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6, 4))
    h_arr = np.asarray(hs, dtype=float)
    for label, errs in series.items():
        e = np.asarray(errs, dtype=float)
        ax.loglog(h_arr, e, "o-", label=label)
        if theoretical_rates and label in theoretical_rates:
            slope = theoretical_rates[label]
            anchor_idx = next(
                (i for i, v in enumerate(e) if np.isfinite(v) and v > 0.0), None
            )
            if anchor_idx is not None:
                e_ref = e[anchor_idx] * (h_arr / h_arr[anchor_idx]) ** slope
                ax.loglog(
                    h_arr, e_ref, "--", alpha=0.5,
                    label=f"O(h^{slope:.0f}) reference",
                )
    ax.set_xlabel("h (mesh size, m)")
    ax.set_ylabel("error norm")
    ax.set_title(title)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def format_table(rows: list[dict], columns: list[str], header: str = "") -> str:
    """Pretty multi-line table for stdout. Numbers use scientific notation."""
    widths = [max(len(c), 12) for c in columns]
    lines: list[str] = []
    if header:
        lines.append(header)
    lines.append("  ".join(c.rjust(w) for c, w in zip(columns, widths, strict=False)))
    for row in rows:
        cells = []
        for c, w in zip(columns, widths, strict=False):
            v = row.get(c, "")
            if isinstance(v, float):
                if math.isnan(v):
                    s = "      ---"
                elif c.startswith("rate"):
                    s = f"{v:8.3f}"
                else:
                    s = f"{v:.3e}"
            else:
                s = str(v)
            cells.append(s.rjust(w))
        lines.append("  ".join(cells))
    return "\n".join(lines)
