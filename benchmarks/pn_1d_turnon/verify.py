"""
Verifier for the pn_1d_turnon benchmark.

Physical basis
--------------
A pn junction is initially at equilibrium (V=0). At t=0+, a forward
bias of V_F = 0.6 V is applied. The anode current J(t) rises as minority
carriers are injected and approach the steady-state distribution. The
time constant of this turn-on transient is related to the minority-carrier
lifetime tau_p.

In the long-diode limit the transient current response is approximately:

    J_anode(t) = J_ss * (1 - exp(-t / tau_eff))

where tau_eff is the effective minority-carrier lifetime extracted from a
least-squares fit to J(t) in the window [tau_p, 4*tau_p]. The benchmark
asserts |tau_eff - tau_p| / tau_p < 0.05 (5 % relative error).

Usage
-----
    python verify.py <result_iv_json>

or import :func:`verify_turnon` and call it directly with the IV data.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from typing import Any


TAU_P_TARGET = 1.0e-7   # SRH hole lifetime configured in config.json, s
TOLERANCE = 0.05        # 5 % relative tolerance on extracted lifetime


def _fit_lifetime(t_arr, J_arr, tau_p: float) -> float:
    """
    Fit J(t) = J_ss * (1 - exp(-t/tau)) over [tau_p, 4*tau_p] using
    a linearised least-squares approach.

    Returns
    -------
    float
        Extracted effective lifetime tau_eff in seconds.
    """
    import numpy as np

    # Window selection
    t_min = tau_p
    t_max = 4.0 * tau_p
    mask = (t_arr >= t_min) & (t_arr <= t_max)
    if mask.sum() < 4:
        raise RuntimeError(
            f"Too few IV points in fit window [{t_min:.2e}, {t_max:.2e}] s; "
            f"got {mask.sum()}. Check t_end and dt in config.json."
        )
    t_fit = t_arr[mask]
    J_fit = J_arr[mask]

    # Estimate J_ss as the value near the end of the fit window
    J_ss_est = float(J_arr[t_arr >= 4.0 * tau_p][0]) if (t_arr >= 4.0 * tau_p).any() else float(J_fit[-1])
    if J_ss_est <= 0.0:
        raise RuntimeError("J_ss estimate is non-positive; check IV data.")

    # Linearise: log(1 - J/J_ss) = -t/tau  => slope = -1/tau
    y = np.log(np.maximum(1.0 - J_fit / J_ss_est, 1.0e-15))
    # Linear regression: y ~ -t/tau
    A = np.column_stack([t_fit, np.ones_like(t_fit)])
    coeffs, _, _, _ = np.linalg.lstsq(A, y, rcond=None)
    slope = coeffs[0]
    if slope >= 0.0:
        raise RuntimeError(
            f"Fitted slope {slope:.3e} is non-negative; current is not rising. "
            f"Check that t_end > 5*tau_p and V_F is well above built-in."
        )
    tau_eff = -1.0 / slope
    return float(tau_eff)


def verify_turnon(
    iv: list[dict[str, Any]],
    tau_p: float = TAU_P_TARGET,
    tol: float = TOLERANCE,
    contact: str = "anode",
    verbose: bool = True,
) -> dict[str, float]:
    """
    Verify the pn diode turn-on transient.

    Parameters
    ----------
    iv : list of dict
        IV records from :class:`TransientResult.iv`. Each entry must
        have ``t``, ``contact``, and ``J`` keys.
    tau_p : float
        Target minority-carrier lifetime (SRH tau_p), seconds.
    tol : float
        Allowed relative error |tau_eff - tau_p| / tau_p.
    contact : str
        Contact name to analyse (default ``"anode"``).
    verbose : bool
        If True, print a summary table.

    Returns
    -------
    dict
        {"tau_eff": float, "tau_p": float, "rel_error": float,
         "passed": bool}
    """
    import numpy as np

    rows = [r for r in iv if r.get("contact") == contact and r["t"] > 0.0]
    if not rows:
        raise RuntimeError(
            f"No IV data found for contact {contact!r} at t > 0."
        )
    t_arr = np.array([r["t"] for r in rows])
    J_arr = np.array([r["J"] for r in rows])

    # Ensure J is positive (forward current out of anode)
    if J_arr.mean() < 0.0:
        J_arr = -J_arr

    tau_eff = _fit_lifetime(t_arr, J_arr, tau_p)
    rel_err = abs(tau_eff - tau_p) / tau_p
    passed = rel_err < tol

    if verbose:
        print(f"pn_1d_turnon verification")
        print(f"  tau_p target : {tau_p:.3e} s")
        print(f"  tau_eff fit  : {tau_eff:.3e} s")
        print(f"  rel error    : {rel_err:.4f} (tolerance {tol:.2f})")
        print(f"  result       : {'PASS' if passed else 'FAIL'}")

    return {
        "tau_eff": tau_eff,
        "tau_p": tau_p,
        "rel_error": rel_err,
        "passed": passed,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description="Verify pn_1d_turnon benchmark.")
    parser.add_argument("iv_json", help="Path to IV JSON file from TransientResult.")
    args = parser.parse_args(argv)

    with open(args.iv_json) as f:
        iv = json.load(f)

    result = verify_turnon(iv, verbose=True)
    if not result["passed"]:
        sys.exit(1)


if __name__ == "__main__":
    main()
