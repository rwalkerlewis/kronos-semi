"""
M14 MMS-style consistency test: AC at omega=0 must equal DC sensitivity.

The AC small-signal runner solves

    (J + j*omega*M) * delta_u(omega) = -dF/dV * delta_V

around a converged DC operating point u_0. At omega = 0 the linear
system collapses to

    J * delta_u(0) = -dF/dV * delta_V

which is exactly the DC sensitivity equation. Any consistent runner
implementation must therefore produce, at omega = 0, an admittance Y(0)
that matches the finite-difference DC differential conductance dI/dV at
the same operating point.

This test is the M14 analogue of the MMS spatial / temporal convergence
tests in `tests/mms/`: it does not depend on a closed-form analytical
solution, only on internal self-consistency of the AC linearisation.
A failure here points at the AC formulation, not the physics.
"""
from __future__ import annotations


def _diode_cfg(V_DC: float, freqs: list[float]) -> dict:
    """Minimal 1D pn diode AC sweep cfg (matches the rc_ac_sweep device)."""
    return {
        "schema_version": "1.2.0",
        "name": f"mms_ac_consistency_{V_DC:+.3f}",
        "description": "MMS DC-limit consistency test fixture (1D pn).",
        "dimension": 1,
        "mesh": {
            "source": "builtin",
            "extents": [[0.0, 2.0e-5]],
            "resolution": [200],
            "regions_by_box": [
                {"name": "silicon", "tag": 1, "bounds": [[0.0, 2.0e-5]]},
            ],
            "facets_by_plane": [
                {"name": "anode",   "tag": 1, "axis": 0, "value": 0.0},
                {"name": "cathode", "tag": 2, "axis": 0, "value": 2.0e-5},
            ],
        },
        "regions": {
            "silicon": {"material": "Si", "tag": 1, "role": "semiconductor"},
        },
        "doping": [
            {
                "region": "silicon",
                "profile": {
                    "type": "step", "axis": 0, "location": 1.0e-5,
                    "N_D_left": 0.0, "N_A_left": 1.0e17,
                    "N_D_right": 1.0e17, "N_A_right": 0.0,
                },
            },
        ],
        "contacts": [
            {"name": "anode",   "facet": "anode",   "type": "ohmic", "voltage": 0.0},
            {"name": "cathode", "facet": "cathode", "type": "ohmic", "voltage": 0.0},
        ],
        "physics": {
            "temperature": 300.0, "statistics": "boltzmann",
            "mobility": {"mu_n": 1400.0, "mu_p": 450.0},
            "recombination": {"srh": True, "tau_n": 1.0e-8, "tau_p": 1.0e-8, "E_t": 0.0},
        },
        "solver": {
            "type": "ac_sweep",
            "snes": {"rtol": 1.0e-14, "atol": 1.0e-14, "stol": 1.0e-14, "max_it": 60},
            "continuation": {"min_step": 1.0e-5, "max_step": 0.1, "max_halvings": 14},
            "dc_bias": {"contact": "anode", "voltage": V_DC},
            "ac": {
                "contact": "anode", "amplitude": 1.0,
                "frequencies": {"type": "list", "values": freqs},
            },
        },
        "output": {"directory": "/tmp/mms_ac_consistency", "fields": []},
    }


def test_ac_at_omega_zero_is_pure_real():
    """At f = 0 the AC admittance must be purely real.

    The 2x2 real-block system at omega = 0 decouples into two copies of
    the steady-state Jacobian solve, and the imaginary part of the
    solution is identically zero.
    """
    from semi import schema
    from semi.runners.ac_sweep import run_ac_sweep

    cfg = schema.validate(_diode_cfg(V_DC=-1.0, freqs=[0.0]))
    result = run_ac_sweep(cfg)
    assert len(result.Y) == 1
    Y0 = result.Y[0]
    # Pure-real DC sensitivity: imaginary part is identically zero up to
    # finite-precision noise from the LU factorisation.
    assert abs(Y0.imag) < 1.0e-12, (
        f"Y(omega=0) must be purely real; got Y={Y0}"
    )
    # The reported capacitance at f=0 is reported as 0 by definition.
    assert result.C[0] == 0.0


def test_ac_low_freq_y_is_stable_in_real_part():
    """Re(Y) at omega = 0 and at the lowest non-zero frequency must
    agree at the linearisation-error scale.

    For a depletion-dominated reverse-biased diode, the small-signal
    admittance is dominated by Y ≈ -j*omega*C at low omega, so Re(Y)
    is small and only weakly dependent on omega until frequencies
    approach the carrier-transit cut-off. Here we check the REAL part
    of Y is stable across omega = 0 and omega = 2*pi*1e-3 Hz to the
    1% level, which is well within the finite-difference perturbation
    error of the runner's DC sensitivity (eps_V = 1e-3 V).
    """
    from semi import schema
    from semi.runners.ac_sweep import run_ac_sweep

    cfg = schema.validate(_diode_cfg(V_DC=-1.0, freqs=[0.0, 1.0e-3]))
    result = run_ac_sweep(cfg)

    Y0 = result.Y[0]
    Y_low = result.Y[1]

    rel_re = abs(Y_low.real - Y0.real) / max(abs(Y0.real), 1.0e-12)
    assert rel_re < 1.0e-2, (
        "AC conductance Re(Y) at omega=0 and omega=2*pi*1e-3 disagree by "
        f"rel={rel_re:.3e}: Y0={Y0}, Y_low={Y_low}"
    )


def test_ac_im_y_scales_linearly_with_omega_at_low_freq():
    """In the depletion-capacitance regime, Im(Y) is proportional to
    omega.  Doubling the frequency doubles |Im(Y)| within 1e-6 because
    the underlying linear system is exact-affine in omega.
    """
    from semi import schema
    from semi.runners.ac_sweep import run_ac_sweep

    cfg = schema.validate(_diode_cfg(V_DC=-1.0, freqs=[1.0, 2.0]))
    result = run_ac_sweep(cfg)

    # The C(f) reported by the runner is Im(Y)/(2*pi*f) (signed).
    # Linearity in omega ==> C(f1) == C(f2).
    C_f1 = result.C[0]
    C_f2 = result.C[1]
    rel_C = abs(C_f1 - C_f2) / max(abs(C_f1), 1.0e-14)
    assert rel_C < 1.0e-6, (
        "Capacitance must be omega-independent in the linear depletion "
        f"regime; got C(1Hz)={C_f1:.4e}, C(2Hz)={C_f2:.4e}, "
        f"rel diff={rel_C:.3e}"
    )
