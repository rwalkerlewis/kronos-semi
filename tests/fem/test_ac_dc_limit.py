"""
M14 acceptance test: low-frequency AC capacitance matches the analytical
depletion capacitance of a 1D pn junction within 5%.

Runs the AC sweep at f = 1 Hz on the same device configured by
`benchmarks/pn_1d_bias_reverse` and compares C(1 Hz) to the closed-form
junction capacitance

    C_dep(V) = sqrt(q * eps_Si * N_eff / (2 * (V_bi - V)))
    N_eff    = N_A * N_D / (N_A + N_D)
    V_bi     = V_t * ln(N_A * N_D / n_i^2)

at three reverse biases: V_DC in {-2.0, -1.0, -0.5} V.

This is the M14 acceptance criterion #2 from
`docs/IMPROVEMENT_GUIDE.md`. A failure here means the AC formulation,
the displacement-current evaluator, or the (J + j*omega*M) assembly is
wrong; in that case the rc_ac_sweep benchmark verifier will also fail
with the same root cause.
"""
from __future__ import annotations


def _diode_cfg(V_DC: float) -> dict:
    return {
        "schema_version": "1.2.0",
        "name": f"ac_dc_limit_{V_DC:+.3f}",
        "description": "AC depletion-C limit test fixture (1D pn).",
        "dimension": 1,
        "mesh": {
            "source": "builtin",
            "extents": [[0.0, 2.0e-5]],
            "resolution": [800],
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
                "frequencies": {"type": "list", "values": [1.0]},
            },
        },
        "output": {"directory": "/tmp/ac_dc_limit", "fields": []},
    }


def _analytical_C_dep(V_DC: float) -> float:
    import numpy as np

    from semi.constants import EPS0, KB, Q, cm3_to_m3
    from semi.materials import get_material

    mat = get_material("Si")
    eps = EPS0 * mat.epsilon_r
    n_i = mat.n_i
    N_A = cm3_to_m3(1.0e17)
    N_D = cm3_to_m3(1.0e17)
    N_eff = N_A * N_D / (N_A + N_D)
    V_t = KB * 300.0 / Q
    V_bi = V_t * np.log(N_A * N_D / (n_i * n_i))
    return float(np.sqrt(Q * eps * N_eff / (2.0 * (V_bi - V_DC))))


def test_ac_dc_limit_pn_depletion_capacitance():
    """C(1 Hz) matches analytical depletion capacitance within 5%."""
    from semi import schema
    from semi.runners.ac_sweep import run_ac_sweep

    rel_tol = 0.05
    for V_DC in (-2.0, -1.0, -0.5):
        cfg = schema.validate(_diode_cfg(V_DC=V_DC))
        result = run_ac_sweep(cfg)
        C_sim = float(result.C[0])
        C_an = _analytical_C_dep(V_DC)
        rel_err = abs(C_sim - C_an) / abs(C_an)
        assert rel_err < rel_tol, (
            f"AC depletion C at V_DC={V_DC:+.2f} V outside {100*rel_tol:.0f}% "
            f"tolerance: C_sim={C_sim:.4e} F/m^2, "
            f"C_an={C_an:.4e} F/m^2, rel_err={100*rel_err:.2f}%"
        )
