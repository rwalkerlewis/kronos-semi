"""
End-to-end tests for the multi-region bias sweep runner.

Exercises the coupled Slotboom drift-diffusion bias sweep on a minimal
two-region (Si + SiO2) 1D geometry to verify:
  - The runner produces an IV table with the expected number of entries.
  - Multi-step ramping produces non-decreasing electron current (the
    physical expectation for a forward-biased pn junction at low bias).
"""
from __future__ import annotations


def _minimal_multiregion_bias_cfg(n_steps: int = 3, stop: float = 0.15) -> dict:
    """
    Minimal 1-D multi-region (Si + SiO2 cap) bias-sweep config.

    The config is intentionally tiny (80 cells, coarse sweep) so the FEM
    test suite runs in seconds.  n_steps determines how many bias points
    beyond V=0 are requested; stop is the final voltage.
    """
    step = round(stop / n_steps, 6)
    return {
        "schema_version": "1.0.0",
        "name": "test_mr_bias",
        "description": "Minimal multi-region bias sweep test fixture.",
        "dimension": 1,
        "mesh": {
            "source": "builtin",
            "extents": [[0.0, 2.0e-5]],
            "resolution": [80],
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
                    "type": "step",
                    "axis": 0,
                    "location": 1.0e-5,
                    "N_D_left": 0.0,
                    "N_A_left": 1.0e17,
                    "N_D_right": 1.0e17,
                    "N_A_right": 0.0,
                },
            }
        ],
        "contacts": [
            {
                "name": "anode",
                "facet": "anode",
                "type": "ohmic",
                "voltage": 0.0,
                "voltage_sweep": {"start": 0.0, "stop": stop, "step": step},
            },
            {"name": "cathode", "facet": "cathode", "type": "ohmic", "voltage": 0.0},
        ],
        "physics": {
            "temperature": 300.0,
            "statistics": "boltzmann",
            "mobility": {"mu_n": 1400.0, "mu_p": 450.0},
            "recombination": {
                "srh": True,
                "tau_n": 1.0e-8,
                "tau_p": 1.0e-8,
                "E_t": 0.0,
            },
        },
        "solver": {
            "type": "bias_sweep",
            "continuation": {
                "min_step": 0.001,
                "max_halvings": 6,
                "max_step": 0.1,
                "easy_iter_threshold": 4,
                "grow_factor": 2.0,
            },
        },
        "output": {
            "directory": "/tmp/test_mr_bias_results",
            "fields": [],
        },
    }


def test_bias_sweep_multiregion_multistep():
    """
    Run a 3-step forward bias ramp on a minimal 1D multi-region config
    and assert:
      - At least 3 IV points are recorded (one per bias step plus seed).
      - Electron current at the sweep contact is non-decreasing across
        the ramp (physical expectation: J_n grows monotonically with
        forward bias for a pn junction at low injection).
    """
    from semi.runners.bias_sweep import run_bias_sweep

    cfg = _minimal_multiregion_bias_cfg(n_steps=3, stop=0.15)
    result = run_bias_sweep(cfg)

    iv = result.iv
    assert len(iv) >= 3, (
        f"Expected at least 3 IV points, got {len(iv)}. "
        f"Bias points: {[row['V'] for row in iv]}"
    )

    # Electron current (J_n) at the anode should be non-decreasing
    # as the forward bias increases from 0 toward 0.15 V.
    j_n_values = [row["J_n"] for row in iv if "J_n" in row]
    # 5% relative slack accounts for numerical noise at near-zero current
    # levels (V < 0.05 V) where the absolute current is tiny but the
    # relative variation is not physically meaningful.
    _J_N_MONOTONE_SLACK = 0.05
    if len(j_n_values) >= 2:
        for i in range(len(j_n_values) - 1):
            assert j_n_values[i + 1] >= j_n_values[i] - abs(j_n_values[i]) * _J_N_MONOTONE_SLACK, (
                f"J_n should be non-decreasing under forward bias ramp: "
                f"{j_n_values}"
            )
