"""
Coverage tests for the M16.7 voltage_t evaluator and its integration
with the transient runner.

The evaluator (`_build_voltage_t_evaluator`) is a pure-Python helper:
the first six tests below exercise it directly without spinning up any
FEM machinery, so they run in the gated `docker-fem-tests` job at zero
mesh-build cost. The last two tests run `run_transient` end to end on
a tiny 1D pn diode (4 cells, 5-10 timesteps) and assert that the
recorded contact voltage in the IV row at each step matches the
expected V(t) within 1e-12 relative.

These tests are independent of audit case 06; they exist so the new
runner branches are covered by the gated suite, not by the audit job
or the benchmark matrix whose coverage is not merged. (M16.7)
"""
from __future__ import annotations

import copy
import math

import numpy as np
import pytest

from semi.runners.transient import (
    _build_voltage_t_evaluator,
    _ramp_target_voltage,
)


def _have_dolfinx() -> bool:
    try:
        import dolfinx  # noqa: F401
    except ImportError:
        return False
    return True


require_dolfinx = pytest.mark.skipif(
    not _have_dolfinx(),
    reason="dolfinx unavailable; FEM-runner tests skipped",
)


def _make_contacts(*specs):
    """Build a minimal `cfg` with one or more ohmic contacts.

    `specs` is a sequence of dicts merged into the contact stub; every
    contact is named by a single character.
    """
    contacts = []
    for i, spec in enumerate(specs):
        c = {"name": f"c{i}", "facet": f"f{i}", "type": "ohmic", "voltage": 0.0}
        c.update(spec)
        contacts.append(c)
    return {"contacts": contacts}


def test_evaluator_step_simple():
    cfg = _make_contacts(
        {"voltage_t": {"type": "step", "t0": 5.0e-9, "v0": 0.0, "v1": 0.6}},
    )
    f = _build_voltage_t_evaluator(cfg)
    assert f(0.0)["c0"] == pytest.approx(0.0, abs=1.0e-15)
    assert f(4.999e-9)["c0"] == pytest.approx(0.0, abs=1.0e-15)
    assert f(5.0e-9)["c0"] == pytest.approx(0.6, abs=1.0e-15)
    assert f(1.0e-8)["c0"] == pytest.approx(0.6, abs=1.0e-15)


def test_evaluator_step_at_t_zero():
    cfg = _make_contacts(
        {"voltage_t": {"type": "step", "t0": 0.0, "v0": 0.0, "v1": 0.6}},
    )
    f = _build_voltage_t_evaluator(cfg)
    assert f(0.0)["c0"] == pytest.approx(0.6, abs=1.0e-15)
    assert f(1.0e-9)["c0"] == pytest.approx(0.6, abs=1.0e-15)


def test_evaluator_table_linear_ramp():
    times = [0.0, 1.0e-9, 2.0e-9, 3.0e-9]
    values = [0.0, 0.2, 0.4, 0.6]
    cfg = _make_contacts(
        {"voltage_t": {"type": "table", "times": times, "values": values}},
    )
    f = _build_voltage_t_evaluator(cfg)
    assert f(0.0)["c0"] == pytest.approx(0.0, abs=1.0e-15)
    assert f(0.5e-9)["c0"] == pytest.approx(0.1, rel=1.0e-12)
    assert f(1.5e-9)["c0"] == pytest.approx(0.3, rel=1.0e-12)
    assert f(2.5e-9)["c0"] == pytest.approx(0.5, rel=1.0e-12)
    assert f(3.0e-9)["c0"] == pytest.approx(0.6, abs=1.0e-15)


def test_evaluator_table_sampled_sinusoid():
    """The audit-test waveform shape: the FFT bin at f recovers the
    fundamental amplitude to within roundoff when the table is sampled
    at an integer multiple of the period."""
    f_hz = 1.0e6
    n_per_cycle = 64
    n_cycles = 4
    n = n_per_cycle * n_cycles
    dt = 1.0 / (f_hz * n_per_cycle)
    times = np.arange(n) * dt
    v_dc = 0.4
    dv = 1.0e-3
    values = v_dc + dv * np.sin(2.0 * math.pi * f_hz * times)
    cfg = _make_contacts({
        "voltage_t": {
            "type": "table", "times": times.tolist(), "values": values.tolist(),
        },
    })
    f_eval = _build_voltage_t_evaluator(cfg)
    # Sampling at the same grid recovers the inputs to floating precision.
    sampled = np.array([f_eval(float(t))["c0"] for t in times])
    np.testing.assert_allclose(sampled, values, rtol=1.0e-12, atol=1.0e-15)
    # FFT bin alignment: the spectral peak is at exactly bin = n_cycles.
    spectrum = np.fft.rfft(values - v_dc)
    bin_amp = np.abs(spectrum)
    assert int(np.argmax(bin_amp)) == n_cycles
    assert bin_amp[n_cycles] == pytest.approx(0.5 * dv * n, rel=1.0e-9)


def test_evaluator_table_nonuniform_times():
    times = [0.0, 1.0e-9, 5.0e-9, 7.0e-9]
    values = [0.0, 0.1, 0.5, 0.7]
    cfg = _make_contacts(
        {"voltage_t": {"type": "table", "times": times, "values": values}},
    )
    f = _build_voltage_t_evaluator(cfg)
    # midpoint between (1e-9, 0.1) and (5e-9, 0.5) is (3e-9, 0.3)
    assert f(3.0e-9)["c0"] == pytest.approx(0.3, rel=1.0e-12)
    # midpoint between (5e-9, 0.5) and (7e-9, 0.7) is (6e-9, 0.6)
    assert f(6.0e-9)["c0"] == pytest.approx(0.6, rel=1.0e-12)


def test_evaluator_table_endpoint_clamp_low():
    times = [0.0, 1.0e-9, 2.0e-9]
    values = [0.1, 0.2, 0.3]
    cfg = _make_contacts(
        {"voltage_t": {"type": "table", "times": times, "values": values}},
    )
    f = _build_voltage_t_evaluator(cfg)
    assert f(-1.0e-9)["c0"] == pytest.approx(0.1, abs=1.0e-15)
    assert f(-1.0)["c0"] == pytest.approx(0.1, abs=1.0e-15)


def test_evaluator_table_endpoint_clamp_high():
    times = [0.0, 1.0e-9, 2.0e-9]
    values = [0.1, 0.2, 0.3]
    cfg = _make_contacts(
        {"voltage_t": {"type": "table", "times": times, "values": values}},
    )
    f = _build_voltage_t_evaluator(cfg)
    assert f(3.0e-9)["c0"] == pytest.approx(0.3, abs=1.0e-15)
    assert f(1.0)["c0"] == pytest.approx(0.3, abs=1.0e-15)


def test_evaluator_two_contacts_independent_voltage_t():
    cfg = _make_contacts(
        {"voltage_t": {"type": "step", "t0": 1.0e-9, "v0": 0.0, "v1": 0.6}},
        {"voltage_t": {
            "type": "table",
            "times": [0.0, 2.0e-9, 4.0e-9],
            "values": [-0.1, -0.2, -0.3],
        }},
    )
    f = _build_voltage_t_evaluator(cfg)
    out_pre = f(0.5e-9)
    out_post = f(3.0e-9)
    assert out_pre["c0"] == pytest.approx(0.0, abs=1.0e-15)
    assert out_pre["c1"] == pytest.approx(-0.125, rel=1.0e-12)
    assert out_post["c0"] == pytest.approx(0.6, abs=1.0e-15)
    assert out_post["c1"] == pytest.approx(-0.25, rel=1.0e-12)


def test_evaluator_no_voltage_t_returns_fixed():
    """Bit-identity guard: when no contact has voltage_t, the
    evaluator returns the same dict at every t (the v0.22.0 path)."""
    cfg = {"contacts": [
        {"name": "anode", "facet": "f", "type": "ohmic", "voltage": 0.4},
        {"name": "cathode", "facet": "g", "type": "ohmic", "voltage": 0.0},
    ]}
    f = _build_voltage_t_evaluator(cfg)
    assert f(0.0) == {"anode": 0.4, "cathode": 0.0}
    assert f(1.0e-6) == {"anode": 0.4, "cathode": 0.0}
    assert f(-1.0) == {"anode": 0.4, "cathode": 0.0}


def test_ramp_target_voltage_table():
    contact = {
        "name": "anode", "facet": "f", "type": "ohmic",
        "voltage_t": {
            "type": "table",
            "times": [0.0, 1.0e-9, 2.0e-9],
            "values": [0.4, 0.6, 0.8],
        },
    }
    assert _ramp_target_voltage(contact) == pytest.approx(0.4)


def test_ramp_target_voltage_step_pre_zero():
    contact = {
        "name": "anode", "facet": "f", "type": "ohmic",
        "voltage_t": {"type": "step", "t0": 5.0e-9, "v0": 0.1, "v1": 0.6},
    }
    assert _ramp_target_voltage(contact) == pytest.approx(0.1)


def test_ramp_target_voltage_step_at_zero():
    contact = {
        "name": "anode", "facet": "f", "type": "ohmic",
        "voltage_t": {"type": "step", "t0": 0.0, "v0": 0.1, "v1": 0.6},
    }
    assert _ramp_target_voltage(contact) == pytest.approx(0.6)


def test_ramp_target_voltage_no_voltage_t():
    contact = {
        "name": "anode", "facet": "f", "type": "ohmic", "voltage": 0.4,
    }
    assert _ramp_target_voltage(contact) == pytest.approx(0.4)


# ---------------------------------------------------------------------
# End-to-end runs against `run_transient`. These exercise the time-loop
# BC build with a real evaluator and the IV recorder. Tiny mesh, short
# t_end, one or two timesteps.
# ---------------------------------------------------------------------


_L = 2.0e-6
_TAU = 1.0e-9


def _make_diode_cfg():
    return {
        "schema_version": "2.7.0",
        "name": "voltage_t_diode",
        "dimension": 1,
        "mesh": {
            "source": "builtin",
            "extents": [[0.0, _L]],
            "resolution": [4],
            "regions_by_box": [
                {"name": "silicon", "tag": 1, "bounds": [[0.0, _L]]},
            ],
            "facets_by_plane": [
                {"name": "anode",   "tag": 1, "axis": 0, "value": 0.0},
                {"name": "cathode", "tag": 2, "axis": 0, "value": _L},
            ],
        },
        "regions": {"silicon": {"material": "Si", "tag": 1,
                                "role": "semiconductor"}},
        "doping": [
            {"region": "silicon",
             "profile": {"type": "step", "axis": 0, "location": _L / 2.0,
                         "N_D_left": 0.0, "N_A_left": 1.0e17,
                         "N_D_right": 1.0e17, "N_A_right": 0.0}}
        ],
        "contacts": [
            {"name": "anode",   "facet": "anode",   "type": "ohmic",
             "voltage": 0.0},
            {"name": "cathode", "facet": "cathode", "type": "ohmic",
             "voltage": 0.0},
        ],
        "physics": {"temperature": 300.0, "statistics": "boltzmann",
                    "mobility": {"mu_n": 1400.0, "mu_p": 450.0},
                    "recombination": {"srh": True, "tau_n": _TAU,
                                      "tau_p": _TAU, "E_t": 0.0}},
        "solver": {
            "type": "transient",
            "t_end": 5.0e-9,
            "dt": 1.0e-9,
            "order": 1,
            "bc_ramp_steps": 0,
            "max_steps": 200,
            "output_every": 1000,
        },
    }


def _iv_v_at(iv_rows, contact, t):
    for r in iv_rows:
        if r["contact"] == contact and abs(r["t"] - t) < 1.0e-15:
            return r["V"]
    raise AssertionError(
        f"no IV row for contact={contact!r} t={t} in {iv_rows!r}"
    )


@require_dolfinx
def test_run_transient_step_records_per_step_voltage():
    """A two-step run: t0 between the two timesteps. The first IV row
    after t=0 has V=v0; the second has V=v1."""
    from semi import schema
    from semi.runners.transient import run_transient

    cfg = _make_diode_cfg()
    cfg["contacts"][0]["voltage_t"] = {
        "type": "step", "t0": 1.5e-9, "v0": 0.0, "v1": 0.4,
    }
    cfg["solver"]["t_end"] = 3.0e-9
    cfg["solver"]["dt"] = 1.0e-9
    cfg = schema.validate(copy.deepcopy(cfg))

    result = run_transient(cfg)

    assert _iv_v_at(result.iv, "anode", 0.0) == pytest.approx(0.0)
    assert _iv_v_at(result.iv, "anode", 1.0e-9) == pytest.approx(0.0)
    # t=2e-9 is past t0=1.5e-9; v1 is in effect.
    assert _iv_v_at(result.iv, "anode", 2.0e-9) == pytest.approx(0.4)
    assert _iv_v_at(result.iv, "anode", 3.0e-9) == pytest.approx(0.4)


@require_dolfinx
def test_run_transient_table_records_interpolated_voltage():
    """A linear ramp: each per-step IV row has V = interp(t)."""
    from semi import schema
    from semi.runners.transient import run_transient

    times = [0.0, 1.0e-9, 2.0e-9, 3.0e-9, 4.0e-9, 5.0e-9]
    values = [0.0, 0.05, 0.1, 0.15, 0.2, 0.25]
    cfg = _make_diode_cfg()
    cfg["contacts"][0]["voltage_t"] = {
        "type": "table", "times": times, "values": values,
    }
    cfg["solver"]["dt"] = 1.0e-9
    cfg["solver"]["t_end"] = 5.0e-9
    cfg = schema.validate(copy.deepcopy(cfg))

    result = run_transient(cfg)
    for t, v in zip(times, values, strict=True):
        assert _iv_v_at(result.iv, "anode", t) == pytest.approx(v, rel=1.0e-12,
                                                                abs=1.0e-15)
