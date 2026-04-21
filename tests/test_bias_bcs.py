"""
Pure-Python tests for bias BC construction helpers in semi.run.

Covers _resolve_sweep (voltage_sweep parsing) and _fmt_tag (PETSc prefix
sanitization). Does not require dolfinx.
"""
from __future__ import annotations

import pytest


def test_fmt_tag_positive():
    from semi.run import _fmt_tag
    assert _fmt_tag(0.0) == "p0d0000"
    assert _fmt_tag(0.3).startswith("p")
    assert "." not in _fmt_tag(0.3)
    assert "+" not in _fmt_tag(0.3)


def test_fmt_tag_negative():
    from semi.run import _fmt_tag
    s = _fmt_tag(-0.25)
    assert s.startswith("m")
    assert "-" not in s
    assert "." not in s


def test_resolve_sweep_ascending():
    from semi.run import _resolve_sweep
    cfg = {
        "contacts": [
            {"name": "anode", "type": "ohmic", "voltage": 0.0},
            {
                "name": "cathode", "type": "ohmic", "voltage": 0.0,
                "voltage_sweep": {"start": 0.0, "stop": 0.2, "step": 0.05},
            },
        ]
    }
    contact, values = _resolve_sweep(cfg)
    assert contact == "cathode"
    assert len(values) == 5
    assert values[0] == pytest.approx(0.0)
    assert values[-1] == pytest.approx(0.2)


def test_resolve_sweep_descending():
    from semi.run import _resolve_sweep
    cfg = {
        "contacts": [
            {"name": "cathode", "type": "ohmic", "voltage": 0.0,
             "voltage_sweep": {"start": 0.0, "stop": -0.1, "step": 0.05}},
        ]
    }
    contact, values = _resolve_sweep(cfg)
    assert contact == "cathode"
    assert values[0] == pytest.approx(0.0)
    assert values[-1] == pytest.approx(-0.1)


def test_resolve_sweep_none():
    from semi.run import _resolve_sweep
    cfg = {
        "contacts": [
            {"name": "a", "type": "ohmic", "voltage": 0.0},
            {"name": "b", "type": "ohmic", "voltage": 0.5},
        ]
    }
    contact, values = _resolve_sweep(cfg)
    assert contact is None
    assert values == []


def test_resolve_sweep_rejects_zero_step():
    from semi.run import _resolve_sweep
    cfg = {
        "contacts": [
            {"name": "c", "type": "ohmic", "voltage": 0.0,
             "voltage_sweep": {"start": 0.0, "stop": 0.2, "step": 0.0}},
        ]
    }
    with pytest.raises(ValueError):
        _resolve_sweep(cfg)
