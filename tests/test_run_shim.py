"""
Pure-Python tests for the `semi.run` dispatcher and its
backward-compat re-exports.

Confirms that:

- `run(cfg)` raises on unknown solver.type.
- The `__getattr__` shim resolves the legacy attribute names
  (`run_equilibrium`, `run_bias_sweep`, `_fmt_tag`, `_resolve_sweep`)
  and rejects unknown attributes.
- `_resolve_sweep` (re-exported from `runners.bias_sweep`) returns
  no-sweep for a config without `voltage_sweep`.
- `runners._common.reference_material` raises when no semiconductor
  region is present.
"""
from __future__ import annotations

import pytest


def test_run_unknown_solver_type_raises():
    import semi.run as semi_run

    with pytest.raises(ValueError, match="Unknown solver.type"):
        semi_run.run({"solver": {"type": "magic"}})


def test_run_unknown_attribute_raises():
    import semi.run as semi_run

    with pytest.raises(AttributeError, match="no attribute"):
        semi_run.no_such_attribute  # noqa: B018


def test_run_shim_exposes_legacy_names():
    from semi.run import _fmt_tag, _resolve_sweep, run_bias_sweep, run_equilibrium

    # All four are callables; exact identity isn't important, only resolution.
    for fn in (_fmt_tag, _resolve_sweep, run_equilibrium, run_bias_sweep):
        assert callable(fn)
    # _fmt_tag is the postprocess.fmt_tag helper.
    assert _fmt_tag(0.0) == "p0d0000"


def test_resolve_sweep_no_sweep_returns_none_and_empty():
    from semi.run import _resolve_sweep

    cfg = {
        "contacts": [
            {"name": "anode",   "type": "ohmic", "voltage": 0.0},
            {"name": "cathode", "type": "ohmic", "voltage": 0.5},
        ],
    }
    contact, values = _resolve_sweep(cfg)
    assert contact is None
    assert values == []


def test_runners_common_reference_material_no_semiconductor_raises():
    from semi.runners._common import reference_material

    cfg = {"regions": {"oxide": {"material": "SiO2"}}}
    with pytest.raises(ValueError, match="No semiconductor region"):
        reference_material(cfg)
