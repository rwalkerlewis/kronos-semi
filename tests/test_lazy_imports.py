"""
Pure-Python tests for the lazy-import ``__getattr__`` in ``semi.run``.

These cover the branches that return the less-frequently-accessed runners
(run_mos_cv, run_mos_cap_ac, run_transient, run_ac_sweep) and the
AttributeError fallthrough, all without loading dolfinx.
"""
import pytest


def test_run_lazy_import_run_mos_cv():
    import semi.run as r
    fn = r.run_mos_cv
    assert callable(fn)


def test_run_lazy_import_run_mos_cap_ac():
    import semi.run as r
    fn = r.run_mos_cap_ac
    assert callable(fn)


def test_run_lazy_import_run_transient():
    import semi.run as r
    fn = r.run_transient
    assert callable(fn)


def test_run_lazy_import_run_ac_sweep():
    import semi.run as r
    fn = r.run_ac_sweep
    assert callable(fn)


def test_run_getattr_raises_attribute_error_for_unknown_name():
    import semi.run as r
    with pytest.raises(AttributeError, match="no attribute"):
        _ = r.no_such_runner_exists
