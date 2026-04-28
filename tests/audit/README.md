# Physics validation audit suite (Phase 1)

This directory holds the cross-runner consistency audit (Phase 1 of
the physics validation suite, see `PLAN.md`). The tests are scaffolded
under `pytest.mark.audit` and **deselected from the default test run**
via `addopts = "-m 'not audit'"` in `pyproject.toml`.

## Running

    docker compose run --rm test pytest tests/audit/ -v -m audit

Each case writes a CSV to `/tmp/audit/<case>.csv` and a markdown
fragment to `/tmp/audit/<case>.md`. The session finalizer in
`conftest.py` aggregates the fragments into `docs/PHYSICS_AUDIT.md`
when `-m audit` is supplied.

## Cases

| # | File | What it compares |
|---|------|------------------|
| 01 | `test_01_bias_vs_transient_steady_state.py` | `run_bias_sweep` vs `run_transient` (deep SS) at multiple V_F |
| 02 | `test_02_ac_omega0_vs_bias_dIdV.py` | `run_ac_sweep` at low omega vs `run_bias_sweep` dI/dV (reverse bias) |
| 03 | `test_03_mos_cv_vs_mos_cap_ac.py` | `run_mos_cv` vs `run_mos_cap_ac` on Q_gate (M14.1 byte-identity claim) |
| 04 | `test_04_equilibrium_vs_bias_sweep_v0.py` | `run_equilibrium` vs `run_bias_sweep` halted at V=0 |
| 05 | `test_05_ac_current_vs_bias_dIdV_forward.py` | AC Re(Y) at low omega vs `run_bias_sweep` dI/dV (forward bias) |
| 06 | `test_06_transient_fft_vs_ac_sweep.py` | Transient FFT vs `run_ac_sweep` Y(omega) at the same frequency |

Case 06 is currently a tracking placeholder: `run_transient` does not
yet accept a time-varying contact voltage, which the FFT comparison
requires. The test is `pytest.skip` with a note recorded in
`docs/PHYSICS_AUDIT.md` so the gap is visible. Closing the gap is
deferred to a follow-up PR.

## Tolerance philosophy

The audit is a discovery exercise, not a regression gate. Tolerances
are intentionally generous on the first pass; tighten in follow-up
PRs once we know what each case reports as a "clean" baseline.

If a case fails on first run, the finding goes into
`docs/PHYSICS_AUDIT.md`, a tracking issue is opened, and a fix is
landed in this PR if it is under 50 lines, or deferred otherwise.
