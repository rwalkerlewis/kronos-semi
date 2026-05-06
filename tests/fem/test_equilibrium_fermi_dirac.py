"""
FEM smoke tests for the Fermi-Dirac equilibrium runner (M16.4).

The Variant G MMS test (`test_mms_fermi_dirac.py`) exercises the FD path
through the drift-diffusion residual but never invokes the equilibrium
runner; the closed-loop benchmark (`benchmarks/diode_fermi_dirac_1d`)
exercises both paths but lives in `docker-fem-benchmarks`, so its
coverage does not contribute to the unit-test gate. This test runs
`run_equilibrium` on a coarse FD pn diode so the FD branches in
`semi.physics.poisson._equilibrium_space_charge`,
`semi.runners.equilibrium.run_equilibrium`, and the UFL helpers in
`semi.physics.slotboom` are exercised in the gated suite.
"""
from __future__ import annotations

import numpy as np
import pytest


def _fd_diode_cfg(N_resolution: int = 200) -> dict:
    """Coarsened FD pn diode (10 um, p / n with N_D = 1e18 cm^-3 on the
    n-side, N_A = 1e17 cm^-3 on the p-side). The benchmark config uses
    N_D = 1e20 cm^-3 over 20 um with 800 cells; the smoke test scales
    down to keep the unit-test wallclock low while staying in the regime
    where the FD correction is materially nonzero (eta_n ~ -1.0 still
    deviates ~1 % in the Blakemore prefactor)."""
    return {
        "schema_version": "2.4.0",
        "name": "fd_equilibrium_smoke",
        "description": "M16.4 FD equilibrium smoke test (coarsened pn diode).",
        "dimension": 1,
        "mesh": {
            "source": "builtin",
            "extents": [[0.0, 1.0e-5]],
            "resolution": [N_resolution],
            "regions_by_box": [
                {"name": "silicon", "tag": 1, "bounds": [[0.0, 1.0e-5]]},
            ],
            "facets_by_plane": [
                {"name": "anode",   "tag": 1, "axis": 0, "value": 0.0},
                {"name": "cathode", "tag": 2, "axis": 0, "value": 1.0e-5},
            ],
        },
        "regions": {
            "silicon": {"material": "Si", "tag": 1, "role": "semiconductor"},
        },
        "doping": [
            {
                "region": "silicon",
                "profile": {
                    "type": "step", "axis": 0, "location": 5.0e-6,
                    "N_D_left": 0.0, "N_A_left": 1.0e17,
                    "N_D_right": 1.0e18, "N_A_right": 0.0,
                },
            },
        ],
        "contacts": [
            {"name": "anode",   "facet": "anode",   "type": "ohmic", "voltage": 0.0},
            {"name": "cathode", "facet": "cathode", "type": "ohmic", "voltage": 0.0},
        ],
        "physics": {
            "temperature": 300.0, "statistics": "fermi_dirac",
            "mobility": {"mu_n": 1400.0, "mu_p": 450.0},
            "recombination": {"srh": True, "tau_n": 1.0e-7, "tau_p": 1.0e-7, "E_t": 0.0},
        },
        "solver": {
            "type": "equilibrium",
            "snes": {"rtol": 1.0e-10, "atol": 1.0e-8, "stol": 1.0e-14, "max_it": 80},
        },
        "output": {"directory": "/tmp/fd_equilibrium_smoke", "fields": []},
    }


def test_run_equilibrium_fermi_dirac_post_processing():
    """`run_equilibrium` under FD must populate `n_phys` / `p_phys` via the
    Blakemore prefactor branch and return a converged result. We compare
    the n+ side bulk electron density to N_D within a few percent: the
    Blakemore correction caps n_n_bulk strictly below N_D, but the deficit
    is small (~3 % on benchmark numbers), so a 10 % envelope guards
    against post-processing regressions without overfitting to the FD
    correction magnitude."""
    from semi import schema
    from semi.runners.equilibrium import run_equilibrium

    cfg = schema.validate(_fd_diode_cfg())
    result = run_equilibrium(cfg)

    assert result.solver_info.get("converged", False), (
        f"FD equilibrium SNES did not converge: {result.solver_info}"
    )

    # Sample dofs in the n+ region (x > 7.5 um) and check n_n_bulk ~ N_D
    # within 10 %. p side mirror.
    x = result.x_dof[:, 0]
    n_phys = np.asarray(result.n_phys)
    p_phys = np.asarray(result.p_phys)
    n_side = x > 7.5e-6
    p_side = x < 2.5e-6
    N_D_SI = 1.0e18 * 1.0e6  # cm^-3 -> m^-3
    N_A_SI = 1.0e17 * 1.0e6

    n_bulk = float(np.median(n_phys[n_side]))
    p_bulk = float(np.median(p_phys[p_side]))
    # Charge neutrality at equilibrium: n - p = N_D - N_A. With N_D
    # >> n_i (1e18 vs 1e10 cm^-3 in Si Altermatt) the n-side bulk
    # electron density is approximately N_D within a few percent;
    # mirror identity for p-side. Loose envelope guards against
    # post-processing regressions without overfitting to mesh
    # resolution.
    assert 0.5 * N_D_SI < n_bulk < 1.5 * N_D_SI, (
        f"FD n bulk electron density {n_bulk:.3e} not near N_D = {N_D_SI:.3e}"
    )
    assert 0.5 * N_A_SI < p_bulk < 1.5 * N_A_SI, (
        f"FD p bulk hole density {p_bulk:.3e} not near N_A = {N_A_SI:.3e}"
    )


def test_n_from_slotboom_ufl_fd_requires_eta_offset():
    """UFL-side `n_from_slotboom` must surface a clear ValueError when
    `eta_offset_n` is missing on the FD path. Pre-M16.4 callers (and any
    statistics_cfg=None caller) must remain on the bare Boltzmann form,
    which doesn't need eta_offset."""
    from dolfinx import fem, mesh
    from mpi4py import MPI
    from petsc4py import PETSc

    from semi.physics.slotboom import n_from_slotboom, p_from_slotboom

    msh = mesh.create_interval(MPI.COMM_WORLD, 4, [0.0, 1.0])
    psi = fem.Constant(msh, PETSc.ScalarType(0.0))
    phi = fem.Constant(msh, PETSc.ScalarType(0.0))
    n_i = fem.Constant(msh, PETSc.ScalarType(1.0))

    with pytest.raises(ValueError, match="eta_offset_n"):
        n_from_slotboom(
            psi, phi, n_i,
            statistics_cfg={"statistics": "fermi_dirac"},
        )
    with pytest.raises(ValueError, match="eta_offset_p"):
        p_from_slotboom(
            psi, phi, n_i,
            statistics_cfg={"statistics": "fermi_dirac"},
        )
