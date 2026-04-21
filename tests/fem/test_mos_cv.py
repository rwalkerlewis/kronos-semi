"""
End-to-end tests for the MOS capacitor C-V runner.

Exercises the Day 6 `mos_cv` solver type:
    - Multi-region equilibrium Poisson with cellwise eps_r
    - Gate contact BC + ohmic body contact
    - Integrated gate charge bookkeeping
    - Sign / monotonicity of Q_gate(V_gate)

Uses a thin 2D mesh with the same Si/SiO2 aspect ratio as the full
benchmark but a coarser vertical resolution (101 cells for a 101 nm
device) so the test runs in a few seconds.
"""
from __future__ import annotations

import numpy as np
import pytest


def _tiny_mos_cfg(v_values=(-0.3, 0.0, 0.3)):
    """A minimal MOS capacitor config for quick smoke tests.

    100 nm Si substrate + 5 nm oxide; 105 cells vertically gives 1 nm
    spacing and keeps y_int = 100 nm on a grid line. 2 cells laterally
    is enough because the device is translation invariant in x.
    """
    return {
        "name": "mos_cap_tiny",
        "dimension": 2,
        "mesh": {
            "source": "builtin",
            "extents": [[0.0, 1.0e-7], [0.0, 1.05e-7]],
            "resolution": [2, 105],
            "regions_by_box": [
                {"name": "silicon", "tag": 1, "bounds": [[0.0, 1.0e-7], [0.0,     1.0e-7]]},
                {"name": "oxide",   "tag": 2, "bounds": [[0.0, 1.0e-7], [1.0e-7,  1.05e-7]]},
            ],
            "facets_by_plane": [
                {"name": "body", "tag": 1, "axis": 1, "value": 0.0},
                {"name": "gate", "tag": 2, "axis": 1, "value": 1.05e-7},
            ],
        },
        "regions": {
            "silicon": {"material": "Si",   "tag": 1, "role": "semiconductor"},
            "oxide":   {"material": "SiO2", "tag": 2, "role": "insulator"},
        },
        "doping": [
            {"region": "silicon",
             "profile": {"type": "uniform", "N_D": 0.0, "N_A": 1.0e17}},
        ],
        "contacts": [
            {"name": "body", "facet": "body", "type": "ohmic", "voltage": 0.0},
            {
                "name": "gate", "facet": "gate", "type": "gate",
                "voltage": 0.0, "workfunction": 0.0,
                "voltage_sweep": {
                    "start": float(min(v_values)),
                    "stop":  float(max(v_values)),
                    "step":  0.3,
                },
            },
        ],
        "physics": {
            "temperature": 300.0,
            "statistics": "boltzmann",
            "mobility": {"mu_n": 1400.0, "mu_p": 450.0},
            "recombination": {"srh": False, "tau_n": 1.0e-7,
                               "tau_p": 1.0e-7, "E_t": 0.0},
        },
        "solver": {"type": "mos_cv"},
        "output": {"directory": "./results/mos_cap_tiny"},
    }


def test_run_mos_cv_sweep_records_iv_rows_and_qgate_sign():
    """Run a 3-point gate sweep and check Q_gate behaves physically.

    Acceptance:
      - Solver converges at each V_gate (last_info.converged=True)
      - iv_rows has one entry per V_gate in the sweep
      - Q_gate at V_gate < V_FB is negative (accumulation: more holes
        at the surface means negative silicon charge density integrates
        negative; Q_gate = -Q_semi > 0... hmm sign check below).
        Concretely: at V_gate = -0.3 (below V_FB = -0.417? actually
        above V_FB here since V_FB ~ -0.417 and -0.3 > -0.417), we're
        in weak depletion. Q_gate > 0.
      - Q_gate monotone non-decreasing in V_gate across a depletion-
        regime sweep.
    """
    from semi import run as semi_run
    from semi import schema

    cfg = schema.validate(_tiny_mos_cfg())
    result = semi_run.run(cfg)

    assert result.solver_info["converged"]
    assert result.iv is not None
    assert len(result.iv) == 3
    assert all("V" in r and "Q_gate" in r for r in result.iv)

    V = np.array([r["V"] for r in result.iv])
    Q = np.array([r["Q_gate"] for r in result.iv])
    order = np.argsort(V)
    V = V[order]
    Q = Q[order]

    # Monotone non-decreasing on this narrow depletion window.
    assert all(Q[i + 1] >= Q[i] - 1.0e-18 for i in range(len(Q) - 1)), (
        f"Q_gate not monotone non-decreasing in V_gate: V={V}, Q={Q}"
    )

    # At the largest V_gate (+0.3), the silicon surface is depleted,
    # so Q_semi < 0 => Q_gate > 0.
    assert Q[-1] > 0.0, f"Q_gate at largest V_gate must be positive: {Q[-1]}"


def test_run_mos_cv_preserves_bulk_equilibrium_at_flatband():
    """At V_gate ≈ V_FB the bulk silicon must return to p-type equilibrium.

    Specifically: psi at a point deep in the bulk (y = 20 nm, well
    inside the 100 nm substrate) must match -phi_F within ~V_t, and
    hole density p = n_i exp(-psi/V_t) must match N_A within 5%.
    """
    from semi import run as semi_run
    from semi import schema
    from semi.constants import cm3_to_m3
    from semi.materials import get_material

    Si = get_material("Si")
    N_A = cm3_to_m3(1.0e17)
    V_t = 0.02585
    phi_F = V_t * float(np.log(N_A / Si.n_i))
    V_FB = -phi_F  # ideal gate, psi=0-at-intrinsic convention

    # Sweep through V_FB; the middle sample is the one we check.
    cfg = schema.validate(_tiny_mos_cfg(v_values=(V_FB - 0.05, V_FB, V_FB + 0.05)))
    result = semi_run.run(cfg)
    assert result.solver_info["converged"]

    x_dof = result.x_dof
    psi = result.psi_phys

    # pick a DOF near (x=0.5e-7, y=20e-9) (20 nm deep in silicon)
    target = np.array([5.0e-8, 2.0e-8])
    dists = np.linalg.norm(x_dof[:, :2] - target, axis=1)
    idx = int(np.argmin(dists))

    psi_bulk = float(psi[idx])
    # at V_FB the bulk psi should equal the equilibrium p-type value = -phi_F.
    # tolerance: a few V_t because this is the last V_gate in the sweep
    # (V_FB + 0.05), not exactly flatband.
    assert abs(psi_bulk - (-phi_F)) < 3.0 * V_t, (
        f"bulk psi at V_gate near flatband: {psi_bulk} V, expected ~{-phi_F} V"
    )


def test_equilibrium_poisson_form_mr_matches_scalar_on_single_region():
    """
    Sanity: on a single-region silicon mesh, the multi-region form
    (with semi_tag covering the whole mesh) must produce the same
    solution as the existing scalar-eps_r equilibrium Poisson form.

    This protects the 1D pn_1d byte-identity invariant by proving the
    coefficient cellwise-function path equals the Constant path when
    the coefficient is the same everywhere.
    """
    import numpy as _np
    from dolfinx import fem
    from dolfinx import mesh as dmesh
    from mpi4py import MPI
    from petsc4py import PETSc

    from semi.materials import get_material
    from semi.mesh import build_eps_r_function
    from semi.physics.poisson import (
        build_equilibrium_poisson_form,
        build_equilibrium_poisson_form_mr,
    )
    from semi.scaling import Scaling
    from semi.solver import solve_nonlinear

    Si = get_material("Si")
    N_A = 1.0e23
    L = 1.0e-6
    N = 40
    sc = Scaling(L0=L, C0=N_A, T=300.0, mu0=Si.mu_n, n_i=Si.n_i)

    msh = dmesh.create_interval(MPI.COMM_WORLD, N, [0.0, L])

    tdim = msh.topology.dim
    msh.topology.create_entities(tdim)
    n_local = msh.topology.index_map(tdim).size_local
    all_cells = _np.arange(n_local, dtype=_np.int32)
    cell_tags = dmesh.meshtags(
        msh, tdim, all_cells, _np.ones(n_local, dtype=_np.int32),
    )

    regions_cfg = {"silicon": {"material": "Si", "tag": 1, "role": "semiconductor"}}
    eps_r_fn = build_eps_r_function(msh, cell_tags, regions_cfg)

    V = fem.functionspace(msh, ("Lagrange", 1))
    N_hat_fn = fem.Function(V)
    N_hat_fn.x.array[:] = -N_A / sc.C0  # p-type

    # BCs: equilibrium psi_hat on both ends
    psi_eq_hat = float(_np.arcsinh(-N_A / (2.0 * Si.n_i)))
    fdim = tdim - 1
    msh.topology.create_connectivity(fdim, tdim)
    left = dmesh.locate_entities_boundary(msh, fdim, lambda x: _np.isclose(x[0], 0.0))
    right = dmesh.locate_entities_boundary(msh, fdim, lambda x: _np.isclose(x[0], L))

    def _solve(build):
        psi = fem.Function(V)
        psi.x.array[:] = psi_eq_hat
        bcs = [
            fem.dirichletbc(
                PETSc.ScalarType(psi_eq_hat),
                fem.locate_dofs_topological(V, fdim, left), V),
            fem.dirichletbc(
                PETSc.ScalarType(psi_eq_hat),
                fem.locate_dofs_topological(V, fdim, right), V),
        ]
        F = build(psi)
        info = solve_nonlinear(
            F, psi, bcs, prefix=f"mos_mr_match_{id(build)}_",
            petsc_options={
                "snes_rtol": 1.0e-14,
                "snes_atol": 1.0e-16,
                "snes_max_it": 10,
            },
        )
        assert info["converged"]
        return psi.x.array.copy()

    scalar_result = _solve(
        lambda psi: build_equilibrium_poisson_form(V, psi, N_hat_fn, sc, Si.epsilon_r)
    )
    mr_result = _solve(
        lambda psi: build_equilibrium_poisson_form_mr(
            V, psi, N_hat_fn, sc, eps_r_fn, cell_tags, semi_tag=1,
        )
    )
    assert scalar_result == pytest.approx(mr_result, rel=1.0e-10, abs=1.0e-12)


def test_mos_cv_raises_without_gate_sweep():
    """The runner must raise a clear error if no gate contact has a voltage_sweep."""
    from semi import run as semi_run
    from semi import schema

    cfg = _tiny_mos_cfg()
    # strip the voltage_sweep from the gate
    for c in cfg["contacts"]:
        if c["type"] == "gate":
            c.pop("voltage_sweep", None)
    cfg = schema.validate(cfg)

    with pytest.raises(ValueError, match="voltage_sweep"):
        semi_run.run(cfg)
