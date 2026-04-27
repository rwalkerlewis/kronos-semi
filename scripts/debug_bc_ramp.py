"""
Diagnostic: verify BC-ramp continuation IC vs bias_sweep, then run the full
transient time loop and report J at t=0, t=mid, and t=end. Reproduces the
M13.1 close-out finding (audit at /tmp/m13.1-integration-blocker.md):

    Step 2 (IC only):      J_tr@t=0  = 2.884213e+00
                           J_ss      = 2.884211e+00
                           rel_err   = 5.7e-7   <-- BC-ramp IC is correct.

    Step 3 (full transient):
                           J at t=0     = 2.884213e+00
                           J at t=mid   = -6.273e+05
                           J at t=end   = -6.273e+05
                           J_ss         = 2.884e+00

    Conclusion: the (n,p) Galerkin time loop drifts to its OWN discrete
    fixed point within ~12 BDF steps; closing the SS-limit test xfail
    needs SG in the time loop, not just the BC-ramp IC.

Kept in scripts/ for the M13.1 follow-up #4 agent. Not run in CI.

Usage:
    docker compose run --rm test python scripts/debug_bc_ramp.py
"""
from __future__ import annotations


def make_cfg(solver_block):
    L = 2.0e-6
    TAU = 1.0e-9
    V_F = 0.3
    return {
        "schema_version": "1.1.0",
        "name": "test_transient_ss_limit",
        "description": "Steady-state limit test for transient solver.",
        "dimension": 1,
        "mesh": {
            "source": "builtin",
            "extents": [[0.0, L]],
            "resolution": [100],
            "regions_by_box": [
                {"name": "silicon", "tag": 1, "bounds": [[0.0, L]]},
            ],
            "facets_by_plane": [
                {"name": "anode",   "tag": 1, "axis": 0, "value": 0.0},
                {"name": "cathode", "tag": 2, "axis": 0, "value": L},
            ],
        },
        "regions": {"silicon": {"material": "Si", "tag": 1, "role": "semiconductor"}},
        "doping": [
            {
                "region": "silicon",
                "profile": {
                    "type": "step",
                    "axis": 0,
                    "location": L / 2.0,
                    "N_D_left": 0.0,
                    "N_A_left": 1.0e17,
                    "N_D_right": 1.0e17,
                    "N_A_right": 0.0,
                },
            }
        ],
        "contacts": [
            {"name": "anode",   "facet": "anode",   "type": "ohmic", "voltage": V_F},
            {"name": "cathode", "facet": "cathode",  "type": "ohmic", "voltage": 0.0},
        ],
        "physics": {
            "temperature": 300.0,
            "statistics": "boltzmann",
            "mobility": {"mu_n": 1400.0, "mu_p": 450.0},
            "recombination": {
                "srh": True,
                "tau_n": TAU,
                "tau_p": TAU,
                "E_t": 0.0,
            },
        },
        "solver": solver_block,
        "output": {"directory": "/tmp/debug_bc_ramp", "fields": []},
    }


def main():
    from semi.runners.bias_sweep import run_bias_sweep
    from semi.runners.transient import run_transient

    V_F = 0.3

    print("=" * 60)
    print("Step 1: bias_sweep ground truth at V_F=0.3")
    print("=" * 60)
    ss_cfg = make_cfg({
        "type": "bias_sweep",
        "snes": {"rtol": 1.0e-14, "atol": 1.0e-14, "stol": 1.0e-14, "max_it": 60},
        "continuation": {
            "min_step": 1.0e-3, "max_halvings": 6,
            "max_step": 0.1, "easy_iter_threshold": 4, "grow_factor": 2.0,
        },
    })
    ss_cfg["contacts"][0]["voltage_sweep"] = {"start": 0.0, "stop": V_F, "step": 0.05}
    ss = run_bias_sweep(ss_cfg)
    rows = [r for r in ss.iv if abs(r.get("V", 0) - V_F) < 1e-4]
    if not rows:
        print(f"available V values: {[r.get('V') for r in ss.iv]}")
        raise SystemExit("no row at V_F")
    J_ss = float(rows[-1]["J"])
    print(f"J_ss (bias_sweep, anode, V=0.3) = {J_ss:.6e} A/m^2")
    print(f"psi[0] = {ss.psi.x.array[0]:.6e}")
    print(f"psi[-1] = {ss.psi.x.array[-1]:.6e}")
    print(f"phi_n[0] = {ss.phi_n.x.array[0]:.6e}")
    print(f"phi_n[-1] = {ss.phi_n.x.array[-1]:.6e}")

    print()
    print("=" * 60)
    print("Step 2: transient with bc_ramp_steps=10, very short t_end")
    print("=" * 60)
    # very short t_end so only 1 step is taken (or zero)
    tr_cfg = make_cfg({
        "type": "transient",
        "t_end": 1.0e-12,  # 1 ps only
        "dt": 5.0e-11,     # bigger than t_end -> 0 steps in time loop
        "order": 1,
        "max_steps": 1,
        "output_every": 1,
        "bc_ramp_steps": 10,
        "snes": {"rtol": 1.0e-10, "atol": 1.0e-7, "stol": 1.0e-14, "max_it": 100},
    })
    tr = run_transient(tr_cfg)
    print(f"# IV rows = {len(tr.iv)}")
    for i, r in enumerate(tr.iv[:6]):
        print(f"  iv[{i}]: t={r['t']:.3e} contact={r['contact']} V={r['V']:.3f} J={r['J']:.6e}")
    anode_first = next(r for r in tr.iv if r.get("contact") == "anode")
    J_tr_t0 = float(anode_first["J"])
    print(f"J_tr at t=0 = {J_tr_t0:.6e} A/m^2")
    print(f"J_ss        = {J_ss:.6e} A/m^2")
    if abs(J_ss) > 1e-30:
        print(f"rel_err     = {abs(J_tr_t0 - J_ss) / abs(J_ss):.3e}")

    print()
    print("=" * 60)
    print("Step 3: full t_end=30ns transient run")
    print("=" * 60)
    tr_cfg2 = make_cfg({
        "type": "transient",
        "t_end": 3.0e-8,
        "dt": 5.0e-11,
        "order": 2,
        "max_steps": 700,
        "output_every": 700,
        "bc_ramp_steps": 10,
        "snes": {"rtol": 1.0e-10, "atol": 1.0e-7, "stol": 1.0e-14, "max_it": 100},
    })
    tr2 = run_transient(tr_cfg2)
    anode_rows = [r for r in tr2.iv if r.get("contact") == "anode"]
    print(f"# anode IV rows = {len(anode_rows)}")
    print(f"J at t=0:    {anode_rows[0]['J']:.6e}")
    print(f"J at t=mid:  {anode_rows[len(anode_rows)//2]['J']:.6e}")
    print(f"J at t=end:  {anode_rows[-1]['J']:.6e}")
    print(f"J_ss:        {J_ss:.6e}")


if __name__ == "__main__":
    main()
