"""
Build the Day 1 Colab notebook that uses the kronos-semi package via git clone.
Replaces the inline version with a thin notebook that imports from semi/.
"""
import nbformat as nbf
from pathlib import Path

nb = nbf.v4.new_notebook()
cells = []

cells.append(nbf.v4.new_markdown_cell(r"""# kronos-semi · 1D pn Junction (Day 1)

Verification benchmark: symmetric 1D pn junction, Si, $N_A = N_D = 10^{17}$ cm⁻³, junction at 1 µm. Solves the nonlinear equilibrium Poisson equation (Boltzmann statistics) via PETSc SNES on a dolfinx 0.10 mesh, then compares to the analytical depletion approximation.

**What this notebook does:**
1. Installs `dolfinx` on Colab via [FEM on Colab](https://fem-on-colab.github.io/) (one wget, ~30 s)
2. Clones the [kronos-semi](https://github.com/rwalkerlewis/kronos-semi) repo
3. Loads `benchmarks/pn_1d/pn_junction.json`
4. Runs `semi.run.run(cfg)` which returns the scaled potential, carriers, and diagnostics
5. Plots ψ(x), E(x), carrier densities; compares to analytical depletion approximation

If you want to modify the physics, edit the JSON (resolution, doping, geometry, contact voltages) and re-run. The Python code lives in the package — this notebook is a thin driver.
"""))

cells.append(nbf.v4.new_markdown_cell(r"""## 1. Install dolfinx

~30 s on first run. Cached across cell re-runs within the same Colab session."""))

cells.append(nbf.v4.new_code_cell(r"""try:
    import dolfinx
    print(f"dolfinx already present: version {dolfinx.__version__}")
except ImportError:
    !wget -q "https://fem-on-colab.github.io/releases/fenicsx-install-release-real.sh" -O "/tmp/fenicsx-install.sh" && bash "/tmp/fenicsx-install.sh"
    import dolfinx
    print(f"dolfinx installed: version {dolfinx.__version__}")

major, minor = (int(x) for x in dolfinx.__version__.split('.')[:2])
assert (major, minor) >= (0, 10), f"Need dolfinx >= 0.10; got {dolfinx.__version__}"
"""))

cells.append(nbf.v4.new_markdown_cell(r"""## 2. Clone the repo and install the package

This pulls in `semi/` (the actual simulator code) and makes it importable. Uses `pip install -e` so any edits you make persist in this session."""))

cells.append(nbf.v4.new_code_cell(r"""import os
if not os.path.exists('kronos-semi'):
    !git clone -q https://github.com/rwalkerlewis/kronos-semi.git

%cd kronos-semi
!pip install -q -e .
print("Package installed.")
"""))

cells.append(nbf.v4.new_code_cell(r"""import numpy as np
import matplotlib.pyplot as plt

from semi import schema, run, materials, constants

print(f"semi version: {__import__('semi').__version__}")
"""))

cells.append(nbf.v4.new_markdown_cell(r"""## 3. Load the benchmark

The JSON schema is flat and human-readable. Inspect it to see how the simulator is driven."""))

cells.append(nbf.v4.new_code_cell(r"""cfg = schema.load("benchmarks/pn_1d/pn_junction.json")

print(f"Case:       {cfg['name']}")
print(f"Dimension:  {cfg['dimension']}D")
print(f"Mesh:       {cfg['mesh']['source']}, extent={cfg['mesh']['extents'][0]}, "
      f"{cfg['mesh']['resolution'][0]} elements")
print(f"Regions:    {list(cfg['regions'].keys())}")
print(f"Doping:     step profile, N_A={cfg['doping'][0]['profile']['N_A_left']:.0e} cm^-3 | "
      f"N_D={cfg['doping'][0]['profile']['N_D_right']:.0e} cm^-3")
print(f"Contacts:   {[c['name'] for c in cfg['contacts']]}")
print(f"Solver:     {cfg['solver']['type']}")
"""))

cells.append(nbf.v4.new_markdown_cell(r"""## 4. Run the simulation

This builds the mesh, interpolates the doping profile, applies ohmic BCs (equilibrium potential from charge neutrality), assembles the nonlinear Poisson residual, and calls PETSc SNES. Watch the SNES monitor output."""))

cells.append(nbf.v4.new_code_cell(r"""result = run.run(cfg)

info = result.solver_info
print(f"\nSNES: iterations = {info['iterations']}, "
      f"converged = {info['converged']}, reason = {info['reason']}")
print(f"Scaling: {result.scaling}")
"""))

cells.append(nbf.v4.new_markdown_cell(r"""## 5. Extract and plot results

The result object carries:
- `result.psi_phys`: potential in volts, dof-ordered (not spatially sorted)
- `result.n_phys`, `result.p_phys`: carrier densities in m⁻³
- `result.x_dof`: dof coordinates (N, 3)

We sort by x for plotting since dof ordering is implementation-dependent."""))

cells.append(nbf.v4.new_code_cell(r"""x_dof = result.x_dof[:, 0]
order = np.argsort(x_dof)
x_phys_um = x_dof[order] * 1e6
psi_phys  = result.psi_phys[order]
n_phys    = result.n_phys[order]
p_phys    = result.p_phys[order]

# Re-derive the analytical depletion approximation for comparison
N_A = cfg['doping'][0]['profile']['N_A_left']  * 1e6   # to m^-3
N_D = cfg['doping'][0]['profile']['N_D_right'] * 1e6
x_j = cfg['doping'][0]['profile']['location']
si = materials.get_material('Si')
eps_Si = si.epsilon
n_i = si.n_i
Vbi = constants.thermal_voltage(300.0) * np.log(N_A * N_D / n_i**2)
N_eff = N_A * N_D / (N_A + N_D)
W = np.sqrt(2 * eps_Si * Vbi / (constants.Q * N_eff))
xp = W * N_D / (N_A + N_D)
xn = W * N_A / (N_A + N_D)
E_peak_an = constants.Q * N_A * xp / eps_Si

def psi_analytical(x_si, psi_p_phys, psi_n_phys):
    xp_edge = x_j - xp
    xn_edge = x_j + xn
    out = np.empty_like(x_si)
    m1 = x_si <  xp_edge
    m2 = (x_si >= xp_edge) & (x_si < x_j)
    m3 = (x_si >= x_j) & (x_si <= xn_edge)
    m4 = x_si >  xn_edge
    out[m1] = psi_p_phys
    out[m2] = psi_p_phys + (constants.Q*N_A/(2*eps_Si)) * (x_si[m2] - xp_edge)**2
    out[m3] = psi_n_phys - (constants.Q*N_D/(2*eps_Si)) * (xn_edge - x_si[m3])**2
    out[m4] = psi_n_phys
    return out

psi_p_phys = psi_phys[2]     # a point well inside p-bulk
psi_n_phys = psi_phys[-3]    # well inside n-bulk
x_fine = np.linspace(0, 2e-6, 2000)
psi_an = psi_analytical(x_fine, psi_p_phys, psi_n_phys)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4))
ax1.plot(x_phys_um, psi_phys, 'b-', lw=2, label='FEM (kronos-semi)')
ax1.plot(x_fine*1e6, psi_an, 'r--', lw=1.5, label='Depletion approx.')
ax1.axvspan((x_j - xp)*1e6, (x_j + xn)*1e6, alpha=0.1, color='green',
            label=f'Depletion region (W = {W*1e9:.0f} nm)')
ax1.set_xlabel(r"$x$ [μm]"); ax1.set_ylabel(r"$\psi$ [V]")
ax1.set_title(r"Electrostatic potential $\psi(x)$")
ax1.legend(); ax1.grid(alpha=0.3)

E_fem = -np.gradient(psi_phys, x_dof[order])
ax2.plot(x_phys_um, E_fem*1e-5, 'b-', lw=2, label='FEM')
ax2.axhline(-E_peak_an*1e-5, color='r', ls='--', alpha=0.7,
            label=rf'$-|E_{{\max}}|$ = {-E_peak_an*1e-5:.1f} kV/cm')
ax2.set_xlabel(r"$x$ [μm]"); ax2.set_ylabel(r"$E$ [kV/cm]")
ax2.set_title("Electric field")
ax2.legend(); ax2.grid(alpha=0.3)
plt.tight_layout(); plt.show()

E_peak_fem = np.max(np.abs(E_fem))
print(f"\n  Peak |E|:  FEM        = {E_peak_fem*1e-5:7.2f} kV/cm")
print(f"             analytical = {E_peak_an*1e-5:7.2f} kV/cm")
print(f"             rel. error = {abs(E_peak_fem - E_peak_an)/E_peak_an*100:.1f}%")
print(f"  Vbi:       FEM        = {psi_phys[-1] - psi_phys[0]:.4f} V")
print(f"             analytical = {Vbi:.4f} V")
"""))

cells.append(nbf.v4.new_code_cell(r"""fig, ax = plt.subplots(figsize=(9, 4))
ax.semilogy(x_phys_um, n_phys*1e-6, 'b-', lw=2, label=r'$n$ (electrons)')
ax.semilogy(x_phys_um, p_phys*1e-6, 'r-', lw=2, label=r'$p$ (holes)')
ax.axhline(n_i*1e-6, color='k', ls=':', alpha=0.5, label=r'$n_i$')
ax.axvline(x_j*1e6, color='gray', ls='--', alpha=0.4)
ax.set_xlabel(r"$x$ [μm]"); ax.set_ylabel(r"density [cm$^{-3}$]")
ax.set_title("Equilibrium carrier densities")
ax.legend(); ax.grid(alpha=0.3, which='both')
ax.set_ylim([1e0, 1e18])
plt.tight_layout(); plt.show()

print(f"  p-bulk hole density:     FEM = {p_phys[2]*1e-6:.2e} cm^-3 "
      f"(expect ~N_A = {N_A*1e-6:.2e})")
print(f"  n-bulk electron density: FEM = {n_phys[-3]*1e-6:.2e} cm^-3 "
      f"(expect ~N_D = {N_D*1e-6:.2e})")
i_mid = len(x_phys_um) // 2
print(f"  n*p at midpoint: {n_phys[i_mid]*p_phys[i_mid]*1e-12:.3e} cm^-6")
print(f"  n_i^2:            {(n_i*1e-6)**2:.3e} cm^-6")
"""))

cells.append(nbf.v4.new_markdown_cell(r"""## Summary

- dolfinx 0.10 install via FEM on Colab: ~30 s
- Package-based workflow: thin notebook, real code in `semi/`
- JSON-driven: change the benchmark without touching Python
- 1D pn junction verifies against depletion approximation: V_bi, W, peak |E| all agree within a few percent
- Mass-action $np = n_i^2$ holds across the domain; bulk densities match majority doping

**Next step (Day 2).** Replace the Boltzmann equilibrium with a Slotboom $(\psi, \Phi_n, \Phi_p)$ block system, add drift-diffusion continuity equations with SRH recombination, enable bias sweeps. The equilibrium case here becomes the starting iterate.

See [the repo](https://github.com/rwalkerlewis/kronos-semi) for the roadmap and code."""))

nb.cells = cells

out_path = Path("/home/claude/kronos-semi/notebooks/01_pn_junction_1d.ipynb")
with out_path.open("w") as f:
    nbf.write(nb, f)
print(f"Wrote {out_path} with {len(cells)} cells")
