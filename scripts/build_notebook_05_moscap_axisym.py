"""Build the Colab walkthrough notebook 03_moscap_axisym_cv.

Generates ``notebooks/03_moscap_axisym_cv.ipynb`` in place. Run from the
repo root::

    python scripts/build_notebook_05_moscap_axisym.py

The notebook drives the M14.2 axisymmetric MOSCAP benchmark and plots
the LF/HF C-V curve (Hu Fig. 5-18 reproduction).
"""
from __future__ import annotations

from pathlib import Path

import nbformat as nbf

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "notebooks" / "03_moscap_axisym_cv.ipynb"

cells = []

cells.append(nbf.v4.new_markdown_cell(r"""# 03 - Axisymmetric 2D MOS capacitor C-V (LF + HF)

This notebook drives the **M14.2** axisymmetric MOSCAP benchmark and
reproduces the Hu Fig. 5-18 LF / HF C-V curves on a cylindrical
(r, z) mesh.

* Body: P-type Si, NA = 1e17 cm^-3
* Oxide: 5 nm SiO2
* Gate radius: 1 µm (axisymmetric circular MOSCAP)
* Sweep: V_g = -2 ... +2.5 V, step 0.05 V

We extract two differential capacitance curves at every bias:

* **LF** (quasi-static): both carriers respond, recovering Cox in
  strong inversion.
* **HF**: minority carriers frozen, saturating at
  Cmin = Cox * Cdep,min / (Cox + Cdep,min).

Wall time on Colab: ~5-8 minutes (a single multi-region equilibrium
solve per bias plus two linearised PDE-sensitivity solves)."""))

cells.append(nbf.v4.new_markdown_cell(r"""## 1. Install dolfinx

~30 s on first run. Cached across cell re-runs within the same Colab session."""))

cells.append(nbf.v4.new_code_cell(r"""try:
    import dolfinx
    print(f"dolfinx already installed: {dolfinx.__version__}")
except ImportError:
    !wget -q "https://fem-on-colab.github.io/releases/fenicsx-install-release-real.sh" -O "/tmp/fenicsx-install.sh" && bash "/tmp/fenicsx-install.sh"
    import dolfinx
    print(f"dolfinx installed: version {dolfinx.__version__}")"""))

cells.append(nbf.v4.new_markdown_cell(r"""## 2. Clone the repo and install the package"""))

cells.append(nbf.v4.new_code_cell(r"""import os, sys, subprocess
if not os.path.isdir("kronos-semi"):
    subprocess.run(["git", "clone", "--depth", "1",
                    "https://github.com/rwalkerlewis/kronos-semi.git"],
                    check=True)
%cd kronos-semi
!pip install -q -e . > /dev/null
print("kronos-semi installed in editable mode")"""))

cells.append(nbf.v4.new_markdown_cell(r"""## 3. Load the benchmark JSON"""))

cells.append(nbf.v4.new_code_cell(r"""import json
from pathlib import Path
from semi import schema

bench_path = Path("benchmarks/moscap_axisym/moscap_axisym.json")
cfg = schema.validate(json.loads(bench_path.read_text()))
print("solver:", cfg["solver"]["type"])
print("axisymmetric:", cfg["mesh"].get("axisymmetric"))
print("N_A:", cfg["doping"][0]["profile"]["N_A"], "cm^-3")
print("t_ox:", (cfg["mesh"]["extents"][1][1] - 1e-6) * 1e9, "nm")"""))

cells.append(nbf.v4.new_markdown_cell(r"""## 4. Run the sweep"""))

cells.append(nbf.v4.new_code_cell(r"""from semi import run as semi_run
import time
t0 = time.time()
result = semi_run.run(cfg)
print(f"sweep wall-clock: {time.time() - t0:.1f} s "
      f"({len(result.iv)} bias points)")"""))

cells.append(nbf.v4.new_markdown_cell(r"""## 5. Analytical reference values

Closed-form Cox / Cmin / V_FB / V_T from Hu Ch. 5 for a P-body with
NA = 1e17 cm^-3, t_ox = 5 nm, phi_ms = -0.977 V."""))

cells.append(nbf.v4.new_code_cell(r"""from semi.physics.cv import analytical_moscap_metrics
from semi.constants import cm3_to_m3

t = analytical_moscap_metrics(
    N_body=cm3_to_m3(1e17), body_type="p", t_ox=5e-9,
    eps_r_si=11.7, eps_r_ox=3.9,
    n_i=cm3_to_m3(1e10), T=300.0, phi_ms=-0.977,
)
for k, v in t.items():
    print(f"  {k:8s} = {v: .4e}")"""))

cells.append(nbf.v4.new_markdown_cell(r"""## 6. LF / HF C-V plot (Hu Fig. 5-18)"""))

cells.append(nbf.v4.new_code_cell(r"""import numpy as np
import matplotlib.pyplot as plt

V    = np.array([r["V"]      for r in result.iv])
C_LF = np.array([r["C_LF"]   for r in result.iv])
C_HF = np.array([r["C_HF"]   for r in result.iv])
Cox, Cmin, V_FB, V_T = t["Cox"], t["Cmin"], t["V_FB"], t["V_T"]

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(V, C_LF / Cox, "C0-",  label="LF (quasi-static)", lw=2)
ax.plot(V, C_HF / Cox, "C3--", label="HF",                lw=2)
ax.axhline(Cmin / Cox, color="grey", ls=":", lw=1,
           label=f"$C_{{min}}/C_{{ox}}={Cmin/Cox:.3f}$")
ax.axvline(V_FB, color="k", ls="--", lw=0.8, alpha=0.6); ax.text(V_FB, 0.05, "$V_{FB}$", ha="center")
ax.axvline(V_T,  color="k", ls="--", lw=0.8, alpha=0.6); ax.text(V_T,  0.05, "$V_T$",    ha="center")
ax.set_xlabel("$V_{gate}$ [V]"); ax.set_ylabel("$C / C_{ox}$")
ax.set_ylim(0, 1.1); ax.set_title("Axisymmetric MOSCAP, NA=1e17 cm$^{-3}$, t$_{ox}$=5 nm")
ax.grid(alpha=0.3); ax.legend(loc="lower right")
plt.show()"""))

cells.append(nbf.v4.new_markdown_cell(r"""## 7. Acceptance check

Quick numerical verification that the simulation hits the Cox / Cmin
plateaus."""))

cells.append(nbf.v4.new_code_cell(r"""acc = V <= V_FB - 1.0
inv = V >= V_T  + 1.5
print(f"C_LF (accumulation) / Cox  = {np.median(C_LF[acc]) / Cox:.3f}  (target ~1.0)")
print(f"C_LF (inversion)    / Cox  = {np.median(C_LF[inv]) / Cox:.3f}  (target ~1.0)")
print(f"C_HF (inversion)    / Cmin = {np.median(C_HF[inv]) / Cmin:.3f}  (target ~1.0)")"""))

nb = nbf.v4.new_notebook()
nb["cells"] = cells
nb["metadata"] = {
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {"name": "python"},
}
OUT.parent.mkdir(parents=True, exist_ok=True)
with OUT.open("w") as fh:
    nbf.write(nb, fh)
print("wrote", OUT)
