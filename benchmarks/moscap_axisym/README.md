# Axisymmetric 2D MOSCAP — LF vs HF C–V (Hu Fig. 5-18)

Cylindrical (r, z) MOS capacitor reproducing the qualitative LF / HF
contrast in **Chenming Hu, Modern Semiconductor Devices for Integrated
Circuits, Ch. 5, Fig. 5-18**.

## Device

| Parameter        | Value                |
| ---------------- | -------------------- |
| Body             | P-type Si, NA=1e17 cm⁻³ |
| Oxide            | SiO₂, 5 nm           |
| Body depth       | 1 µm                 |
| Gate radius      | 1 µm                 |
| Outer radius     | 1 µm (no fringing)   |
| Gate workfn (φms)| −0.977 V (≈ flat-band)|
| T                | 300 K                |

The (r, z) computational domain is the half-cross-section; the full 3D
solid is recovered by rotation about r = 0. The symmetry-axis natural
condition is automatic via the r-weighted volume measure.

## Run

```
python -m semi.run benchmarks/moscap_axisym/moscap_axisym.json
```

Outputs the manifest, IV CSV, and field artifacts under
`./results/moscap_axisym/`. Each IV row contains `{V, Q_gate, J, C_LF,
C_HF}`; capacitances are in F/m² (per unit gate area).

## Expected curves

* `C_LF` returns to `Cox` in both deep accumulation and strong
  inversion — minority carriers track the AC signal in the LF limit.
* `C_HF` returns to `Cox` only in deep accumulation; in inversion it
  saturates at `C_min = Cox * Cdep,min / (Cox + Cdep,min)`.
* The two curves coincide through accumulation, flat band, and depletion;
  they diverge only past `V_T`.

## Reference

Hu §5.6 and Fig. 5-18:
<https://www.chu.berkeley.edu/wp-content/uploads/2020/01/Chenming-Hu_ch5-1.pdf>
