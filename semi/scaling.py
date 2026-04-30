"""
Nondimensional scaling for Poisson and drift-diffusion systems.

Why it matters. The raw semiconductor device equations have wildly
different scales: potentials ~ 1 V, carrier densities ~ 10^22 m^-3,
permittivity ~ 10^-11 F/m, charge ~ 10^-19 C. Feeding those directly
to Newton gives a Jacobian with condition number > 10^30 and the
solver fails. Scaling pulls everything to O(1).

Scales used:
    L0 = device characteristic length (m)
    V0 = thermal voltage kT/q (V)
    C0 = characteristic density, typically max|doping| (m^-3)
    mu0 = reference mobility (m^2/(V s))
    D0 = V0 * mu0 (Einstein diffusivity)
    t0 = L0^2 / D0
    J0 = q D0 C0 / L0

The dimensionless Poisson equation becomes

    -lambda^2 eps_r nabla^2 psi_hat = p_hat - n_hat + N_D_hat - N_A_hat

with lambda^2 = eps_0 V0 / (q C0 L0^2).
"""
from __future__ import annotations

from dataclasses import dataclass

from .constants import EPS0, Q, thermal_voltage


@dataclass
class Scaling:
    """Scale factors for nondimensionalizing the device equations."""

    L0: float      # characteristic length, m
    C0: float      # characteristic density, m^-3
    T: float       # temperature, K
    mu0: float     # reference mobility, m^2/(V s)
    n_i: float     # intrinsic density of reference material, m^-3

    @property
    def V0(self) -> float:
        """Thermal voltage, V."""
        return thermal_voltage(self.T)

    @property
    def D0(self) -> float:
        """Einstein diffusivity for reference mobility, m^2/s."""
        return self.V0 * self.mu0

    @property
    def t0(self) -> float:
        """Diffusion time over L0, s."""
        return self.L0 ** 2 / self.D0

    @property
    def J0(self) -> float:
        """Current density scale, A/m^2."""
        return Q * self.D0 * self.C0 / self.L0

    @property
    def lambda2(self) -> float:
        """
        Bare (eps_0-only) squared Debye-length-to-device-length ratio.
        Multiply by eps_r in the UFL form for the physical coefficient.
        Small values (<<1) mean the Debye length is much shorter than
        the device, which is the space-charge-dominant regime.
        """
        return EPS0 * self.V0 / (Q * self.C0 * self.L0 ** 2)

    def debye_length(self, n_ref: float | None = None) -> float:
        """
        Debye length in meters for a reference carrier density n_ref (m^-3).
        Uses C0 if n_ref not given.
        """
        if n_ref is None:
            n_ref = self.C0
        import math
        return math.sqrt(EPS0 * self.V0 / (Q * n_ref))

    def __repr__(self) -> str:
        return (
            f"Scaling(L0={self.L0:.2e} m, C0={self.C0:.2e} m^-3, T={self.T} K, "
            f"V0={self.V0:.4f} V, lambda2={self.lambda2:.3e})"
        )


def make_scaling_from_config(cfg: dict, reference_material) -> Scaling:
    """
    Construct a Scaling from a validated config dict and reference material.

    L0 is inferred from mesh extents (builtin) or defaults to 1 um (file-based).
    C0 is the max absolute doping across all profiles, floored at 1e16 cm^-3.
    """
    L0 = _infer_length(cfg)
    C0 = _infer_density(cfg)
    T = cfg["physics"]["temperature"]
    return Scaling(
        L0=L0, C0=C0, T=T,
        mu0=reference_material.mu_n,
        n_i=reference_material.n_i,
    )


def _infer_length(cfg: dict) -> float:
    mesh = cfg["mesh"]
    if mesh["source"] == "builtin":
        extents = mesh["extents"]
        return max(b[1] - b[0] for b in extents)
    if mesh["source"] == "builtin_axi":
        extents = mesh["extents"]
        r_range = extents["r"][1] - extents["r"][0]
        z_range = extents["z"][1] - extents["z"][0]
        return max(r_range, z_range)
    return 1.0e-6  # fallback for file-based meshes


def _infer_density(cfg: dict) -> float:
    from .constants import cm3_to_m3
    C_max = cm3_to_m3(1.0e16)  # floor
    for entry in cfg["doping"]:
        p = entry["profile"]
        if p["type"] == "uniform":
            C_max = max(C_max, cm3_to_m3(p["N_D"]), cm3_to_m3(p["N_A"]))
        elif p["type"] == "step":
            for k in ("N_D_left", "N_A_left", "N_D_right", "N_A_right"):
                C_max = max(C_max, cm3_to_m3(p[k]))
        elif p["type"] == "gaussian":
            C_max = max(C_max, cm3_to_m3(p["peak"]))
    return C_max
