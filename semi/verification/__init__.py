"""
Verification & Validation subpackage.

This module hosts Method-of-Manufactured-Solutions (MMS) studies, mesh
convergence sweeps, and conservation checks. Unlike the pure-Python
core (constants, materials, scaling, doping, schema, continuation,
diode_analytical) this subpackage is allowed to import dolfinx because
MMS necessarily exercises the full FEM stack.

Public submodules:
    mms_poisson       MMS for the equilibrium Poisson equation
"""
