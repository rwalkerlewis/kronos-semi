// M7: 3D doped resistor gmsh fixture: 1 um x 200 nm x 200 nm rectangular bar.
//
// Physical groups:
//   - Volume "silicon"       (tag 1) covers the whole bar.
//   - Surface "contact_left" (tag 10) is the x = 0 face.
//   - Surface "contact_right"(tag 20) is the x = L face.
// The four remaining side faces are intentionally left without a
// physical group; they remain natural (insulating) boundaries for
// both Poisson and the continuity equations, matching the builtin
// create_box configuration.
//
// Regenerate with:
//   gmsh -3 -format msh22 box.geo -o box.msh
// The msh22 format keeps the file small (~2-5 KB for this resolution)
// and is supported by dolfinx.io.gmsh.read_from_msh.

SetFactory("OpenCASCADE");

L  = 1.0e-6;   // x-length (1 um)
W  = 2.0e-7;   // y,z cross-section side (200 nm)
lc = 5.0e-8;   // target mesh size (50 nm)

// Build a box [0, L] x [0, W] x [0, W].
Box(1) = {0, 0, 0, L, W, W};

// Coarsen/refine via the single characteristic length at all points.
MeshSize{ PointsOf{ Volume{1}; } } = lc;

// Physical groups for dolfinx.io.gmsh.read_from_msh.
Physical Volume("silicon", 1) = {1};

// Default OCC tagging:
//   Surface 1: x = 0 face   (contact_left)
//   Surface 2: x = L face   (contact_right)
// Remaining surfaces (3..6) are the y=0, y=W, z=0, z=W faces and
// are left untagged so the solver sees them as natural boundaries.
Physical Surface("contact_left", 10)  = {1};
Physical Surface("contact_right", 20) = {2};
