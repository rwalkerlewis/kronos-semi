// M19: 3D n-channel MOSFET capstone geometry.
//
// A planar bulk NMOS: a p-type silicon body with a thin SiO2 gate
// oxide grown over the full top footprint. Source and drain are the
// two x-end faces of the silicon (n+ counter-doped by Gaussian
// implants declared in the JSON); the body contact is the silicon
// bottom face; the gate contact is the oxide top face.
//
//   x : source-to-drain (channel) direction, [0, L_total]
//   y : channel-width direction,             [0, W]
//   z : depth / gate-stack direction,        [0, H_si + t_ox]
//
// Physical volumes:
//   silicon (tag 1) : z in [0, H_si]
//   oxide   (tag 2) : z in [H_si, H_si + t_ox]
//
// Physical surfaces:
//   source (tag 10) : silicon face x = 0
//   drain  (tag 11) : silicon face x = L_total
//   gate   (tag 12) : oxide top face z = H_si + t_ox
//   body   (tag 13) : silicon bottom face z = 0
//   all remaining faces are left untagged (natural insulating,
//   tag 0), matching the resistor_3d / builtin-box convention.
//
// UNITS: the geometry is built in MICROMETERS because OpenCASCADE's
// default linear tolerance (~1e-7) chokes on nanometre-scale solids
// (a 5 nm oxide box in metre units fails the boolean fragment). The
// final mesh is rescaled to METERS on output via Mesh.ScalingFactor
// = 1e-6, so the written node coordinates are in meters per ADR 0002
// / PLAN invariant 3. Every length below is therefore in um.
//
// Regenerate the shipped ~200k-DOF mesh with:
//   gmsh -3 mosfet3d.geo -o mosfet3d.msh
// or, to sweep the characteristic length (smoke / GPU meshes):
//   gmsh -3 -setnumber cl 0.05 mosfet3d.geo -o mosfet3d_coarse.msh
// The bundled generator `benchmarks/mosfet_3d/generate_mesh.py`
// overrides `cl` via the gmsh Python API for reproducible builds.

SetFactory("OpenCASCADE");

// Device dimensions in micrometers.
L_total = 1.0;      // body length (source to drain)  = 1 um
W       = 1.0;      // body width                     = 1 um
H_si    = 0.200;    // silicon body thickness         = 200 nm
t_ox    = 0.005;    // gate oxide thickness           = 5 nm

// Characteristic mesh length in um. 0.02 um (20 nm) yields ~200k DOFs
// on the coupled (psi, phi_n, phi_p) Slotboom system. Overridable via
// -setnumber cl <value_in_um>.
DefineConstant[ cl = 0.02 ];

eps = 1.0e-4;  // bounding-box selection tolerance (0.1 nm in um)

// Emit node coordinates in meters (um geometry -> m mesh).
Mesh.ScalingFactor = 1.0e-6;

// Silicon body and oxide slab as two stacked boxes.
Box(1) = {0, 0, 0,    L_total, W, H_si};   // silicon
Box(2) = {0, 0, H_si, L_total, W, t_ox};   // oxide

// Glue the stack so the Si/SiO2 interface is a single shared internal
// surface (conformal mesh across the junction).
BooleanFragments{ Volume{1}; Delete; }{ Volume{2}; Delete; }

// Uniform target size at every geometry point.
MeshSize{ PointsOf{ Volume{:}; } } = cl;

// --- Physical volumes (select by bounding box; fragment renumbers) ---
sil_vol() = Volume In BoundingBox
  { -eps, -eps, -eps,  L_total + eps, W + eps, H_si + eps };
ox_vol()  = Volume In BoundingBox
  { -eps, -eps, H_si - eps,  L_total + eps, W + eps, H_si + t_ox + eps };
Physical Volume("silicon", 1) = { sil_vol() };
Physical Volume("oxide", 2)   = { ox_vol() };

// --- Physical surfaces (select the device faces by bounding box) ---
src_surf()  = Surface In BoundingBox
  { -eps, -eps, -eps,  eps, W + eps, H_si + eps };
drn_surf()  = Surface In BoundingBox
  { L_total - eps, -eps, -eps,  L_total + eps, W + eps, H_si + eps };
gate_surf() = Surface In BoundingBox
  { -eps, -eps, H_si + t_ox - eps,  L_total + eps, W + eps, H_si + t_ox + eps };
body_surf() = Surface In BoundingBox
  { -eps, -eps, -eps,  L_total + eps, W + eps, eps };

Physical Surface("source", 10) = { src_surf() };
Physical Surface("drain", 11)  = { drn_surf() };
Physical Surface("gate", 12)   = { gate_surf() };
Physical Surface("body", 13)   = { body_surf() };
