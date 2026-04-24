// 2D long-channel n-MOSFET cross-section (M12).
//
//   y (m)
//     |        gate contact (tag 13)
//  2.005e-6    +--------+
//              | oxide  |  (tag 2)
//  2.0e-6 +----+--------+----+
//         |src | Si body | dr|   y=2e-6 top of silicon
//         |    |         |   |
//         |    | (tag 1) |   |
//         |    |         |   |
//         +----+---------+---+   y=0  body contact (tag 12, full width)
//         0   1.5e-6   3.5e-6  5e-6  (x)
//
// Geometry
//   Body silicon:  x in [0, 5e-6], y in [0, 2e-6]
//   Gate oxide:    x in [1.5e-6, 3.5e-6], y in [2e-6, 2.005e-6]
//   Source contact line: y = 2e-6, x in [0, 1e-6]       (tag 10)
//   Drain contact line:  y = 2e-6, x in [4e-6, 5e-6]    (tag 11)
//   Body contact line:   y = 0,    x in [0, 5e-6]       (tag 12)
//   Gate contact line:   y = 2.005e-6, x in [1.5e-6, 3.5e-6]  (tag 13)
//
// Sizing is driven per-point. The top-level -clmax from the engine
// command line caps all of these uniformly; the convergence test uses
// -clmax = {1e-7, 5e-8, 2.5e-8}.

cl_bulk    = 2.0e-7;  // far from the channel
cl_contact = 7.5e-8;  // source/drain leading edges at x = {0, 5e-6}
cl_channel = 3.0e-8;  // near the channel under the oxide
cl_gate    = 5.0e-9;  // through-oxide (oxide is only 5 nm thick)

Lx   = 5.0e-6;
Ly   = 2.0e-6;
tox  = 5.0e-9;

xs1  = 1.0e-6;    // source contact right edge (start of ungated top)
xg1  = 1.5e-6;    // gate-oxide left edge
xg2  = 3.5e-6;    // gate-oxide right edge
xd1  = 4.0e-6;    // drain contact left edge

// Silicon body corners
Point(1) = {0.0,  0.0, 0.0, cl_contact};
Point(2) = {Lx,   0.0, 0.0, cl_contact};
Point(3) = {Lx,   Ly,  0.0, cl_contact};
Point(4) = {0.0,  Ly,  0.0, cl_contact};

// Top-of-silicon break points
Point(5)  = {xs1, Ly, 0.0, cl_contact};   // source end
Point(6)  = {xg1, Ly, 0.0, cl_channel};   // gate-oxide start
Point(7)  = {xg2, Ly, 0.0, cl_channel};   // gate-oxide end
Point(8)  = {xd1, Ly, 0.0, cl_contact};   // drain start

// Oxide top corners
Point(9)  = {xg1, Ly + tox, 0.0, cl_gate};
Point(10) = {xg2, Ly + tox, 0.0, cl_gate};

// Silicon-body edges (counter-clockwise starting at bottom-left)
Line(1) = {1, 2};        // body contact (bottom)
Line(2) = {2, 3};        // right edge
Line(3) = {3, 8};        // top right, drain -> right corner
Line(4) = {8, 7};        // drain contact, right-gate-corner -> drain-start
Line(5) = {7, 6};        // under-gate silicon top
Line(6) = {6, 5};        // source-right to gate-left on silicon top
Line(7) = {5, 4};        // source contact, left-corner -> source end
Line(8) = {4, 1};        // left edge

// Oxide edges (counter-clockwise: bottom on silicon, up right, across top, down left)
Line(9)  = {6, 7};       // oxide bottom (shares geometry with line 5, but oriented opposite)
Line(10) = {7, 10};      // oxide right
Line(11) = {10, 9};      // oxide top (gate contact)
Line(12) = {9, 6};       // oxide left

// Surfaces
Line Loop(100) = {1, 2, 3, 4, 5, 6, 7, 8};
Plane Surface(1) = {100};

Line Loop(200) = {9, 10, 11, 12};
Plane Surface(2) = {200};

// Physical groups -- names and tags must match the input JSON's
// geometry.physical_groups map.
Physical Surface("silicon", 1) = {1};
Physical Surface("oxide",   2) = {2};

// Source contact line: x in [0, xs1] at y=Ly, which is Line 7 traversed
// from point 5 (xs1) to point 4 (0). Note Line(7) = {5, 4}.
Physical Line("source", 10) = {7};

// Drain contact line: x in [xd1, Lx] at y=Ly, which is Line 3 = {3, 8}
// from point 3 (Lx, Ly) to point 8 (xd1, Ly).
Physical Line("drain", 11) = {3};

// Body contact line: full bottom, Line 1.
Physical Line("body", 12) = {1};

// Gate contact line: top of the oxide, Line 11.
Physical Line("gate", 13) = {11};

// Make sure gmsh writes physical entities only (no orphans).
Mesh.SaveAll = 0;
Mesh.MshFileVersion = 2.2;
