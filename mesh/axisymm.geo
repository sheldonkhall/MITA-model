L = 1e0;
M = 1e1;
eps = 1e-5;

ncl = 0.00002 / L;
fcl = 0.00125 / L;
dcl = 0.000015 / L;
ccl = 0.00002 / L;
scl = 0.000015 / L;

// Liver Tissue
Point(1) = { (0.+eps), 0. / M, 0., ncl };
Point(2) = { (0.+eps), 0.009105 / M, 0., ncl };
Point(3) = { 0.000895 / L, 0.01 / M, 0., ncl };
Point(4) = { 0.000895 / L, 0.016 / M, 0., ncl };
Point(5) = { 0.000895 / L, 0.017 / M, 0., ncl };
Point(6) = { 0.000895 / L, 0.08 / M, 0., ncl };
Point(7) = { 0.03 / L, 0.08 / M, 0., fcl };
Point(8) = { 0.03 / L, 0.00 / M, 0., fcl };

// Slot
Point(11) = { 0.000470 / L, 0.016 / M, 0., scl };
Point(12) = { 0.000470 / L, 0.017 / M, 0., scl };
Point(13) = { 0.000595 / L, 0.017 / M, 0., scl };
Point(14) = { 0.000595 / L, 0.016 / M, 0., scl };

// Dielectric
Point(21) = { 0.000470 / L, 0.08 / M, 0., dcl };
Point(22) = { 0.000470 / L, 0.010135 / M, 0., dcl };
Point(23) = { 0.000135 / L, 0.010135 / M, 0., dcl };
Point(24) = { 0.000135 / L, 0.08 / M, 0., dcl };

// Catheter
Point(31) = { 0.000595 / L, 0.08 / M, 0., ccl };
Point(32) = { 0.000595 / L, 0.01 / M, 0., ccl };
Point(33) = { (0.0+eps) / L, 0.01 / M, 0., ccl };

// Liver Tissue
Line(1) = { 1, 2 };
Line(2) = { 2, 3 };
Line(3) = { 3, 4 };
Line(4) = { 4, 5 };
Line(5) = { 5, 6 };
Line(6) = { 6, 7 };
Line(7) = { 7, 8 };
Line(8) = { 8, 1 };

// Slot
Line(11) = { 11, 12 };
Line(12) = { 12, 13 };
Line(13) = { 13, 14 };
Line(14) = { 14, 11 };

// Dielectric
Line(21) = { 21, 12 };
Line(22) = { 11, 22 };
Line(23) = { 22, 23 };
Line(24) = { 23, 24 };
Line(25) = { 24, 21 };

// Catheter
Line(31) = { 31, 13 };
Line(32) = { 14, 32 };
Line(33) = { 32, 33 };
Line(34) = { 33, 2 };
Line(35) = { 6, 31 };

// Liver Tissue
Line Loop(1) = { 1, 2, 3, 4, 5, 6, 7, 8 };

// Slot
Line Loop(11) = { 11, 12, 13, 14 };

// Dielectric
Line Loop(21) = { 21, -11, 22, 23, 24, 25 };

// Catheter
Line Loop(31) = { 31, 13, 32, 33, 34, 2, 3, 4, 5, 35 };

/* Subdomains */

// Liver Tissue
Plane Surface(1) = {1};

// Slot
Plane Surface(2) = {11};

// Dielectric
Plane Surface(3) = {21};

// Catheter
Plane Surface(4) = {31};

Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};
Physical Surface(4) = {4};

/* Boundaries */

// Scattering Boundary Condition
Physical Line(1) = { -35, 6, 7, 8 };

// Conductor Boundary Condition (Upper + Lower)
Physical Line(2) = { 31, -12, -21,
                     -33, -32, 14, 22, 23, 24 };

// Axisymmetric Boundary Condition
Physical Line(3) = { 1, -34 };

// Port Boundary Condition
Physical Line(4) = { 25 };
