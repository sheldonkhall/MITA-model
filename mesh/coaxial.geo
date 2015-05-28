L = 1;
M = 1e2;

ncl = 0.00002 / L;
fcl = 0.00125 / L;
dcl = 0.00002 / L;
dncl = 0.000005 / L;
ccl = 0.00005 / L;
scl = 0.000025 / L;

// Dielectric
Point(21) = { 0.000470/L, 0.08/M, 0., dcl };
Point(22) = { 0.000470/L, 0.0/M, 0., dcl };
Point(23) = { 0.000135/L, 0.0/M, 0., dncl };
Point(24) = { 0.000135/L, 0.08/M, 0., dncl };

// Dielectric
Line(21) = { 21, 22 };
Line(23) = { 22, 23 };
Line(24) = { 23, 24 };
Line(25) = { 24, 21 };

// Dielectric
Line Loop(21) = { 21, 23, 24, 25 };

/* Subdomains */

// Dielectric
Plane Surface(3) = {21};

Physical Surface(3) = {3};

/* Boundaries */

// Scattering Boundary Condition
Physical Line(1) = { 23 };

// Conductor Boundary Condition (Upper + Lower)
Physical Line(2) = { 21, 24 };

// Port Boundary Condition
Physical Line(4) = { 25 };
