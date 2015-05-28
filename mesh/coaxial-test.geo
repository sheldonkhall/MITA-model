dcl = 0.00001;
dncl = 0.00001;

// Dielectric
Point(21) = { 0.000470, 0.005, 0., dcl };
Point(22) = { 0.000470, 0.0, 0., dcl };
Point(23) = { 0.000135, 0.0, 0., dncl };
Point(24) = { 0.000135, 0.005, 0., dncl };

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

Physical Surface(1) = {3};

/* Boundaries */

// Scattering Boundary Condition
Physical Line(2) = { 25 };

// Conductor Boundary Condition (Upper + Lower)
// Physical Line(2) = { 21, 24 };

// Port Boundary Condition
Physical Line(4) = { 23 };