// Gmsh project created on Wed Jan  8 16:44:32 2014

// this is the air case (figure 2)
// material matters as in paper k*a = 2.0
//                              b/a = 2.5
// infer that the length of the coaxial needs to be some multiplier of lambda
//     lc = n*2*pi/k
a = 0.0389;
b = 0.0974;
lc = 0.245;
ld = 1;
ms = 0.01;
as = 0.001;

// (figure 6a)
// material matters as in paper k*a = 0.1
//                              b/a = 1.5
// infer that the length of the coaxial needs to be some multiplier of lambda
//     lc = n*2*pi/k
//a = 0.00195;
//b = 0.00292;
//lc = 0.245;
//ld = 0.04;
//ms = 0.001;
//as = 0.0001;

Point(1) = {0, 0, 0, as};
Point(2) = {a, 0, 0, as};
Point(4) = {b, 0, 0, as};
Point(5) = {ld, 0, 0, ms};
Point(6) = {a, -lc, 0, as};
Point(7) = {b, -lc, 0, as};
Point(8) = {0, ld, 0, ms};
Point(9) = {ld, ld, 0, ms};
Line(1) = {1, 2};
Line(2) = {2, 4};
Line(3) = {4, 5};
Line(4) = {5, 9};
Line(5) = {9, 8};
Line(6) = {8, 1};
Line(7) = {2, 6};
Line(8) = {6, 7};
Line(9) = {7, 4};
Line Loop(10) = {6, 1, 2, 3, 4, 5};
Plane Surface(11) = {10};
Line Loop(12) = {7, 8, 9, -2};
Plane Surface(13) = {12};

// symmetry bc
Physical Line(1) = {6};
// PEC 
// Physical Line(15) = {1, 7, 9, 3};
// WPBC
Physical Line(4) = {8};
// scattering
Physical Line(2) = {4, 5};

// coax
Physical Surface(1) = {13};
// lossy tissue
Physical Surface(2) = {11};
