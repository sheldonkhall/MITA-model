as = 0.00003; // antenna mesh size
r1 = 0.000135; // inner conductor
r2 = 0.000470; // inside outer conductor
r3 = 0.000595; // outside outer conductor
r4 = 0.000895; // catheter
ms = 0.005; // coarse mesh size
Point(1) = {0, 0.010, 0, as};
Point(2) = {0, 0.0095, 0, as};
Point(3) = {r1, 0.025, 0, as};
Point(4) = {r1, 0.01+(r3-r2), 0, as};
Point(5) = {r2, 0.025, 0, as};
Point(6) = {r2, 0.0165, 0, as};
Point(7) = {r2, 0.0155, 0, as};
Point(8) = {r2, 0.01+(r3-r2), 0, as};
Point(9) = {r3, 0.05, 0, as};
Point(10) = {r3, 0.0165, 0, as};
Point(11) = {r3, 0.0155, 0, as};
Point(12) = {r3, 0.01, 0, as};
Point(13) = {r4, 0.05, 0, as};
Point(14) = {r4, 0.01, 0, as};
//Point(15) = {0.05, 0.05, 0, ms};
//Point(16) = {0.05, -0.03, 0, ms};
Point(17) = {0, -0.01, 0, ms};
Point(18) = {0, 0.02, 0, as};
Point(19) = {0.03, 0.02, 0, ms};
Line(1) = {3, 4};
Line(2) = {4, 8};
Line(3) = {8, 7};
Line(4) = {7, 11};
Line(5) = {11, 12};
Line(6) = {12, 1};
Line(7) = {1, 2};
Line(8) = {2, 14};
Line(9) = {14, 13};
Line(10) = {13, 9};
Line(11) = {9, 10};
Line(12) = {10, 6};
Line(13) = {6, 5};
Line(14) = {5, 3};
Line(15) = {6, 7};
Line(16) = {10, 11};
//Line(17) = {13, 15};
//Line(18) = {15, 16};
//Line(19) = {16, 17};
Line(20) = {17, 2};
Line(21) = {12, 14};
Line Loop(21) = {1, 2, 3, -15, 13, 14};
Plane Surface(22) = {21};
Line Loop(23) = {15, 4, -16, 12};
Plane Surface(24) = {23};
//Line Loop(27) = {17, 18, 19, 20, 8, 9};
//Plane Surface(28) = {27};
Line Loop(29) = {7, 8, -21, 6};
Plane Surface(30) = {29};
Line Loop(31) = {5, 21, 9, 10, 11, 16};
Plane Surface(32) = {31};
Physical Line(1) = {20, 7}; // symmetry
//Physical Line(2) = {19, 18, 17, 10}; // abc
Physical Line(4) = {14}; // wpbc
//Physical Line(5) = {19, 18, 17}; // dirichlet temp bc
//Physical Surface(1) = {28}; // liver
Physical Surface(2) = {22}; // dielectric (coaxial)
Physical Surface(3) = {24}; // air gap
Physical Surface(4) = {30,32}; // catheter


Circle(33) = {13, 18, 19};
Circle(34) = {19, 18, 17};
Line Loop(35) = {33, 34, 20, 8, 9};
Plane Surface(36) = {35};
Physical Line(5) = {33, 34}; // dirichlet temp bc
Physical Surface(1) = {36}; // liver
Physical Line(2) = {33, 34, 10}; // abc