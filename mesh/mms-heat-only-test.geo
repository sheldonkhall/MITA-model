// Gmsh project created on Thu Jan 16 08:45:44 2014
R = 1;
L = 10;
ms = 0.1;
Point(1) = {0, -L, 0, ms};
Point(2) = {0, L, 0, ms};
Point(3) = {R, L, 0, ms};
Point(4) = {R, -L, 0, ms};
Line(1) = {1, 4};
Line(2) = {4, 3};
Line(3) = {3, 2};
Line(4) = {2, 1};
Line Loop(5) = {2, 3, 4, 1};
Plane Surface(6) = {5};
Physical Surface(1) = {6}; // tissue
Physical Line(5) = {3, 2, 1}; // external boundary
Physical Line(1) = {4}; // symmetry boundary
