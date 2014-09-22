// Gmsh project created on Wed Jun 11 11:38:25 2014
ms = 0.0005;
Point(1) = {0, -0.02, 0, ms};
Point(2) = {0.02, -0.02, 0, ms};
Point(3) = {0.02, 0.02, 0, ms};
Point(4) = {0, 0.02, 0, ms};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {2, 3, 4, 1};
Plane Surface(6) = {5};
Physical Line(1) = {4};
Physical Line(5) = {3, 2, 1};
Physical Surface(1) = {6};
