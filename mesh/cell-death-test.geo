// Gmsh project created on Sun Jun  8 16:19:58 2014
ms = 0.1;
Point(1) = {0, 0, 0, ms};
Point(2) = {1, 0, 0, ms};
Point(3) = {1, 1, 0, ms};
Point(4) = {0, 1, 0, ms};
Line(1) = {1, 2};
Line(3) = {2, 3};
Line(4) = {3, 4};
Line(5) = {4, 1};
Line Loop(6) = {4, 5, 1, 3};
Plane Surface(7) = {6};
Physical Line(1) = {5};
Physical Surface(1) = {7};
