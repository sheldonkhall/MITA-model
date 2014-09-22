// Gmsh project created on Thu Sep  4 14:41:39 2014
Point(1) = {0.0, 0.0, 0, .03};
Point(2) = {0.0, 3.141592, 0, .03};
Point(3) = {1.5707963, 3.141592, 0, .03};
Point(4) = {1.5707963, 0.0, 0, .03};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {2, 3, 4, 1};
Plane Surface(6) = {5};
Physical Surface(1) = {6};
Physical Line(5) = {2, 3, 4};
Physical Line(1) = {1};
