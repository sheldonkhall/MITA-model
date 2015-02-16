// Gmsh project created on Wed Aug  6 16:50:26 2014
am = 0.0006; // antenna mesh
cm = 0.01; // coarse tissue mesh
rd = 0.09; // radius of tissue domain
Point(1) = {0.0015, 0, 0, am};
Point(2) = {0.0015, Sqrt(rd^2-0.0015^2)-0.015, 0, cm};
Point(3) = {0, -0.03, 0, am};
Point(4) = {0.0015, -0.0285, 0, am};
Point(5) = {rd, -0.015, 0, cm};
Point(6) = {0.0, -0.015-rd, 0, cm};
Point(7) = {0.0, -0.015, 0.0, am};
Line(1) = {3, 4};
Line(2) = {4, 1};
Line(3) = {1, 2};
Line(4) = {3, 6};
Circle(5) = {2, 7, 5};
Circle(6) = {5, 7, 6};
Line Loop(7) = {5, 6, -4, 1, 2, 3};
Plane Surface(8) = {7};
Physical Line(4) = {1, 2}; // active tip
Physical Line(3) = {3}; // insulating probe surface
Physical Line(2) = {5, 6}; //bulk tissue boundary
Physical Line(1) = {4}; // line of cylindrical symmetry
Physical Surface(1) = {8}; // tissue