// Gmsh project created on Wed Aug  6 16:50:26 2014
am = 0.0005;
cm = 0.01;
at = 0.03;
Point(1) = {0.0015, 0, 0, am};
Point(2) = {0.0015, 0.06, 0, cm};
Point(3) = {0.06, 0.06, 0, cm};
Point(4) = {0.06, -0.1, 0, cm};
Point(5) = {0, -0.1, 0, cm};
Point(6) = {0, -0.03, 0, am};
Point(7) = {0.0015, -0.0285, 0, am};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 1};
Line Loop(8) = {2, 3, 4, 5, 6, 7, 1};
Plane Surface(9) = {8};
Physical Line(4) = {6, 7};// active tip
Physical Line(3) = {1}; // insulating probe surface
Physical Line(2) = {2, 3, 4}; //bulk tissue boundary
Physical Line(1) = {5}; // line of cylindrical symmetry
Physical Surface(1) = {9}; // tissue
//Physical Line(5) = {2, 3, 4}; //bulk tissue boundary (thermal)
