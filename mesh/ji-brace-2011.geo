// Gmsh project created on Wed May 28 15:01:15 2014
am = 0.00003;
lm = 0.005;
wp = 0.01;
ds = 0.05;
Point(1) = {0, -0.013, 0, am};
Point(2) = {0.000255, -0.013, 0, am};
Point(3) = {0.000255, wp, 0, am};
Point(4) = {0.00084, wp, 0, am};
Point(5) = {0.00084, 0, 0, am};
Point(6) = {0.0011, 0, 0, am};
Point(7) = {0.0011, ds, 0, lm};
Point(8) = {0.05, ds, 0, lm};
Point(9) = {0.05, -ds, 0, lm};
Point(10) = {0, -ds, 0, lm};
Point(11) = {0.000255, 0, 0, am};
Line(1) = {1, 2};
Line(2) = {2, 11};
Line(3) = {11, 3};
Line(4) = {3, 4};
Line(5) = {4, 5};
Line(6) = {5, 11};
Line(7) = {5, 6};
Line(8) = {6, 7};
Line(9) = {7, 8};
Line(10) = {8, 9};
Line(11) = {9, 10};
Line(12) = {10, 1};
Line Loop(13) = {12, 1, 2, -6, 7, 8, 9, 10, 11};
Plane Surface(14) = {13};
Line Loop(15) = {3, 4, 5, 6};
Plane Surface(16) = {15};
Physical Line(1) = {12};
Physical Line(2) = {9, 10, 11};
Physical Line(4) = {4};
Physical Line(5) = {9, 10, 11};
Physical Surface(2) = {16};
Physical Surface(1) = {14};
