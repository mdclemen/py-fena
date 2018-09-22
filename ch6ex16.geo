//+
SetFactory("OpenCASCADE");
//+
lc = 0.3;
//+
Point(1) = {0, 0, 0, lc};
//+
Point(2) = {1, 0, 0, lc};
//+
Point(3) = {1, .5, 0, lc};
//+
Point(4) = {.5, 1, 0, lc};
//+
Point(5) = {0, 1, 0, lc};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 1};
//+
Curve Loop(1) = {4, 5, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
//Characteristic Length {1, 5, 4, 3, 2} = lc;
//+
Physical Point("points") = {1, 2, 3, 4, 5};
//+
Physical Curve("lines") = {1, 2, 3, 4, 5};
//+
Physical Surface("surface") = {1};
//+
Compound Curve {1, 2, 3, 4, 5};
