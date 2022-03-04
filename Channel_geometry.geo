// Gmsh project created on Fri Mar 04 07:49:03 2022
SetFactory("OpenCASCADE");
//+
Point(1) = {0, -0.5, 0, 1.0};
//+
Point(2) = {0, 0.5, 0, 1.0};
//+
Point(3) = {20, 0.5, 0, 1.0};
//+
Point(4) = {20, -0.5, 0, 1.0};
//+
Point(5) = {0, 0, 0, 1.0};
//+
Point(6) = {20, 0, 0, 1.0};
//+
Line(1) = {1, 5};
//+
Line(2) = {2, 5};
//+
Line(3) = {2, 3};
//+
Line(4) = {3, 6};
//+
Line(5) = {4, 6};
//+
Line(6) = {4, 1};
//+
Line(7) = {5, 6};
//+
Curve Loop(1) = {2, 3, 4, 7};
//+
Curve Loop(2) = {1, 7, 5, 6};
//+
Plane Surface(1) = {1};
//+
Plane Surface(2) = {2};
//+
Transfinite Curve {1, 2} = 1000 Using Bump 1;
//+
Physical Curve("Inlet", 8) = {1, 2};
//+
Physical Curve("Outlet", 9) = {5, 4};
//+
Physical Curve("Top_Wall", 10) = {3};
//+
Physical Curve("Bottom_Wall", 11) = {6};
//+
Physical Surface("Fluid", 12) = {1, 2};
//+
Transfinite Curve {3, 6, 7} = 500 Using Progression 1;
//+
Transfinite Curve {2, 1, 4, 5} = 100 Using Progression 1.05;
//+
Transfinite Surface {1, 2};
//+
Recombine Surface {1, 2};

