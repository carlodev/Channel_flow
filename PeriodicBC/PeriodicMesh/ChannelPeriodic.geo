// Gmsh project created on Tue Mar 22 08:30:46 2022
SetFactory("OpenCASCADE");
//+
pi = DefineNumber[ 3.14159, Name "Parameters/pi" ];
//+
N = DefineNumber[ 32, Name "Parameters/N" ];
//+
Lx = DefineNumber[ 2*pi, Name "Parameters/Lx" ];
//+
Ly = DefineNumber[ 2, Name "Parameters/Ly" ];
//+
Lz = DefineNumber[ 2/3*pi, Name "Parameters/Lz" ];
//+
Point(1) = {0, -Ly/2, -Lz/2, 1.0};
//+
Point(2) = {0, Ly/2, -Lz/2, 1.0};
//+
Point(3) = {Lx, -Ly/2, -Lz/2, 1.0};
//+
Point(4) = {Lx, Ly/2, -Lz/2, 1.0};
//+
Point(5) = {0, -Ly/2, Lz/2, 1.0};
//+
Point(6) = {0, Ly/2, Lz/2, 1.0};
//+
Point(7) = {Lx, -Ly/2, Lz/2, 1.0};
//+
Point(8) = {Lx, Ly/2, Lz/2, 1.0};
//+
Line(1) = {1, 5};
//+
Line(2) = {5, 7};
//+
Line(3) = {7, 8};
//+
Line(4) = {8, 6};
//+
Line(5) = {6, 5};
//+
Line(6) = {6, 2};
//+
Line(7) = {8, 4};
//+
Line(8) = {7, 3};
//+
Line(9) = {1, 2};
//+
Line(10) = {2, 4};
//+
Line(11) = {4, 3};
//+
Line(12) = {3, 1};
//+
Curve Loop(1) = {10, 11, 12, 9};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {8, 12, 1, 2};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {6, -9, 1, -5};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {3, 7, 11, -8};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {6, 10, -7, 4};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {5, 2, 3, 4};
//+
Plane Surface(6) = {6};
//+
Physical Surface("Inlet", 13) = {3};
//+
Physical Surface("Outlet", 14) = {4};
//+
Physical Surface("top_wall", 15) = {5};
//+
Physical Surface("bottom_wall", 16) = {2};
//+
Physical Surface("zmin", 17) = {1};
//+
Physical Surface("zmax", 18) = {6};
//+
Surface Loop(1) = {5, 3, 1, 4, 6, 2};
//+
Volume(1) = {1};
//+
Physical Volume("Fluid_V", 19) = {1};
//+
Transfinite Curve {1,2,3,4,5,6,7,8,9,10,11,12} = N Using Progression 1;
//+
Transfinite Surface {1,2,3,4,5,6};
//+
Recombine Surface {1,2,3,4,5,6};
//+
Transfinite Volume{1};
//+ 
Periodic Surface {{13}} = {{14}} Translate{1, 0, 0};
//+ 
Periodic Surface {{17}} = {{18} Translate{0, 0, 1};
