using Gridap
using GridapGmsh
using GridapDistributed
using PartitionedArrays
using GridapPETSc
using LineSearches: BackTracking

import Gridap: ∇
import GridapODEs.TransientFETools: ∂t

"""
2D channel flow
laminar - vorticity ω extracted
I assume is stable because is 2D
In 3D, same (nu) -> strange results (no parabolic profile for u)
Periodic BC in x-normal faces
? define e CartesianDiscreteModel with no equally spaced nodes; In y direction custom function (Nodey) by Bazilev
?Import GMSH mesh (there we have better control of the mesh) - but how to set periodic BC in msh file
"""


Lx=2*pi;
Ly=2;
Lz=2/3*pi;
nx = 32;
ny = 32 ;
nz = 32;
domain2d = (0,Lx,-Ly/2,Ly/2)
partition2d = (nx,ny)
domain3d = (0,Lx,-Ly/2,Ly/2,-Lz/2,Lz/2)
partition3d = (nx,ny,nz)




model = CartesianDiscreteModel(domain2d,partition2d, isperiodic=(true,false))

#model = CartesianDiscreteModel(domain3d,partition3d, isperiodic=(false,false,true))

#model = GmshDiscreteModel("Channel_geometry.msh")
#writevtk(model,"model")


labels = get_face_labeling(model);

D = 2
order = 2
reffeᵤ = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
#V = TestFESpace(model,reffeᵤ,conformity=:H1,dirichlet_tags=["tag_6", "tag_5", "tag_7", "tag_8"]) #tag6 top wall; tag7 laft wall. ; tag5 bottom wall; tag8 right wall
V = TestFESpace(model,reffeᵤ,conformity=:H1,dirichlet_tags=["tag_6", "tag_5"])
#V = TestFESpace(model,reffeᵤ,conformity=:H1,labels=labels,dirichlet_tags=["Inlet", "Bottom_Wall","Top_Wall"])


#reffeₚ = ReferenceFE(lagrangian,Float64,order-1;space=:P)


reffeₚ = ReferenceFE(lagrangian,Float64,order-1;space=:P)
Q = TestFESpace(model,reffeₚ,conformity=:L2, constraint=:zeromean)
h = Ly/2;
Re = 10
G = -0.01;
nu = 0.0001472;
u_t = Re*nu/h;
u_0=0;
#u_i(x) = VectorValue(u_0*x[2]/h-(h^2)*(1-(x[2]/h)^2)*G, 0)
#u_i(x) = VectorValue(u_0, 0)

u_top = VectorValue(u_0,0)
u_bottom = VectorValue(-u_0,0)

u_wall(x,t) = VectorValue(0.0, 0.0)
u_wall(t::Real) = x -> u_wall(x,t)

U = TransientTrialFESpace(V, [u_wall, u_wall])
#U = TrialFESpace(V,[u_top, u_bottom])
h1(x) = -0.0033*x[1];
P = TransientTrialFESpace(Q)
#P = TrialFESpace(Q)

#P = TrialFESpace(Q,[1000,-1000])


Y = TransientMultiFieldFESpace([V, Q])
X = TransientMultiFieldFESpace([U, P])

degree = order
Ωₕ = Triangulation(model)
dΩ = Measure(Ωₕ,degree)



hf(x)=VectorValue(0.0033, 0);

conv(u,∇u) =(∇u')⋅u
a((u,p),(v,q)) = ∫( nu*∇(v)⊙∇(u)  - (∇⋅v)*p + q*(∇⋅u) )dΩ  - ∫( v⋅hf )*dΩ
#a((u,p),(v,q)) = ∫( (1/Re)*∇(v)⊙∇(u) + q*(∇⋅u) )dΩ +∫( v*hf )*dΩ

c(u,v) = ∫( v⊙(conv∘(u,∇(u))) )dΩ

res(u,p),(v,q)) = a((u,p),(v,q)) + c(u,v)

op = TransientFEOperator(res,X,Y)


U0 = U(0.0)
P0 = P(0.0)
X0 = X(0.0)

uh0 = interpolate_everywhere(VectorValue(0.0, 0.0), U0)
ph0 = interpolate_everywhere(0.0, P0)
xh0 = interpolate_everywhere([uh0, ph0], X0)


t0 = 0.0
tF = 2
dt = 0.2
θ = 1
print("Nls")

using LineSearches: BackTracking
nls = NLSolver(
  show_trace=true, method=:newton, linesearch=BackTracking())

ode_solver = ThetaMethod(nls,dt,θ)
print("starting solving")

sol_t = solve(ode_solver,op,xh0,t0,tF)


ω = ∇×uh;
_t_nn = t0
createpvd("Channel2d_td") do pvd
  for (xh_tn, tn) in sol_t
    global _t_nn
    _t_nn += dt
    uh_tn = xh_tn[1]
    ph_tn = xh_tn[2]
    pvd[tn] = createvtk(Ω,"Channel_2d_$_t_nn"*".vtu",cellfields=["uh"=>uh_tn,"ph"=>ph_tn])
  end

end