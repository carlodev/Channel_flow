using Gridap
using GridapGmsh
using GridapDistributed
using PartitionedArrays
using GridapPETSc
using LineSearches: BackTracking


"""
For stability I imposed nu=100; 
"""


Lx=2*pi;
Ly=2;
Lz=2/3*pi;
nx = 20;
ny = 10 ;
nz = 12;
domain2d = (0,Lx,-Ly/2,Ly/2)
partition2d = (nx,ny)
domain3d = (0,Lx,-Ly/2,Ly/2,-Lz/2,Lz/2)
partition3d = (nx,ny,nz)
#model = CartesianDiscreteModel(domain2d,partition2d, isperiodic=(false,false))

model = CartesianDiscreteModel(domain3d,partition3d, isperiodic=(true,false,true))

#model = GmshDiscreteModel("Channel_geometry.msh")
#writevtk(model,"model3d")


labels = get_face_labeling(model);

D = 2
order = 2
reffeᵤ = ReferenceFE(lagrangian,VectorValue{3,Float64},order)
#V = TestFESpace(model,reffeᵤ,conformity=:H1,dirichlet_tags=["tag_6", "tag_5", "tag_7", "tag_8"]) #tag6 top wall; tag7 laft wall. ; tag5 bottom wall; tag8 right wall
V = TestFESpace(model,reffeᵤ,conformity=:H1,dirichlet_tags=["tag_24", "tag_23"]) #tag24 top wall;  tag23 bottom wall; tag25 inlet
#V = TestFESpace(model,reffeᵤ,conformity=:H1,labels=labels,dirichlet_tags=["Inlet", "Bottom_Wall","Top_Wall"])


#reffeₚ = ReferenceFE(lagrangian,Float64,order-1;space=:P)


reffeₚ = ReferenceFE(lagrangian,Float64,order-1;space=:P)
Q = TestFESpace(model,reffeₚ,conformity=:L2, constraint=:zeromean)
h = Ly/2;
Re = 10
G = -0.01;
nu = 0.0001472;
u_0 = Re*nu/h;
u_0=0;
#u_i(x) = VectorValue(u_0*x[2]/h-(h^2)*(1-(x[2]/h)^2)*G, 0,0)
#u_i(x) = VectorValue(u_0, 0)

u_top = VectorValue(0,0,0)
u_bottom = VectorValue(0,0,0)
u_walls=VectorValue(0,0)
U = TrialFESpace(V,[u_top, u_bottom])
#U = TrialFESpace(V,[u_top, u_bottom])

P = TrialFESpace(Q)
#P = TrialFESpace(Q)

#P = TrialFESpace(Q,[1000,-1000])


Y = MultiFieldFESpace([V, Q])
X = MultiFieldFESpace([U, P])

degree = order
Ωₕ = Triangulation(model)
dΩ = Measure(Ωₕ,degree)



hf(x)=VectorValue(0.0033, 0, 0);

conv(u,∇u) =(∇u')⋅u
a((u,p),(v,q)) = ∫( 100*∇(v)⊙∇(u) - (∇⋅v)*p + q*(∇⋅u) )dΩ- ∫( v⋅hf )*dΩ
#a((u,p),(v,q)) = ∫( nu*∇(v)⊙∇(u) + q*(∇⋅u) )dΩ +∫( v*hf )*dΩ

c(u,v) = ∫( v⊙(conv∘(u,∇(u))) )dΩ

res((u,p),(v,q)) = a((u,p),(v,q)) + c(u,v)

op = FEOperator(res,X,Y)

  

nls = NLSolver(
show_trace=true, method=:newton, linesearch=BackTracking())

solver = FESolver(nls)
uh, ph = solve(solver,op)

ω = ∇×uh;
writevtk(Ωₕ,"results-ch3d-dpdx",cellfields=["uh"=>uh,"ph"=>ph, "ω"=>ω])
