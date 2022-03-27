using Gridap
using GridapGmsh
using GridapDistributed
using PartitionedArrays
using GridapPETSc
using LineSearches: BackTracking
using Revise

"""
2D and 3D channel flow
laminar - vorticity ω extracted
I assume is stable because is 2D
In 3D, same (nu) -> strange results (no parabolic profile for u)
Periodic BC in x-normal faces
For D=3, 3d case, wrong velocity results - instability?
?Import GMSH mesh (there we have better control of the mesh) - but how to set periodic BC in msh file
"""

include("Channel_Mesh.jl")

D=2; #dimensions number 2 or 3
N=32; #cells per dimensions

model=mesh_channel(D=D, N=N, printmodel=true)


#model = CartesianDiscreteModel(domain2d,partition2d, isperiodic=(true,false))

#model = CartesianDiscreteModel(domain3d,partition3d, isperiodic=(false,false,true))

#model = GmshDiscreteModel("Channel_geometry.msh")
#writevtk(model,"model")


order = 2
reffeᵤ = ReferenceFE(lagrangian,VectorValue{D,Float64},order)

if D==2
    wall_tags=["tag_6", "tag_5"]
    u_walls=VectorValue(0,0)
    nu = 0.0001472
    hf(x)=VectorValue(0.0033, 0)
elseif D==3
    wall_tags=["tag_24", "tag_23"]
    u_walls=VectorValue(0,0,0)
    nu=100
    hf(x)=VectorValue(0.0033, 0, 0);
end

V = TestFESpace(model,reffeᵤ,conformity=:H1,dirichlet_tags=wall_tags)


reffeₚ = ReferenceFE(lagrangian,Float64,order-1;space=:P)
Q = TestFESpace(model,reffeₚ,conformity=:L2, constraint=:zeromean)
h = 1;
Re = 10



u_t = Re*nu/h;
u_0=0;
#u_i(x) = VectorValue(u_0*x[2]/h-(h^2)*(1-(x[2]/h)^2)*G, 0)

U = TrialFESpace(V,[u_walls, u_walls])

P = TrialFESpace(Q)


Y = MultiFieldFESpace([V, Q])
X = MultiFieldFESpace([U, P])

degree = order
Ω = Triangulation(model)
dΩ = Measure(Ω,degree)



conv(u,∇u) =(∇u')⋅u
a((u,p),(v,q)) = ∫( nu*∇(v)⊙∇(u)   - (∇⋅v)*p + q*(∇⋅u))dΩ  - ∫( v⋅ hf )*dΩ

c(u,v) = ∫( v⊙(conv∘(u,∇(u))) )dΩ

res((u,p),(v,q)) = a((u,p),(v,q)) + c(u,v)

op = FEOperator(res,X,Y)

  

nls = NLSolver(
show_trace=true, method=:newton, linesearch=BackTracking())

solver = FESolver(nls)
uh, ph = solve(solver,op)

ω = ∇×uh;


writevtk(Ω,"results-channel-d$D",cellfields=["uh"=>uh,"ph"=>ph, "ω"=>ω])
