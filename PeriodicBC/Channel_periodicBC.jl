using Revise
using Gridap
using GridapGmsh
using GridapDistributed
using PartitionedArrays
using GridapPETSc
using LineSearches: BackTracking

Lx=2*pi;
Ly=2;
Lz=2/3*pi;
nx = 12;
ny = 12;
nz = 12;
domain = (0,Lx,-Ly/2,Ly/2,-Lz/2,Lz/2,)
partition = (nx,ny,nz)
model = CartesianDiscreteModel(domain,partition, isperiodic=(false,false,true))
#model = GmshDiscreteModel("ChannelPeriodic5.msh")
writevtk(model,"model")
D = 2
order = 2
reffeᵤ = ReferenceFE(lagrangian,VectorValue{3,Float64},order)
#V = TestFESpace(model,reffeᵤ,conformity=:H1,dirichlet_tags=["top_wall", "bottom_wall", "Inlet"])
V = TestFESpace(model,reffeᵤ,conformity=:H1,dirichlet_tags=["tag_23", "tag_24", "tag_25"])
reffeₚ = ReferenceFE(lagrangian,Float64,order)
Q = TestFESpace(model,reffeₚ,conformity=:H1, constraint=:zeromean, dirichlet_tags=["tag_25", "tag_26"])

mu = 8.9*10^(-4)
rho = 10
nu = mu/rho

u_0 = 10
Re = 10

#u_free = VectorValue(u_0,0,0)
u_free(x) = VectorValue(1.5 * u_0 * x[2] * ( Ly - x[2] ) / ( (Ly/2)^2 ), 0, 0)
u_walls=VectorValue(0,0,0)
U = TrialFESpace(V,[u_walls, u_walls, u_free])
U
P = TrialFESpace(Q,[0, -1/6.14])
#P = TrialFESpace(Q,[1000,-1000])


Y = MultiFieldFESpace([V, Q])
X = MultiFieldFESpace([U, P])

degree = order
Ωₕ = Triangulation(model)
dΩ = Measure(Ωₕ,degree)


conv(u,∇u) =(∇u')⋅u
a((u,p),(v,q)) = ∫( (1/Re)*∇(v)⊙∇(u) - (∇⋅v)*p + q*(∇⋅u) )dΩ
c(u,v) = ∫( v⊙(conv∘(u,∇(u))) )dΩ

res((u,p),(v,q)) = a((u,p),(v,q)) + c(u,v)

op = FEOperator(res,X,Y)

  

nls = NLSolver(
show_trace=true, method=:newton, linesearch=BackTracking())

solver = FESolver(nls)
uh, ph = solve(solver,op)
writevtk(Ωₕ,"ins-results-1",cellfields=["uh"=>uh,"ph"=>ph])
