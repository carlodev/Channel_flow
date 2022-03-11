using Gridap
using GridapGmsh
using GridapDistributed
using PartitionedArrays
using GridapPETSc
using LineSearches: BackTracking


model = GmshDiscreteModel("Channel_geometry.msh")
#writevtk(model,"model")


labels = get_face_labeling(model);

D = 2
order = 2
reffeᵤ = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
V = TestFESpace(model,reffeᵤ,conformity=:H1,labels=labels,dirichlet_tags=["Inlet", "Bottom_Wall","Top_Wall"])
  #V = TestFESpace(model,reffeᵤ,conformity=:H1,labels=labels,dirichlet_tags=["Bottom_Wall","Top_Wall"])

  #reffeₚ = ReferenceFE(lagrangian,Float64,order-1;space=:P)
reffeₚ = ReferenceFE(lagrangian,Float64,order)
Q = TestFESpace(model,reffeₚ,conformity=:H1, constraint=:zeromean)

mu = 8.9*10^(-4)
rho = 10
nu = mu/rho

u_0 = 10
Re = 10

u_free = VectorValue(u_0,0)
u_walls=VectorValue(0,0)
U = TrialFESpace(V,[u_free, u_walls, u_walls])
U
P = TrialFESpace(Q)
#P = TrialFESpace(Q,[1000,-1000])


Y = MultiFieldFESpace([V, Q])
X = MultiFieldFESpace([U, P])

degree = order
Ωₕ = Triangulation(model)
dΩ = Measure(Ωₕ,degree)


conv(u,∇u) =(∇u')⋅u
a((u,p),(v,q)) = ∫( (1/Re)*∇(v)⊙∇(u) - (∇⋅v)*p/rho + q*(∇⋅u) )dΩ
c(u,v) = ∫( v⊙(conv∘(u,∇(u))) )dΩ

res((u,p),(v,q)) = a((u,p),(v,q)) + c(u,v)

op = FEOperator(res,X,Y)

  

nls = NLSolver(
show_trace=true, method=:newton, linesearch=BackTracking())

solver = FESolver(nls)
uh, ph = solve(solver,op)
writevtk(Ωₕ,"ins-results-1",cellfields=["uh"=>uh,"ph"=>ph])
