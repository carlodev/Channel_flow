using Gridap
using GridapGmsh
using GridapDistributed
using PartitionedArrays
using GridapPETSc
using LineSearches: BackTracking
#using Testing
"""
2D and 3D channel flow
laminar - vorticity ω extracted
I assume is stable because is 2D
In 3D, same (nu) -> strange results (no parabolic profile for u)
Periodic BC in x-normal faces
For D=3, 3d case, wrong velocity results - instability?
high nu-> less strong covective term, stable
?Import GMSH mesh (there we have better control of the mesh) - but how to set periodic BC in msh file
"""

include("../Channel_Mesh.jl")

D = 2; #dimensions number 2 or 3
N = 10; #cells per dimensions

model = mesh_channel(D=D, N=N, printmodel=false)


#model = CartesianDiscreteModel(domain2d,partition2d, isperiodic=(true,false))

#model = CartesianDiscreteModel(domain3d,partition3d, isperiodic=(false,false,true))

#model = GmshDiscreteModel("Channel_geometry.msh")
#writevtk(model,"model")


order = 2
reffeᵤ = ReferenceFE(lagrangian, VectorValue{D,Float64}, order)

if D == 2
  wall_tags = ["tag_6", "tag_5"]
  u_walls(x, t::Real) = VectorValue(0, 0)
  nu = 0.0001472
  #nu=100
  hf(x, t::Real) = VectorValue(0.0033, 0)
elseif D == 3
  wall_tags = ["tag_24", "tag_23"]
  u_walls(x, t::Real) = VectorValue(0, 0, 0)
  nu = 0.0001472
  hf(x, t::Real) = VectorValue(0.0033, 0, 0)
end

u_walls(t::Real) = x -> u_walls(x, t)
hf(t::Real) = x -> hf(x, t)

V = TestFESpace(model, reffeᵤ, conformity=:H1, dirichlet_tags=wall_tags)


reffeₚ = ReferenceFE(lagrangian, Float64, order - 1; space=:P)
Q = TestFESpace(model, reffeₚ, conformity=:L2, constraint=:zeromean)
h = 1;
Re = 10



u_t = Re * nu / h;
u_0 = 0;
#u_i(x) = VectorValue(u_0*x[2]/h-(h^2)*(1-(x[2]/h)^2)*G, 0)

U = TransientTrialFESpace(V, [u_walls, u_walls])

P = TransientTrialFESpace(Q)


Y = MultiFieldFESpace([V, Q])
X = TransientMultiFieldFESpace([U, P])

degree = order
Ω = Triangulation(model)
dΩ = Measure(Ω, degree)



conv(u, ∇u) = (∇u') ⋅ u
m(ut, v) = ∫(ut ⋅ v)dΩ

a(t, (u, p), (v, q)) = ∫(nu * ∇(v) ⊙ ∇(u) - (∇ ⋅ v) * p + q * (∇ ⋅ u))dΩ

c(u, v) = ∫(v ⊙ (conv ∘ (u, ∇(u))))dΩ

res(t, (u, p), (v, q)) = a(t, (u, p), (v, q)) + c(u, v) + m(∂t(u), v) - ∫(v ⋅ hf(t)) * dΩ

op = TransientFEOperator(res, X, Y)


U0 = U(0.0)
P0 = P(0.0)
X0 = X(0.0)
if D == 2
  uh0 = interpolate_everywhere(VectorValue(0.0, 0.0), U0)
elseif D == 3
  uh0 = interpolate_everywhere(VectorValue(0.0, 0.0, 00), U0)
end
ph0 = interpolate_everywhere(0.0, P0)
xh0 = interpolate_everywhere([uh0, ph0], X0)


t0 = 0.0
tF = 1
dt = 0.01
θ = 0.5

nls = NLSolver(
  show_trace=true, method=:newton, linesearch=BackTracking())

solver = FESolver(nls)
ode_solver = ThetaMethod(nls, dt, θ)


sol_t = solve(ode_solver, op, xh0, t0, tF)


_t_nn = t0
createpvd("Channel2d_td") do pvd
  for (xh_tn, tn) in sol_t
    global _t_nn
    _t_nn += dt
    uh_tn = xh_tn[1]
    ω = ∇ × uh_tn
    ph_tn = xh_tn[2]
    pvd[tn] = createvtk(Ω, "Results/Channel_2d_$_t_nn" * ".vtu", cellfields=["uh" => uh_tn, "ω" => ω, "ph" => ph_tn])
  end

end