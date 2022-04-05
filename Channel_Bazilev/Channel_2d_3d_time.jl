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

D = 3; #dimensions number 2 or 3
N = 32; #cells per dimensions

h = 1;
Re_tau = 10

u_0 = 0;

#Time settings
t0 = 0.0
tF = 1
dt = 0.1
θ = 0.5







model = mesh_channel(D=D, N=N, printmodel=false)



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


u_tau = Re_tau * nu / h;


V = TestFESpace(model, reffeᵤ, conformity=:H1, dirichlet_tags=wall_tags)


reffeₚ = ReferenceFE(lagrangian, Float64, order - 1; space=:P)
Q = TestFESpace(model, reffeₚ, conformity=:L2, constraint=:zeromean)


U = TransientTrialFESpace(V, [u_walls, u_walls])

P = TransientTrialFESpace(Q)


Y = MultiFieldFESpace([V, Q])
X = TransientMultiFieldFESpace([U, P])

degree = order
Ω = Triangulation(model)
dΩ = Measure(Ω, degree)

conv(u, ∇u) = (∇u') ⋅ u

rc(u) = ∇⋅u
rm(u,p,t) = ∂t(u) + conv(u, ∇(u)) + ∇(p)  - hf(t)  #- nu*Δ(u)

c(u, v) = ∫(v ⊙ (conv ∘ (u, ∇(u))))dΩ


t_m = (4/dt^2)^(-1/2)
t_c = (t_m)^(-1)
T1(t, (u, p), (v, q)) = ∫((u⋅∇(v)+∇(q))⋅(t_m*rm(u,p,t)))dΩ
T2((u, p), (v, q)) = ∫((∇⋅(v))⋅(t_c*rc(u)))dΩ
T3(t, (u, p), (v, q)) = ∫((u⋅∇(v)')⋅(t_m*rm(u,p,t)))dΩ
T4(t, (u, p), (v, q)) = (t_m^2)*∫((∇(v))⋅(outer(rm(u,p,t),rm(u,p,t))))dΩ



B_G(t, (u, p), (v, q)) = ∫(∂t(u) ⋅ v)dΩ + c(u, v)+ ∫(nu * ∇(v) ⊙ ∇(u) - (∇ ⋅ v) * p + q * (∇ ⋅ u))dΩ
L_MS(t, v) = ∫(v ⋅ hf(t)) * dΩ




#res(t, (u, p), (v, q)) = a(t, (u, p), (v, q)) + c(u, v) + m(∂t(u), v) - ∫(v ⋅ hf(t)) * dΩ
res(t, (u, p), (v, q)) = B_G(t, (u, p), (v, q))+T1(t, (u, p), (v, q)) + T1(t, (u, p), (v, q)) + T2((u, p), (v, q))+T3(t, (u, p), (v, q))-L_MS(t, v)#-T4(t, (u, p), (v, q))
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

