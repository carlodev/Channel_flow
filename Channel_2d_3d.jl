using Gridap
using LineSearches: BackTracking, Static, MoreThuente

"""
2D and 3D channel flow
laminar - vorticity ω extracted
I assume is stable because is 2D
In 3D, same (ν) -> strange results (no parabolic profile for u)
Periodic BC in x-normal faces
For D=3, 3d case, wrong velocity results - instability?
?Import GMSH mesh (there we have better control of the mesh) - but how to set periodic BC in msh file
"""

# Settings
periodic = false # If set to false, will put a uniform velocity u_in at the inlet
u_in = 1.0
ν = 0.0005 # Kinematic vicosity
const D=2; #dimensions number 2 or 3
const N=32; #cells per dimensions
order = 1
u_0 = u_in

include("Channel_Mesh.jl")

model=mesh_channel(;D=D, N=N, printmodel=true, periodic)

@static if D==2
    top = "tag_5"
    bottom = "tag_6"
    inlet = "tag_7"
    outlet = "tag_8"
    outlet_top = "tag_2" # top right corner, not adding corners results in a "jump" at the outlet
    outlet_bottom = "tag_4" # bottom right corner
    inlet_top = "tag_1"
    inlet_bottom = "tag_3"
    u_diri_tags=[top, bottom]
    u_walls=VectorValue(0,0)
    u_in_v = VectorValue(u_in,0)
    u_diri_values = [u_walls, u_walls]
    p_diri_tags=String[]
    p_diri_values=Float64[]
    if !periodic
        append!(u_diri_tags, [inlet, inlet_top, inlet_bottom, outlet_top, outlet_bottom])
        append!(p_diri_tags, [outlet,outlet_top,outlet_bottom])
        append!(u_diri_values, [u_in_v, u_in_v, u_in_v, u_walls, u_walls])
        append!(p_diri_values, [0,0,0])
    end
    
    body_force = periodic ? 0.033 : 0.0
    hf(x)=VectorValue(body_force, 0)
elseif D==3
    u_diri_tags=["tag_24", "tag_23"]
    u_walls=VectorValue(0,0,0)
    hf(x)=VectorValue(0.0033, 0, 0);
end

reffeᵤ = ReferenceFE(lagrangian,VectorValue{D,Float64},order)
V = TestFESpace(model,reffeᵤ,conformity=:H1,dirichlet_tags=u_diri_tags)

reffeₚ = ReferenceFE(lagrangian,Float64,order)
if periodic
    Q = TestFESpace(model, reffeₚ, conformity=:H1, constraint=:zeromean)
else
    Q = TestFESpace(model, reffeₚ, conformity=:H1, dirichlet_tags=p_diri_tags)
end

U = TrialFESpace(V,u_diri_values)
if periodic
    P = TrialFESpace(Q)
else
    P = TrialFESpace(Q,p_diri_values)
end

Y = MultiFieldFESpace([V, Q])
X = MultiFieldFESpace([U, P])

degree = order
Ω = Triangulation(model)
dΩ = Measure(Ω,degree)


# Momentum residual, without the viscous term
Rm(u,p) = u⋅∇(u) + ∇(p) - hf

# Continuity residual
Rc(u) = ∇⋅u

h = lazy_map(h->h^(1/D),get_cell_measure(Ω))

function τ(u,h)
    τ₂ = h^2/(4*ν)
    val(x) = x
    val(x::Gridap.Fields.ForwardDiff.Dual) = x.value
    u = val(norm(u))
    if iszero(u)
        return τ₂
    end
    τ₁ = h/(2*u)
    return 1/(1/τ₁ + 1/τ₂)
end

τb(u,h) = (u⋅u)*τ(u,h)

res((u,p),(v,q)) = ∫(
    ν*∇(v)⊙∇(u) # Viscous term
    + v⊙Rm(u,p) # Other momentum terms
    + q*Rc(u) # Continuity
    + (τ∘(u,h)*(u⋅∇(v) + ∇(q)))⊙Rm(u,p) # First term: SUPG, second term: PSPG
    + τb∘(u,h)*(∇⋅v)⊙Rc(u) # Bulk viscosity. Try commenting out both stabilization terms to see what happens in periodic and non-periodic cases
)dΩ

op = FEOperator(res,X,Y)
nls = NLSolver(show_trace=true, method=:newton, linesearch=MoreThuente(), iterations=30)
solver = FESolver(nls)

uh0 = interpolate_everywhere(VectorValue(u_0,0.0), U)
ph0 = interpolate_everywhere(0.0, P)
xh0 = interpolate_everywhere([uh0, ph0],MultiFieldFESpace([U,P]))
@time (uh, ph), _ = solve!(xh0,solver,op)

println("solve complete")

ω = ∇×uh;


writevtk(Ω,"results-channel-d$D",cellfields=["uh"=>uh,"ph"=>ph, "ω"=>ω])
