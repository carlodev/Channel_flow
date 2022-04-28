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
ν = 0.0001472 # Kinematic vicosity
D=2; #add const, dimensions number 2 or 3
N=32; #add const, cells per dimensions
order = 1
u_0 = u_in

include("Channel_Mesh.jl")

model=mesh_channel(;D=D, N=N, printmodel=false, periodic)
body_force = periodic ? 0.00337204 : 0.0

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
    u_walls=VectorValue(0.0,0.0)
    u_in_v = VectorValue(u_in, 0.0)
    u_diri_values = [u_walls, u_walls]
    p_diri_tags=String[]
    p_diri_values=Float64[]
    if !periodic
        append!(u_diri_tags, [inlet, inlet_top, inlet_bottom, outlet_top, outlet_bottom])
        append!(p_diri_tags, [outlet,outlet_top,outlet_bottom])
        append!(u_diri_values, [u_in_v, u_in_v, u_in_v, u_walls, u_walls])
        append!(p_diri_values, [0,0,0])
    end
    
    
    hf(x)=VectorValue(body_force, 0.0)
    
elseif D==3
    top = ["tag_23", "tag_09", "tag_11"] # face right left
    bottom = ["tag_24", "tag_10", "tag_12"]

    inlet = "tag_05"
    inlet_corners = ["tag_01", "tag_03", "tag_05", "tag_07"] #topleft bottomleft topright bottomright
    inlet_sides = ["tag_13", "tag_15", "tag_17", "tag_19"] #left right top bottom

    outlet = "tag_26"
    outlet_corners = ["tag_02", "tag_04", "tag_06", "tag_08"] #topleft bottomleft topright bottomright
    outlet_sides = ["tag_14", "tag_16", "tag_18", "tag_20"]

    sides = ["tag_21", "tag_22"] # left right

    u_diri_tags=append!(top, bottom)
    u_walls=VectorValue(0, 0, 0)
    u_in_v= VectorValue(u_in, 0, 0)

    u_diri_values = [u_walls, u_walls, u_walls, u_walls, u_walls, u_walls]
    p_diri_tags=String[]
    p_diri_values=Float64[]
    if !periodic
        append!(u_diri_tags, [inlet], inlet_corners, inlet_sides, outlet_corners, outlet_sides)
        append!(p_diri_tags, [outlet], outlet_corners, outlet_sides)
        append!(u_diri_values, [u_in_v, 
                                u_in_v, u_in_v, u_in_v, u_in_v,
                                u_in_v, u_in_v, u_in_v, u_in_v, 
                                u_walls, u_walls, u_walls, u_walls,
                                u_walls, u_walls, u_walls, u_walls])
        append!(p_diri_values, [0,0,0])
    end


    hf(x)=VectorValue(body_force, 0, 0)
end

reffeᵤ = ReferenceFE(lagrangian,VectorValue{D,Float64},order)
V = TestFESpace(model,reffeᵤ,conformity=:H1,dirichlet_tags=u_diri_tags)
reffeₚ = ReferenceFE(lagrangian,Float64, order)

if periodic
    Q = TestFESpace(model, reffeₚ, conformity=:H1, constraint=:zeromean)

else
    
    reffeₚ = ReferenceFE(lagrangian,Float64,order) 
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

degree = 2 
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

var_equations((u,p),(v,q)) = ∫(
    ν*∇(v)⊙∇(u) # Viscous term
    + v⊙Rm(u,p) # Other momentum terms
    + q*Rc(u) )dΩ
    # Continuity
stab_equations((u,p),(v,q)) = ∫((τ∘(u,h)*(u⋅∇(v) + ∇(q)))⊙Rm(u,p) # First term: SUPG, second term: PSPG
    + τb∘(u,h)*(∇⋅v)⊙Rc(u) # Bulk viscosity. Try commenting out both stabilization terms to see what happens in periodic and non-periodic cases
)dΩ

res((u,p),(v,q)) = var_equations((u,p),(v,q)) + stab_equations((u,p),(v,q))
op = FEOperator(res,X,Y)
nls = NLSolver(show_trace=true, method=:newton, linesearch=MoreThuente(), iterations=30)
solver = FESolver(nls)

uh0 = interpolate_everywhere(u_in_v, U)
ph0 = interpolate_everywhere(0.0, P)
xh0 = interpolate_everywhere([uh0, ph0],MultiFieldFESpace([U,P]))
@time (uh, ph), _ = solve!(xh0,solver,op)

println("solve complete")

ω = ∇×uh;


writevtk(Ω,"results-channel-d$D",cellfields=["uh"=>uh,"ph"=>ph, "ω"=>ω])


using Plots

nn = 1000
ev_nodes = LinRange(-0.999, 0.999, nn)
y =1
vt(i) = VectorValue(y,i,0)
ve(v) = v[1] 

vm = evaluate(uh, vt.(ev_nodes))
using Trapz
Um = 0.5 * trapz(ev_nodes, ve.(vm))


plot(ve.(vm), ev_nodes,  seriestype = :scatter)



