using Gridap
using LineSearches: BackTracking, Static, MoreThuente
using Statistics

include("MeshChannel.jl")
using .MeshChannel: mesh_channel, h_cell
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
periodic = true # If set to false, will put a uniform velocity u_in at the inlet
u_in = 1.0
ν = 0.0001472 # Kinematic vicosity
Reτ = 395
D = 3; #add const, dimensions number 2 or 3
N = 64; #add const, cells per dimensions
order = 1
u_0 = u_in




model = mesh_channel(; D=D, N=N, printmodel=true, periodic)
body_force = periodic ? 0.00337204 : 0.0

@static if D == 2
    top = "tag_5"
    bottom = "tag_6"
    inlet = "tag_7"
    outlet = "tag_8"
    outlet_top = "tag_2" # top right corner, not adding corners results in a "jump" at the outlet
    outlet_bottom = "tag_4" # bottom right corner
    inlet_top = "tag_1"
    inlet_bottom = "tag_3"
    u_diri_tags = [top, bottom]
    u_walls(x, t::Real) = VectorValue(0.0, 0.0)
    u_in_v(x, t::Real) = VectorValue(u_in, 0.0)
    u_walls(t::Real) = x -> u_walls(x, t)
    u_in_v(t::Real) = x -> u_in_v(x, t)
    p_diri_tags = String[]
    p_diri_values = Float64[]

    if periodic
        u_diri_values = [u_walls, u_walls]
    else

        append!(u_diri_tags, [inlet, inlet_top, inlet_bottom, outlet_top, outlet_bottom])
        append!(p_diri_tags, [outlet, outlet_top, outlet_bottom])
        u_diri_values = [u_walls, u_walls, u_in_v, u_in_v, u_in_v, u_walls, u_walls]
        append!(p_diri_values, [0, 0, 0])
    end


    hf(x, t::Real) = VectorValue(body_force, 0.0)

elseif D == 3
    top = ["tag_23", "tag_09", "tag_11"] # face right left
    bottom = ["tag_24", "tag_10", "tag_12"]

    inlet = "tag_05"
    inlet_corners = ["tag_01", "tag_03", "tag_05", "tag_07"] #topleft bottomleft topright bottomright
    inlet_sides = ["tag_13", "tag_15", "tag_17", "tag_19"] #left right top bottom

    outlet = "tag_26"
    outlet_corners = ["tag_02", "tag_04", "tag_06", "tag_08"] #topleft bottomleft topright bottomright
    outlet_sides = ["tag_14", "tag_16", "tag_18", "tag_20"]

    sides = ["tag_21", "tag_22"] # left right

    u_diri_tags = append!(top, bottom)
    u_walls(x, t) = VectorValue(0, 0, 0)
    u_in_v(x, t::Real) = VectorValue(u_in, 0, 0)
    u_walls(t::Real) = x -> u_walls(x, t)
    u_in_v(t::Real) = x -> u_in_v(x, t)

    p_diri_tags = String[]
    p_diri_values = Float64[]
    if periodic
        u_diri_values = [u_walls, u_walls, u_walls, u_walls, u_walls, u_walls]

    else
        append!(u_diri_tags, [inlet], inlet_corners, inlet_sides, outlet_corners, outlet_sides)
        append!(p_diri_tags, [outlet], outlet_corners, outlet_sides)
        u_diri_values = [u_walls, u_walls, u_walls, u_walls, u_walls, u_walls,
            u_in_v,
            u_in_v, u_in_v, u_in_v, u_in_v,
            u_in_v, u_in_v, u_in_v, u_in_v,
            u_walls, u_walls, u_walls, u_walls,
            u_walls, u_walls, u_walls, u_walls]
        append!(p_diri_values, [0, 0, 0])
    end


    hf(x, t::Real) = VectorValue(body_force, 0, 0)
end


hf(t::Real) = x -> hf(x, t)


reffeᵤ = ReferenceFE(lagrangian, VectorValue{D,Float64}, order)
V = TestFESpace(model, reffeᵤ, conformity=:H1, dirichlet_tags=u_diri_tags)
reffeₚ = ReferenceFE(lagrangian, Float64, order)

if periodic
    Q = TestFESpace(model, reffeₚ, conformity=:H1, constraint=:zeromean)

else

    reffeₚ = ReferenceFE(lagrangian, Float64, order)
    Q = TestFESpace(model, reffeₚ, conformity=:H1, dirichlet_tags=p_diri_tags)
end

U = TransientTrialFESpace(V, u_diri_values)
if periodic
    P = TrialFESpace(Q)
else
    P = TrialFESpace(Q, p_diri_values)
end

Y = MultiFieldFESpace([V, Q])
X = TransientMultiFieldFESpace([U, P])

degree = 2
Ω = Triangulation(model)
dΩ = Measure(Ω, degree)



# Momentum residual, without the viscous term
Rm(t, (u, p)) = ∂t(u) + u ⋅ ∇(u) + ∇(p) - hf(t)

# Continuity residual
Rc(u) = ∇ ⋅ u

h = lazy_map(h -> h^(1 / D), get_cell_measure(Ω))
h  = h_cell(N,D)
function τ(u, h)

    τ₂ = h^2 / (4 * ν)
    val(x) = x
    val(x::Gridap.Fields.ForwardDiff.Dual) = x.value

    u = val(norm(u))
    if iszero(u)
        return τ₂
    end
    τ₁ = h / (2 * u)
    τ₃ = dt / 2
    return 1 / (1 / τ₁ + 1 / τ₂ + 1 / τ₃)

end

τb(u, h) = (u ⋅ u) * τ(u, h)


var_equations(t, (u, p), (v, q)) = ∫(
    ν * ∇(v) ⊙ ∇(u) # Viscous term
    + v ⊙ Rm(t, (u, p)) # Other momentum terms
    + q * Rc(u)
)dΩ # Continuity
stab_equations(t, (u, p), (v, q)) = ∫((τ ∘ (u, h) * (u ⋅ ∇(v) + ∇(q))) ⊙ Rm(t, (u, p)) # First term: SUPG, second term: PSPG
                                      +
                                      τb ∘ (u, h) * (∇ ⋅ v) ⊙ Rc(u) # Bulk viscosity. Try commenting out both stabilization terms to see what happens in periodic and non-periodic cases
)dΩ

res(t, (u, p), (v, q)) = var_equations(t, (u, p), (v, q)) + stab_equations(t, (u, p), (v, q))
op = TransientFEOperator(res, X, Y)

nls = NLSolver(show_trace=true, method=:newton, linesearch=MoreThuente(), iterations=30)
solver = FESolver(nls)

U0 = U(0.0)
P0 = P(0.0)
X0 = X(0.0)

uh0 = interpolate_everywhere(u_in_v(0), U0)
ph0 = interpolate_everywhere(0.0, P0)
xh0 = interpolate_everywhere([uh0, ph0], X0)

t0 = 0.0
dt = 0.03 #Colmes 
Nt = 25000
Ntc = 5000
tF = 100*dt #(Nt + Ntc) * dt

θ = 0.5

ode_solver = ThetaMethod(nls, dt, θ)


sol_t = solve(ode_solver, op, xh0, t0, tF)



#Extract Node in for y=const
model_nodes = DiscreteModel(Polytope{0},model)
y_coord = zeros(N+1)
#2D
if D == 2
    num_nodes_plane = (N + 1)
    planes_nodes = model_nodes.grid.node_coordinates[1:N+1]
    for i = 2:1:N+1
        node_start = (i - 1) * N + i
        node_end = node_start + N
        planey = model_nodes.grid.node_coordinates[node_start:node_end]
        y_coord[i] = model_nodes.grid.node_coordinates[node_start][2]
        planes_nodes = hcat(planes_nodes, planey)
    end
end

#3D
if D == 3
    planes_nodes = model_nodes.grid.node_coordinates[1:N+1]
    num_nodes_plane = (N + 1) * (N + 1)
    y_coord[1] = model_nodes.grid.node_coordinates[1][2]
    for k = 2:1:N+1
        node_start2 = (k - 1) * num_nodes_plane + 1
        node_end2 = node_start2 + N
        planey2 = model_nodes.grid.node_coordinates[node_start2:node_end2]
        planes_nodes = vcat(planes_nodes, planey2)
    end

    for i = 2:1:N+1
        node_start1 = (i - 1) * N + i
        node_end1 = node_start1 + N
        planes_nodes1 = model_nodes.grid.node_coordinates[node_start1:node_end1]

        for k = 2:1:N+1
            node_start2 = (k - 1) * num_nodes_plane + (i - 1) * N + i
            node_end2 = node_start2 + N
            planey2 = model_nodes.grid.node_coordinates[node_start2:node_end2]
            planes_nodes1 = vcat(planes_nodes1, planey2)
        end


        y_coord[i] = model_nodes.grid.node_coordinates[node_start1][2]
        planes_nodes = hcat(planes_nodes, planes_nodes1)
    end

end



#Testing that extraction is correct
planes_nodes
y_coord
planes_nodes[:, 5]
"""
for i = 1:1:N+1
    for j = 2:1:num_nodes_plane
        print("i = $i, j = $j\n")
        @test ≈(planes_nodes[1, i][2], planes_nodes[j, i][2], atol=1 / (N + 1) / 100)

    end
end
"""



V_mean = zeros(Float64, N+ 1, D) #V_mean in each plane
Δc = 1/100 #1 / Ntc
nstep = 0;
_t_nn = t0 
if D==2
    ve(v) = [v[1], v[2]]
elseif D==3
    ve(v) = [v[1], v[2], v[3]]
end

createpvd("Channel3d_td") do pvd

    for (xh_tn, tn) in sol_t
        global _t_nn
        _t_nn += dt
        global nstep
        nstep += 1
        print("iteration = $nstep\n")
        
        uh_tn = xh_tn[1]
        ω = ∇ × uh_tn
        ph_tn = xh_tn[2]
        
        #if mod(nstep, 100) < 1
            pvd[tn] = createvtk(Ω, "Results/Channel_2d_$_t_nn" * ".vtu", cellfields=["uh" => uh_tn, "ω" => ω, "ph" => ph_tn])
        #end
        
        #if nstep > Nt #Channel turbulence devoped
        if D==2
            for p = 2:1:N
                pl_value = evaluate(uh_tn, planes_nodes[:, p])
                M = zeros(Float64, num_nodes_plane, D)
                for i = 1:1:num_nodes_plane
                    M[i, :] = ve(pl_value[i])
                end
                V_mean[p, :] = V_mean[p, :]' + Δc .* mean!([1.0 1.0], M)

            end
        elseif D==3
            for p = 2:1:N
                pl_value = evaluate(uh_tn, planes_nodes[:, p])
                M = zeros(Float64, num_nodes_plane, D)
                for i = 1:1:num_nodes_plane
                    M[i, :] = ve(pl_value[i])
                end
                V_mean[p, :] = V_mean[p, :]' + Δc .* mean!([1.0 1.0 1.0], M)

            end

        end
        #end

    end
end


V_mean
println("solve complete")

using JLD2
@save "V_mean.jld2" V_mean y_coord

