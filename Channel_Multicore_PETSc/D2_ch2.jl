module Channel

using Gridap
using Gridap.FESpaces
using Gridap.ReferenceFEs
using Gridap.Arrays
using Gridap.Geometry
using Gridap.Fields
using Gridap.CellData
using GridapPETSc
using FillArrays
using Test
using InteractiveUtils
using GridapDistributed
using PartitionedArrays

include("MeshChannel.jl")
using .MeshChannel: mesh_channel

periodic = true # If set to false, will put a uniform velocity u_in at the inlet
u_in = 1.0
ν = 0.0001472 # Kinematic vicosity
D=3; #add const, dimensions number 2 or 3
N=32; #add const, cells per dimensions
order = 1
np = 4

function main(N::Int,order::Int,np::Int)
    prun(mpi,(1,np,1)) do parts
      options = "-snes_type newtonls -snes_linesearch_type basic  -snes_linesearch_damping 1.0 -snes_rtol 1.0e-12 -snes_atol 0.0 -snes_monitor -ksp_error_if_not_converged true -ksp_converged_reason -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps"
      GridapPETSc.with(args=split(options)) do
        main_ch(parts,N,order)
      end
    end
  end


function main_ch(parts, N::Int, order::Int)
 
    u_0 = u_in

    model = mesh_channel(; D=D, N=N, parts=parts, printmodel=true, periodic=periodic)
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
    u_walls(x, t::Real)=VectorValue(0.0, 0.0)
    u_in_v(x, t::Real) = VectorValue(u_in, 0)

    u_diri_values = [u_walls, u_walls]
    p_diri_tags=String[]
    p_diri_values=Float64[]
    if !periodic
        append!(u_diri_tags, [inlet, inlet_top, inlet_bottom, outlet_top, outlet_bottom])
        append!(p_diri_tags, [outlet,outlet_top,outlet_bottom])
        append!(u_diri_values, [u_in_v, u_in_v, u_in_v, u_walls, u_walls])
        append!(p_diri_values, [0, 0, 0])
    end
    
    
    hf(x, t::Real)=VectorValue(body_force, 0.0)
    
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
    u_walls(x, t) = VectorValue(0, 0, 0)
    
    u_in_v(x, t::Real) = VectorValue(u_in, 0, 0)

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


    hf(x, t::Real)=VectorValue(body_force, 0, 0)
end

u_walls(t::Real) = x -> u_walls(x, t)
u_in_v(t::Real) = x -> u_in_v(x, t)
hf(t::Real) = x -> hf(x, t)


reffeᵤ = ReferenceFE(lagrangian,VectorValue{D,Float64},order)
V = TestFESpace(model,reffeᵤ,conformity=:H1,dirichlet_tags=u_diri_tags)
reffeₚ = ReferenceFE(lagrangian,Float64, order)

if periodic
    Q = TestFESpace(model, reffeₚ, conformity=:H1, constraint=:zeromean)

else
    
    reffeₚ = ReferenceFE(lagrangian,Float64,order) 
    Q = TestFESpace(model, reffeₚ, conformity=:H1, dirichlet_tags=p_diri_tags)
end

U = TransientTrialFESpace(V,u_diri_values)
if periodic
    P = TrialFESpace(Q)
else
    P = TrialFESpace(Q, p_diri_values)
end

Y = MultiFieldFESpace([V, Q])
X = TransientMultiFieldFESpace([U, P])

degree = 2 
Ω = Triangulation(model)
dΩ = Measure(Ω,degree)



# Momentum residual, without the viscous term
Rm(t,(u,p)) = ∂t(u) + u⋅∇(u) + ∇(p) - hf(t)

# Continuity residual
Rc(u) = ∇⋅u

h = lazy_map(h->h^(1/2),get_cell_measure(Ω))

function τ(u,h)
    
    τ₂ = h^2/(4*ν)
    val(x) = x
    val(x::Gridap.Fields.ForwardDiff.Dual) = x.value

    u = val(norm(u))
    if iszero(u)
        return τ₂
    end
    τ₁ = h/(2*u)
    τ₃ = dt/2
    return 1/(1/τ₁ + 1/τ₂ + 1/τ₃)
  
end

τb(u,h) = (u⋅u)*τ(u,h)


var_equations(t,(u,p),(v,q)) = ∫(
    ν*∇(v)⊙∇(u) # Viscous term
    + v⊙Rm(t,(u,p)) # Other momentum terms
    + q*Rc(u)
 )dΩ # Continuity
stab_equations(t,(u,p),(v,q)) = ∫(  (τ∘(u,h)*(u⋅∇(v) + ∇(q)))⊙Rm(t,(u,p)) # First term: SUPG, second term: PSPG
    + τb∘(u,h)*(∇⋅v)⊙Rc(u) # Bulk viscosity. Try commenting out both stabilization terms to see what happens in periodic and non-periodic cases
)dΩ

res(t,(u,p),(v,q)) = var_equations(t,(u,p),(v,q)) + stab_equations(t,(u,p),(v,q))
op = TransientFEOperator(res,X,Y)

ls = PETScLinearSolver()
nls = NLSolver(ls,show_trace=true,method=:newton,iterations=10)

solver = FESolver(nls)

U0 = U(0.0)
P0 = P(0.0)
X0 = X(0.0)

uh0 = interpolate_everywhere(u_in_v(0), U0)
ph0 = interpolate_everywhere(0.0, P0)
xh0 = interpolate_everywhere([uh0, ph0],X0)
t0 = 0.0
tF = 0.6
dt = 0.2
θ = 0.5

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
end

end