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

include("../Channel_Mesh.jl")

function main_ex1(part)

    #options = "-ksp_type cg -pc_type gamg -ksp_monitor"
    options = "-ksp_type gmres -pc_type gamg -ksp_monitor"
    #options = "-ksp_type gmres  --ksp_gmres_haptol 6.993e+06 -ksp_monitor"
    #options = "-snes_type ngmres -pc_type asm -snes_max_it 100 -snes_monitor"
    #options = "-snes_type ncg -snes_monitor"
    #options = "-ksp_type cg  -ksp_rtol 1.0e-12  -ksp_atol 0.0  -pc_type jacobi"
    #options = "-ksp_type cg  -ksp_rtol 1  -ksp_atol 0.0  -pc_type jacobi"

    #options = "-snes_type ngmres 	-snes_ngmres_select_type linesearch -snes_ngmres_restart_it 10000 -snes_ngmres_gammaC 1	-snes_ngmres_epsilonB 1 snes_ngmres_deltaB 1 -snes_ngmres_monitor"
    #options = "-snes_type newtonls -snes_linesearch_type basic  -snes_linesearch_damping 1.0 -snes_rtol 1.0e-14 -snes_atol 0.0 -snes_monitor -pc_type jacobi -ksp_type gmres -snes_converged_reason"
    GridapPETSc.with(args=split(options)) do
        D = 2 #dimensions number 2 or 3
        N = 7 #cells per dimensions

        model = mesh_channel(D=D, N=N, printmodel=true)


        #model = CartesianDiscreteModel(domain2d,partition2d, isperiodic=(true,false))

        #model = CartesianDiscreteModel(domain3d,partition3d, isperiodic=(false,false,true))

        #model = GmshDiscreteModel("Channel_geometry.msh")
        #writevtk(model,"model")


        order = 2
        reffeᵤ = ReferenceFE(lagrangian, VectorValue{D,Float64}, order)

        if D == 2
            wall_tags = ["tag_6", "tag_5"]
            u_walls = VectorValue(0, 0)
            nu = 0.0001472
            hf(x) = VectorValue(0.0033, 0)
        elseif D == 3
            wall_tags = ["tag_24", "tag_23"]
            u_walls = VectorValue(0, 0, 0)
            nu = 100
            hf(x) = VectorValue(0.0033, 0, 0)
        end

        V = TestFESpace(model, reffeᵤ, conformity=:H1, dirichlet_tags=wall_tags)


        reffeₚ = ReferenceFE(lagrangian, Float64, order - 1; space=:P)
        Q = TestFESpace(model, reffeₚ, conformity=:L2, constraint=:zeromean)
        h = 1
        Re = 10



        u_t = Re * nu / h
        u_0 = 0
        #u_i(x) = VectorValue(u_0*x[2]/h-(h^2)*(1-(x[2]/h)^2)*G, 0)

        U = TrialFESpace(V, [u_walls, u_walls])

        P = TrialFESpace(Q)


        Y = MultiFieldFESpace([V, Q])
        X = MultiFieldFESpace([U, P])

        degree = order
        Ω = Triangulation(model)
        dΩ = Measure(Ω, degree)



        conv(u, ∇u) = (∇u') ⋅ u
        a((u, p), (v, q)) = ∫(nu * ∇(v) ⊙ ∇(u) - (∇ ⋅ v) * p + q * (∇ ⋅ u))dΩ - ∫(v*hf) * dΩ

        c(u, v) = ∫(v ⊙ (conv ∘ (u, ∇(u))))dΩ

        res((u, p), (v, q)) = a((u, p), (v, q)) + c(u, v)

        op = FEOperator(res, X, Y)



        nls = NLSolver(PETScLinearSolver(); show_trace=true, method=:newton, linesearch=BackTracking())

        solver = FESolver(nls)
        uh, ph = solve(solver, op)

        ω = ∇ × uh


        writevtk(Ω, "results-channel-d$D", cellfields=["uh" => uh, "ph" => ph, "ω" => ω])
    end
end

partition = (1, 1)
prun(main_ex1, mpi, partition)
