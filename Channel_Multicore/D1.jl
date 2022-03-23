using Gridap
using GridapDistributed
using PartitionedArrays
using GridapGmsh
using GridapPETSc
using MPI
using LineSearches: BackTracking

function main_ex1(parts)
    #options = "-ksp_type cg -pc_type gamg -ksp_monitor"
    options = "-ksp_type gmres -pc_type gamg -ksp_monitor"
    GridapPETSc.with(args=split(options)) do
    domain = (0,1,0,1)
    mesh_partition = (4,4)
    model = CartesianDiscreteModel(parts,domain,mesh_partition) #with this line is not workig, in splitting the mesh on the cores; It seems it want the jacobian

    #model = CartesianDiscreteModel(domain,mesh_partition) # with this line on, and the other off it works
    order = 2
    u((x,y)) = (x+y)^order  
    f(x) = -Δ(u,x)
    reffe = ReferenceFE(lagrangian,Float64,order)
    V = TestFESpace(model,reffe,dirichlet_tags="boundary")
    U = TrialFESpace(u,V)
    Ω = Triangulation(model)
    dΩ = Measure(Ω,2*order)
    a(u,v) = ∫( ∇(v)⋅∇(u) )dΩ
    l(v) = ∫( v*f )dΩ
    #op = AffineFEOperator(a,l,U,V)
    res(u,v)=a(u,v)-l(v)

    op = FEOperator(res,U,V)

  

    nls = NLSolver(PETScLinearSolver(); 
    show_trace=true, method=:newton, linesearch=BackTracking()) #It is a non linear solver for a linear problem, it should work anyway, just for test, it should be just non efficient

    solver = FESolver(nls)
    uh = solve(solver, op)
    writevtk(Ω,"results_ex1",cellfields=["uh"=>uh,"grad_uh"=>∇(uh)])
    end
end

partition = (2,2)
nparts = 2;
prun(main_ex1, mpi, partition)