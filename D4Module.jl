module D4Module
using Gridap
using GridapDistributed
using PartitionedArrays
using GridapGmsh
using GridapPETSc
using MPI

export main_ex4

  function main_ex4(parts)
    options = "-ksp_type cg -pc_type gamg -ksp_monitor"
    #options = "-ksp_type gmres -pc_type gamg -ksp_monitor"
    GridapPETSc.with(args=split(options)) do
      #model = GmshDiscreteModel(parts, "demo.msh")
      model = GmshDiscreteModel(parts, "demo.msh")
      order = 2
      u((x,y)) = (x+y)^order
      f(x) = -Δ(u,x)
      reffe = ReferenceFE(lagrangian,Float64,order)
      
      V = TestFESpace(model,reffe,dirichlet_tags=["boundary1","boundary2"])
      print("1\n")
      U = TrialFESpace(u,V)
      Ω = Triangulation(model)
      dΩ = Measure(Ω,2*order)
      a(u,v) = ∫( ∇(v)⋅∇(u) )dΩ
      l(v) = ∫( v*f )dΩ
      op = AffineFEOperator(a,l,U,V)
      solver = PETScLinearSolver()
      uh = solve(solver,op)
      writevtk(Ω,"results_ex4",cellfields=["uh"=>uh,"grad_uh"=>∇(uh)])
    end
  end

end


