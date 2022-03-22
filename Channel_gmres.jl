using Gridap
using GridapGmsh
using GridapDistributed
using PartitionedArrays
using GridapPETSc
using LineSearches: BackTracking

function main_ex4(parts)
  #options = "-ksp_type cg -pc_type gamg -ksp_monitor"
  #options = "-ksp_type gmres -pc_type gamg -ksp_monitor"
  #options = "-ksp_type gmres  --ksp_gmres_haptol 6.993e+06 -ksp_monitor"
  #options = "-snes_type ngmres -pc_type asm -snes_max_it 100 -snes_monitor"
  #options = "-snes_type ncg -snes_monitor"
  #options = "-ksp_type cg  -ksp_rtol 1.0e-12  -ksp_atol 0.0  -pc_type jacobi"
  #options = "-ksp_type cg  -ksp_rtol 1  -ksp_atol 0.0  -pc_type jacobi"

  options = "-snes_type ngmres 	-snes_ngmres_select_type linesearch -snes_ngmres_restart_it 10000 -snes_ngmres_gammaC 1	-snes_ngmres_epsilonB 1 snes_ngmres_deltaB 1 -snes_ngmres_monitor"
  #options = "-snes_type newtonls -snes_linesearch_type basic  -snes_linesearch_damping 1.0 -snes_rtol 1.0e-14 -snes_atol 0.0 -snes_monitor -pc_type jacobi -ksp_type gmres -snes_converged_reason"
  GridapPETSc.with(args=split(options)) do

  model = GmshDiscreteModel("Channel_geometry.msh")
  #writevtk(model,"model")
  
  labels = get_face_labeling(model);

  D = 2
  order = 2
  print("1\n")

  reffeᵤ = ReferenceFE(lagrangian, VectorValue{2,Float64}, order)
  print("1.1\n")

  V = TestFESpace(model, reffeᵤ, dirichlet_tags=["Inlet", "Bottom_Wall","Top_Wall"])

  #V = TestFESpace(model, reffeᵤ, conformity=:H1, labels=labels,dirichlet_tags=["Inlet", "Bottom_Wall","Top_Wall"])
  print("1.2\n")

  #V = TestFESpace(model,reffeᵤ,conformity=:H1,labels=labels,dirichlet_tags=["Bottom_Wall","Top_Wall"])

  #reffeₚ = ReferenceFE(lagrangian,Float64,order-1;space=:P)
  reffeₚ = ReferenceFE(lagrangian,Float64,order)
  Q = TestFESpace(model,reffeₚ,conformity=:H1, constraint=:zeromean)
  print("2\n")

  mu = 8.9*10^(-4)
  rho = 10
  nu = mu/rho

  u_0 = 10
  Re = 10
  u_free = VectorValue(u_0,0)
  u_walls=VectorValue(0,0)
  U = TrialFESpace(V,[u_free, u_walls, u_walls])
  
  P = TrialFESpace(Q)
  #P = TrialFESpace(Q,[1000,-1000])
  print("3\n")


  Y = MultiFieldFESpace([V, Q])
  X = MultiFieldFESpace([U, P])
  print("4\n")

  degree = order
  Ωₕ = Triangulation(model)
  dΩ = Measure(Ωₕ,degree)


  conv(u,∇u) =(∇u')⋅u
  dconv(du,∇du,u,∇u) = conv(u,∇du)+conv(du,∇u)

  a((u,p),(v,q)) = ∫( (1/Re)*∇(v)⊙∇(u) - (∇⋅v)*p/rho + q*(∇⋅u) )dΩ

  c(u,v) = ∫( v⊙(conv∘(u,∇(u))) )dΩ
  

  res((u,p),(v,q)) = a((u,p),(v,q)) + c(u,v)
  
  print("5\n")

  op = FEOperator(res,X,Y)
  print("6\n")

  
  
  nls = NLSolver(PETScLinearSolver(); show_trace=true, method=:newton, linesearch=BackTracking())

  solver = FESolver(nls)
 
 print("PreSolver ok\n")
"""
 solver = FESolver(PETScNonlinearSolver())
"""
 print("AfterSolver ok\n")

  uh, ph = solve(solver,op)
  writevtk(Ωₕ,"ins-results-1",cellfields=["uh"=>uh,"ph"=>ph])
  end
end

nparts = 1
prun(main_ex4, mpi, nparts)