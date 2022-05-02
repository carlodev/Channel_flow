using Gridap
using GridapDistributed
using PartitionedArrays
using GridapGmsh
using GridapPETSc
using MPI
using LineSearches: BackTracking

function main_ex1(parts)
    #options = "-ksp_type cg -pc_type gamg -ksp_monitor"
    #options = "-ksp_type gmres -pc_type gamg -ksp_monitor"
    #GridapPETSc.with(args=split(options)) do
        𝒯 = CartesianDiscreteModel(parts,(0,1,0,1),(20,20))
        #𝒯 = CartesianDiscreteModel((0,1,0,1),(20,20))

        Ω = Interior(𝒯)
        dΩ = Measure(Ω,2)
        
        # ## FE space
        
        # In this tutorial we will use linear Lagrangian Finite Elements.
        refFE = ReferenceFE(lagrangian,Float64,1)
        
        # The space of test functions is constant in time and is defined in steady problems:
        V = TestFESpace(𝒯,refFE,dirichlet_tags="boundary")
        
        # The trial space is now a `TransientTrialFESpace`, wich is constructed from a `TestFESpace` and a function (or vector of functions) for the Dirichlet boundary condition/s. In that case, the boundary condition function is a time-independent constant, but it could also be a time-dependent field depending on the coordinates $x$ and time $t$.
        function g(x,t::Real)
            0.0
        end
        function g(t::Real)
          x -> g(x,t)
        end
        U =   TransientTrialFESpace(V,g)
        
        # ## Weak form
        
        # The weak form of the problem follows the same structure as other `Gridap` tutorials, where we define the bilinear and linear forms to define the `FEOperator`. In this case we need to deal with time-dependent quantities and with the presence of time derivatives. The former is handled by passing the time, $t$, as an additional argument to the form, i.e. $a(t,u,v)$. The latter is defined using the time derivative operator `∂t`.
        
        # The most general way of constructing a transient FE operator is by using the `TransientFEOperator` function, which receives a residual, a jacobian with respect to the unknown and a jacobian with respect to the time derivative.
        κ(t) = 1.0 + 0.95*sin(2π*t)
        f(t) = sin(π*t)
        res(t,u,v) = ∫( ∂t(u)*v + κ(t)*(∇(u)⋅∇(v)) - f(t)*v )dΩ
        jac(t,u,du,v) = ∫( κ(t)*(∇(du)⋅∇(v)) )dΩ
        jac_t(t,u,duₜ,v) = ∫( duₜ*v )dΩ
        op = TransientFEOperator(res,jac,jac_t,U,V)
        
        # We can also take advantage of automatic differentiation techniques to compute both Jacobians and use the `TransientFEOperator` function sending just the residual.
        op_AD = TransientFEOperator(res,U,V)
        
        # Alternatively, we can exploit the fact that the problem is linear and use the transient Affine FE operator signature `TransientAffineFEOperator`. In that case, we send a form for the mass contribution, $m$, a form for the stiffness contribution, $a$, and the forcing term, $b$.
        m(t,u,v) = ∫( u*v )dΩ
        a(t,u,v) = ∫( κ(t)*(∇(u)⋅∇(v)) )dΩ
        b(t,v) = ∫( f(t)*v )dΩ
        op_Af = TransientAffineFEOperator(m,a,b,U,V)
        
        # ### Alternative FE operator definitions
        
        # For time-dependent problems with constant coefficients, which is not the case of this tutorial, one could use the optimized operator `TransientConstantMatrixFEOperator`, which assumes that the matrix contributions ($m$ and $a$) are time-independent. That is:
        m₀(u,v) = ∫( u*v )dΩ
        a₀(u,v) = ∫( κ(0.0)*(∇(u)⋅∇(v)) )dΩ
        op_CM = TransientConstantMatrixFEOperator(m,a,b,U,V)
        
        # Going further, if we had a problem with constant forcing term, i.e. constant force and constant boundary conditions, we could have used the `TransientConstantFEOperator`. In that case the linear form is also time-independent.
        b₀(v) = ∫( f(0.0)*v )dΩ
        op_C = TransientConstantFEOperator(m,a,b,U,V)
        
        # ## Transient solver
        
        # Once we have the FE operator defined, we proceed with the definition of the transient solver. First, we define a linear solver to be used at each time step. Here we use the `LUSolver`, but other choices are possible.
        linear_solver = LUSolver()
        
        # Then, we define the ODE solver. That is, the scheme that will be used for the time integration. In this tutorial we use the `ThetaMethod` with $\theta = 0.5$, resulting in a 2nd order scheme. The `ThetaMethod` function receives the linear solver, the time step size $\Delta t$ (constant) and the value of $\theta $.
        Δt = 0.05
        θ = 0.5
        ode_solver = ThetaMethod(linear_solver,Δt,θ)
        
        # Finally, we define the solution using the `solve` function, giving the ODE solver, the FE operator, an initial solution, an initial time and a final time. To construct the initial condition we interpolate the initial value (in that case a constant value of 0.0) into the FE space $U(t)$ at $t=0.0$.
        u₀ = interpolate_everywhere(0.0,U(0.0))
        t₀ = 0.0
        T = 10.0
        uₕₜ = solve(ode_solver,op,u₀,t₀,T)
        
        # ## Postprocessing
        
        # We should highlight that `uₕₜ` is just an _iterable_ function and the results at each time steps are only computed when iterating over it, i.e., lazily. We can post-process the results and generate the corresponding `vtk` files using the `createpvd` and `createvtk` functions. The former will create a `.pvd` file with the collection of `.vtu` files saved at each time step by `createvtk`. The computation of the problem solutions will be triggered in the following loop:
        createpvd("poisson_transient_solution") do pvd
          for (uₕ,t) in uₕₜ
            pvd[t] = createvtk(Ω,"poisson_transient_solution_$t"*".vtu",cellfields=["u"=>uₕ])
          end
        end
        
    #end
end

partition = (2,2)
prun(main_ex1, mpi, partition)

#mpiexecjl --project=. -n 4 julia D1_transient.jl