# crank_nicolson.jl
# Crank-Nicolson implicit solver for PDEs

export crank_nicolson

# push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
# push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/solver/euler"))
# push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/NonlinearSolvers"))
# push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Debugging"))
# push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Utils"))


@doc """
crank_nicolson

  This function performs Crank-Nicolson implicit time solution, using a function 
  of the form du/dt = f(u, t)
  

"""->
function crank_nicolson(f::Function, h::AbstractFloat, t_max::AbstractFloat,
                        q_vec::AbstractVector, res_vec::AbstractVector, pre_func,
                        post_func, ctx, opts)

  myrank = MPI.Comm_rank(MPI.COMM_WORLD)
  fstdout = BufferedIO(STDOUT)
  if myrank == 0
    println(fstdout, "\nEntered Crank-Nicolson")
    println(fstdout, "res_tol = ", res_tol)
  end

  flush(fstdout)
  for i = 2:(t_steps + 1)

    # TODO: output_freq
    @mpi_master if i % output_freq == 0
       println(fstdout,"\ntimestep ", i)
       if i % 5*output_freq == 0
         flush(fstdout)
        end
    end

    # MASTER TODO:
    # assign q_i-1 and q_i
    # form R vector, then enter newton
    # Newton somehow needs Jac and r' inv
    # use newton to drive norm R below newton_tol

    # NOTE:
    # majorIterationCallback: called before every step of Newton's method
    # majorIterationCallback(itr, mesh::AbstractMesh, sbp::AbstractSBP, eqn::AbstractEulerData, opts)

    # TODO:
    # All the parenthesis stuff fits in eqn.q?
    # need to build up the func argument to Newton, from evalEuler components

    # TODO: tear Jac alloc out of newton so it doesn't need to be called every time iteration 
    #   (instead: only one alloc at first time step, then future time steps use that alloc)

    @time newton(cnResidual, mesh, sbp, eqn, opts, pmesh, itermax=opts["itermax"])


  end

  #returns t?

end

function cnResidual(mesh::AbstractMesh, sbp::AbstractSBP, eqn::EulerData, opts, t=0.0)

  # u_(n+1) - 0.5*dt* (del dot G_(n+1)) 0 u_n - 0.5*dt* (del dot G_n)
  # u is q_vec (ref: rk4.jl)

  # allocate u_(n+1)

  #q_np1 = eqn.q
  eqn_nextstep = copy(eqn)

  # evalEuler n+1 args
  residual = eqn_nextstep.q - 0.5*dt*evalEuler(mesh, sbp, eqn_nextstep, opts, t=0.0) - 
              eqn.q - 0.5*dt*evalEuler(mesh, sbp, eqn, opts, t=0.0)
  
  return residual


end


