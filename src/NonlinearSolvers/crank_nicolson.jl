# crank_nicolson.jl
# Crank-Nicolson implicit solver for PDEs

export crank_nicolson, cnResidual

push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/solver/euler"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/NonlinearSolvers"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Debugging"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Utils"))


@doc """
crank_nicolson

  This function performs Crank-Nicolson implicit time solution, using a function 
  of the form du/dt = f(u, t)
  
  Arguments:
    * f  : function evaluation, must have signature (ctx..., opts, t)
    * h  : time step size
    * t_max: time value to stop time stepping (time starts at 0)
    * q_vec: vector of the u values
    * res_vec: vector of du/dt values (the output of the function f)
    * pre_func: function to to be called after the new u values are put into
                q_vec but before the function f is evaluated.  Must have
                signature: pre_func(ctx..., opts)
    * post_func: function called immediately after f is called.  The function
                 must have the signature res_norm = post_func(ctx..., opts, 
                 calc_norm=true),
                 where res_norm is a norm of res_vec, and calc_norm determines
                 whether or not to calculate the norm.
    * ctx: a tuple (or any iterable container) of the objects needed by
           f, pre_func, and post func.  The tuple is splatted before being
           passed to the functions.
    * opts : options dictionary

    Keyword Arguments:
    * majorIterationCallback: a callback function called after the first
                              stage, useful to do output or logging
    * res_tol : keyword arg, residual topping tolerance
    * real_time : do actual time marching, not pseudo-time marching

   The eqn.q_vec should hold the whichever variables (conservative or
   entropy) that the simulation should use.

   The idea behind pre_func, post func, and ctx is that they abstract what kind
   of system rk4 is timestepping.  rk4 only needs to know about q_vec and
   res_vec.

   For physics modules, ctx should be (mesh, sbp, eqn) and q_vec and res_vec 
   should be eqn.q_vec and eqn.res_vec.
"""->
function crank_nicolson(f::Function, h::AbstractFloat, t_max::AbstractFloat,
                        mesh::AbstractMesh, sbp::AbstractSBP, eqn::AbstractSolutionData,
                        opts, res_tol=-1.0, real_time=true)

  myrank = MPI.Comm_rank(MPI.COMM_WORLD)
  fstdout = BufferedIO(STDOUT)
  if myrank == 0
    println(fstdout, "\nEntered Crank-Nicolson")
    println(fstdout, "res_tol = ", res_tol)
  end

  output_freq = opts["output_freq"]::Int
  write_vis = opts["write_vis"]::Bool
  use_itermax = opts["use_itermax"]::Bool
  if use_itermax
    itermax = opts["itermax"]
  end

  if myrank == 0
    _f1 = open("convergence.dat", "a+")
    f1 = BufferedIO(_f1)
  end

  t = 0.0
  t_steps = round(Int, t_max/h)

  flush(fstdout)
  for i = 2:(t_steps + 1)

    println("CN: entered time-stepping loop")

    # Allow for some kind of stage loop

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

    println("mark1")
    eqn_nextstep = deepcopy(eqn)
    # TODO: don't copy the entire AbstractSolutionData

    # TODO: pre_func & post_func?
#     pre_func(cts..., opt)
#     if real_time treal = t end
#     f( ctx..., opts, treal)
#     sol_norm = post_func(ctx..., opts)

    println("mark2")
    if use_itermax && i > itermax
      if myrank == 0
        println(fstdout, "breaking due to itermax")
        close(f1)
        flush(fstdout)
      end
      break
    end
    println("mark3")

    # pmesh is preconditioning mesh, just passing mesh for arg after opts
#     @time newton(f, cnResidual, mesh, sbp, eqn, opts, mesh, itermax=opts["itermax"])

    # NOTE: Must include a comma in the ctx tuple to indicate tuple
    @time newton(cnResidual, mesh, sbp, eqn, opts, mesh, itermax=opts["itermax"], ctx=(eqn_nextstep,))
    println("mark4")

    # TODO: disassembleSolution?  

  end   # end of t step loop

  #returns t?

end

# TODO: does cnResidual need to have (mesh, sbp, eqn, opts) signature?
function cnResidual(f, mesh::AbstractMesh, sbp::AbstractSBP, eqn::AbstractSolutionData, eqn_nextstep::AbstractSolutionData,
                    opts, t=0.0)

  # u_(n+1) - 0.5*dt* (del dot G_(n+1)) 0 u_n - 0.5*dt* (del dot G_n)
  # u is q_vec (ref: rk4.jl)

  # allocate u_(n+1)

  #q_np1 = eqn.q

  @mpi_master println(fstdout, "entered cnResidual")
  flush(fstdout)

  # evalEuler n+1 args
#   residual = eqn_nextstep.q - 0.5*delta_t*evalEuler(mesh, sbp, eqn_nextstep, opts, t=0.0) - 
#               eqn.q - 0.5*delta_t*evalEuler(mesh, sbp, eqn, opts, t=0.0)

  # TODO:
  #   update only the eqn.q, then eval residual, then replace eqn.q 

  # NOTE: changed evalEuler to f for generic usage
  residual = eqn_nextstep.q - 0.5*delta_t*f(mesh, sbp, eqn_nextstep, opts, t) - 
              eqn.q - 0.5*delta_t*f(mesh, sbp, eqn, opts, t)
  
  return residual


end

function odeConstantResidual(f, mesh::AbstractMesh, sbp::AbstractSBP, eqn::AbstractSolutionData, 
                             eqn_nextstep::AbstractSolutionData, opts, t=0.0)

  return 6

end

@doc """
### NonlinearSolvers.pde_pre_func

  The pre-function for solving partial differential equations with a physics
  module.  The only operation it performs is disassembling eqn.q_vec into
  eqn.q

  Inputs:
    mesh
    sbp
    eqn
    opts
"""->
function pde_pre_func(mesh, sbp, eqn, opts)
  
  eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
end


@doc """
### NonlinearSolvers.pde_post_func

  The post-function for solving partial differential equations with a physics
  module.  This function multiplies by A0inv, assembles eqn.res into
  eqn.res_vec, multiplies by the inverse mass matrix, and calculates
  the SBP approximation to the integral L2 norm

  Inputs:
    mesh
    sbp
    eqn
    opts

"""->
function pde_post_func(mesh, sbp, eqn, opts; calc_norm=true)
  eqn.multiplyA0inv(mesh, sbp, eqn, opts, eqn.res)
  eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
  for j=1:length(eqn.res_vec) eqn.res_vec[j] = eqn.Minv[j]*eqn.res_vec[j] end
  if calc_norm
    local_norm = calcNorm(eqn, eqn.res_vec)
    eqn.params.time.t_allreduce += @elapsed global_norm = MPI.Allreduce(local_norm*local_norm, MPI.SUM, mesh.comm)
    return sqrt(global_norm)
  end

   return nothing
end

