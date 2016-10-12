# rk4.jl
# Runge Kutta 4th order solver for ODEs
# Anthony Ashley


export rk4

# base RK4 method:
# dxdt = f(t,x)



# Inputs:
#   f:      function, that accepts input: (scalar t, vector x_old, vector x_new)
#   h:      delta t
#   x_ic:   initial condition for x
#   t_max:  length of time to step through
# Outputs:
#   x:      solved x at t_max

@doc """
rk4

  This function does 4th order Runge Kutta time stepping, using a function of
  the form du/dt = f(u, t)

  Arguments:
    * f  : function evaluation, must have signature (ctx..., opts, t), must have signature (ctx..., opts, t)
    * h  : time step size
    * t_max : time value to stop time stepping (time starts at 0)
    * q_vec vector of the u values
    * res_vec: vector of du/dt values (the output of the function f)
    * pre_func: function to to be called after the new u values are put into
                q_vec but before the function f is evaluated.  Mut have
                signature: post_func(ctx..., opts)
    * post_func: function called immediately after f is called.  The function
                 must have the signature res_norm = post_func(ctx..., opts, 
                 calc_norm=true),
                 where res_norm is a norm of res_vec, and calc_norm determins
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
function rk4(f::Function, h::AbstractFloat, t_max::AbstractFloat, 
             q_vec::AbstractVector, res_vec::AbstractVector, pre_func, 
             post_func, ctx, opts; majorIterationCallback=((a...) -> (a...)), 
             res_tol = -1.0, real_time=false)
#function rk4(f::Function, h::AbstractFloat, t_max::AbstractFloat, mesh::AbstractMesh, sbp, eqn::AbstractSolutionData, opts; res_tol = -1.0, real_time=false) 
#function rk4(f, h, x_new, x_ic, t_max, extra_args)

  myrank = MPI.Comm_rank(MPI.COMM_WORLD)
  fstdout = BufferedIO(STDOUT)
  if myrank == 0
    println(fstdout, "\nEntered rk4")
    println(fstdout, "res_tol = ", res_tol)
  end
# res_tol is alternative stopping criteria


  # unpack options
  output_freq = opts["output_freq"]::Int
  write_vis = opts["write_vis"]::Bool
  use_itermax = opts["use_itermax"]::Bool
  if use_itermax
    itermax = opts["itermax"]
  end

  t = 0.0  # timestepper time
  treal = 0.0  # real time (as opposed to pseudo-time)
  t_steps = round(Int, t_max/h)
  println(fstdout, "t_steps: ",t_steps)
  println(fstdout, "delta_t = ", h)

  (m,) = size(q_vec)

  if myrank == 0
    _f1 = open("convergence.dat", "a+")
    f1 = BufferedIO(_f1)
  end


  x_old = copy(q_vec)
  k1 = zeros(x_old)
  k2 = zeros(x_old)
  k3 = zeros(x_old)
  k4 = zeros(x_old)


  flush(fstdout)
  for i=2:(t_steps + 1)

    @mpi_master if i % output_freq == 0
       println(fstdout, "\ntimestep ",i)
       if i % 5*output_freq == 0
         flush(fstdout)
       end
    end

    pre_func(ctx..., opts)
    if real_time treal = t end
    f( ctx..., opts, treal)
    eqn = ctx[3]
    sol_norm = post_func(ctx..., opts)
    
    majorIterationCallback(i, ctx..., opts, fstdout)
    for j=1:m
      k1[j] = res_vec[j]
      q_vec[j] = x_old[j] + (h/2)*k1[j]
    end

   
    @mpi_master if i % 1 == 0
      println(f1, i, " ", sol_norm)
    end
    
    @mpi_master if i % output_freq == 0
      println(fstdout, "flushing convergence.dat to disk")
      flush(f1)
    end

    # check stopping conditions
    if (sol_norm < res_tol)
      if myrank == 0
        println(fstdout, "breaking due to res_tol")
        close(f1)
        flush(fstdout)
      end
      break
    end

    if use_itermax && i > itermax
      if myrank == 0
        println(fstdout, "breaking due to itermax")
        close(f1)
        flush(fstdout)
      end
      break
    end

    # stage 2
    pre_func(ctx..., opts) 
    if real_time  treal = t + h/2 end
    f( ctx..., opts, treal)
    post_func(ctx..., opts, calc_norm=false)
    for j=1:m
      k2[j] = res_vec[j]
      q_vec[j] = x_old[j] + (h/2)*k2[j]
    end

    # stage 3
    pre_func(ctx..., opts)
    if real_time treal= t + h/2 end
    f( ctx..., opts, treal)
    post_func(ctx..., opts, calc_norm=false)

    for j=1:m
      k3[j] = res_vec[j]
      q_vec[j] = x_old[j] + h*k3[j]
    end

    # stage 4
    pre_func(ctx..., opts)
    if real_time treal = t + h end
    f( ctx..., opts, treal)
    post_func(ctx..., opts, calc_norm=false)
    for j=1:m
      k4[j] = res_vec[j]
    end
#=
    println("k1 = \n", k1)
    println("k2 = \n", k2)
    println("k3 = \n", k3)
    println("k4 = \n", k4)

    println("q_old = \n", x_old)
=#
    # update
    for j=1:m
      x_old[j] = x_old[j] + (h/6)*(k1[j] + 2*k2[j] + 2*k3[j] + k4[j])
      q_vec[j] = x_old[j]
    end


    fill!(k1, 0.0)
    fill!(k2, 0.0)
    fill!(k3, 0.0)
    fill!(k4, 0.0)


    t = t + h

  end

  if myrank == 0
    close(f1)
  end

  return t

end

# this is the version for solving PDEs
# it uses the pde_pre_func and pde_post_func below
@doc """
### NonlinearSolvers.rk4

  This method of rk4 takes in the ctx, but not the pre_func and post_func, 
  using pde_pre_func and pde_post_func instead.

  All argument names are the same as for the main rk4 method

  Inputs:
    f: 
    h
    t_max
    q_vec
    res_vec
    ctx
    opts
    majorIterationCallback
    res_tol
    real_time
"""->
function rk4(f::Function, h::AbstractFloat, t_max::AbstractFloat, 
             q_vec::AbstractVector, res_vec::AbstractVector, ctx, opts; 
             majorIterationCallback=((a...) -> (a...)), res_tol=-1.0, 
             real_time=false)
    rk4(f::Function, h::AbstractFloat, t_max::AbstractFloat, q_vec::AbstractVector, res_vec::AbstractVector, pde_pre_func, pde_post_func, ctx, opts; majorIterationCallback=majorIterationCallback, res_tol =res_tol, real_time=real_time)

end

@doc """
### NonlinearSolvers.rk4

  This is the original (non-general) interface for rk4.

  The argument names are the same as for the main rk4 method, unless noted
  otherwise

  Inputs:
    f
    h
    t_max
    mesh:  mesh object
    sbp: sbp object
    eqn: equation object
    opts

  Keyword Arguments:
    res_tol
    real_time

  eqn.q_vec for q_vec, eqn.res_vec for res_vec, pde_pre_func and pde_post func
  for the pre and post functions, eqn.majorIterationCallback for the 
  majorIterationCallback, and (mesh, sbp, eqn) as the ctx
"""->
function rk4(f::Function, h::AbstractFloat, t_max::AbstractFloat, mesh, sbp, eqn, opts; res_tol=-1.0, real_time=false)

  rk4(f, h, t_max, eqn.q_vec, eqn.res_vec, pde_pre_func, pde_post_func, (mesh, sbp, eqn), opts; majorIterationCallback=eqn.majorIterationCallback, res_tol=res_tol, real_time=real_time)
end


@doc """
### NonlinearSolvers.pde_pre_func

  The pre-function for solve partial differential equations with a physics
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

  The post-function for solver partial differential equations with a physics
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


