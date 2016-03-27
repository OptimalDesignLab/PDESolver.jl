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
                 must have the signature res_norm = post_func(ctx..., opts),
                 where res_norm is a norm of res_vec
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

  println("\nEntered rk4")
# res_tol is alternative stopping criteria

  # unpack options
  output_freq = opts["output_freq"]::Int
  write_vis = opts["write_vis"]::Bool

  t = 0.0  # timestepper time
  treal = 0.0  # real time (as opposed to pseudo-time)
  t_steps = round(Int, t_max/h)
  println("t_steps: ",t_steps)

  (m,) = size(q_vec)

  f1 = open("convergence.dat", "a+")

  x_old = zeros(q_vec)
  x_old[:] = q_vec
  k1 = zeros(x_old)
  k2 = zeros(x_old)
  k3 = zeros(x_old)
  k4 = zeros(x_old)

  # the algorithm could be implemented without explicitly storing these
  x2 = zeros(x_old)
  x3 = zeros(x_old)
  x4 = zeros(x_old)


  for i=2:(t_steps + 1)

    if i % output_freq == 0
       println("\ntimestep ",i)
    end

    pre_func(ctx..., opts)
    if real_time treal = t end
    f( ctx..., opts, treal)
    sol_norm = post_func(ctx..., opts)
    
    for j=1:m
      k1[j] = res_vec[j]
    end
    # TODO: make this faster
    x2[:] = x_old + (h/2)*k1


     majorIterationCallback(i, ctx..., opts)
   
    if i % 1 == 0
      println(f1, i, " ", sol_norm)
    end
    
    if i % output_freq == 0
      println("flushing convergence.dat to disk")
      flush(f1)
    end



    if (sol_norm < res_tol)
      println("breaking due to res_tol")
      flush(f1)
      break
    end

    # stage 2
    q_vec[:] = x2
    pre_func(ctx..., opts) 
    if real_time  treal = t + h/2 end
    f( ctx..., opts, treal)
    post_func(ctx..., opts)
    for j=1:m
      k2[j] = res_vec[j]
    end

    # TODO: make this faster
    x3[:] = x_old + (h/2)*k2

    # stage 3
    q_vec[:] = x3
    pre_func(ctx..., opts)
    if real_time treal= t + h/2 end
    f( ctx..., opts, treal)
    post_func(ctx..., opts)

    for j=1:m
      k3[j] = res_vec[j]
    end
    x4[:] = x_old + h*k3

    # stage 4
    q_vec[:] = x4
    pre_func(ctx..., opts)
    if real_time treal = t + h end
    f( ctx..., opts, treal)
    post_func(ctx..., opts)
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
    # TODO: make this faster
    x_old[:] = x_old + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
    q_vec[:] = x_old


    fill!(k1, 0.0)
    fill!(k2, 0.0)
    fill!(k3, 0.0)
    fill!(k4, 0.0)


    t = t + h

  end

  close(f1)

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
function pde_post_func(mesh, sbp, eqn, opts)
  eqn.multiplyA0inv(mesh, sbp, eqn, opts, eqn.res)
  eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
  for j=1:length(eqn.res_vec) eqn.res_vec[j] = eqn.Minv[j]*eqn.res_vec[j] end
  return calcNorm(eqn, eqn.res_vec)
end


