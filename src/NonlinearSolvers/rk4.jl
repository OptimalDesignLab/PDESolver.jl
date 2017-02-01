# rk4.jl
# Runge Kutta 4th order solver for ODEs

export rk4

global const REVOLVELIB = "librevolvewrap"
global const REVOLVEINFO = 3

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
    * timing: a Timing object, a new one will be created if not provided

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

  For searching: RK4_1
"""->
function rk4(f::Function, h::AbstractFloat, t_start::AbstractFloat, t_end::AbstractFloat,
             q_vec::AbstractVector, res_vec::AbstractVector, pre_func, 
             post_func, ctx, opts, timing::Timings=Timings(); majorIterationCallback=((a...) -> (a...)), 
             res_tol = -1.0, real_time=false)
# TODO cleanup
# TODO revolve structure
# TODO revolve if-else

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

  if myrank == 0 _f1 = open("convergence.dat", "a+")
    f1 = BufferedIO(_f1)
  end


  x_old = copy(q_vec)
  k1 = zeros(x_old)
  k2 = zeros(x_old)
  k3 = zeros(x_old)
  k4 = zeros(x_old)

  # Note: q_vec_old_DEBUG is a tool for showing the change in q between timesteps for comparison with CN (for ex)
#   q_vec_old_DEBUG = zeros(q_vec)

  flush(fstdout)
  #-----------------------------------------------------
  ### Main timestepping loop ###
  # beginning of RK4 time stepping loop

  timing.t_timemarch += @elapsed for i=2:(t_steps + 1)

#     q_vec_old_DEBUG = deepcopy(q_vec)

    @mpi_master if i % output_freq == 0
      println(fstdout, "\ntimestep ",i)
      if i % 5*output_freq == 0
        flush(fstdout)
      end
    end

    pre_func(ctx..., opts)
    if real_time treal = t end
    timing.t_func += @elapsed f( ctx..., opts, treal)
    sol_norm = post_func(ctx..., opts)
    
    timing.t_callback += @elapsed majorIterationCallback(i, ctx..., opts, fstdout)
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
    timing.t_func += @elapsed f( ctx..., opts, treal)
    post_func(ctx..., opts, calc_norm=false)
    for j=1:m
      k2[j] = res_vec[j]
      q_vec[j] = x_old[j] + (h/2)*k2[j]
    end

    # stage 3
    pre_func(ctx..., opts)
    if real_time treal= t + h/2 end
    timing.t_func += @elapsed f( ctx..., opts, treal)
    post_func(ctx..., opts, calc_norm=false)

    for j=1:m
      k3[j] = res_vec[j]
      q_vec[j] = x_old[j] + h*k3[j]
    end

    # stage 4
    pre_func(ctx..., opts)
    if real_time treal = t + h end
    timing.t_func += @elapsed f( ctx..., opts, treal)
    post_func(ctx..., opts, calc_norm=false)
    for j=1:m
      k4[j] = res_vec[j]
    end

#     println("k1 = \n", k1)
#     println("k2 = \n", k2)
#     println("k3 = \n", k3)
#     println("k4 = \n", k4)

#     println("q_old = \n", x_old)

    # update
    for j=1:m
      x_old[j] = x_old[j] + (h/6)*(k1[j] + 2*k2[j] + 2*k3[j] + k4[j])
      q_vec[j] = x_old[j]
    end

#     println("q_vec = \n", q_vec)

    fill!(k1, 0.0)
    fill!(k2, 0.0)
    fill!(k3, 0.0)
    fill!(k4, 0.0)

    t = t + h

  end   # end of RK4 time stepping loop

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

    For searching: RK4_2
"""->
function rk4(f::Function, h::AbstractFloat, t_max::AbstractFloat, 
             q_vec::AbstractVector, res_vec::AbstractVector, ctx, opts, timing::Timings=Timings(); 
             majorIterationCallback=((a...) -> (a...)), res_tol=-1.0, 
             real_time=false)

    rk4(f::Function, h::AbstractFloat, t_max::AbstractFloat, q_vec::AbstractVector, 
        res_vec::AbstractVector, pde_pre_func, pde_post_func, ctx, opts, timing::Timings=Timings(); 
        majorIterationCallback=majorIterationCallback, res_tol =res_tol, real_time=real_time)

end

@doc """
### NonlinearSolvers.rk4

  This method of rk4 is the 'original' rk4, before revolve/adjoint implementation. 
  It passes in only t_max, and sets t_start to equal 0. This preserves readability 
  and reverse compatibility with non-revolve/adjoint rk4 calls.

  For searching: RK4_3
"""->
function rk4(f::Function, h::AbstractFloat, t_max::AbstractFloat, 
             q_vec::AbstractVector, res_vec::AbstractVector, pre_func, 
             post_func, ctx, opts, timing::Timings=Timings(); majorIterationCallback=((a...) -> (a...)), 
             res_tol = -1.0, real_time=false)
             
  t_start = 0.0
  t_end = t_max
  rk4(f::Function, h::AbstractFloat, t_start::AbstractFloat, t_end::AbstractFloat,
      q_vec::AbstractVector, res_vec::AbstractVector, pre_func, 
      post_func, ctx, opts, timing::Timings=Timings(); majorIterationCallback=((a...) -> (a...)), 
      res_tol = -1.0, real_time=false)


end

# TODO TODO
@doc """
### NonlinearSolvers.rk4_revolve

  Document this

  For searching: RK4_4
"""->
function rk4_revolve(f::Function, h::AbstractFloat, t_max::AbstractFloat, 
             q_vec::AbstractVector, res_vec::AbstractVector, pre_func, 
             post_func, ctx, opts, timing::Timings=Timings(); majorIterationCallback=((a...) -> (a...)), 
             res_tol = -1.0, real_time=false)

  revolve_type = 1
  steps = round(Integer, t_max/h)
  snaps = 10
  snaps_in_ram = 4
  info = 3

  RevolveClassPtr = ccall((:revolve_init, REVOLVELIB), Ptr{Void}, (Cint, Cint, Cint, Cint, Cint), 
                          revolve_type, steps, snaps, snaps_in_ram, info)
  
  t = 0.0

  println(" entering rk4_revolve do-while")
  whilestop = false

  while (whilestop == false)

    whatodo = ccall((:revolve, REVOLVELIB), Cint, (Ptr{Void},), RevolveClassPtr)
    # whatodo enum:
    #   1:  ACTION::takeshot
    #   2:  ACTION::advance
    #   3:  ACTION::firsturn
    #   4:  ACTION::youturn
    #   5:  ACTION::restore
    #   98: ACTION::terminate
    #   99: ACTION::error

    # Revolve function outputs/relevant variables
    #   r_getcheck: index of checkpoint
    #   TODO r_getcapo & r_getoldcapo comments

    if (whatodo == 1)   # takeshot
      r_getcheck = ccall((:getcheck, REVOLVELIB), Cint, (Ptr{Void},), RevolveClassPtr)
      r_getcapo = ccall((:getcapo, REVOLVELIB), Cint, (Ptr{Void},), RevolveClassPtr)
      ### store(F, F_Check, t, r_getcheck)
      # store function should store all state vars (i.e. q_vec) & t
      revolve_store(q_vec, t, r_getcheck)
      if (REVOLVEINFO > 1) 
        println(" takeshot at ", r_getcapo, " in CP ", r_getcheck)
      end
    elseif (whatodo == 2)   # advance
      r_getoldcapo = ccall((:getoldcapo, REVOLVELIB), Cint, (Ptr{Void},), RevolveClassPtr)
      r_getcapo = ccall((:getcapo, REVOLVELIB), Cint, (Ptr{Void},), RevolveClassPtr)
      for j = r_getoldcapo:(r_getcapo-1)    # cpp: j=r->getoldcapo();j<r->getcapo();j++
        #### advance(F, F_H, t, h)
        # advance is intended to perform one time integration step
        t_start = t
        t_end = t+h
        rk4(f, h, t_start, t_end, q_vec, res_vec, pre_func, post_func, ctx, opts, 
            timing::Timings=Timings(); majorIterationCallback=((a...) -> (a...)), 
            res_tol = -1.0, real_time=false)
        t = t+h
      end
      if (REVOLVEINFO > 1) 
        println(" advance to ", r_getcapo)
      end
    elseif (whatodo == 3)   # firsturn
      ### advance(F_final, F_H, t, h)
      t_start = t
      t_end = t+h
      rk4(f, h, t_start, t_end, q_vec, res_vec, pre_func, post_func, ctx, opts, 
          timing::Timings=Timings(); majorIterationCallback=((a...) -> (a...)), 
          res_tol = -1.0, real_time=false)

      t = 1.0 - h
      ### t[1] = 1.0 - h        # caution, must assign to t[1], not t, or else t changes to Float64 type
      ### adjoint(L_H, F_H, L, t, h)
      revolve_adjoint(q_vec, t, h)
      r_getcapo = ccall((:getcapo, REVOLVELIB), Cint, (Ptr{Void},), RevolveClassPtr)
      if (REVOLVEINFO > 1) 
        println(" firsturn at ", r_getcapo)
      end
    elseif (whatodo == 4)   # youturn
      ### adjoint(L_H, F, L, t, h)
      revolve_adjoint(q_vec, t, h)
      t = t - h
      r_getcapo = ccall((:getcapo, REVOLVELIB), Cint, (Ptr{Void},), RevolveClassPtr)
      if (REVOLVEINFO > 1) 
        println(" youturn at ", r_getcapo)
      end
    elseif (whatodo == 5)   # restore
      r_getcheck = ccall((:getcheck, REVOLVELIB), Cint, (Ptr{Void},), RevolveClassPtr)
      ### restore(F, F_Check, t, r_getcheck)
      (q_vec, t) = revolve_restore(r_getcheck)
      r_getcapo = ccall((:getcapo, REVOLVELIB), Cint, (Ptr{Void},), RevolveClassPtr)
      if (REVOLVEINFO > 1) 
        println(" restore at ", r_getcapo, " in CP ", r_getcheck)
      end
    elseif (whatodo == 99)  # error
      r_getinfo = ccall((:getinfo, REVOLVELIB), Cint, (Ptr{Void},), RevolveClassPtr)
      if (r_getinfo == 10)
        println(" number of checkpoints stored exceeds checkup,")
        println(" increase constant 'checkup' and recompile")
        error()
      elseif (r_getinfo == 11)
        r_getcheck = ccall((:getcheck, REVOLVELIB), Cint, (Ptr{Void},), RevolveClassPtr)
        println(" number of checkpoints stored = ", r_getcheck+1, " exceeds snaps = ", snaps)
        println(" ensure 'snaps' > 0 and increase initial 'fine'")
        error()
      elseif (r_getinfo == 12)
        println(" error occurs in numforw")
        error()
      elseif (r_getinfo == 13)
        println(" enhancement of 'fine', 'snaps' checkpoints stored, increase 'snaps'")
        error()
      elseif (r_getinfo == 14)
        println(" number of snaps exceeds snapsup, increase constant 'snapsup' and recompile")
        error()
      elseif (r_getinfo == 15)
        println(" number of reps exceeds repsup, increase constant 'repsup' and recompile")
        error()
      end   # end of r_getinfo checking if-block

    end   # end of whatodo checking if-block

    if ((whatodo == 98) || (whatodo == 99))   # 98: terminate, 99: error
      whilestop = true
    end
  end   # end of do-while


  # this is here just as an example of function signature
#   rk4(f::Function, h::AbstractFloat, t_start::AbstractFloat, t_end::AbstractFloat,
#       q_vec::AbstractVector, res_vec::AbstractVector, pre_func, 
#       post_func, ctx, opts, timing::Timings=Timings(); majorIterationCallback=((a...) -> (a...)), 
#       res_tol = -1.0, real_time=false)

end

# TODO put these in RevolveCheckpointing
# TODO comment
function revolve_store(q_vec, t, ix)

  filename = string("revolve_chkpt_", ix, ".dat")
  # TODO file existence check, fail if already exists
  data = vcat(q_vec, t)
  writedlm(filename, data)
    
  # TODO disassemble/assemble?
end

# distant TODO: ram storage

# TODO comment
function revolve_restore(ix)

  filename = string("revolve_chkpt_", ix, ".dat")
  data = readdlm(filename)
  q_vec = data[1:end-1]
  t = data[end]

  return (q_vec, t)

  # TODO disassemble/assemble?
end

# TODO
function revolve_adjoint(eqn, t, ix)

  # need objective

  # need dRdy
  # since linear, inverse of dRdy is transpose of dRdy
  adj = transpose(eqn.res)


  # what is the jacobian produced by newton



end

@doc """
### NonlinearSolvers.rk4


  This is the original (non-general) interface for rk4.
  It is called by call_nlsolver within startup_func.

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

  For searching: RK4_5
"""->
function rk4(f::Function, h::AbstractFloat, t_max::AbstractFloat, mesh, sbp, eqn, opts;
             res_tol=-1.0, real_time=false)

  if opts["revolve"] == true
    rk4_revolve(f, h, t_max, eqn.q_vec, eqn.res_vec, pde_pre_func, pde_post_func,
                (mesh, sbp, eqn), opts, eqn.params.time;
                majorIterationCallback=eqn.majorIterationCallback, res_tol=res_tol, real_time=real_time)
  else
    rk4(f, h, t_max, eqn.q_vec, eqn.res_vec, pde_pre_func, pde_post_func,
        (mesh, sbp, eqn), opts, eqn.params.time;
        majorIterationCallback=eqn.majorIterationCallback, res_tol=res_tol, real_time=real_time)
  end

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


