# rk4.jl
# Runge Kutta 4th order solver for ODEs

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

"""
  This type stores all the data RK4 needs to restart.  It is a subtype of
  [`AbstractCheckpointData`](@ref Utils.AbstractCheckpointData).

  **Fields**

   * i: the current time step

"""
type RK4CheckpointData <: AbstractCheckpointData
  i::Int  # current time step
end

"""
  This constructor loads a RK4CheckpointData from the most recently saved
  checkpoint.

  **Inputs**

   * chkpointer: a [`Checkpointer`](@ref Utils.Checkpointer)
   * comm_rank: MPI rank of this process

  **Outputs**

   * chkpoint_data: a RK4CheckpointData object
"""
function RK4CheckpointData(chkpointer::Checkpointer, comm_rank::Integer)

  chkpoint_data = readLastCheckpointData(chkpointer, comm_rank)

  return chkpoint_data::RK4CheckpointData
end

"""
  This function assists in setting up checkpoinging related things for
  self-starting explicit time marching methods

  **Inputs**

   * opts: the options dictionary
   * myrank: MPI rank of this process

  **Outputs*

   * chkpointer: a Checkpointer fully initialized
   * chkpointdata: a RK4CheckpoinntData object, fully initialized
   * skip_checkpoint: a bool indicating if the next checkpoint write should be
                      skipped
"""
function explicit_checkpoint_setup(opts, myrank)
  is_restart = opts["is_restart"]
  ncheckpoints = opts["ncheckpoints"]

  if !is_restart
    # this is a new simulation, create all the stuff needed to checkpoint
    # note that having zero checkpoints is valid
    istart = 2
    chkpointdata = RK4CheckpointData(istart)
    chkpointer = Checkpointer(myrank, ncheckpoints)
    skip_checkpoint = false
  else  # this is a restart, load existing data
    # using default value of 0 checkpoints is ok
    chkpointer = Checkpointer(opts, myrank)
    chkpointdata = RK4CheckpointData(chkpointer, myrank)
    skip_checkpoint = true  # when restarting, don't immediately write a
                            # checkpoint
                            # doing so it not strictly incorrect, but not useful
  end

  return chkpointer, chkpointdata, skip_checkpoint
end



@doc """
rk4

  This function does 4th order Runge Kutta time stepping, using a function of
  the form du/dt = f(u, t)

  Arguments:
    * f  : function evaluation, must have signature (ctx..., opts, t)
    * h  : time step size
    * t_max: time value to stop time stepping (time starts at 0)
    * q_vec: vector of the u values, must be eqn.q_vec
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
   should be eqn.q_vec and eqn.res_vec.  For physics modules, pre_func should
   take the values from q_vec and put them in q, and post_func should take
   the values in res and put them in res_vec.  Thus pre_func and post_func
   provide the link between the way the rk4 represents the data and the 
   way the physics modules represent the data.

   Options Keys

   Implementation Notes
     sol_norm check is only performed in real_time mode
"""->
function rk4(f::Function, h::AbstractFloat, t_max::AbstractFloat, 
             q_vec::AbstractVector, res_vec::AbstractVector, pre_func, 
             post_func, ctx, opts, timing::Timings=Timings(); 
             majorIterationCallback=((a...) -> (a...)), 
             res_tol = -1.0, real_time=false)

  myrank = MPI.Comm_rank(MPI.COMM_WORLD)  #???
  if myrank == 0
    println(BSTDOUT, "\nEntered rk4")
    println(BSTDOUT, "res_tol = ", res_tol)
  end
# res_tol is alternative stopping criteria


  # unpack options
  output_freq = opts["output_freq"]::Int
  use_itermax = opts["use_itermax"]::Bool
  if use_itermax
    itermax = opts["itermax"]
  end

  use_checkpointing = opts["use_checkpointing"]::Bool
  chkpoint_freq = opts["checkpoint_freq"]::Int
  ncheckpoints = opts["ncheckpoints"]::Int

  t = 0.0  # timestepper time
  treal = 0.0  # real time (as opposed to pseudo-time)
  t_steps = round(Int, t_max/h)
  @mpi_master println(BSTDOUT, "t_steps: ",t_steps)
  @mpi_master println(BSTDOUT, "delta_t = ", h)

  (m,) = size(q_vec)

  if myrank == 0
    _f1 = open("convergence.dat", "a")
    f1 = BufferedIO(_f1)
  end

  x_old = copy(q_vec)
  k1 = zeros(x_old)
  k2 = zeros(x_old)
  k3 = zeros(x_old)
  k4 = zeros(x_old)

  # Note: q_vec_old_DEBUG is a tool for showing the change in q between timesteps for comparison with CN (for ex)
#   q_vec_old_DEBUG = zeros(q_vec)

  # setup all the checkpointing related data
  chkpointer, chkpointdata, skip_checkpoint = explicit_checkpoint_setup(opts, myrank)
  istart = chkpointdata.i

  flush(BSTDOUT)
  #-----------------------------------------------------
  ### Main timestepping loop ###
  # beginning of RK4 time stepping loop
  #TODO: make this loop 1-based
  # this loop is 2:(t_steps + 1) when not restarting
  timing.t_timemarch += @elapsed for i=istart:(t_steps + 1)

    # compute time value from time step
    t = (i - 2)*h

#     q_vec_old_DEBUG = deepcopy(q_vec)

    @mpi_master if i % output_freq == 0
       println(BSTDOUT, "\ntimestep ",i)
    end


    if use_checkpointing && i % chkpoint_freq == 0 && !skip_checkpoint
      @mpi_master println(BSTDOUT, "Saving checkpoint at timestep ", i)
      skip_checkpoint = false
      # save all needed variables to the chkpointdata
      chkpointdata.i = i

      if countFreeCheckpoints(chkpointer) == 0
        freeOldestCheckpoint(chkpointer)  # make room for a new checkpoint
      end
      
      # save the checkpoint
      saveNextFreeCheckpoint(chkpointer, ctx..., opts, chkpointdata)
    end

    # flush after all printing
    if i % output_freq == 0
      flush(BSTDOUT)
    end

		# TODO TODO TODO: remove when done debugging viscous par
		(mesh, sbp, eqn) = ctx
    # k3d = zeros(eqn.res)    # DEBUGAA

    # DEBUGAA
    #=
    if mesh.myrank == 0
      for j = 1:8
        println(" ")
      end
      println(" RK4RK4RK4RK4RK4RK4RK4RK4 RK4 top: i = $i")
      for j = 1:4
        println(" ")
      end
    end
    =#

    #=
    if mesh.myrank == 0
      println(">>>>> in RK4, starting stage 1 call <<<<<")
      println(" size(eqn.vecflux_faceL): ", size(eqn.vecflux_faceL))
      println(" size(eqn.vecflux_faceR): ", size(eqn.vecflux_faceR))
      if mesh.commsize > 1
        println(" size(eqn.vecflux_faceL_shared): ", size(eqn.vecflux_faceL_shared))
        println(" size(eqn.vecflux_faceL_shared[1]): ", size(eqn.vecflux_faceL_shared[1]))
      end
		end
    =#


    # stage 1
    pre_func(ctx..., opts)
    if real_time treal = t end
    timing.t_func += @elapsed f( ctx..., opts, treal)
    sol_norm = post_func(ctx..., opts)

    timing.t_callback += @elapsed majorIterationCallback(i, ctx..., opts, BSTDOUT)
    for j=1:m
      k1[j] = res_vec[j]
      q_vec[j] = x_old[j] + (h/2)*k1[j]
    end

    # DEBUGAA
    # disassembleSolution(mesh, sbp, eqn, opts, k3d, k1)
    # print_all_q_res_coords(mesh, eqn, k3d, "after stage 1")
    # print_all_vecflux(mesh, eqn)

    # logging
    @mpi_master if i % 1 == 0
      println(f1, i, " ", sol_norm)
    end
    
    @mpi_master if i % output_freq == 0
      println(BSTDOUT, "flushing convergence.dat to disk")
      flush(f1)
    end

    # check stopping conditions
    if (sol_norm < res_tol) && !real_time
      if myrank == 0
        println(BSTDOUT, "breaking due to res_tol, res norm = $sol_norm")
        close(f1)
        flush(BSTDOUT)
      end
      break
    end

    if use_itermax && i > itermax
      if myrank == 0
        println(BSTDOUT, "breaking due to itermax")
        close(f1)
        flush(BSTDOUT)
      end
      break
    end

    # DEBUGAA
    #=
    if mesh.myrank == 0
      println(">>>>> in RK4, starting stage 2 call <<<<<")
    end
    =#

    # stage 2
    pre_func(ctx..., opts) 
    if real_time  treal = t + h/2 end
    timing.t_func += @elapsed f( ctx..., opts, treal)
    post_func(ctx..., opts, calc_norm=false)
    for j=1:m
      k2[j] = res_vec[j]
      q_vec[j] = x_old[j] + (h/2)*k2[j]
    end

    # DEBUGAA
    # disassembleSolution(mesh, sbp, eqn, opts, k3d, k2)
    # print_all_q_res_coords(mesh, eqn, k3d, "after stage 2")
    # print_all_vecflux(mesh, eqn)

    # DEBUGAA
    #=
    if mesh.myrank == 0
      println(">>>>> in RK4, starting stage 3 call <<<<<")
    end
    =#

    # stage 3
    pre_func(ctx..., opts)
    if real_time treal= t + h/2 end
    timing.t_func += @elapsed f( ctx..., opts, treal)
    post_func(ctx..., opts, calc_norm=false)

    for j=1:m
      k3[j] = res_vec[j]
      q_vec[j] = x_old[j] + h*k3[j]
    end

    # DEBUGAA
    #=
    disassembleSolution(mesh, sbp, eqn, opts, k3d, k3)
    print_all_q_res_coords(mesh, eqn, k3d, "after stage 3")
    print_all_vecflux(mesh, eqn)
    =#

    #=
    if mesh.myrank == 0
      println(">>>>> in RK4, starting stage 4 call <<<<<")
    end
    =#

    # stage 4
    pre_func(ctx..., opts)
    if real_time treal = t + h end
    timing.t_func += @elapsed f( ctx..., opts, treal)
    post_func(ctx..., opts, calc_norm=false)
    for j=1:m
      k4[j] = res_vec[j]
    end

    # DEBUGAA
    #=
    if mesh.myrank == 0
      println(">>>>> in RK4, end <<<<<")
		end
    =#

    # update
    for j=1:m
      x_old[j] = x_old[j] + (h/6)*(k1[j] + 2*k2[j] + 2*k3[j] + k4[j])
      q_vec[j] = x_old[j]
    end

    # DEBUGAA
    #=
    disassembleSolution(mesh, sbp, eqn, opts, k3d, k4)
    print_all_q_res_coords(mesh, eqn, k3d, "after stage 4")
    print_all_vecflux(mesh, eqn)
    =#


    #TODO: is this necessary?
   # fill!(k1, 0.0)
   # fill!(k2, 0.0)
   # fill!(k3, 0.0)
   # fill!(k4, 0.0)

   # t = t + h

  end   # end of RK4 time stepping loop

  t += h  # final time step

  if myrank == 0
    close(f1)
  end

  flush(BSTDOUT)

  # should this be treal?
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
             q_vec::AbstractVector, res_vec::AbstractVector, ctx, opts, timing::Timings=Timings(); 
             majorIterationCallback=((a...) -> (a...)), res_tol=-1.0, 
             real_time=false)

    rk4(f::Function, h::AbstractFloat, t_max::AbstractFloat, q_vec::AbstractVector, 
        res_vec::AbstractVector, pde_pre_func, pde_post_func, ctx, opts; 
        majorIterationCallback=majorIterationCallback, res_tol =res_tol, real_time=real_time)

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

  # DEBUGAA
  #=
  println(eqn.params.f, " >>>>>>>> Start of RK4 loop <<<<<<<<< ")
  print_qvec_coords(mesh, eqn, filename=eqn.params.f)
  =#

  t = rk4(f, h, t_max, eqn.q_vec, eqn.res_vec, pde_pre_func, pde_post_func,
      (mesh, sbp, eqn), opts, eqn.params.time;
      majorIterationCallback=eqn.majorIterationCallback, res_tol=res_tol, real_time=real_time)

  # DEBUGAA
  #=
  println(eqn.params.f, " >>>>>>>> End of RK4 loop <<<<<<<<< ")
  print_qvec_coords(mesh, eqn, filename=eqn.params.f)
  =#

  return t

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
  
  disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
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
  assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
  for j=1:length(eqn.res_vec) eqn.res_vec[j] = eqn.Minv[j]*eqn.res_vec[j] end
  if calc_norm
    local_norm = calcNorm(eqn, eqn.res_vec)
    eqn.params.time.t_allreduce += @elapsed global_norm = MPI.Allreduce(local_norm*local_norm, MPI.SUM, mesh.comm)
    return sqrt(global_norm)
  end

  return nothing
end


#DEBUGGING

function globalNorm(vec)

  local_norm = norm(vec)
  global_norm = MPI.Allreduce(local_norm*local_norm, MPI.SUM, MPI.COMM_WORLD)
  return sqrt(global_norm)
end

# another debugging function
function print_all_q_res_coords(mesh, eqn, k3d, phrase)
  if mesh.commsize == 1
    max_el_num = 8
  elseif mesh.commsize == 2
    max_el_num = 4
  elseif mesh.commsize == 4
    max_el_num = 2
  end
  if mesh.myrank == 0
    println("'''''''''''''''''''''' in RK4: ", phrase, " '''''''''''''''''''''''''")
    for el_ix = 1:max_el_num
      for node_ix = 1:3
        # println(" q[:, $node_ix, $el_ix]: ", eqn.q[:, node_ix, el_ix])
        # println(" res[:, $node_ix, $el_ix]: ", eqn.res[:, node_ix, el_ix])
        # println(" mesh.coords[:, $node_ix, $el_ix]: ", round(mesh.coords[:, node_ix, el_ix], 2))
        println(" k3d[:, $node_ix, $el_ix]: ", k3d[:, node_ix, el_ix])
      end
      println(" ")
    end
  end
end

function print_all_vecflux(mesh, eqn)
  if mesh.myrank == 0
    println(" ")
    println(" (in rk4)")
    for f_ix = 1:size(eqn.vecflux_faceL,4)
      println(" eqn.vecflux_faceL[:, :, :, $f_ix]: ", eqn.vecflux_faceL[:, :, :, f_ix])
    end
    println(" ")
    if mesh.commsize == 1
      for f_ix = 1:size(eqn.vecflux_faceR,4)
        println(" eqn.vecflux_faceR[:, :, :, $f_ix]: ", eqn.vecflux_faceR[:, :, :, f_ix])
      end
    elseif mesh.commsize > 1
      peeridx = 1
      # println(" size(eqn.vecflux_faceL_shared): ", size(eqn.vecflux_faceL_shared))
      # println(" peeridx: ", peeridx)
      # println(" size(eqn.vecflux_faceL_shared[peeridx]): ", size(eqn.vecflux_faceL_shared[peeridx]))
        for f_ix = 1:size(eqn.vecflux_faceL_shared[peeridx], 4)
          println(" eqn.vecflux_faceL_shared[peeridx][:, :, :, $f_ix]: ", eqn.vecflux_faceL_shared[peeridx][:, :, :, f_ix])
        end
        # TODO: print out the mesh coords of these fluxes. use peeridx and f_ix to index mesh.shared_interfaces or shared coords or whatever
  #=
  =#
    end
    println(" ")
  end
end
