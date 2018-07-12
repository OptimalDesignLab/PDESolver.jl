
# this code borrows much of the infrastructure (pde_pre_func, pde_post_func) from
# LSERK, and therefore the classical rk4 file
# Also borrows the RK$CheckpointData


"""
  This function implements the explicit Euler. 

  Arguments:
    f: a function that evalutes du/dt = f(q, t)
    delta_t: the time step
    t_max: the maximum time value
    q_vec: vector (of length numDof) containing initial solution.  Will contain final solution 
           at exit
    res_vec: vector to store the residual in during evaluation of f.  The contents of this vector
             at exit is undefined
    pre_func: function to call after new values are written into q_vec but before f is called
    post_func: function to call after f is called but before res_vec is accessed
    ctx: tuple arguments of f (ie. f = f(ctx...))
    opts: options dictionary
    timing: a Timings object
    
  Keyword Arguments:
    majorIterationCallback: function to call after first function evaluation of each time step, ie. 
                            when q_vec and res_vec have been updated.  Useful for logging.  Defaults
                            to no-op
    res_tol: stopping tolerance for residual (useful for pseudo-timestepping), default -1.0
    real_time: whether or not to advance time (ie. pseudo timestepping or not) default faulse

  See the documentation for rk4.
"""

function explicit_euler(f::Function, delta_t::AbstractFloat, t_max::AbstractFloat, 
             q_vec::AbstractVector, res_vec::AbstractVector, pre_func, 
             post_func, ctx, opts, timing::Timings=Timings(); 
             majorIterationCallback=((a...) -> (a...)), 
             res_tol = -1.0, real_time=false)
  myrank = MPI.Comm_rank(MPI.COMM_WORLD)
#  MPI.Barrier(MPI.COMM_WORLD)
  if myrank == 0
    println(BSTDOUT, "\nEntered explicit_euler")
    println(BSTDOUT, "res_tol = ", res_tol)
  end
#  flush(BSTDOUT)
#  MPI.Barrier(MPI.COMM_WORLD)
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
  t_steps = round(Int, t_max/delta_t)
  println(BSTDOUT, "t_steps: ",t_steps)
  println(BSTDOUT, "delta_t = ", delta_t)

  (m,) = size(q_vec)

  # allocate storage
  # this is actually a 3N scheme because the function f cannot overwrite its
  # input
  dq_vec = zeros(q_vec)

  if myrank == 0
    _f1 = open("convergence.dat", "a")
    f1 = BufferedIO(_f1)
  end

  # setup all the checkpointing related data
  chkpointer, chkpointdata, skip_checkpoint = explicit_checkpoint_setup(opts, myrank)
  istart = chkpointdata.i

  #------------------------------------------------------------------------------
  # direct sensitivity of Cd wrt M : setup
  if opts["perturb_Ma"]
    term23 = 0.0
    Ma_pert_mag = opts["perturb_Ma_magnitude"]
    Ma_pert = complex(0, Ma_pert_mag)

  end   # end if opts["perturb_Ma"]
  # (mesh, sbp, eqn) = ctx...
  mesh = ctx[1]   # fastest way to grab mesh from ctx?
  sbp = ctx[2]   # fastest way to grab mesh from ctx?
  eqn = ctx[3]   # fastest way to grab mesh from ctx?

  # Initialize quadrature weight for trapezoidal rule
  #   This will be adjusted within the loop for the first & final time steps

  #------------------------------------------------------------------------------
  # capture direct sensitivity at the IC
  # v is the direct sensitivity, du/dM
  # Ma has been perturbed during setup, in types.jl when eqn.params is initialized
  # objective = EulerEquationMod.createObjectiveFunctionalData(mesh, sbp, eqn, opts)
  objective = EulerEquationMod.createFunctional(mesh, sbp, eqn, opts, 1)    # 1 is the functional num
  drag = real(evalFunctional(mesh, sbp, eqn, opts, objective))
  @mpi_master f_drag = eqn.file_dict[opts["write_drag_fname"]]
  @mpi_master println(f_drag, 1, " ", drag)
  @mpi_master flush(f_drag)
  if opts["perturb_Ma"]

    # this is the IC, so it gets the first time step's quad_weight
    i = 1       # note that timestep loop below starts at i = 2
    finaliter = calcFinalIter(t_steps, itermax)
    quad_weight = calcQuadWeight(i, delta_t, finaliter)

    v_vec = zeros(q_vec)      # direct sensitivity vector
    for v_ix = 1:length(v_vec)
      v_vec[v_ix] = imag(q_vec[v_ix])/imag(Ma_pert)
    end
    term2 = zeros(eqn.q)
    # evalFunctional calls disassembleSolution, which puts q_vec into q
    # should be calling evalFunctional, not calcFunctional. disassemble isn't getting called. but it doesn't seem to matter?
    EulerEquationMod.evalFunctionalDeriv(mesh, sbp, eqn, opts, objective, term2)    # term2 is func_deriv_arr

    # do the dot product of the two terms, and save
    term2_vec = zeros(Complex128, mesh.numDofPerNode * mesh.numNodesPerElement * mesh.numEl,)
    assembleSolution(mesh, sbp, eqn, opts, term2, term2_vec)      # term2 -> term2_vec

    new_contrib = 0.0
    for v_ix = 1:length(v_vec)
      # this accumulation occurs across all dofs and all time steps.
      new_contrib = quad_weight * term2_vec[v_ix] * v_vec[v_ix]     
      term23 += new_contrib
    end

  end   # end if opts["perturb_Ma"]

  flush(BSTDOUT)
  #------------------------------------------------------------------------------
  # Main timestepping loop
  finaliter = 0
  println("---- Ma @ EE start: ", eqn.params.Ma, " ----")
  timing.t_timemarch += @elapsed for i=istart:(t_steps + 1)

    if opts["perturb_Ma"]

      finaliter = calcFinalIter(t_steps, itermax)
      quad_weight = calcQuadWeight(i, delta_t, finaliter)

    end   # end if opts["perturb_Ma"]

    t = (i-2)*delta_t

    @mpi_master if i % output_freq == 0
       println(BSTDOUT, "\ntimestep ",i)
       if i % output_freq == 0
         flush(BSTDOUT)
       end
    end

    if use_checkpointing && i % chkpoint_freq == 0
      if skip_checkpoint    # skip only the first checkpoint
        skip_checkpoint = false
      else
        @mpi_master println(BSTDOUT, "Saving checkpoint at timestep ", i)
        skip_checkpoint = false
        # save all needed variables to the chkpointdata
        chkpointdata.i = i

        if countFreeCheckpoints(chkpointer) == 0
          freeOldestCheckpoint(chkpointer)  # make room for a new checkpoint
        end

        # save the checkpoint
        saveNextFreeCheckpoint(chkpointer, ctx..., opts, chkpointdata)
      end   # end of if skip_checkpoint check
    end   # end of if use_checkpointing check


    #--------------------------------------------------------------------------
    # single stage: q_{n+1} = q_{n} + dt*f(q_{n})

    # evalResidual form: f(params, u, F_vals, t_i)

    pre_func(ctx..., opts)
    if real_time treal = t end
    timing.t_func += @elapsed f(ctx..., opts, treal)            # evalResidual call
    sol_norm = post_func(ctx..., opts)

    for j=1:length(q_vec)
      dq_vec[j] = delta_t*res_vec[j]
      q_vec[j] += dq_vec[j]
    end

    #--------------------------------------------------------------------------
    # callback and logging
    timing.t_callback += @elapsed majorIterationCallback(i, ctx..., opts, BSTDOUT) # dirsens note: here is where drag is written


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

    #------------------------------------------------------------------------------
    # direct sensitivity of Cd wrt M : calculation each time step
    if opts["perturb_Ma"]

      # v is the direct sensitivity, du/dM
      # Ma has been perturbed during setup, in types.jl when eqn.params is initialized
      v_vec = zeros(q_vec)      # direct sensitivity vector
      for v_ix = 1:length(v_vec)
        v_vec[v_ix] = imag(q_vec[v_ix])/imag(Ma_pert)
      end

      # term2 is the partial deriv of the functional wrt the state: dCd/du
      term2 = zeros(eqn.q)
      # evalFunctional calls disassembleSolution, which puts q_vec into q
      # should be calling evalFunctional, not calcFunctional. disassemble isn't getting called. but it doesn't seem to matter?
      # objective = EulerEquationMod.createObjectiveFunctionalData(mesh, sbp, eqn, opts)
      objective = EulerEquationMod.createFunctional(mesh, sbp, eqn, opts, 1)    # 1 is the functional num
      EulerEquationMod.evalFunctionalDeriv(mesh, sbp, eqn, opts, objective, term2)    # term2 is func_deriv_arr

      # do the dot product of the two terms, and save
      term2_vec = zeros(Complex128, mesh.numDofPerNode * mesh.numNodesPerElement * mesh.numEl,)
      assembleSolution(mesh, sbp, eqn, opts, term2, term2_vec)      # term2 -> term2_vec

      new_contrib = 0.0
      for v_ix = 1:length(v_vec)
        # this accumulation occurs across all dofs and all time steps.
        new_contrib = quad_weight * term2_vec[v_ix] * v_vec[v_ix]     
        term23 += new_contrib
      end

    end   # end if opts["perturb_Ma"]

    # 201805
    # moved after the q_vec update part of lserk - needed to handle itermax == 1 case. 
    # -------------->>>>> move back after.
    # needs to go after this check: println(BSTDOUT, "breaking due to res_tol, res norm = $sol_norm")
    if use_itermax && i > itermax
      if myrank == 0
        println(BSTDOUT, "breaking due to itermax")
        close(f1)
        flush(BSTDOUT)
      end
      break
    end

  end  # end loop over timesteps


  #------------------------------------------------------------------------------
  # EE end of time step stuff

  # final update
  t += delta_t

  @mpi_master println("------------------------------------------------------------")
  @mpi_master println("   explicit_euler: final time step reached. t = $t")
  @mpi_master println("------------------------------------------------------------")

  if opts["perturb_Ma"]

    @mpi_master println(" eqn.params.Ma: ", eqn.params.Ma)
    @mpi_master println(" Ma_pert: ", Ma_pert)
    eqn.params.Ma -= Ma_pert      # need to remove perturbation now
    @mpi_master println(" pert removed from Ma")
    @mpi_master println(" eqn.params.Ma: ", eqn.params.Ma)

    finaliter = calcFinalIter(t_steps, itermax)
    Cd, dCddM = calcDragTimeAverage(mesh, sbp, eqn, opts, delta_t, finaliter)   # will use eqn.params.Ma
    term23 = term23 * 1.0/t     # final step of time average: divide by total time
    global_term23 = MPI.Allreduce(term23, MPI.SUM, mesh.comm)
    total_dCddM = dCddM + global_term23

    # Cd calculations
    @mpi_master f_total_dCddM = open("total_dCddM.dat", "w")
    @mpi_master println(f_total_dCddM, " dCd/dM: ", dCddM)
    @mpi_master println(f_total_dCddM, " global_term23: ", global_term23)
    @mpi_master println(f_total_dCddM, " total dCd/dM: ", total_dCddM)
    @mpi_master flush(f_total_dCddM)
    @mpi_master close(f_total_dCddM)
    @mpi_master println(" dCd/dM: ", dCddM)
    @mpi_master println(" global_term23: ", global_term23)
    @mpi_master println(" total dCd/dM: ", total_dCddM)

  end   # end if opts["perturb_Ma"]
  @mpi_master f_Ma = open("Ma.dat", "w")
  @mpi_master println(f_Ma, eqn.params.Ma)
  @mpi_master close(f_Ma)
  @mpi_master f_dt = open("delta_t.dat", "w")
  @mpi_master println(f_dt, delta_t)
  @mpi_master close(f_dt)

  @mpi_master println(" ")
  @mpi_master println(" run parameters that were used:")
  if opts["perturb_Ma"] 
    @mpi_master println("    Ma: ", eqn.params.Ma + Ma_pert)
  else
    @mpi_master println("    Ma: ", eqn.params.Ma)
  end
  @mpi_master println("    delta_t: ", delta_t)
  @mpi_master println("    a_inf: ", eqn.params.a_free)
  @mpi_master println("    rho_inf: ", eqn.params.rho_free)
  @mpi_master println("    c: ", 1.0)
  @mpi_master println("    mesh.coord_order: ", mesh.coord_order)
  @mpi_master println(" ")


  if myrank == 0
    close(f1)
  end

  flush(BSTDOUT)

  return t
end  # end explicit_euler

"""
  See rk4 method with same signature
"""
function explicit_euler(f::Function, h::AbstractFloat, t_max::AbstractFloat, 
             q_vec::AbstractVector, res_vec::AbstractVector, ctx, opts,
             timing::Timings=Timings(); 
             majorIterationCallback=((a...) -> (a...)), res_tol=-1.0, 
             real_time=false)

  t = explicit_euler(f::Function, h::AbstractFloat, t_max::AbstractFloat,
              q_vec::AbstractVector, res_vec::AbstractVector,
              pde_pre_func, pde_post_func, ctx, opts, timing; 
              majorIterationCallback=majorIterationCallback, res_tol=res_tol,
              real_time=real_time)

        return t
end
