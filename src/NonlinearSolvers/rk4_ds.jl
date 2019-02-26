# rk4_ds.jl
# Runge Kutta 4th order solver for ODEs

export rk4_ds

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
function rk4_ds(f::Function, h::AbstractFloat, t_max::AbstractFloat, 
             q_vec::AbstractVector, res_vec::AbstractVector, pre_func, 
             post_func, ctx, opts, timing::Timings=Timings(); 
             majorIterationCallback=((a...) -> (a...)), 
             res_tol = -1.0, real_time=false)

  dt = h      # for clarity, I use this everywhere

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

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Start NEW
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #------------------------------------------------------------------------------
  # splat ctx to get mesh, sbp, eqn
  # (mesh, sbp, eqn) = ctx...
  mesh = ctx[1]   # fastest way to grab mesh from ctx?
  sbp = ctx[2]   # fastest way to grab mesh from ctx?
  eqn = ctx[3]   # fastest way to grab mesh from ctx?

  #------------------------------------------------------------------------------
  # direct sensitivity of Cd wrt M : setup
  if opts["perturb_Ma"]
    term23 = zero(eqn.params.Ma)      # type stable version of 'term23 = 0.0'
    Ma_pert_mag = opts["perturb_Ma_magnitude"]
    Ma_pert = complex(0, Ma_pert_mag)

    term2 = zeros(eqn.q)
    term2_vec = zeros(Complex128, mesh.numDofPerNode * mesh.numNodesPerElement * mesh.numEl,)
    @mpi_master println("Direct sensitivity setup done.")
  end   # end if opts["perturb_Ma"]

  #------------------------------------------------------------------------------
  # capture direct sensitivity at the IC
  # v is the direct sensitivity, du/dM
  # Ma has been perturbed during setup, in types.jl when eqn.params is initialized
  if opts["write_drag"]
    objective = createFunctional(mesh, sbp, eqn, opts, 1)    # 1 is the functional num
    drag = real(evalFunctional(mesh, sbp, eqn, opts, objective))
    # note about drag writing: file_dict populated and file opened in src/solver/euler/types.jl
    @mpi_master f_drag = eqn.file_dict[opts["write_drag_fname"]]
    @mpi_master println(f_drag, 1, " ", drag)
    println("i: ", 1, "  myrank: ", myrank,"  drag: ", drag)
    @mpi_master flush(f_drag)
  end
  if opts["write_L2vnorm"]
    @mpi_master f_L2vnorm = eqn.file_dict[opts["write_L2vnorm_fname"]]
  end

  if opts["write_L2vnorm"]
    # for visualization of element level DS energy
    old_q_vec = zeros(q_vec)
    R_stab = zeros(q_vec)
    v_energy = zeros(q_vec)

    @mpi_master f_v_energy_stage1 = open("v_energy_data_stage1.dat", "w")
    @mpi_master f_v_energy_stageall = open("v_energy_data_stageall.dat", "w")
  end
  if opts["perturb_Ma"]

    # this is the IC, so it gets the first time step's quad_weight
    i = 1       # note that timestep loop below starts at i = 2
    finaliter = calcFinalIter(t_steps, itermax)
    quad_weight = calcQuadWeight(i, dt, finaliter)

    #------------------------------------------------------------------------------
    # allocation of objects for stabilization routine
    if opts["stabilize_v"]

      stab_A = DiagJac(Complex128, mesh.numDofPerNode*mesh.numNodesPerElement, mesh.numEl)
      stab_assembler = AssembleDiagJacData(mesh, sbp, eqn, opts, stab_A)
      clipJacData = ClipJacData(mesh.numDofPerNode*mesh.numNodesPerElement)

      # Bv = zeros(Float64, length(q_vec), )
      Bv = zeros(Complex128, length(q_vec))
      tmp_imag = zeros(Float64, length(eqn.q_vec))
      dqimag_vec = zeros(Bv)

      @mpi_master f_stabilize_v = open("stabilize_v_updates.dat", "w")        # TODO: buffered IO

    end

    # Note: no stabilization of q_vec at the IC
    #   no evalResidual yet, and hasn't entered the time-stepper yet
    #   so no appropriate scaling factors like delta_t or fac or anything

    v_vec = zeros(q_vec)      # direct sensitivity vector
    for v_ix = 1:length(v_vec)
      v_vec[v_ix] = imag(q_vec[v_ix])/Ma_pert_mag
    end
    term2 = zeros(eqn.q)      # First allocation of term2. fill! used below, during timestep loop
    # evalFunctional calls disassembleSolution, which puts q_vec into q
    # should be calling evalFunctional, not calcFunctional.
    evalFunctionalDeriv(mesh, sbp, eqn, opts, objective, term2)    # term2 is func_deriv_arr
    println(" >>>> i: ", i, "  quad_weight: ", quad_weight, "  term2: ", vecnorm(term2), "  v_vec: ", vecnorm(v_vec))

    # do the dot product of the two terms, and save
    # this dot product is: dJdu*dudM
    term2_vec = zeros(Complex128, mesh.numDofPerNode * mesh.numNodesPerElement * mesh.numEl,)
    # assembleSolution(mesh, sbp, eqn, opts, term2, term2_vec)      # term2 -> term2_vec
    array3DTo1D(mesh, sbp, eqn, opts, term2, term2_vec)      # term2 -> term2_vec

    for v_ix = 1:length(v_vec)
      # this accumulation occurs across all dofs and all time steps.
      term23 += quad_weight * term2_vec[v_ix] * v_vec[v_ix]
    end

    if opts["write_L2vnorm"]
      L2_v_norm = calcNorm(eqn, v_vec)
      @mpi_master println(f_L2vnorm, i, "  ", L2_v_norm)
    end

  end   # end if opts["perturb_Ma"]

  flush(BSTDOUT)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # End NEW
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
  # this loop is 2:(t_steps + 1) when not restarting
  timing.t_timemarch += @elapsed for i=istart:(t_steps + 1)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if opts["perturb_Ma"]
      quad_weight = calcQuadWeight(i, dt, finaliter)
    end   # end if opts["perturb_Ma"]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # compute time value from time step
    t = (i - 2)*h

#     q_vec_old_DEBUG = deepcopy(q_vec)

    @mpi_master if i % output_freq == 0
       println(BSTDOUT, "\ntimestep ",i)
    end

    if use_checkpointing && i % chkpoint_freq == 0
      if skip_checkpoint    # skips only the first checkpoint
        skip_checkpoint = false
      else

        @mpi_master println(BSTDOUT, "Saving checkpoint at timestep ", i)
        # save all needed variables to the chkpointdata
        chkpointdata.i = i

        if countFreeCheckpoints(chkpointer) == 0
          freeOldestCheckpoint(chkpointer)  # make room for a new checkpoint
        end

        # save the checkpoint
        saveNextFreeCheckpoint(chkpointer, ctx..., opts, chkpointdata)

      end   # end of if skip_checkpoint check
    end   # end of if use_checkpointing check

    # flush after all printing
    if i % output_freq == 0
      flush(BSTDOUT)
    end

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # save q_vec from the last time step in old_q_vec: used for v_energy calcs
    if opts["write_L2vnorm"]
      # for visualization of element level DS energy
      fill!(old_q_vec, 0.0)
      for j = 1:length(q_vec)
        old_q_vec[j] = q_vec[j]
      end
    end
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FAC = 1.0       # for making stabilization +=    (now that explicit fix has been implemented: no fac needed)
    # FAC = -1.0       # for making stabilization -=  (needed for eig clipping?)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    #------------------------------------------------------------------------------
    # Stage 1
    pre_func(ctx..., opts)
    if real_time treal = t end
    timing.t_func += @elapsed f( ctx..., opts, treal)     # evalResidual, stage 1
    sol_norm = post_func(ctx..., opts)
    timing.t_callback += @elapsed majorIterationCallback(i, ctx..., opts, BSTDOUT)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # stabilize q_vec: needs to be before q_vec update (NO, it needs to be after, according to Lorenz LSERK)
    if opts["stabilize_v"] && i != 2
      # Stage 1: get B*v
      calcStabilizedQUpdate!(mesh, sbp, eqn, opts, 
                            stab_A, stab_assembler, clipJacData,
                            treal, Bv, tmp_imag)   # q_vec now obtained from eqn.q_vec
      # for j=1:length(q_vec)
        # q_vec[j] += FAC*h*Bv[j]*im     # needs to be -=, done w/ FAC
      # end
    end

    for j = 1:m
      k1[j] = h*res_vec[j]
      if opts["stabilize_v"] && i != 2
        k1[j] += FAC*im*dt*Bv[j]      # needs to be -=, done w/ FAC
      end
      q_vec[j] = x_old[j] + 0.5*k1[j]
    end
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # for visualization of element level DS energy
    if opts["write_L2vnorm"]
      # special case for calculating v_energy:
      #   Need to calculate v_vec after the first stage for use in the first-stage-only v_energy calculation.
      # TODO: why only first-stage?
      # TODO: stage parameter divide
      # TODO: store old_q_vec
      #   Previously, this would only be done at the end of all the stages
      for v_ix = 1:length(v_vec)
        v_vec[v_ix] = imag(q_vec[v_ix])/Ma_pert_mag         # v_vec alloc'd outside timestep loop
      end
      fill!(R_stab, 0.0)
      for j = 1:length(q_vec)
        # R_stab[j] = (q_vec[j] - old_q_vec[j])/(fac*delta_t)     # This was LSERK's
        R_stab[j] = (q_vec[j] - old_q_vec[j])/(dt)
      end
      fill!(v_energy, 0.0)
      for j = 1:length(q_vec)
        v_energy[j] = v_vec[j]*eqn.M[j]*imag(R_stab[j])/Ma_pert_mag
      end

      if (i % output_freq) == 0
        saveSolutionToMesh(mesh, v_energy)
        fname = string("v_energy_stage1_", i)
        writeVisFiles(mesh, fname)
      end

      # i  calcNorm    avg   min   max
      v_energy_mean_local = mean(v_energy)
      # v_energy_min_local = minimum(v_energy)      # not doing abs on purpose
      # v_energy_max_local = maximum(v_energy)

      v_energy_mean = MPI.Allreduce(v_energy_mean_local, MPI.SUM, mesh.comm)
      # v_energy_min = MPI.Allreduce(v_energy_min_local, MPI.SUM, mesh.comm)
      # v_energy_max = MPI.Allreduce(v_energy_max_local, MPI.SUM, mesh.comm)

      v_energy_norm = calcNorm(eqn, v_energy)

      # @mpi_master println(f_v_energy, i, "  ", real(v_energy_norm), "  ", real(v_energy_mean), "  ", real(v_energy_min), "  ", real(v_energy_max))
      @mpi_master println(f_v_energy_stage1, i, "  ", real(v_energy_norm))
      if (i % 500) == 0
        @mpi_master flush(f_v_energy_stage1)
      end

    end
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    #------------------------------------------------------------------------------
    # Stage 2: q_vec now contains the updated state after stage 1
    pre_func(ctx..., opts) 
    if real_time  treal = t + h/2 end
    timing.t_func += @elapsed f( ctx..., opts, treal)       # evalResidual, stage 2
    post_func(ctx..., opts, calc_norm=false)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if opts["stabilize_v"] && i != 2
      # Stage 2: get B*v
      calcStabilizedQUpdate!(mesh, sbp, eqn, opts, 
                            stab_A, stab_assembler, clipJacData,
                            treal, Bv, tmp_imag)   # q_vec now obtained from eqn.q_vec
    end

    for j=1:m
      k2[j] = dt*res_vec[j]
      if opts["stabilize_v"] && i != 2
        k2[j] += FAC*im*dt*Bv[j]      # needs to be -=, done w/ FAC
      end
      q_vec[j] = x_old[j] + 0.5*k2[j]
    end
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    #------------------------------------------------------------------------------
    # Stage 3: q_vec now contains the updated state after stage 2
    pre_func(ctx..., opts)
    if real_time treal= t + h/2 end
    timing.t_func += @elapsed f( ctx..., opts, treal)       # evalResidual, stage 3
    post_func(ctx..., opts, calc_norm=false)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if opts["stabilize_v"] && i != 2
      # Stage 3: get B*v
      calcStabilizedQUpdate!(mesh, sbp, eqn, opts, 
                            stab_A, stab_assembler, clipJacData,
                            treal, Bv, tmp_imag)   # q_vec now obtained from eqn.q_vec
    end

    for j=1:m
      k3[j] = dt*res_vec[j]
      if opts["stabilize_v"] && i != 2
        k3[j] += FAC*im*dt*Bv[j]      # needs to be -=, done w/ FAC
      end
      q_vec[j] = x_old[j] + k3[j]
    end

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    #------------------------------------------------------------------------------
    # Stage 4
    pre_func(ctx..., opts)
    if real_time treal = t + h end
    timing.t_func += @elapsed f( ctx..., opts, treal)       # evalResidual, stage 4
    post_func(ctx..., opts, calc_norm=false)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if opts["stabilize_v"] && i != 2
      # Stage 4: get B*v
      calcStabilizedQUpdate!(mesh, sbp, eqn, opts, 
                            stab_A, stab_assembler, clipJacData,
                            treal, Bv, tmp_imag)   # q_vec now obtained from eqn.q_vec
    end

    for j=1:m
      k4[j] = dt*res_vec[j]
      if opts["stabilize_v"] && i != 2
        k4[j] += FAC*im*dt*Bv[j]      # needs to be -=, done w/ FAC
      end
    end
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    # update
    for j=1:m
      # x_old[j] = x_old[j] + (h/6)*(k1[j] + 2*k2[j] + 2*k3[j] + k4[j])
      # Note: factor in front of the k's needs to be 1/6: 
      #       This DS code has factors of h in the k's, so we don't need to multiply by (h/6)
      x_old[j] = x_old[j] + (1/6)*(k1[j] + 2*k2[j] + 2*k3[j] + k4[j])   
      q_vec[j] = x_old[j]
    end


    #-----------------------------------------------------------------------------------------------
    # End of RK4's stages
    #-----------------------------------------------------------------------------------------------

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if opts["write_drag"]
      drag = real(evalFunctional(mesh, sbp, eqn, opts, objective))
      @mpi_master f_drag = eqn.file_dict[opts["write_drag_fname"]]
      @mpi_master println(f_drag, i, " ", drag)
      @mpi_master if (i % opts["output_freq"]) == 0
        flush(f_drag)
      end
    end

    #------------------------------------------------------------------------------
    # direct sensitivity of Cd wrt M : calculation each time step
    if opts["perturb_Ma"]

      # v is the direct sensitivity, du/dM
      # Ma has been perturbed during setup, in types.jl when eqn.params is initialized
      for v_ix = 1:length(v_vec)
        v_vec[v_ix] = imag(q_vec[v_ix])/Ma_pert_mag         # v_vec alloc'd outside timestep loop
      end

      # term2 is the partial deriv of the functional wrt the state: dCd/du
      fill!(term2, 0.0)     # initialized before timestepping loop
      # evalFunctional calls disassembleSolution, which puts q_vec into q
      # should be calling evalFunctional, not calcFunctional.
      #     disassemble isn't getting called. but it shouldn't matter b/c DG
      # EulerEquationMod.evalFunctionalDeriv(mesh, sbp, eqn, opts, objective, term2)    # term2 is func_deriv_arr
      evalFunctionalDeriv(mesh, sbp, eqn, opts, objective, term2)    # term2 is func_deriv_arr
      println(" >>>> i: ", i, "  quad_weight: ", quad_weight, "  term2: ", vecnorm(term2), "  v_vec: ", vecnorm(v_vec))

      # do the dot product of the two terms, and save
      fill!(term2_vec, 0.0)     # not sure this is necessary
      array3DTo1D(mesh, sbp, eqn, opts, term2, term2_vec)      # term2 -> term2_vec

      old_term23 = term23
      for v_ix = 1:length(v_vec)
        # this accumulation occurs across all dofs and all time steps.
        term23 += quad_weight * term2_vec[v_ix] * v_vec[v_ix]
      end
      println("-------------- term23 debugging ---------------")
      println(" i: ", i)
      println(" quad_weight: ", quad_weight)
      println(" vecnorm(term2_vec) (term2 is dDdu in CN): ", vecnorm(term2_vec))
      println(" vecnorm(v_vec): ", vecnorm(v_vec))
      println(" term23: ", term23)
      # println("  term23 change: ", (term23 - old_term23)*1.0/dt)
      println("-----------------------------------------------")

      #------------------------------------------------------------------------------
      # here is where we should be calculating the 'energy' to show that it is increasing over time
      #   'energy' = L2 norm of the solution
      # JEH: So, at each time step, evaluate: sum_{i,j,k} q[i,j,k]*q[i,j,k]*sbp.w[j]*jac[j,k]
      #      (here I assume jac is proportional to the element volume)
      if opts["write_L2vnorm"]
        L2_v_norm = calcNorm(eqn, v_vec)
        @mpi_master println(f_L2vnorm, i, "  ", L2_v_norm)

        if (i % opts["output_freq"]) == 0
          @mpi_master flush(f_L2vnorm)
        end
      end

    end   # end if opts["perturb_Ma"]

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # for visualization of element level DS energy: end of all stages
    if opts["write_L2vnorm"]
      # special case for calculating v_energy:
      #   Need to calculate v_vec after the first stage for use in the first-stage-only v_energy calculation.
      # TODO: why only first-stage?
      # TODO: stage parameter divide
      # TODO: store old_q_vec
      #   Previously, this would only be done at the end of all the stages
      # for v_ix = 1:length(v_vec)
        # v_vec[v_ix] = imag(q_vec[v_ix])/Ma_pert_mag         # v_vec alloc'd outside timestep loop
      # end
      # RK4: the v_vec assignment is commented out because it should be set inside perturb_Ma just above
      fill!(R_stab, 0.0)
      for j = 1:length(q_vec)
        # R_stab[j] = (q_vec[j] - old_q_vec[j])/(fac*delta_t)     # This was LSERK's
        R_stab[j] = (q_vec[j] - old_q_vec[j])/(dt)
      end
      fill!(v_energy, 0.0)
      for j = 1:length(q_vec)
        v_energy[j] = v_vec[j]*eqn.M[j]*imag(R_stab[j])/Ma_pert_mag
      end

      if (i % output_freq) == 0
        saveSolutionToMesh(mesh, v_energy)
        fname = string("v_energy_stageall_", i)
        writeVisFiles(mesh, fname)
      end

      # i  calcNorm    avg   min   max
      v_energy_mean_local = mean(v_energy)
      # v_energy_min_local = minimum(v_energy)      # not doing abs on purpose
      # v_energy_max_local = maximum(v_energy)

      v_energy_mean = MPI.Allreduce(v_energy_mean_local, MPI.SUM, mesh.comm)
      # v_energy_min = MPI.Allreduce(v_energy_min_local, MPI.SUM, mesh.comm)
      # v_energy_max = MPI.Allreduce(v_energy_max_local, MPI.SUM, mesh.comm)

      v_energy_norm = calcNorm(eqn, v_energy)

      # @mpi_master println(f_v_energy, i, "  ", real(v_energy_norm), "  ", real(v_energy_mean), "  ", real(v_energy_min), "  ", real(v_energy_max))
      @mpi_master println(f_v_energy_stageall, i, "  ", real(v_energy_norm))
      if (i % 500) == 0
        @mpi_master flush(f_v_energy_stageall)
      end

    end
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if use_itermax && i > itermax
      if myrank == 0
        println(BSTDOUT, "breaking due to itermax")
        close(f1)
        flush(BSTDOUT)
      end
      break
    end

  end   # end of RK4 time stepping loop

  t += h  # final time step

  if myrank == 0
    close(f1)
  end

  flush(BSTDOUT)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  @mpi_master begin
    println("---------------------------------------------")
    println("   RK4: final time step reached. t = $t")
    println("---------------------------------------------")
  end

  if opts["perturb_Ma"]

    @mpi_master if opts["stabilize_v"]
      close(f_stabilize_v)
    end

    @mpi_master close(f_drag)

    @mpi_master println(" eqn.params.Ma: ", eqn.params.Ma)
    @mpi_master println(" Ma_pert: ", Ma_pert)
    eqn.params.Ma -= Ma_pert      # need to remove perturbation now
    @mpi_master println(" pert removed from Ma")
    @mpi_master println(" eqn.params.Ma: ", eqn.params.Ma)

    finaliter = calcFinalIter(t_steps, itermax)
    Cd, dCddM = calcDragTimeAverage(mesh, sbp, eqn, opts, dt, finaliter)   # will use eqn.params.Ma
    term23 = term23 * 1.0/t     # final step of time average: divide by total time
    global_term23 = MPI.Allreduce(term23, MPI.SUM, mesh.comm)
    total_dCddM = dCddM + global_term23

    # Cd calculations
    @mpi_master begin
      f_total_dCddM = open("total_dCddM.dat", "w")
      println(f_total_dCddM, " Cd: ", Cd)
      println(f_total_dCddM, " dCd/dM: ", dCddM)
      println(f_total_dCddM, " global_term23: ", global_term23)
      println(f_total_dCddM, " total dCd/dM: ", total_dCddM)
      flush(f_total_dCddM)
      close(f_total_dCddM)
      println(" Cd: ", Cd)
      println(" dCd/dM: ", dCddM)
      println(" global_term23: ", global_term23)
      println(" total dCd/dM: ", total_dCddM)
    end

  end   # end if opts["perturb_Ma"]

  @mpi_master begin
    f_Ma = open("Ma.dat", "w")
    println(f_Ma, eqn.params.Ma)
    close(f_Ma)
    f_dt = open("delta_t.dat", "w")
    println(f_dt, dt)
    close(f_dt)

    println(" ")
    println(" run parameters that were used:")
    if opts["perturb_Ma"]
      println("    Ma: ", eqn.params.Ma + Ma_pert)
    else
      println("    Ma: ", eqn.params.Ma)
    end
    println("    aoa: ", eqn.params.aoa)
    println("    dt: ", dt)
    println("    a_inf: ", eqn.params.a_free)
    println("    rho_inf: ", eqn.params.rho_free)
    println("    c: ", 1.0)
    println("    mesh.coord_order: ", mesh.coord_order)
    println(" ")
    println("    opts[stabilization_method]: ", opts["stabilization_method"])
    println("    opts[output_freq]: ", opts["output_freq"])
    println("    opts[use_itermax]: ", opts["use_itermax"])
    println("    opts[itermax]: ", opts["itermax"])
    println("    opts[use_checkpointing]: ", opts["use_checkpointing"])
    println("    opts[checkpoint_freq]: ", opts["checkpoint_freq"])
    println("    opts[ncheckpoints]: ", opts["ncheckpoints"])
    println(" ")
  end

  if opts["write_L2vnorm"]
    @mpi_master close(f_L2vnorm)
    # @mpi_master close(f_v_energy)
    @mpi_master close(f_v_energy_stage1)
    @mpi_master close(f_v_energy_stageall)
  end
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


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
function rk4_ds(f::Function, h::AbstractFloat, t_max::AbstractFloat, 
             q_vec::AbstractVector, res_vec::AbstractVector, ctx, opts, timing::Timings=Timings(); 
             majorIterationCallback=((a...) -> (a...)), res_tol=-1.0, 
             real_time=false)

    rk4_ds(f::Function, h::AbstractFloat, t_max::AbstractFloat, q_vec::AbstractVector, 
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
function rk4_ds(f::Function, h::AbstractFloat, t_max::AbstractFloat, mesh, sbp, eqn, opts; res_tol=-1.0, real_time=false)

  rk4_ds(f, h, t_max, eqn.q_vec, eqn.res_vec, pde_pre_func, pde_post_func,
      (mesh, sbp, eqn), opts, eqn.params.time;
      majorIterationCallback=eqn.majorIterationCallback, res_tol=res_tol, real_time=real_time)

end


