# implementation of 5 stage, 4th order Low Storage Explicit Runge Kutta
# from Carpenter and Kennedy's Fourth-Order 2N-Storage Runge Kutta Schemes
# NASA Technical Memorandum 109112

# this code borrows some infrastructure (pde_pre_func, pde_post_func) from
# the classical rk4 file
# Also borrows the RK$CheckpointData

"""
  This function implements the 5 stage, 4th order Low Storage Explicit Runge Kutta scheme of
  Carpenter and Kennedy

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
function lserk54_ds(f::Function, delta_t::AbstractFloat, t_max::AbstractFloat, 
             q_vec::AbstractVector, res_vec::AbstractVector, pre_func, 
             post_func, ctx, opts, timing::Timings=Timings(); 
             majorIterationCallback=((a...) -> (a...)), 
             res_tol = -1.0, real_time=false)

  # LSERK coefficients
  const a_coeffs = [0; 
                    -567301805773.0/1357537059087.0; 
                    -2404267990393.0/2016746695238.0; 
                    -3550918686646.0/2091501179385.0; 
                    -1275806237668.0/842570457699.0]

  const b_coeffs = [1432997174477.0/9575080441755.0; 
                    5161836677717.0/13612068292357.0; 
                    1720146321549.0/2090206949498.0; 
                    3134564353537.0/4481467310338.0; 
                    2277821191437.0/14882151754819.0]

  const c_coeffs = [0; 
                    1432997174477/9575080441755; 
                    2526269341429/6820363962896; 
                    2006345519317/3224310063776; 
                    2802321613138/2924317926251]

  myrank = MPI.Comm_rank(MPI.COMM_WORLD)
  # MPI.Barrier(MPI.COMM_WORLD)
  if myrank == 0
    println(BSTDOUT, "\nEntered lserk54")
    println(BSTDOUT, "res_tol = ", res_tol)
  end
  # flush(BSTDOUT)
  # MPI.Barrier(MPI.COMM_WORLD)
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
  @mpi_master println(BSTDOUT, "t_steps: ", t_steps)
  @mpi_master println(BSTDOUT, "delta_t = ", delta_t)

  (m,) = size(q_vec)    # TODO: what is this

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
    # Complex128(
    # Tsol(     # TODO

    term2 = zeros(eqn.q)      # TODO: change var name to dJdu
    term2_vec = zeros(Complex128, mesh.numDofPerNode * mesh.numNodesPerElement * mesh.numEl,)     # TODO: Tsol
    @mpi_master println("Direct sensitivity setup done.")     # TODO: BSTDOUT
  end   # end if opts["perturb_Ma"]

  #------------------------------------------------------------------------------
  # capture direct sensitivity at the IC
  # v is the direct sensitivity, du/dM
  # Ma has been perturbed during setup, in types.jl when eqn.params is initialized
  if opts["write_drag"]
    # objective = EulerEquationMod.createObjectiveFunctionalData(mesh, sbp, eqn, opts)
    # objective = EulerEquationMod.createFunctional(mesh, sbp, eqn, opts, 1)    # 1 is the functional num
    objective = createFunctional(mesh, sbp, eqn, opts, 1)    # 1 is the functional num
    drag = real(evalFunctional(mesh, sbp, eqn, opts, objective))
    # note about drag writing: file_dict populated and file opened in src/solver/euler/types.jl
    @mpi_master f_drag = eqn.file_dict[opts["write_drag_fname"]]
    @mpi_master println(f_drag, 1, " ", drag)
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

    @mpi_master f_v_energy = open("v_energy_data.dat", "w")
  end
  if opts["perturb_Ma"]

    # this is the IC, so it gets the first time step's quad_weight
    i = 1       # note that timestep loop below starts at i = 2
    finaliter = calcFinalIter(t_steps, itermax)
    quad_weight = calcQuadWeight(i, delta_t, finaliter)

    @mpi_master f_quadweights = open("quadweights.dat", "w")
    @mpi_master println(f_quadweights, "finaliter: ", finaliter)
    @mpi_master println(f_quadweights, "i    quad_weight")
    @mpi_master println(f_quadweights, i, "    ", quad_weight)

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

  #------------------------------------------------------------------------------
  # Main timestepping loop
  @mpi_master println(BSTDOUT, "---- Ma @ LSERK start: ", eqn.params.Ma, " ----")
  timing.t_timemarch += @elapsed for i=istart:(t_steps + 1)

    if opts["perturb_Ma"]
      quad_weight = calcQuadWeight(i, delta_t, finaliter)
      @mpi_master println(f_quadweights, i, "    ", quad_weight)
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
    # stage 1
    # f(params, u, F_vals, t_i)

    pre_func(ctx..., opts)
    if real_time treal = t end
    timing.t_func += @elapsed f(ctx..., opts, treal)            # evalResidual call
    sol_norm = post_func(ctx..., opts)

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


    if opts["write_L2vnorm"]
      # for visualization of element level DS energy
      fill!(old_q_vec, 0.0)
      for j = 1:length(q_vec)
        old_q_vec[j] = q_vec[j]
      end
    end


    # Stage 1 update
    fac = b_coeffs[1]
    for j=1:length(q_vec)
      dq_vec[j] = delta_t*res_vec[j]
      q_vec[j] += fac*dq_vec[j]
    end

    #------------------------------------------------------------------------------
    # stabilize q_vec: needs to be before q_vec update (NO, it needs to be after, according to Lorenz LSERK)
    if opts["stabilize_v"]        # need to have i != 2??? TODO
      # Stage 1: get B*v
      calcStabilizedQUpdate!(mesh, sbp, eqn, opts, 
                             stab_A, stab_assembler, clipJacData,
                             treal, Bv, tmp_imag)   # q_vec now obtained from eqn.q_vec
      for j=1:length(q_vec)
        q_vec[j] -= fac*delta_t*Bv[j]*im     # needs to be -=, or else v grows exponentially even faster
      end
    end


    # div by (fac*delta_t) ---- put v_energy calc here

    #=
    if opts["stabilize_v"]

      # Stage 1: stabilize q_vec (this only affects the imaginary part of q_vec)
      # The below is doing this:
      #     imag(q) -= fac*delta_t*Bv
      #   or
      #     imag(q) -= fac*delta_t*B*imag(q)
      # TODO: WHY MINUS
      update_tmp = 0.0
      for j = 1:length(q_vec)
        dqimag_vec[j] = delta_t*Bv[j]
        # q_vec[j] = complex(real(q_vec[j]), imag(q_vec[j]) - fac*dqimag_vec[j]) 
        update_tmp += fac*dqimag_vec[j]
        q_vec[j] = complex(real(q_vec[j]), real(imag(q_vec[j]) - fac*dqimag_vec[j]))
      end
      @mpi_master println(f_stabilize_v, "i: $i   stage: 1   sum of stab update:", update_tmp)
    end
    =#

    # for visualization of element level DS energy
    if opts["write_L2vnorm"]
      # special case for calculating v_energy:
      #   Need to calculate v_vec after the first stage for use in the first-stage-only v_energy calculation.
      #   Previously, this would only be done at the end of all the stages
      for v_ix = 1:length(v_vec)
        v_vec[v_ix] = imag(q_vec[v_ix])/Ma_pert_mag         # v_vec alloc'd outside timestep loop
      end
      fill!(R_stab, 0.0)
      for j = 1:length(q_vec)
        R_stab[j] = (q_vec[j] - old_q_vec[j])/(fac*delta_t)
      end
      fill!(v_energy, 0.0)
      for j = 1:length(q_vec)
        v_energy[j] = v_vec[j]*eqn.M[j]*imag(R_stab[j])/Ma_pert_mag
      end

      if (i % output_freq) == 0
        saveSolutionToMesh(mesh, v_energy)
        fname = string("v_energy_", i)
        writeVisFiles(mesh, fname)
      end

      # i  calcNorm    avg   min   max
      v_energy_mean_local = mean(v_energy)
      v_energy_min_local = minimum(v_energy)      # not doing abs on purpose
      v_energy_max_local = maximum(v_energy)

      v_energy_mean = MPI.Allreduce(v_energy_mean_local, MPI.SUM, mesh.comm)
      v_energy_min = MPI.Allreduce(v_energy_min_local, MPI.SUM, mesh.comm)
      v_energy_max = MPI.Allreduce(v_energy_max_local, MPI.SUM, mesh.comm)

      v_energy_norm = calcNorm(eqn, v_energy)

      @mpi_master println(f_v_energy, i, "  ", real(v_energy_norm), "  ", real(v_energy_mean), "  ", real(v_energy_min), "  ", real(v_energy_max))
      if (i % 500) == 0
        @mpi_master flush(f_v_energy)
      end

    end

    #--------------------------------------------------------------------------
    # loop over remaining stages
    for stage=2:5
      pre_func(ctx..., opts) 
      if real_time
        treal = t + c_coeffs[stage]*delta_t 
      end
      timing.t_func += @elapsed f( ctx..., opts, treal)           # evalResidual call
      post_func(ctx..., opts, calc_norm=false)

      # LSERK solution update
      fac = a_coeffs[stage]
      fac2 = b_coeffs[stage]
      for j=1:length(q_vec)
        dq_vec[j] = fac*dq_vec[j] + delta_t*res_vec[j]
        q_vec[j] += fac2*dq_vec[j]
      end

      # call to calcStabilizedQUpdate should be HERE, before the update to q_vec at this stage
      # NO, it needs to be after according to Lorenz LSERK
      if opts["stabilize_v"]        # need to have i != 2??? TODO
        calcStabilizedQUpdate!(mesh, sbp, eqn, opts, 
                              stab_A, stab_assembler, clipJacData,
                              treal, Bv, tmp_imag)   # q_vec now obtained from eqn.q_vec
        for j=1:length(q_vec)
          q_vec[j] -= fac2*delta_t*Bv[j]*im     # needs to be -=, or else v grows exponentially even faster
        end
      end


      #=
      if opts["stabilize_v"]
        # Stages 2-5: get B*v

        # Stages 2-5: stabilize q_vec (this only affects the imaginary part of q_vec)
        # This is doing a -= on the imaginary part of q_vec.
        # The steps for doing this are the same as the full update on q_vec above as part
        #   of LSERK. The difference here is that the update is with (B*v) instead of 
        #   res_vec, and with a separate holding vector for the previous stage update.
        update_tmp = 0.0
        for j = 1:length(q_vec)
          dqimag_vec[j] = fac*dqimag_vec[j] + delta_t*Bv[j]
          # q_vec[j] = complex(real(q_vec[j]), imag(q_vec[j]) - fac2*dqimag_vec[j])

          update_tmp += fac2*dqimag_vec[j]
          q_vec[j] = complex(real(q_vec[j]), real(imag(q_vec[j]) - fac2*dqimag_vec[j]))
        end
        @mpi_master println(f_stabilize_v, "i: $i   stage: $stage   sum of stab update:", update_tmp)
      end
      =#



    end  # end loop over stages
    # -----> after this point, majorIterCallback needs to be called (or at least where drag needs to be written)

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

      # do the dot product of the two terms, and save
      fill!(term2_vec, 0.0)     # not sure this is necessary
      # assembleSolution(mesh, sbp, eqn, opts, term2, term2_vec)      # term2 -> term2_vec
      array3DTo1D(mesh, sbp, eqn, opts, term2, term2_vec)      # term2 -> term2_vec

      for v_ix = 1:length(v_vec)
        # this accumulation occurs across all dofs and all time steps.
        term23 += quad_weight * term2_vec[v_ix] * v_vec[v_ix]
      end

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
  # LSERK end of time step stuff

  # final update
  t += delta_t

  @mpi_master begin
    println("---------------------------------------------")
    println("   LSERK: final time step reached. t = $t")
    println("---------------------------------------------")
  end

  if opts["perturb_Ma"]

    @mpi_master if opts["stabilize_v"]
      close(f_stabilize_v)
    end

    @mpi_master close(f_drag)
    @mpi_master close(f_quadweights)

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
    println(f_dt, delta_t)
    close(f_dt)

    println(" ")
    println(" run parameters that were used:")
    if opts["perturb_Ma"] 
      println("    Ma: ", eqn.params.Ma + Ma_pert)
    else
      println("    Ma: ", eqn.params.Ma)
    end
    println("    aoa: ", eqn.params.aoa)
    println("    delta_t: ", delta_t)
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
    @mpi_master close(f_v_energy)
  end


  if myrank == 0
    close(f1)
  end

  flush(BSTDOUT)


  global STAB_ctr_usmallandnotstabing
  global STAB_ctr_prodgt0andnotstabing
  global STAB_ctr_prodle0andstabing

  # println(" STAB_ctr_usmallandnotstabing: ", STAB_ctr_usmallandnotstabing)
  # println(" STAB_ctr_prodgt0andnotstabing: ", STAB_ctr_prodgt0andnotstabing)
  # println(" STAB_ctr_prodle0andstabing: ", STAB_ctr_prodle0andstabing)

  return t
end  # end lserk54

"""
  See rk4 method with same signature
"""
function lserk54_ds(f::Function, h::AbstractFloat, t_max::AbstractFloat, 
             q_vec::AbstractVector, res_vec::AbstractVector, ctx, opts,
             timing::Timings=Timings(); 
             majorIterationCallback=((a...) -> (a...)), res_tol=-1.0, 
             real_time=false)

  t = lserk54_ds(f::Function, h::AbstractFloat, t_max::AbstractFloat,
              q_vec::AbstractVector, res_vec::AbstractVector,
              pde_pre_func, pde_post_func, ctx, opts, timing; 
              majorIterationCallback=majorIterationCallback, res_tol=res_tol,
              real_time=real_time)

        return t
end




