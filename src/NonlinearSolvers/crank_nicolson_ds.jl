# crank_nicolson_ds.jl
# Crank-Nicolson implicit solver for PDEs

export crank_nicolson_ds, cnResidual

#TODO: stop doing this
push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/solver/euler"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/NonlinearSolvers"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Debugging"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Utils"))

@doc """
crank_nicolson_ds

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

   This function supported jacobian/preconditioner freezing with the prefix
   "CN", with the default setting to never recalculate either.  newtonInner
   will use its recalculation policy to recalculate the PC and jacobian.
"""->
function crank_nicolson_ds(f::Function, h::AbstractFloat, t_max::AbstractFloat,
                          mesh::AbstractMesh, sbp::AbstractSBP, eqn::AbstractSolutionData,
                          opts, res_tol=-1.0, real_time=true)

  delta_t = h
  dt = h

  myrank = eqn.myrank
  if myrank == 0
    println(BSTDOUT, "Entered Crank-Nicolson")
    println(BSTDOUT, "res_tol = ", res_tol)
  end
  flush(BSTDOUT)

  output_freq = opts["output_freq"]::Int
  write_vis = opts["write_vis"]::Bool
  use_itermax = opts["use_itermax"]::Bool
  jac_type = opts["jac_type"]
  if use_itermax
    itermax = opts["itermax"]
  end

  use_checkpointing = opts["use_checkpointing"]::Bool
  chkpoint_freq = opts["checkpoint_freq"]::Int
  ncheckpoints = opts["ncheckpoints"]::Int


#  if jac_type == 4
#    throw(ErrorException("CN not implemented for matrix-free ops. (jac_type cannot be 4)"))
#  end

  @mpi_master if myrank == 0
    _f1 = open("convergence.dat", "a")
    f1 = BufferedIO(_f1)
  end

  t = 0.0
  t_steps = round(Int, t_max/h)
  finaliter = calcFinalIter(t_steps, itermax)
  println(BSTDOUT, "finaliter: ", finaliter)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Start direct sensitivity setup    1111
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if opts["perturb_Ma_CN"]
      # term23 = zero(eqn.params.Ma)
      dDdu = zeros(eqn.res)         # dDdu is of size eqn.q or eqn.res (dDdu is term2 in rk4_ds.jl)
      dDdu_vec = zeros(Complex128, mesh.numDofPerNode * mesh.numNodesPerElement * mesh.numEl,)
      v_vec = zeros(Float64, length(eqn.res_vec))   # needs to be Float64, even if res_vec is cplx

      ###### CSR: This section is the only section that contains things that are 
      #           new to the CN version of the complex step method (CSR: 'complex step residual')
      # TODO: make sure we don't need old/new_q_vec_Ma_imag
      # old_q_vec_Maimag = zeros(eqn.q_vec)    # corresponds to eqn
      # new_q_vec_Maimag = zeros(eqn.q_vec)    # corresponds to eqn_nextstep
      res_hat_vec = zeros(eqn.res_vec)    # unsteady residual
      res_hat_vec_cnRhs_CS = zeros(eqn.res_vec)    # unsteady residual
      res_hat_vec_nopert = zeros(eqn.res_vec)    # unsteady residual
      res_hat_vec_pospert = zeros(eqn.res_vec)    # unsteady residual
      # TODO: ensure all of these are being used when debugging complete; otherwise, cleanup

      beforeDS_eqn_q_vec = zeros(eqn.q_vec)
      beforeDS_eqn_nextstep_q_vec = zeros(eqn.q_vec)
      beforeDS_eqn_res_vec = zeros(eqn.res_vec)
      beforeDS_eqn_nextstep_res_vec = zeros(eqn.res_vec)

      dRdM_vec = zeros(eqn.res_vec)
      dRdM_vec_FD = zeros(eqn.res_vec)
      b_vec = zeros(Float64, length(eqn.res_vec))   # needs to be Float64, even if res_vec is cplx
      TEST_b_vec = zeros(b_vec)
      TEST_dRdq_vn_prod = zeros(b_vec)
      dRdq_vn_prod = zeros(eqn.res_vec)

      #------------------------------------------------------------------------------
      # DS linear solver
      pc_ds, lo_ds = getCNDSPCandLO(mesh, sbp, eqn, opts)
      ls_ds = StandardLinearSolver(pc_ds, lo_ds, eqn.comm, opts)
      println(" typeof(pc_ds): ", typeof(pc_ds))
      println(" typeof(lo_ds): ", typeof(lo_ds))
      println(" typeof(ls_ds): ", typeof(ls_ds))

      # these are needed for direct manipulation of the Petsc Jac later
      lo_ds_innermost = getBaseLO(lo_ds)     # makes sure to return the innermost LO object (basically returning the Jac)
      # not needed?
      println(" typeof(lo_ds_innermost): ", typeof(lo_ds_innermost))    # returns LinearSolvers.PetscMatLO
      println(" typeof(lo_ds_innermost.A): ", typeof(lo_ds_innermost.A))  # returns PETSc2.PetscMat

      newton_data_ds, rhs_vec_ds = setupNewton(mesh, mesh, sbp, eqn, opts, ls_ds)
      newton_data_ds.itermax = 30
      recalc_policy_ds = getRecalculationPolicy(opts, "CN")

    end   # end if opts["perturb_Ma_CN"]


    # setup all the checkpointing related data
    chkpointer, chkpointdata, skip_checkpoint = CNDS_checkpoint_setup(mesh, opts, myrank, finaliter)
    istart = chkpointdata.i
    i_test = chkpointdata.i_test
    v_vec = chkpointdata.v_vec
    drag_array = chkpointdata.drag_array
    term23 = chkpointdata.term23
    println(BSTDOUT, "\n >>>> Loaded checkpoint")
    println(BSTDOUT, " istart: ", istart)
    println(BSTDOUT, " i_test: ", i_test)
    println(BSTDOUT, " vecnorm(v_vec): ", vecnorm(v_vec))
    println(BSTDOUT, " vecnorm(drag_array): ", vecnorm(drag_array))
    println(BSTDOUT, " term23: ", term23)

    #------------------------------------------------------------------------------
    # capture direct sensitivity at the IC
    # v is the direct sensitivity, du/dM
    # Ma has been perturbed during setup, in types.jl when eqn.params is initialized (old method)
    if opts["write_drag"]
      objective = createFunctional(mesh, sbp, eqn, opts, 1)    # 1 is the functional num
      drag = real(evalFunctional(mesh, sbp, eqn, opts, objective))
      # note about drag writing: file_dict populated and file opened in src/solver/euler/types.jl
      @mpi_master f_drag = eqn.file_dict[opts["write_drag_fname"]]
    end
    # only perform this at truly the first tstep, not if restarting
    if opts["write_drag"] && ! opts["is_restart"]
      @mpi_master println(f_drag, 1, " ", drag)
      println("i: ", 1, "  myrank: ", myrank,"  drag: ", drag)
      drag_array = zeros(Float64, finaliter)
      drag_array[1] = real(drag)
      @mpi_master flush(f_drag)
    end
    if opts["write_L2vnorm"]

      @mpi_master f_L2vnorm = eqn.file_dict[opts["write_L2vnorm_fname"]]

      # for visualization of element level DS energy
      old_q_vec = zeros(eqn.q_vec)
      R_stab = zeros(eqn.q_vec)
      v_energy = zeros(eqn.q_vec)

      @mpi_master f_v_energy_norm = open("v_energy_norm.dat", "w")
      @mpi_master f_v_vec_norm = open("v_vec_norm.dat", "w")
      @mpi_master f_i_test = open("i_test.dat", "w")
      # @mpi_master f_check1 = open("check1.dat", "w")
    end
    if opts["perturb_Ma_CN"]


      #------------------------------------------------------------------------------
      # allocation of objects for stabilization routine
      if opts["stabilize_v"]

        blocksize = mesh.numDofPerNode*mesh.numNodesPerElement
        stab_A = DiagJac(Complex128, blocksize, mesh.numEl)
        stab_assembler = AssembleDiagJacData(mesh, sbp, eqn, opts, stab_A)
        clipJacData = ClipJacData(mesh.numDofPerNode*mesh.numNodesPerElement)

        # Bv = zeros(Float64, length(q_vec), )
        Bv = zeros(Complex128, length(eqn.q_vec))
        tmp_imag = zeros(Float64, length(eqn.q_vec))
        dqimag_vec = zeros(Bv)

        @mpi_master f_stabilize_v = open("stabilize_v_updates.dat", "w")        # TODO: buffered IO

      end

    end

    # only perform this at truly the first tstep, not if restarting
    if opts["perturb_Ma_CN"] && ! opts["is_restart"]

      # this is the IC, so it gets the first time step's quad_weight
      i = 1       # note that timestep loop below starts at i = 2
      quad_weight = calcQuadWeight(i, dt, finaliter)

      i_test = 10

      # Note: no stabilization of q_vec at the IC
      #   no evalResidual yet, and hasn't entered the time-stepper yet
      #   so no appropriate scaling factors like delta_t or fac or anything

      # direct sensitivity vector       This is at the IC
      fill!(v_vec, 0.0)     # We should do this instead of zeros, so it's not re-allocated every time. 
                            # This also permits v_vec to be Complex{Float64} or Float64 independent 
                            # of the type of q_vec or res_vec, and doesn't change it.
      # for v_ix = 1:length(v_vec)
        # v_vec[v_ix] = imag(eqn.q_vec[v_ix])/Ma_pert_mag
      # end

      # Note: CSR comment on all lines new to CS'ing R method

      dDdu = zeros(eqn.q)      # First allocation of dDdu. fill! used below, during timestep loop
      # evalFunctional calls disassembleSolution, which puts q_vec into q
      # should be calling evalFunctional, not calcFunctional.
      evalFunctionalDeriv(mesh, sbp, eqn, opts, objective, dDdu)    # dDdu is func_deriv_arr
      println(" >>>> i: ", i, "  quad_weight: ", quad_weight, "  dDdu: ", vecnorm(dDdu), "  v_vec: ", vecnorm(v_vec))

      # do the dot product of the two terms, and save
      # this dot product is: dJdu*dudM
      dDdu_vec = zeros(Complex128, mesh.numDofPerNode * mesh.numNodesPerElement * mesh.numEl,)
      array3DTo1D(mesh, sbp, eqn, opts, dDdu, dDdu_vec)      # dDdu -> dDdu_vec

      for v_ix = 1:length(v_vec)
        # this accumulation occurs across all dofs and all time steps.
        term23 += quad_weight * dDdu_vec[v_ix] * v_vec[v_ix]
      end

      if opts["write_L2vnorm"]
        L2_v_norm = calcNorm(eqn, v_vec)
        @mpi_master println(f_L2vnorm, i, "  ", L2_v_norm)
      end

    end   # end if opts["perturb_Ma_CN"]

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # End direct sensitivity setup: needs to be before eqn_deepcopy call to setup eqn_nextstep
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  flush(BSTDOUT)

  # eqn_nextstep = deepcopy(eqn)
  eqn_nextstep = eqn_deepcopy(mesh, sbp, eqn, opts)

  # TODO: copyForMultistage! does not give correct values.
  #     deepcopy works for now, but uses more memory than copyForMultistage!, if it worked
  # eqn_nextstep = copyForMultistage!(eqn)
  eqn_nextstep.q = reshape(eqn_nextstep.q_vec, mesh.numDofPerNode,
                           mesh.numNodesPerElement, mesh.numEl)
  eqn_nextstep.res = reshape(eqn_nextstep.res_vec, mesh.numDofPerNode,
                             mesh.numNodesPerElement, mesh.numEl)

  @debug1 println("============ In CN ============")

  #-------------------------------------------------------------------------------
  # allocate Jac outside of time-stepping loop

  # NOTE: eqn_nextstep changed to eqn 20161013
  pc, lo = getCNPCandLO(mesh, sbp, eqn, opts)
  ls = StandardLinearSolver(pc, lo, eqn.comm, opts)
  newton_data, rhs_vec = setupNewton(mesh, mesh, sbp, eqn, opts, ls)
  newton_data.itermax = 30
  recalc_policy = getRecalculationPolicy(opts, "CN")

  #------------------------------------------------------------------------------
  # ### Main timestepping loop ###
  #   this loop is 2:(t_steps+1) when not restarting
  for i = istart:(t_steps + 1) ##################################################################################################

    println(BSTDOUT, " ----- i = $i -----")
    finaliter = calcFinalIter(t_steps, itermax)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if opts["perturb_Ma_CN"]
      quad_weight = calcQuadWeight(i, dt, finaliter)
    end
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    t = (i-2)*h
    if (i % opts["output_freq"]) == 0
      @mpi_master println(BSTDOUT, "\n============== i = ", i, ", t = ", t, " ==============\n")
    end

    if (i % output_freq) == 0
      @mpi_master flush(BSTDOUT)
      @mpi_master if opts["perturb_Ma_CN"]
        flush(f_v_energy_norm)
        flush(f_v_vec_norm)
      end
      @mpi_master flush(f_i_test)
      @mpi_master if opts["write_drag"]
        flush(f_drag)
      end
    end

    if use_checkpointing && i % chkpoint_freq == 0
      if skip_checkpoint
        skip_checkpoint = false
      else

        @mpi_master println(BSTDOUT, "Saving checkpoint at timestep ", i)
        skip_checkpoint = false
        # save all needed variables to the chkpointdata
        chkpointdata.i = i
        chkpointdata.i_test = i_test
        chkpointdata.v_vec = v_vec
        chkpointdata.drag_array = drag_array
        chkpointdata.term23 = term23
        println(BSTDOUT, " i: ", i)
        println(BSTDOUT, " i_test: ", i_test)
        println(BSTDOUT, " vecnorm(v_vec): ", vecnorm(v_vec))
        println(BSTDOUT, " vecnorm(drag_array): ", vecnorm(drag_array))
        println(BSTDOUT, " term23: ", term23, "\n")

        if countFreeCheckpoints(chkpointer) == 0
          freeOldestCheckpoint(chkpointer)  # make room for a new checkpoint
        end

        # save the checkpoint
        saveNextFreeCheckpoint(chkpointer, mesh, sbp, eqn, opts, chkpointdata)
      end   # end of if skip_checkpoint check
    end   # end of if use_checkpointing check

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FAC = 1.0       # for making stabilization +=    (now that explicit fix has been implemented: no fac needed)
    # FAC = -1.0       # for making stabilization -=  (needed for eig clipping?)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    t_nextstep = t + h

    # allow for user to select CN's internal Newton's method. Only supports dense FD Jacs, so only for debugging
    if opts["cleansheet_CN_newton"]
      cnNewton(mesh, sbp, opts, h, f, eqn, eqn_nextstep, t_nextstep)
    else

      ctx_residual = (f, eqn, h, newton_data)

      # recalculate PC and LO if needed
      doRecalculation(recalc_policy, i,
                    ls, mesh, sbp, eqn_nextstep, opts, ctx_residual, t_nextstep)

      newtonInner(newton_data, mesh, sbp, eqn_nextstep, opts, cnRhs, ls, 
                  rhs_vec, ctx_residual, t_nextstep)

                  # The calcResidual call is buried inside newtonInner.
                  #   newtonInner -> rhs_func -> physicsRhs (residual_evaluation.jl)
                  # TODO: need to save the complex part of R inside physicsRhs
    end

    # do the callback using the current eqn object at time t
    eqn.majorIterationCallback(i, mesh, sbp, eqn, opts, BSTDOUT)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # TODO: stabilization in newtonInner
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # need to assemble solution into res_vec?
    res_norm = calcNorm(eqn, eqn.res_vec)
    # logging
    @mpi_master if i % 1 == 0
      println(f1, i, " ", res_norm)
    end
    
    @mpi_master if i % output_freq == 0
      println(BSTDOUT, " flushing convergence.dat to disk")
      flush(f1)
    end


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if opts["write_drag"]
      drag = real(evalFunctional(mesh, sbp, eqn, opts, objective))
      @mpi_master f_drag = eqn.file_dict[opts["write_drag_fname"]]
      @mpi_master println(f_drag, i, " ", drag)
      @mpi_master drag_array[i] = real(drag)
      @mpi_master if (i % opts["output_freq"]) == 0
        flush(f_drag)
      end
    end

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Start direct sensitivity calc's for each time step      2222
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if opts["perturb_Ma_CN"]

      # ctx_residual: f, eqn, h, newton_data

      if i != 1     # can't do F(q^(n+1)) + F(q^(n)) until we have both

        # Store both q_vec's & both res_vec's. Will be put back after all DS calcs.
        # TODO: do we need all of these?
        for ix_dof = 1:mesh.numDof
          beforeDS_eqn_q_vec[ix_dof] = eqn.q_vec[ix_dof]
          beforeDS_eqn_nextstep_q_vec[ix_dof] = eqn_nextstep.q_vec[ix_dof]

          beforeDS_eqn_res_vec[ix_dof] = eqn.res_vec[ix_dof]
          beforeDS_eqn_nextstep_res_vec[ix_dof] = eqn_nextstep.res_vec[ix_dof]
        end
        println(BSTDOUT, " +++ Storing beforeDS vecs +++")
        println(BSTDOUT, " vecnorm(beforeDS_eqn_q_vec): ", vecnorm(beforeDS_eqn_q_vec))
        println(BSTDOUT, " vecnorm(beforeDS_eqn_nextstep_q_vec): ", vecnorm(beforeDS_eqn_nextstep_q_vec))
        println(BSTDOUT, " vecnorm(beforeDS_eqn_res_vec): ", vecnorm(beforeDS_eqn_res_vec))
        println(BSTDOUT, " vecnorm(beforeDS_eqn_nextstep_res_vec): ", vecnorm(beforeDS_eqn_nextstep_res_vec), "\n")
        # println(BSTDOUT, " q_vec and res_vec, both eqn & eqn_nextstep saved")


        #------------------------------------------------------------------------------
        # This section calculates dR/dM
        # println(BSTDOUT, "> Now solving dR/dM")

        Ma_pert_mag = opts["perturb_Ma_magnitude"]
        pert = complex(0, Ma_pert_mag)

        ### using cnRhs, complex-step
        eqn_nextstep.params.Ma += pert
        eqn.params.Ma += pert
        ctx = (f, eqn, h)

        # Form unsteady residual (res_hat) with F(q^(n)) and F(q^(n+1))
        #   res_hat_vec = q^(n+1) - q^(n) - 0.5*Minv*dt* (F(q^(n+1)) - F(q^(n)))
        cnRhs(mesh, sbp, eqn_nextstep, opts, res_hat_vec, ctx, t)
        eqn_nextstep.params.Ma -= pert
        eqn.params.Ma -= pert

        for ix_dof = 1:mesh.numDof
          dRdM_vec[ix_dof] = imag(res_hat_vec[ix_dof])/Ma_pert_mag
        end
        ### end cnRhs, complex-step

        ### FD check of dRdM
        #=
        FD_pert = 1e-7      # looks like minimum of the v is 1e-7
        eqn_nextstep.params.Ma += FD_pert
        eqn.params.Ma += FD_pert
        ctx = (f, eqn, h)

        cnRhs(mesh, sbp, eqn_nextstep, opts, res_hat_vec_pospert, ctx, t)

        eqn.params.Ma -= FD_pert
        eqn_nextstep.params.Ma -= FD_pert

        cnRhs(mesh, sbp, eqn_nextstep, opts, res_hat_vec_nopert, ctx, t)
        for ix_dof = 1:mesh.numDof
          dRdM_vec_FD[ix_dof] = (res_hat_vec_pospert[ix_dof] - res_hat_vec_nopert[ix_dof])/FD_pert
        end
        println(BSTDOUT, "  vecnorm(dRdM_vec_FD): ", vecnorm(dRdM_vec_FD))
        flush(BSTDOUT)
        check1 = vecnorm(dRdM_vec_FD - dRdM_vec)
        print(BSTDOUT, "  >>> dRdM verify: vecnorm(dRdM_vec_FD - dRdM_vec): ", check1)
        flush(BSTDOUT)
        if check1 < 10*FD_pert
          println(BSTDOUT, "   PASS")
        else
          println(BSTDOUT, "   FAIL")
        end
        flush(BSTDOUT)
        =#
        dRdM_norm_global = calcNorm(eqn, dRdM_vec)
        println(BSTDOUT, " +++ dRdM_norm_global: ", dRdM_norm_global)
        flush(BSTDOUT)
        ### end check

        # should I be collecting into q?

        #------------------------------------------------------------------------------

        #------------------------------------------------------------------------------
        # This section calculates dR/dq * v^(n)
        # println(BSTDOUT, "> Now solving dR/dq * v^(n)")

        # v_vec currently holds v at timestep n: v^(n)

        # need to add epsilon*v_vec*im (evi) to q_vec, then re-evaluate res_hat at the updated q
        for ix_dof = 1:mesh.numDof
          eqn.q_vec[ix_dof]          += Ma_pert_mag*im*v_vec[ix_dof]
          # TODO: commented this eqn_nextstep application out. Verify with formulation.
          #       Do I need to save and restore eqn_nextstep??? Could save cycles
          # eqn_nextstep.q_vec[ix_dof]          += Ma_pert_mag*im*v_vec[ix_dof]
        end

        if opts["parallel_type"] == 2 && mesh.npeers > 0
          startSolutionExchange(mesh, sbp, eqn, opts)
        end
        f(mesh, sbp, eqn, opts)             # F(q^(n) + evi) now in eqn.res_vec (confirmed w/ random test)
        array3DTo1D(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)

        for ix_dof = 1:mesh.numDof
          
          # This calculation:
          #   dR_hat^(n)/dq^(n) * v^(n) = -v^(n) - 0.5*Minv*dt* Im[F(q^(n) + epsilon*v^(n)*im)]/epsilon
          # Note that Minv is applied inside evalResidual already.
          dRdq_vn_prod[ix_dof] = - v_vec[ix_dof] - 0.5*dt*imag(eqn.res_vec[ix_dof])/Ma_pert_mag

          # combine (-dRdM - dRdq * v^(n)) into b
          b_vec[ix_dof] = - dRdq_vn_prod[ix_dof] - dRdM_vec[ix_dof] 

        end
        #=
        # output of norms of quantities every time step (debugging only)
        =#
        dRdq_vn_prod_norm_global = calcNorm(eqn, dRdq_vn_prod)
        println(BSTDOUT, " +++ dRdq_vn_prod_norm_global: ", dRdq_vn_prod_norm_global)
        b_vec_norm_global = calcNorm(eqn, b_vec)
        println(BSTDOUT, " +++ b_vec_norm_global: ", b_vec_norm_global)
        flush(BSTDOUT)

        ### Only for serial julia sparse.
        ### If serial Petsc Jac, A_mul_B is very slow (bc of mixing PetscMat & Julia vecs, improper method called)
        ### If parallel, get an "only local values currently supported"
        #=
        modifyCNJacForMatFreeCheck(lo_ds, mesh, sbp, eqn, opts, ctx_residual, t)
        if opts["jac_type"] == 3 error("jac_type is 3 but you're doing check 2") end
        A_mul_B!(TEST_dRdq_vn_prod, lo_ds_innermost.A, v_vec)
        modifyCNJacForMatFreeCheck_reverse(lo_ds, mesh, sbp, eqn, opts, ctx_residual, t)
        check2 = vecnorm(dRdq_vn_prod - TEST_dRdq_vn_prod)
        print(BSTDOUT, "  >>> mat-vec product verify: vecnorm(dRdq_vn_prod - TEST_dRdq_vn_prod): ", check2)
        if check2 < 1e-15
          println(BSTDOUT, "   PASS")
        else
          println(BSTDOUT, "   FAIL")
        end
        =#

        #------------------------------------------------------------------------------
        # Now the calculation of v_ix at n+1
        # Solve: dR/dq^(n+1) v^(n+1) = b
        #--------
        # println(BSTDOUT, "> Now solving for v^(n+1)")

        #------------------------------------------------------------------------------
        # Restore q_vec & res_vec for both eqn & eqn_nextstep, as they were
        #   before the DS calcs
        for ix_dof = 1:mesh.numDof
          eqn.q_vec[ix_dof] = beforeDS_eqn_q_vec[ix_dof]
          eqn_nextstep.q_vec[ix_dof] = beforeDS_eqn_nextstep_q_vec[ix_dof]

          eqn.res_vec[ix_dof] = beforeDS_eqn_res_vec[ix_dof]
          eqn_nextstep.res_vec[ix_dof] = beforeDS_eqn_nextstep_res_vec[ix_dof]
        end

        # Contents of ctx_residual: f, eqn, h, newton_data
        if opts["stabilize_v"]
          ctx_residual = (f, eqn, h, newton_data, stab_A, stab_assembler, clipJacData, v_vec)
        end
        
        # Update linear operator:
        #   The Jacobian ∂R_hat/∂q^(n+1) is lo_ds_innermost.A
        calcLinearOperator(ls_ds, mesh, sbp, eqn_nextstep, opts, ctx_residual, t)
        # Note: this is properly modifying the Jac for CN.
        flush(BSTDOUT)

        println(BSTDOUT, " Sleeping for 5s...")
        flush(BSTDOUT)
        run(`sleep 5`)
        lo_ds_innermost_A_norm_global = calcNorm(eqn, lo_ds_innermost.A)
        println(BSTDOUT, " +++ vecnorm(lo_ds_innermost.A): ", vecnorm(lo_ds_innermost.A))
        println(BSTDOUT, " +++ lo_ds_innermost_A_norm_global (after cLO & stab): ", lo_ds_innermost_A_norm_global)

        fill!(v_vec, 0.0)

        eqn_q_vec_norm_global = calcNorm(eqn, eqn.q_vec)
        eqnnextstep_q_vec_norm_global = calcNorm(eqn, eqn_nextstep.q_vec)
        eqn_res_vec_norm_global = calcNorm(eqn, eqn.res_vec)
        eqnnextstep_q_vec_norm_global = calcNorm(eqn, eqn_nextstep.res_vec)

        println(BSTDOUT, " eqn_q_vec_norm_global: ", eqn_q_vec_norm_global)
        println(BSTDOUT, " eqnnextstep_q_vec_norm_global: ", eqnnextstep_q_vec_norm_global)
        println(BSTDOUT, " eqn_res_vec_norm_global: ", eqn_res_vec_norm_global)
        println(BSTDOUT, " eqnnextstep_q_vec_norm_global: ", eqnnextstep_q_vec_norm_global)
        # linearSolve: solves Ax=b for x. 
        #   ls::StandardLinearSolver, b::AbstractVector (RHS), x::AbstractVector  (what is solved for)
        linearSolve(ls_ds, b_vec, v_vec)

        v_vec_norm_global = calcNorm(eqn, v_vec)
        println(BSTDOUT, " +++ v_vec_norm_global: ", v_vec_norm_global)

        ### Only for serial julia sparse.
        ### If serial Petsc Jac, A_mul_B is very slow (bc of mixing PetscMat & Julia vecs, improper method called)
        ### If parallel, get an "only local values currently supported"
        #=
        if opts["jac_type"] == 3 error("jac_type is 3 but you're doing check 3") end
        A_mul_B!(TEST_b_vec, lo_ds_innermost.A, v_vec)
        # println(BSTDOUT, " cond(full(lo_ds_innermost.A)): ", cond(full(lo_ds_innermost.A)))
        println(BSTDOUT, " vecnorm(b_vec): ", vecnorm(b_vec))
        println(BSTDOUT, " vecnorm(TEST_b_vec): ", vecnorm(TEST_b_vec))
        check3 = vecnorm(TEST_b_vec - b_vec)
        print(BSTDOUT, "  >>> b_vec verify: ", check3)
        if check3 < 1e-15
          println(BSTDOUT, "   PASS")
        else
          println(BSTDOUT, "   FAIL")
        end
        flush(BSTDOUT)
        =#

        #### The above is all new CSR code

        #------------------------------------------------------------------------------
        # getting from v_vec to term23
        fill!(dDdu, 0.0)
        evalFunctionalDeriv(mesh, sbp, eqn, opts, objective, dDdu)    # dDdu is func_deriv_arr
        fill!(dDdu_vec, 0.0)     # not sure this is necessary
        array3DTo1D(mesh, sbp, eqn, opts, dDdu, dDdu_vec)      # dDdu -> dDdu_vec

        old_term23 = term23
        for v_ix = 1:length(v_vec)
          # this accumulation occurs across all dofs and all time steps.
          term23 += quad_weight * dDdu_vec[v_ix] * v_vec[v_ix]
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



      end     # end 'if i != 1'

      # for visualization of element level DS energy: end of all stages
      if opts["write_L2vnorm"]
        for ix_dof = 1:mesh.numDof
          v_energy[ix_dof] = v_vec[ix_dof] * eqn.M[ix_dof] * v_vec[ix_dof]
        end

        if (i % output_freq) == 0
          # saveSolutionToMesh(mesh, v_energy)
          # fname = string("v_energy_", i)
          # writeVisFiles(mesh, fname)

          saveSolutionToMesh(mesh, v_vec)
          fname = string("v_vec_", i)
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
        v_vec_norm = calcNorm(eqn, v_vec)

      end   # end of if opts["write_L2vnorm"]
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    end   # end of if opts["perturb_Ma_CN"]

    # This needs to be above the checkpoint write, in case the checkpoint is written before the 
    #   files are flushed. This would cause a gap in the data files.
    @mpi_master println(f_v_energy_norm, i, "  ", real(v_energy_norm))
    @mpi_master println(f_v_vec_norm, i, "  ", real(v_vec_norm))
    @mpi_master println(f_i_test, i, "  ", i_test, "  ", t)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # End direct sensitivity calc's for each time step
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    #--------------------------------------------------------------------------------------------------
    # Return to original CN code: moving eqn_nextstep into eqn
  
    # This allows the solution to be updated from _nextstep without a deepcopy.
    # There are two memory locations used by eqn & eqn_nextstep, 
    #   and this flips the location of eqn & eqn_nextstep every time step
    eqn_temp = eqn
    eqn = eqn_nextstep
    eqn_nextstep = eqn_temp

    # Note: we now need to copy the updated q over for the initial newton guess
    for j = 1:mesh.numDof
      eqn_nextstep.q_vec[j] = eqn.q_vec[j]
    end

    array1DTo3D(mesh, sbp, eqn_nextstep, opts, eqn_nextstep.q_vec, eqn_nextstep.q)


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # direct-sensitivity: moved these breaks here. Before they were before the eqn/eqn_nextstep/eqn_temp flip
    # check stopping conditions

    # commenting this out - was not allowing run w/ steady-generated IC to proceed. 
    # Why is it here anyway?
    #=
    if (res_norm < res_tol)  # ???
      if myrank == 0
        println(BSTDOUT, "breaking due to res_tol, res norm = $res_norm")
        close(f1)
        flush(BSTDOUT)
      end
      break       # EXIT condition
    end
    =#


    if use_itermax && i > itermax
      if myrank == 0
        println(BSTDOUT, "breaking due to itermax")
        close(f1)
        flush(BSTDOUT)
      end
      break       # EXIT condition
    end

    i_test = i_test+10

  end   # end of t step loop

  # final time update
  t += h

  # depending on how many timesteps we do, this may or may not be necessary
  #   usage: copyForMultistage!(dest, src)
  copyForMultistage!(eqn, eqn_nextstep)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  @mpi_master begin
    println("---------------------------------------------")
    println("   CN: final time step reached. t = $t")
    println("---------------------------------------------")
  end

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Start direct sensitivity calc's after time step loop    3333
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      if opts["perturb_Ma_CN"]

        @mpi_master if opts["stabilize_v"]
          close(f_stabilize_v)
        end

        @mpi_master close(f_drag)
        @mpi_master writedlm("drag_array.dat", drag_array)

        println(BSTDOUT, " size(drag_array): ", size(drag_array))

        Cd, dCddM = calcDragTimeAverage(mesh, sbp, eqn, opts, dt, finaliter, useArray=false)   # will use eqn.params.Ma
        Cd_file, dCddM_file = calcDragTimeAverage(mesh, sbp, eqn, opts, dt, finaliter, useArray=true, drag_array=drag_array)   # will use eqn.params.Ma
        # println(BSTDOUT, " Cd_file - Cd: ", Cd_file - Cd)
        # println(BSTDOUT, " dCddM_file - dCddM: ", dCddM_file - dCddM)
        term23 = term23 * 1.0/t     # final step of time average: divide by total time
        global_term23 = MPI.Allreduce(term23, MPI.SUM, mesh.comm)
        total_dCddM = dCddM + global_term23

        # Cd calculations
        @mpi_master begin
          println(BSTDOUT, " Drag calculations made from array, not file.")
          f_total_dCddM = open("total_dCddM.dat", "w")
          println(f_total_dCddM, " Cd: ", Cd)
          println(f_total_dCddM, " dCd/dM: ", dCddM)
          println(f_total_dCddM, " global_term23: ", global_term23)
          println(f_total_dCddM, " total dCd/dM: ", total_dCddM)
          flush(f_total_dCddM)
          close(f_total_dCddM)
          println(BSTDOUT, " Cd: ", Cd)
          println(BSTDOUT, " dCd/dM: ", dCddM)
          println(BSTDOUT, " global_term23: ", global_term23)
          println(BSTDOUT, " total dCd/dM: ", total_dCddM)
        end

      end   # end if opts["perturb_Ma_CN"]

      @mpi_master begin
        f_Ma = open("Ma.dat", "w")
        println(f_Ma, eqn.params.Ma)
        close(f_Ma)
        f_dt = open("delta_t.dat", "w")
        println(f_dt, dt)
        close(f_dt)

        println(BSTDOUT, " ")
        println(BSTDOUT, " run parameters that were used:")
        # if opts["perturb_Ma_CN"]
          # println("    Ma: ", eqn.params.Ma + Ma_pert_mag)
        # else
          println(BSTDOUT, "    Ma: ", eqn.params.Ma)
        # end
        println(BSTDOUT, "    aoa: ", eqn.params.aoa)
        println(BSTDOUT, "    dt: ", dt)
        println(BSTDOUT, "    a_inf: ", eqn.params.a_free)
        println(BSTDOUT, "    rho_inf: ", eqn.params.rho_free)
        println(BSTDOUT, "    c: ", 1.0)
        println(BSTDOUT, "    mesh.coord_order: ", mesh.coord_order)
        println(BSTDOUT, " ")
        println(BSTDOUT, "    opts[stabilization_method]: ", opts["stabilization_method"])
        println(BSTDOUT, "    opts[output_freq]: ", opts["output_freq"])
        println(BSTDOUT, "    opts[use_itermax]: ", opts["use_itermax"])
        println(BSTDOUT, "    opts[itermax]: ", opts["itermax"])
        println(BSTDOUT, "    opts[use_checkpointing]: ", opts["use_checkpointing"])
        println(BSTDOUT, "    opts[checkpoint_freq]: ", opts["checkpoint_freq"])
        println(BSTDOUT, "    opts[ncheckpoints]: ", opts["ncheckpoints"])
        println(BSTDOUT, " ")
        println(BSTDOUT, " sbp.w: ", sbp.w)
      end

      if opts["write_L2vnorm"]
        @mpi_master close(f_L2vnorm)
        @mpi_master close(f_v_energy_norm)
        @mpi_master close(f_v_vec_norm)
        @mpi_master close(f_i_test)
        # @mpi_master close(f_check1)
      end
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # End direct sensitivity calc's after time step loop
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  flush(BSTDOUT)

  #TODO: return the NewtonData?
  free(newton_data)

  @mpi_master println("============= end of CN: t = $t ===============")
  return t        # EXIT condition

  flush(BSTDOUT)

end   # end of crank_nicolson_ds function

DS_LO_file = string(Pkg.dir("PDESolver"),"/src/NonlinearSolvers/","crank_nicolson_ds-LO.jl")
include(DS_LO_file)

DS_functions_file = string(Pkg.dir("PDESolver"),"/src/NonlinearSolvers/","crank_nicolson_ds-functions.jl")
include(DS_functions_file)

DS_stab_file = string(Pkg.dir("PDESolver"),"/src/NonlinearSolvers/","crank_nicolson_ds-stabilize.jl")
include(DS_stab_file)

