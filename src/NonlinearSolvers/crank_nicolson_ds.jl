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

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Start direct sensitivity setup    1111
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if opts["perturb_Ma_CN"]
      term23 = zero(eqn.params.Ma)
      dDdu = zeros(eqn.res)         # dDdu is of size eqn.q or eqn.res (dDdu is term2 in rk4_ds.jl)
      dDdu_vec = zeros(Complex128, mesh.numDofPerNode * mesh.numNodesPerElement * mesh.numEl,)
      v_vec = zeros(Float64, length(eqn.res_vec))   # needs to be Float64, even if res_vec is cplx

      ###### CSR: This section is the only section that contains things that are 
      #           new to the CN version of the complex step method (CSR: 'complex step residual')
      old_res_vec_Maimag = zeros(eqn.res_vec)    # corresponds to eqn
      new_res_vec_Maimag = zeros(eqn.res_vec)    # corresponds to eqn_nextstep
      old_q_vec_Maimag = zeros(eqn.q_vec)    # corresponds to eqn
      new_q_vec_Maimag = zeros(eqn.q_vec)    # corresponds to eqn_nextstep
      res_hat_vec = zeros(eqn.res_vec)    # unsteady residual
      # TODO: ensure all of these are being used when debugging complete; otherwise, cleanup

      beforeDS_eqn_q_vec = zeros(eqn.q_vec)
      beforeDS_eqn_nextstep_q_vec = zeros(eqn.q_vec)
      beforeDS_eqn_res_vec = zeros(eqn.res_vec)
      beforeDS_eqn_nextstep_res_vec = zeros(eqn.res_vec)

      dRdM_vec = zeros(eqn.res_vec)
      b_vec = zeros(Float64, length(eqn.res_vec))   # needs to be Float64, even if res_vec is cplx
      TEST_b_vec = zeros(b_vec)
      dRdq_vn_prod = zeros(eqn.res_vec)

      #------------------------------------------------------------------------------
      # DS linear solver
      pc_ds, lo_ds = getCNDSPCandLO(mesh, sbp, eqn, opts)
      ls_ds = StandardLinearSolver(pc_ds, lo_ds, eqn.comm, opts)
      println(" typeof(pc_ds): ", typeof(pc_ds))
      println(" typeof(lo_ds): ", typeof(lo_ds))
      println(" typeof(ls_ds): ", typeof(ls_ds))

      #------------------------------------------------------------------------------
      # debug linear solver
      # pc_debug, lo_debug = getCNDSPCandLO(mesh, sbp, eqn, opts)
      # ls_debug = StandardLinearSolver(pc_debug, lo_debug, eqn.comm, opts)


      # these are needed for direct manipulation of the Petsc Jac later
      lo_ds_innermost = getBaseLO(lo_ds)     # makes sure to return the innermost LO object (basically returning the Jac)
      # not needed?
      println(" typeof(lo_ds_innermost): ", typeof(lo_ds_innermost))    # returns LinearSolvers.PetscMatLO
      println(" typeof(lo_ds_innermost.A): ", typeof(lo_ds_innermost.A))  # returns PETSc2.PetscMat

      newton_data_ds, rhs_vec_ds = setupNewton(mesh, mesh, sbp, eqn, opts, ls_ds)
      newton_data_ds.itermax = 30
      recalc_policy_ds = getRecalculationPolicy(opts, "CN")

    end   # end if opts["perturb_Ma_CN"]


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
      old_q_vec = zeros(eqn.q_vec)
      R_stab = zeros(eqn.q_vec)
      v_energy = zeros(eqn.q_vec)

      @mpi_master f_v_energy_stageall = open("v_energy_data_stageall.dat", "w")
    end
    if opts["perturb_Ma_CN"]

      # this is the IC, so it gets the first time step's quad_weight
      i = 1       # note that timestep loop below starts at i = 2
      finaliter = calcFinalIter(t_steps, itermax)
      quad_weight = calcQuadWeight(i, dt, finaliter)

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
      println(" >>>> i: ", i, "  quad_weight: ", quad_weight, "  dDdu: ", vecnorm(dDdu), "  v_vec: ", vecnorm(v_vec))    # matches rk4_ds to 1e-5 or so, consistent with drag variation

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

   # setup all the checkpointing related data
  chkpointer, chkpointdata, skip_checkpoint = explicit_checkpoint_setup(opts, myrank)
  istart = chkpointdata.i

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
  for i = istart:(t_steps + 1)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if opts["perturb_Ma_CN"]
      quad_weight = calcQuadWeight(i, dt, finaliter)

      # for j = 1:length(eqn.q_vec)                 # store imaginary part of q_vec
        # q_vec_save_imag[i] = imag(eqn.q_vec[i])
        # Shouldn't need to remove the imaginary component here - PETSc call only takes in the real part
      # end

    end   # end if opts["perturb_Ma_CN"]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    t = (i-2)*h
    @mpi_master println(BSTDOUT, " ")
    @mpi_master println(BSTDOUT, "============== i = ", i, ", t = ", t, " ==============\n")
    @debug1 println(eqn.params.f, "====== CN: at the top of time-stepping loop, t = $t, i = $i")
    @debug1 flush(eqn.params.f)

    if use_checkpointing && i % chkpoint_freq == 0
      if skip_checkpoint
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
        saveNextFreeCheckpoint(chkpointer, mesh, sbp, eqn, opts, chkpointdata)
      end   # end of if skip_checkpoint check
    end   # end of if use_checkpointing check

#=
    #----------------------------
    # zero out Jac
    #   this works for both PETSc and Julia matrices.
    #   when jac is a Julia matrix, this is effectively wrapping: fill!(jac, 0.0)
    if jac_type != 4
      MatZeroEntries(jac)
    end
=#
    # NOTE: Must include a comma in the ctx tuple to indicate tuple
    # f is the physics function, like evalEuler

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # save q_vec from the last time step in old_q_vec: used for v_energy calcs
    if opts["write_L2vnorm"]
      # for visualization of element level DS energy
      fill!(old_q_vec, 0.0)
      for j = 1:length(eqn.q_vec)
        old_q_vec[j] = eqn.q_vec[j]
      end
    end
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

      # println(BSTDOUT, " Calling newtonInner from CN.")
      flush(BSTDOUT)
      newtonInner(newton_data, mesh, sbp, eqn_nextstep, opts, cnRhs, ls, 
                  rhs_vec, ctx_residual, t_nextstep)

                  # The calcResidual call is buried inside newtonInner.
                  #   newtonInner -> rhs_func -> physicsRhs (residual_evaluation.jl)
                  # TODO: need to save the complex part of R inside physicsRhs
    end
    # print(BSTDOUT, " Newton solve on primal complete.")

    # do the callback using the current eqn object at time t
    eqn.majorIterationCallback(i, mesh, sbp, eqn, opts, BSTDOUT)
    # print(BSTDOUT, " majorIterationCallback called.")
    println(BSTDOUT, " ")
    flush(BSTDOUT)

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
      @mpi_master if (i % opts["output_freq"]) == 0
        flush(f_drag)
      end
    end

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Start direct sensitivity calc's for each time step      2222
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if opts["perturb_Ma_CN"]

      # This entire code block is called with an include().
      # It is wrapped in:  if opts["perturb_Ma_CN"]
      # Such that it should only be included/evaluated when opts["perturb_Ma_CN"] is set

      print(BSTDOUT, " Now entering perturb_Ma_CN block.")
      # ctx_residual: f, eqn, h, newton_data
      flush(BSTDOUT)

      if i != 1     # can't do F(q^(n+1)) + F(q^(n)) until we have both

        # Store both q_vec's & both res_vec's. Will be put back after all DS calcs.
        for ix_dof = 1:mesh.numDof
          beforeDS_eqn_q_vec[ix_dof] = eqn.q_vec[ix_dof]
          beforeDS_eqn_nextstep_q_vec[ix_dof] = eqn_nextstep.q_vec[ix_dof]

          beforeDS_eqn_res_vec[ix_dof] = eqn.res_vec[ix_dof]
          beforeDS_eqn_nextstep_res_vec[ix_dof] = eqn_nextstep.res_vec[ix_dof]
        end
        println(BSTDOUT, " q_vec and res_vec, both eqn & eqn_nextstep saved")


        #------------------------------------------------------------------------------
        # This section calculates dR/dM

        println(BSTDOUT, "> Now solving dR/dM")

        Ma_pert_mag = opts["perturb_Ma_magnitude"]
        pert = complex(0, Ma_pert_mag)
        # eqn.params.Ma += pert
        eqn_nextstep.params.Ma += pert
        println(BSTDOUT, " Perturbing Ma. eqn_nextstep.params.Ma: ", eqn_nextstep.params.Ma)

        # now evalResidual to store into F(q^(n+1))
        f(mesh, sbp, eqn_nextstep, opts)      # F(q^(n+1)) -- with Ma perturbation, now in eqn_nextstep.res_vec

        for ix_dof = 1:mesh.numDof
          # store F(q^(n)): it was calculated as F(q^(n+1)) during the last time step
          old_res_vec_Maimag[ix_dof] = new_res_vec_Maimag[ix_dof]         # F(q^(n)) -- with Ma perturbation
          old_q_vec_Maimag[ix_dof] = new_q_vec_Maimag[ix_dof]             # q^(n) -- with Ma perturbation

          new_res_vec_Maimag[ix_dof] = eqn_nextstep.res_vec[ix_dof]       # F(q^(n+1)) -- with Ma perturbation
          new_q_vec_Maimag[ix_dof] = eqn_nextstep.q_vec[ix_dof]           # q^(n+1) -- with Ma perturbation

          # Form unsteady residual (res_hat) with F(q^(n)) and F(q^(n+1))
          # Note: q^(n+1) and q^(n) should have no imaginary component. Therefore,
          #       we do not need to include them in the calculation of res_hat_vec,
          #       as it is only formed to later consider only its imaginary component.

          # R_hat = q^(n+1) - q^(n) - 0.5*Minv*dt* (F(q^(n+1)) - F(q^(n)))
          # res_hat_vec[ix_dof] = new_q_vec_Maimag[ix_dof] - old_q_vec_Maimag[ix_dof] 
          #                       - 0.5*eqn.Minv[ix_dof]*dt * (new_res_vec_Maimag[ix_dof] + old_res_vec_Maimag[ix_dof])
          res_hat_vec[ix_dof] = -0.5*eqn.Minv[ix_dof]*dt * (new_res_vec_Maimag[ix_dof] + old_res_vec_Maimag[ix_dof])

        end
        # DUPEDEBUG
        #=
        println(BSTDOUT, " typeof(res_hat_vec): ", typeof(res_hat_vec))
        println(BSTDOUT, " vecnorm(new_q_vec_Maimag): ", vecnorm(new_q_vec_Maimag))
        println(BSTDOUT, " vecnorm(old_q_vec_Maimag): ", vecnorm(old_q_vec_Maimag))
        if vecnorm(imag(new_q_vec_Maimag)) > 1e-15
          println(BSTDOUT, " NON-ZERO imag(new_q_vec_Maimag)!")
          println(BSTDOUT, " vecnorm(imag(new_q_vec_Maimag)): ", vecnorm(imag(new_q_vec_Maimag)))
        end
        if vecnorm(imag(old_q_vec_Maimag)) > 1e-15
          println(BSTDOUT, " NON-ZERO imag(old_q_vec_Maimag)!")
          println(BSTDOUT, " vecnorm(imag(old_q_vec_Maimag)): ", vecnorm(imag(old_q_vec_Maimag)))
        end
        println(BSTDOUT, " vecnorm(new_res_vec_Maimag): ", vecnorm(new_res_vec_Maimag))
        println(BSTDOUT, " vecnorm(imag(new_res_vec_Maimag)): ", vecnorm(imag(new_res_vec_Maimag)))
        println(BSTDOUT, " vecnorm(old_res_vec_Maimag): ", vecnorm(old_res_vec_Maimag))
        println(BSTDOUT, " vecnorm(imag(old_res_vec_Maimag)): ", vecnorm(imag(old_res_vec_Maimag)))
        println(BSTDOUT, " vecnorm(res_hat_vec): ", vecnorm(res_hat_vec))
        println(BSTDOUT, " vecnorm(imag(res_hat_vec)): ", vecnorm(imag(res_hat_vec)))
        =#
        # should I be collecting into q?

        # obtain dR/dM using the complex step method
        for ix_dof = 1:mesh.numDof
          dRdM_vec[ix_dof] = imag(res_hat_vec[ix_dof])/Ma_pert_mag     # should this be res_hat_vec??
        end
        println(BSTDOUT, " vecnorm(dRdM_vec): ", vecnorm(dRdM_vec))

        # eqn.params.Ma -= pert
        eqn_nextstep.params.Ma -= pert
        println(BSTDOUT, " Removing Ma perturbation. eqn_nextstep.params.Ma: ", eqn_nextstep.params.Ma)
        #------------------------------------------------------------------------------

        #------------------------------------------------------------------------------
        # This section calculates dR/dq * v^(n)
        println(BSTDOUT, "> Now solving dR/dq * v^(n)")

        # v_vec currently holds v at timestep n: v^(n)

        # need to add epsilon*v_vec*im (evi) to q_vec, then re-evaluate res_hat at the updated q
        # so: need to add evi to old_q_vec too. Because to form res_hat, need F(q^(n)) and F(q^(n+1))
        for ix_dof = 1:mesh.numDof
          eqn.q_vec[ix_dof]          += Ma_pert_mag*im*v_vec[ix_dof]
          # TODO: commented this eqn_nextstep application out. Verify with formulation.
          #       Do I need to save and restore eqn_nextstep??? Could save cycles
          # eqn_nextstep.q_vec[ix_dof]          += Ma_pert_mag*im*v_vec[ix_dof]
        end

        f(mesh, sbp, eqn, opts)             # F(q^(n) + evi) now in eqn.res_vec (confirmed w/ random test)

        for ix_dof = 1:mesh.numDof
          
          # This calculation:
          #   dR_hat^(n)/dq^(n) * v^(n) = 
          #     -v^(n) - 0.5*Minv*dt* Im[F(q^(n) + epsilon*v^(n)*im)]/epsilon
          # Note that Minv is applied inside evalResidual already.
          dRdq_vn_prod[ix_dof] = - v_vec[ix_dof] - 0.5*dt*imag(eqn.res_vec[ix_dof])/Ma_pert_mag

          # combine (-dRdM - dRdq * v^(n)) into b
          b_vec[ix_dof] = - dRdq_vn_prod[ix_dof] - dRdM_vec[ix_dof] 

        end

        # Here was scratch content 1

        #------------------------------------------------------------------------------
        # Now the calculation of v_ix at n+1
        # Solve:
        #   dR/dq^(n+1) v^(n+1) = b
        #--------
        println(BSTDOUT, "> Now solving for v^(n+1)")

        #------------------------------------------------------------------------------
        # Restore q_vec & res_vec for both eqn & eqn_nextstep, as they were
        #   before the DS calcs
        for ix_dof = 1:mesh.numDof
          eqn.q_vec[ix_dof] = beforeDS_eqn_q_vec[ix_dof]
          eqn_nextstep.q_vec[ix_dof] = beforeDS_eqn_nextstep_q_vec[ix_dof]

          eqn.res_vec[ix_dof] = beforeDS_eqn_res_vec[ix_dof]
          eqn_nextstep.res_vec[ix_dof] = beforeDS_eqn_nextstep_res_vec[ix_dof]
        end

        # Update linear operator:
        #   The Jacobian ∂R_hat/∂q^(n+1) is lo_ds_innermost.A
        println(BSTDOUT, "  Now updating linear operator.")
        flush(BSTDOUT)
        # typeof(ls_ds): LinearSolvers.StandardLinearSolver{NonlinearSolvers.CNMatPC,NonlinearSolvers.CNPetscMatLO}
        # typeof(ls_ds.lo): NonlinearSolvers.CNPetscMatLO
        # typeof(ls_ds.lo.lo_inner): NonlinearSolvers.NewtonPetscMatLO
        # fieldnames(ls_ds): Symbol[:pc, :lo, :shared_mat, :comm, :myrank, 
        #                           :commsize, :ksp, :is_finalized, :reltol, :abstol, :dtol, :itermax]
        calcLinearOperator(ls_ds, mesh, sbp, eqn_nextstep, opts, ctx_residual, t)
        # this is calling: 
        # - calcLinearOperator(ls::LinearSolvers.StandardLinearSolver, 
        #                      mesh::ODLCommonTools.AbstractMesh, 
        #                      sbp::SummationByParts.AbstractSBP, 
        #                      eqn::ODLCommonTools.AbstractSolutionData, 
        #                      opts::Dict, ctx_residual, t; start_comm) 
        #   in LinearSolvers at /users/ashlea/.julia/v0.6/PDESolver/src/linearsolvers/ls_standard.jl:76
        # which has
        #   calcLinearOperator(ls.lo, mesh, sbp, eqn, opts, ctx_residual, t)
        # which is calling
        # - calcLinearOperator(lo::Union{NonlinearSolvers.CNDenseLO, 
        #                                NonlinearSolvers.CNPetscMatLO, NonlinearSolvers.CNSparseDirectLO}, 
        #                      mesh::ODLCommonTools.AbstractMesh, sbp::SummationByParts.AbstractSBP, 
        #                      eqn::ODLCommonTools.AbstractSolutionData, 
        #                      opts::Dict, ctx_residual, t) 
        #   in NonlinearSolvers at /users/ashlea/.julia/v0.6/PDESolver/src/NonlinearSolvers/crank_nicolson.jl:412
        # which has
        #   calcLinearOperator(lo.lo_inner, mesh, sbp, eqn, opts, ctx_residual, t)
        #   modifyJacCN(lo, mesh, sbp, eqn, opts, ctx_residual, t)
        # which is calling
        # - calcLinearOperator(lo::Union{NonlinearSolvers.NewtonDenseLO, NonlinearSolvers.NewtonPetscMatLO, 
        #                                NonlinearSolvers.NewtonSparseDirectLO}, 
        #                      mesh::ODLCommonTools.AbstractMesh, sbp::SummationByParts.AbstractSBP, 
        #                      eqn::ODLCommonTools.AbstractSolutionData, opts::Dict, ctx_residual, t) 
        #   in NonlinearSolvers at /users/ashlea/.julia/v0.6/PDESolver/src/NonlinearSolvers/newton_setup.jl:516
        # which has 
        #   lo2 = getBaseLO(lo)
        #   physicsJac(mesh, sbp, eqn, opts, lo2.A, ctx_residual, t)
        #
        # SO! this is modifying the Jac for CN.
        println(BSTDOUT, "  ---> Linear operator updated. i = ", i, ", t = ", t, " ==============\n")

        # Here was scratch content 2

        fill!(v_vec, 0.0)

        ### Same test as below, but only for non-Petsc
        # ccc = ones(Float64, mesh.numDof)
        # Accc = zeros(Float64, mesh.numDof)
        # A_mul_B!(Accc, lo_ds_innermost.A, ccc)
        # println(BSTDOUT, " vecnorm(Accc): ", vecnorm(Accc))
        # flush(BSTDOUT)

        # println(BSTDOUT, " cond(lo_ds_innermost.A): ", cond(full(lo_ds_innermost.A)))
        # flush(BSTDOUT)
        #=
        # This section is intended to assess the magnitude of A at every iteration 
        #   by contracting it with a vector of ones
        ccc = ones(Float64, mesh.numDof)

        for ix_dof = 1:mesh.numDof
          PETSC_insert[ix_dof] = ccc[ix_dof]
          PETSC_index[ix_dof] = ix_dof
        end
        # set values of ccc_petsc to ccc
        set_values1!(ccc_petsc, PETSC_index, PETSC_insert)
        A_mul_B!(Accc_petsc, lo_ds_innermost.A, ccc_petsc)
        println(BSTDOUT, " norm(ccc_petsc): ", norm(ccc_petsc))
        println(BSTDOUT, " norm(Accc_petsc): ", norm(Accc_petsc))
        =#


        if opts["stabilize_v"]

          # Recalculate dRdq
          filterDiagJac(mesh, opts, real(tmp_imag), clipJacData, stab_A, eigs_to_remove="neg")
          # filterDiagJac(mesh, opts, real(tmp_imag), clipJacData, stab_A, eigs_to_remove="pos")

          # TODO TODO
          # evalJacobianStrong(mesh, sbp, eqn, opts, stab_assembler, t)
          # Need to modify strong jacobian like we modify the CN Jac:
          # Modify

          # loop over blocks
          # blocksize is set above (during DiagJac init) as mesh.numDofPerNode*mesh.numNodesPerElement
          nblocks = size(stab_A.A, 3)       # third dimension of our block diag Jac is the block index
          ix_petsc_row = zeros(PetscInt, blocksize)
          ix_petsc_col = zeros(PetscInt, blocksize)
          block_to_add = zeros(PetscScalar, blocksize, blocksize)
          for block_ix = 1:nblocks

            for row_ix = 1:length(ix_petsc_row)
              # set the row indicies that we will insert into
              ix_petsc_row[row_ix] = blocksize*(block_ix-1)+row_ix
            end
            for col_ix = 1:length(ix_petsc_col)
              ix_petsc_col[col_ix] = blocksize*(block_ix-1)+col_ix
            end

            for row_ix = 1:length(ix_petsc_row)
              for col_ix = 1:length(ix_petsc_col)
                block_to_add[row_ix, col_ix] = stab_A.A[row_ix, col_ix, block_ix]
              end
            end

            # We should be subtracting, so we should scale block_to_add by -1.0
            scale!(block_to_add, -1.0)

            # now subtract the filtered DiagJac to the actual Jacobian, which will remove the positive eigenvalues of
            #   the strong Jacobian from the actual Jacobian
        
            # Add the negated block to the existing Jac inside the ls_ds LO object
            set_values1!(lo_ds_innermost.A, ix_petsc_row, ix_petsc_col, block_to_add, ADD_VALUES)

          end

          MatAssemblyBegin(lo_ds_innermost.A, MAT_FINAL_ASSEMBLY)
          MatAssemblyEnd(lo_ds_innermost.A, MAT_FINAL_ASSEMBLY)

        end   # end if opts["stabilize_v"]

        # linearSolve: solves Ax=b for x
        #   ls::StandardLinearSolver
        #   b::AbstractVector     -> RHS
        #   x::AbstractVector     -> what is solved for
        #   verbose=5
        #--------
        # linearSolve(ls_ds, b_vec, v_vec, verbose=5)
        linearSolve(ls_ds, b_vec, v_vec)

        #------------------------------------------------------------------------------
        # Checking linearSolve

        ### Only for julia sparse?
        # A_mul_B!(TEST_b_vec, lo_ds_innermost.A, v_vec)
        # println(BSTDOUT, " cond(full(lo_ds_innermost.A)): ", cond(full(lo_ds_innermost.A)))
        # println(BSTDOUT, " vecnorm(b_vec): ", vecnorm(b_vec))
        # println(BSTDOUT, " vecnorm(TEST_b_vec): ", vecnorm(TEST_b_vec))
        # println(BSTDOUT, "  >>> b_vec verify: ", vecnorm(TEST_b_vec - b_vec))
        # flush(BSTDOUT)

        #------------------------------------------------------------------------------


        #=
        # This section is extremely slow for Petsc. Don't use.
        println(BSTDOUT, " lo_ds_innermost.A:  vecnorm: ", vecnorm(lo_ds_innermost.A), 
                         "  max: ", maximum(lo_ds_innermost.A), 
                         "  min: ", minimum(lo_ds_innermost.A), "  minabs: ", minimum(abs.(lo_ds_innermost.A)))
        println(BSTDOUT, " b_vec:  vecnorm: ", vecnorm(b_vec), "  max: ", maximum(b_vec), 
                         "  min: ", minimum(b_vec), "  minabs: ", minimum(abs.(b_vec)))
        println(BSTDOUT, " v_vec:  vecnorm: ", vecnorm(v_vec), "  max: ", maximum(v_vec), 
                         "  min: ", minimum(v_vec), "  minabs: ", minimum(abs.(v_vec)))
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

      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # start of old DS code (the above is mostly new CSR stuff (except for the dDdu & term23 stuff))
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

        #=
        fill!(R_stab, 0.0)
        for j = 1:length(eqn.q_vec)
          # R_stab[j] = (q_vec[j] - old_q_vec[j])/(fac*delta_t)     # This was LSERK's
          R_stab[j] = (eqn.q_vec[j] - old_q_vec[j])/(dt)
        end
        fill!(v_energy, 0.0)
        for j = 1:length(eqn.q_vec)
          v_energy[j] = v_vec[j]*eqn.M[j]*imag(R_stab[j])/Ma_pert_mag
        end
        =#
        for ix_dof = 1:mesh.numDof

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

        # @mpi_master println(f_v_energy, i, "  ", real(v_energy_norm), "  ", real(v_energy_mean), "  ", real(v_energy_min), "  ", real(v_energy_max))
        @mpi_master println(f_v_energy_stageall, i, "  ", real(v_energy_norm))
        if (i % 500) == 0
          @mpi_master flush(f_v_energy_stageall)
        end

      end
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    end   # end of opts["perturb_Ma_CN"]

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

    flush(BSTDOUT)


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

        #=
        @mpi_master println(" eqn.params.Ma: ", eqn.params.Ma)
        @mpi_master println(" Ma_pert: ", Ma_pert_mag)
        eqn.params.Ma -= Ma_pert_mag      # need to remove perturbation now
        @mpi_master println(" pert removed from Ma")
        @mpi_master println(" eqn.params.Ma: ", eqn.params.Ma)
        =#

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
        # @mpi_master close(f_v_energy)
        @mpi_master close(f_v_energy_stageall)
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

DS_functions_file = string(Pkg.dir("PDESolver"),"/src/NonlinearSolvers/","crank_nicolson_ds-DS_functions.jl")
include(DS_functions_file)

