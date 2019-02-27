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
    println(BSTDOUT, "\nEntered Crank-Nicolson")
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
      v_vec = zeros(eqn.res_vec)

      ###### CSR: This section is the only section that contains things that are 
      #           new to the CN version of the complex step method (CSR: 'complex step residual')
      old_res_vec_Maimag = zeros(eqn.res_vec)    # corresponds to eqn
      new_res_vec_Maimag = zeros(eqn.res_vec)    # corresponds to eqn_nextstep
      res_hat_vec = zeros(eqn.res_vec)    # unsteady residual

      beforeDS_eqn_q_vec = zeros(eqn.q_vec)
      beforeDS_eqn_nextstep_q_vec = zeros(eqn.q_vec)
      beforeDS_eqn_res_vec = zeros(eqn.res_vec)
      beforeDS_eqn_nextstep_res_vec = zeros(eqn.res_vec)

      dRdM_vec = zeros(eqn.res_vec)
      b_vec = zeros(eqn.res_vec)
      dRdq_vn_prod = zeros(eqn.res_vec)

      #------------------------------------------------------------------------------
      # DS linear solver
      pc_ds, lo_ds = getCNDSPCandLO(mesh, sbp, eqn, opts)
      ls_ds = StandardLinearSolver(pc_ds, lo_ds, eqn.comm, opts)
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

        stab_A = DiagJac(Complex128, mesh.numDofPerNode*mesh.numNodesPerElement, mesh.numEl)
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

      v_vec = zeros(eqn.q_vec)      # direct sensitivity vector       This is at the IC
      # for v_ix = 1:length(v_vec)
        # v_vec[v_ix] = imag(eqn.q_vec[v_ix])/Ma_pert_mag
      # end
      fill!(v_vec, 0.0)       # v at IC is all 0's  ++++ new CS'ing R method (CSR) <- CSR comment on all lines new to CS'ing R method

      # q_vec_save_imag = Array{Float64}(length(eqn.q_vec))     # set up array for saving imag component of q (CSR)

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
    @mpi_master println(BSTDOUT, "\ni = ", i, ", t = ", t)
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
      println(BSTDOUT, "flushing convergence.dat to disk")
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

      if i != 1     # can't do F(q^(n+1)) + F(q^(n)) until we have both

        # Store both q_vec's & both res_vec's. Will be put back after all DS calcs.
        for ix_dof = 1:mesh.numDof
          beforeDS_eqn_q_vec[ix_dof] = eqn.q_vec[ix_dof]
          beforeDS_eqn_nextstep_q_vec[ix_dof] = eqn_nextstep.q_vec[ix_dof]

          beforeDS_eqn_res_vec[ix_dof] = eqn.res_vec[ix_dof]
          beforeDS_eqn_nextstep_res_vec[ix_dof] = eqn_nextstep.res_vec[ix_dof]
        end


        #------------------------------------------------------------------------------
        # This section calculates dR/dM

        Ma_pert_mag = opts["perturb_Ma_magnitude"]
        pert = complex(0, Ma_pert_mag)
        # eqn.params.Ma += pert
        eqn_nextstep.params.Ma += pert
        #=
        # DEBUG_CNTHES
        println("-------------- term23 debugging 0 -------------")
        println(" i: ", i)
        println(" eqn.params.Ma: ", eqn.params.Ma)
        println(" eqn_nextstep.params.Ma: ", eqn.params.Ma)
        println(" f: ", f)
        println("-----------------------------------------------")
        =#

        # now evalResidual to store into F(q^(n+1))
        f(mesh, sbp, eqn_nextstep, opts)      # F(q^(n+1)) -- with Ma perturbation, now in eqn_nextstep.res_vec

        for ix_dof = 1:mesh.numDof
          old_res_vec_Maimag[ix_dof] = new_res_vec_Maimag[ix_dof]        # F(q^(n)) -- with Ma perturbation
          new_res_vec_Maimag[ix_dof] = eqn_nextstep.res_vec[ix_dof]      # F(q^(n+1)) -- with Ma perturbation

          # form unsteady residual (res_hat) with F(q^(n)) and F(q^(n+1))
          res_hat_vec[ix_dof] = -0.5*eqn.Minv[ix_dof]*dt * (new_res_vec_Maimag[ix_dof] + old_res_vec_Maimag[ix_dof])
        end

        # obtain dR/dM using the complex step method
        for ix_dof = 1:mesh.numDof
          # dRdM_vec[ix_dof] = imag(dRdM_vec[ix_dof])/Ma_pert_mag     # TODO TODO, what am I doing here
          dRdM_vec[ix_dof] = imag(res_hat_vec[ix_dof])/Ma_pert_mag     # should this be res_hat_vec??
        end

        # eqn.params.Ma -= pert
        eqn_nextstep.params.Ma -= pert
        #------------------------------------------------------------------------------

        #------------------------------------------------------------------------------
        # This section calculates dR/dq * v^(n)

        # v_vec currently holds v at timestep n: v^(n)

        # need to add epsilon*v_vec*im (evi) to q_vec, then re-evaluate res_hat at the updated q
        # so: need to add evi to old_q_vec too. Because to form res_hat, need F(q^(n)) and F(q^(n+1))
        for ix_dof = 1:mesh.numDof
          eqn.q_vec[ix_dof]          += Ma_pert_mag*im*v_vec[ix_dof]
          eqn_nextstep.q_vec[ix_dof] += Ma_pert_mag*im*v_vec[ix_dof]
        end
        #=
        # DEBUG_CNTHES
        println("-------------- term23 debugging 1 -------------")
        println(" i: ", i)
        println(" eqn.params.Ma: ", eqn.params.Ma)
        println(" vecnorm(real(eqn.q_vec)): ", vecnorm(real(eqn.q_vec)))
        println(" vecnorm(imag(eqn.q_vec)): ", vecnorm(imag(eqn.q_vec)))
        println(" vecnorm(real(eqn_nextstep.q_vec)): ", vecnorm(real(eqn_nextstep.q_vec)))
        println(" vecnorm(imag(eqn_nextstep.q_vec)): ", vecnorm(imag(eqn_nextstep.q_vec)))
        println("-----------------------------------------------")
        =#


        f(mesh, sbp, eqn, opts)             # F(q^(n) + evi) now in eqn.res_vec
        f(mesh, sbp, eqn_nextstep, opts)    # F(q^(n+1) + evi) now in eqn_nextstep.res_vec

        for ix_dof = 1:mesh.numDof
          
          # form unsteady residual (res_hat) with F(q^(n) + evi) and F(q^(n+1) + evi)
          res_hat_vec[ix_dof] = -0.5*eqn.Minv[ix_dof]*dt * (eqn_nextstep.res_vec[ix_dof] + eqn.res_vec[ix_dof])

          # calc dRdq * v^(n) by doing matrix-free complex step
          dRdq_vn_prod[ix_dof] = imag(res_hat_vec[ix_dof])/Ma_pert_mag

          # remove the imaginary component from q used for matrix_dof-free product    # TODO: necessary?
          # eqn.q_vec[ix_dof] = complex(real(eqn.q_vec[ix_dof]))
          # eqn_nextstep.q_vec[ix_dof] = complex(real(eqn_nextstep.q_vec[ix_dof]))

          # combine (-dRdM - dRdq * v^(n)) into b
          b_vec[ix_dof] = - dRdM_vec[ix_dof] - dRdq_vn_prod[ix_dof]

        end

        #------------------------------------------------------------------------------
        # Now the calculation of v_ix at n+1
        # Solve:
        #   dR/dq^(n+1) v^(n+1) = b
        #--------

        fill!(v_vec, 0.0)

        # Recalculate dRdq
        calcLinearOperator(ls_ds, mesh, sbp, eqn, opts, ctx_residual, t)

        # HERE is where we stabilize?
        # evalJacobianStrong(mesh, sbp, eqn, opts, stab_assembler, t)
        # filterDiagJac(mesh, opts, real(tmp_imag), clipJacData, stab_A)

        # linearSolve: solves Ax=b for x
        #   ls::StandardLinearSolver
        #   b::AbstractVector     -> RHS
        #   x::AbstractVector     -> what is solved for
        #   verbose=5
        linearSolve(ls_ds, b_vec, v_vec)      
        #--------

        #------------------------------------------------------------------------------
        # Restore q_vec & res_vec for both eqn & eqn_nextstep, as they were
        #   before the DS calcs
        for ix_dof = 1:mesh.numDof
          eqn.q_vec[ix_dof] = beforeDS_eqn_q_vec[ix_dof]
          eqn_nextstep.q_vec[ix_dof] = beforeDS_eqn_nextstep_q_vec[ix_dof]

          eqn.res_vec[ix_dof] = beforeDS_eqn_res_vec[ix_dof]
          eqn_nextstep.res_vec[ix_dof] = beforeDS_eqn_nextstep_res_vec[ix_dof]
        end

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
        #=  
        # DEBUG_CNTHES
        println("-------------- term23 debugging 2 -------------")
        println(" i: ", i)
        println(" vecnorm(b_vec): ", vecnorm(b_vec))
        println(" vecnorm(dRdM_vec): ", vecnorm(dRdM_vec))
        println(" vecnorm(dRdq_vn_prod): ", vecnorm(dRdq_vn_prod))
        println(" vecnorm(real(res_hat_vec)): ", vecnorm(real(res_hat_vec)))
        println(" vecnorm(imag(res_hat_vec)): ", vecnorm(imag(res_hat_vec)))
        println(" ")
        println(" quad_weight: ", quad_weight)
        println(" vecnorm(dDdu_vec): ", vecnorm(dDdu_vec))
        println(" vecnorm(v_vec): ", vecnorm(v_vec))
        println(" term23: ", term23)
        # println("  term23 change: ", (term23 - old_term23)*1.0/dt)
        println("-----------------------------------------------")
        =#

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
        fill!(R_stab, 0.0)
        for j = 1:length(eqn.q_vec)
          # R_stab[j] = (q_vec[j] - old_q_vec[j])/(fac*delta_t)     # This was LSERK's
          R_stab[j] = (eqn.q_vec[j] - old_q_vec[j])/(dt)
        end
        fill!(v_energy, 0.0)
        for j = 1:length(eqn.q_vec)
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
    if (res_norm < res_tol)  # ???
      if myrank == 0
        println(BSTDOUT, "breaking due to res_tol, res norm = $res_norm")
        close(f1)
        flush(BSTDOUT)
      end
      break       # EXIT condition
    end


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

