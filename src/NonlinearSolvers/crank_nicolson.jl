# crank_nicolson.jl
# Crank-Nicolson implicit solver for PDEs

export crank_nicolson, cnResidual, negativeZeroCheck

push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/solver/euler"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/NonlinearSolvers"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Debugging"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Utils"))

"""
crank_nicolson

  This function performs Crank-Nicolson implicit time solution, using a function 
  of the form du/dt = f(u, t)
  
  Arguments:
    * physics_func  : function evaluation, must have signature (ctx..., opts, t)
    * h  : time step size
    * t_max: time value to stop time stepping (time starts at 0)
    * q_vec: vector of the u values
    * res_vec: vector of du/dt values (the output of the function physics_func)
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
    * neg_time : step through time in the negative direction,
                 starting at t_max, stepping with h, and ending at 0.0.


   The eqn.q_vec should hold the whichever variables (conservative or
   entropy) that the simulation should use.

   The idea behind pre_func, post func, and ctx is that they abstract what kind
   of system CN is timestepping. CN only needs to know about q_vec and
   res_vec.

   For physics modules, ctx should be (mesh, sbp, eqn) and q_vec and res_vec 
   should be eqn.q_vec and eqn.res_vec.

   TODO: fully document eqn/eqn_nextstep
"""
# function crank_nicolson{Tmsh, Tsol}(physics_func::Function, h::AbstractFloat, t_max::AbstractFloat,
                        # mesh::AbstractMesh{Tmsh}, sbp::AbstractSBP, eqn::AbstractSolutionData{Tsol},
                        # opts, res_tol=-1.0; neg_time=false, obj_fn=obj_zero, store_u_to_disk=false)
function crank_nicolson{Tmsh, Tsol}(physics_func::Function, h::AbstractFloat, t_max::AbstractFloat,
                        mesh::AbstractMesh{Tmsh}, sbp::AbstractSBP, eqn::AbstractSolutionData{Tsol},
                        opts, 
                        WWW, ZZZ, dRdu_global_fwd, dRdu_global_rev, dRdu_global_rev_PM,
                        res_tol=-1.0; neg_time=false, obj_fn=obj_zero, store_u_to_disk=false)
                        # NEWNEW: neg_time, obj_fn, store_u_to_disk
  #----------------------------------------------------------------------
  if opts["uadj_global"]
    println(" GLOBAL, CN start: size(dRdu_global_fwd): ", size(dRdu_global_fwd))
    println(" GLOBAL, CN start: size(dRdu_global_rev): ", size(dRdu_global_rev))
    println(" GLOBAL, CN start: size(WWW): ", size(WWW))
    println(" GLOBAL, CN start: size(ZZZ): ", size(ZZZ))
  end

  myrank = MPI.Comm_rank(MPI.COMM_WORLD)
  fstdout = BufferedIO(STDOUT)
  if myrank == 0
    println(fstdout, "\nEntered Crank-Nicolson")
    println(fstdout, "res_tol = ", res_tol)
  end

  output_freq = opts["output_freq"]::Int
  write_vis = opts["write_vis"]::Bool
  use_itermax = opts["use_itermax"]::Bool
  jac_type = opts["jac_type"]
  if use_itermax
    itermax = opts["itermax"]
  end

  if jac_type == 4
    error("CN not implemented for matrix-free ops. (jac_type cannot be 4)")
  end

  if myrank == 0
    _f1 = open("convergence.dat", "a+")
    f1 = BufferedIO(_f1)
  end
 
  # calculate t_steps, the number of time steps that CN will take
  t_steps = floor(Int, t_max/h)   # this allows t_max to be any value up to h greater the final time step
  # TODO: Ideally should be t_steps*h for clarity, issue #92. Should be fixed for RK4 also.
  time_of_final_step = (t_steps-1)*h

  println("=========== t_steps: ", t_steps, " =========")

  if neg_time == false    # negative time is for unsteady adjoint
    t = 0.0     # start time at 0.0
  else    
    t = time_of_final_step
  end

  # NOTE: WWW, ZZZ, dRdu_global_fwd & rev formed in initialization.jl

  if neg_time == false
    # make a copy of the eqn object for storage of t_(n+1) information
    eqn_nextstep = eqn_deepcopy(eqn, mesh, sbp, opts)
    # TODO: copyForMultistage does not give correct values.
    #     deepcopy works for now, but uses more memory than copyForMultistage, if it worked
    # eqn_nextstep = copyForMultistage(eqn)
    # eqn_nextstep.q = reshape(eqn_nextstep.q_vec, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)
    # eqn_nextstep.res = reshape(eqn_nextstep.res_vec, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)
  else
    # Initialize adjoint pde eqn objects
    adj = eqn_deepcopy(eqn, mesh, sbp, opts)
    # adj.q = reshape(adj.q_vec, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)
    # adj.res = reshape(adj.res_vec, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)

    adj_nextstep = eqn_deepcopy(eqn, mesh, sbp, opts)
    # adj_nextstep.q = reshape(adj_nextstep.q_vec, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)
    # adj_nextstep.res = reshape(adj_nextstep.res_vec, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)
  end

  println("============ Starting CN ============")

  #-------------------------------------------------------------------------------
  # allocate Jac outside of time-stepping loop
  # these jacs are for full CN or CN adj jac
  println("================== neg_time: ", neg_time, " .... t = ", t, " .... t_max = ", t_max, "  ===================")
  if neg_time == false
    newton_data, jac, rhs_vec = setupNewton(mesh, mesh, sbp, eqn, opts, physics_func)
  else
    newton_data, jac, rhs_vec = setupNewton(mesh, mesh, sbp, adj, opts, physics_func)
  end

  # Setting IC for reverse sweep
  if neg_time == true

    #----------------
    # this section:
    #   1) reads the checkpointed q_vec at the last time step of the forward sweep (n'th time step)
    #   2) uses calcJacobianComplex to calculate dRdu at time step n
    i_fwd = t_steps + 2  # index of the last time step's checkpoint. 
                         #  It was saved after the end of the time stepping loop during the forward sweep

    # load checkpoint to calculate dRdu at this time step
    println("Setting IC for reverse sweep, i_fwd (forward sweep time step index): ", i_fwd)
    # eqn_fwd = cnAdjLoadChkpt(mesh, sbp, opts, adj, physics_func, i_fwd, t)
    eqn_fwd = cnAdjLoadChkpt(mesh, sbp, opts, adj, physics_func, i_fwd)
    check_q_qvec_consistency(mesh, sbp, eqn_fwd, opts)

    jac = cnAdjCalcdRdu(mesh, sbp, opts, eqn_fwd, physics_func, i_fwd, t)
    dRdu_n = jac      # TODO: check transpose
    println(" size of dRdu: ", size(dRdu_n))
    #----------------

    dJdu_CS = calcdJdu_CS(mesh, sbp, eqn_fwd, opts)  # obtain dJdu at time step n
    dJdu_FD = calcdJdu_FD(mesh, sbp, eqn_fwd, opts)  # obtain dJdu at time step n
    dJdu = dJdu_CS
    dJdu_analytical = calcObjectiveFn(mesh, sbp, eqn_fwd, opts, isDeriv=true)

    # println(" dJdu_CS: ", dJdu)
    # println(" dJdu_CS - dJdu_analytical: ", dJdu)
    writedlm("dJdu_IC_CS.dat", dJdu_CS)
    writedlm("dJdu_IC_analytical.dat", reshape(dJdu_analytical, (mesh.numDof, 1)))

    # now that dRdu and dJdu at time step n has been obtained, we can now set the IC for the adjoint eqn
    I = eye(length(eqn_fwd.q_vec))
    B = (I - (h/2) * (dRdu_n))
    psi = transpose(B)\(-dJdu)
    adj.q_vec = copy(psi)
    disassembleSolution(mesh, sbp, adj, opts, adj.q, adj.q_vec)

  end

  #--------------------------------------------------------------------------------------------------------------------------
  # Start of CN loop
  #--------------------------------------------------------------------------------------------------------------------------
  # TODO: Ideally should be t_steps + 2 for clarity, issue #92. Should be fixed for RK4 also.
  t_steps_end = t_steps + 1
  for i = 2:t_steps_end

    println(" ")
    println(" -------------- start of this CN iter. t = $t, i = $i, neg_time = ", neg_time, " --------------")
    @debug1 flush(eqn.params.f)

    println("   eqn.params.alpha_x: ",eqn.params.alpha_x)
    println("   eqn.params.alpha_y: ",eqn.params.alpha_y)

    # Write checkpoint data at start of time step
    if neg_time == false
      # for adjoint_straight option: stores every time step's q to disk
      # Note: cannot store full eqn object without extending one of the julia write methods
      if store_u_to_disk == true
        filename = string("qvec_for_adj-", i, ".dat")
        writedlm(filename, eqn.q_vec)
        vis_filename = string("solution_storedtodisk_i-", i)
        saveSolutionToMesh(mesh, real(eqn.q_vec))
        writeVisFiles(mesh, vis_filename)

        time_filename = string("t-",i,".dat")
        writedlm(time_filename, t)              # make sure that this time value, whether t or t_nextstep,
                                                #   correctly corresponds to the q_vec stored just above
      end
    else
      # save every time step's adjoint to disk
      if opts["adjoint_saveall"]
        filename = string("adj-", i, ".dat")
        writedlm(filename, adj.q_vec)

        vis_filename = string("adj_i-", i)
        saveSolutionToMesh(mesh, real(adj.q_vec))
        writeVisFiles(mesh, vis_filename)
      end
    end

    if neg_time == false
      println(" -------------- eqn.q_vec start of this CN iter. t = $t, i = $i --------------")
      # print_qvec_coords(mesh, sbp, eqn, opts)
      println(" -------------- J for this eqn.q_vec --------------")
      J_arr = calcObjectiveFn(mesh, sbp, eqn, opts)
      # println(" -------------- eqn.q_bndry for this eqn.q_vec --------------")
      # print_qvec_coords(mesh, sbp, eqn, opts; bndry=true)
      J = J_arr[1]
      println("  J: ", J)
      println(" -------------- Now using eqn.q_vec to compute eqn_nextstep.q_vec --------------")
    else
      println(" -------------- adj.q_vec start of this CN iter. t = $t, i = $i --------------")
      # print_qvec_coords(mesh, sbp, adj, opts)
      println(" -------------- Now using adj.q_vec to compute adj_nextstep.q_vec --------------")
    end
    println(" ")

    #----------------------------
    # zero out Jac
    #   this works for both PETSc and Julia matrices.
    #   when jac is a Julia matrix, this is effectively wrapping: fill!(jac, 0.0)
    PetscMatZeroEntries(jac)

    # NOTE:
    # majorIterationCallback: called before every step of Newton's method. signature: 
    #   majorIterationCallback(itr, mesh::AbstractMesh, sbp::AbstractSBP, eqn::AbstractEulerData, opts)

    # NOTE: Must include a comma in the ctx tuple to indicate tuple
    # f is the physics function, like evalEuler

    #-----------------------------------------------------------------------
    # time step update: h is passed in as argument to crank_nicolson
    if neg_time == false
      # need to add h in the forward time usage
      t_nextstep = t + h
    else
      # need to subtract h in the reverse time usage
      t_nextstep = t - h
      t_nextstep = negativeZeroCheck(t_nextstep)   # ensure negative zero is changed to zero
    end

    #-------------
    # objective function section
    # 1. read option to indicate which obj fun
    # 2. call it, complex step it, and store it in dJdu
    if neg_time == true
      dJdu = zeros(Tsol, length(eqn.q_vec))
      dJdu_CS = calcdJdu_CS(mesh, sbp, eqn_fwd, opts)  # obtain dJdu at time step n
      dJdu_FD = calcdJdu_FD(mesh, sbp, eqn_fwd, opts)  # obtain dJdu at time step n
      dJdu = dJdu_CS
      dJdu_analytical = calcObjectiveFn(mesh, sbp, eqn_fwd, opts, isDeriv=true)

      filename = string("dJdu_",i,"_CS.dat")
      writedlm(filename, dJdu_CS)
      filename = string("dJdu_",i,"_analytical.dat")
      writedlm(filename, reshape(dJdu_analytical, (mesh.numDof, 1)))

      println("       checking direct method: size(dJdu): ", size(dJdu))

      # VV is the algebraic v, which is dudA, calculated for the advection adjoint test

      # TODO: figure out if calcVV is to be called with t or t_nextstep. I think it's t_nextstep. 
      #       It is for sure t_nextstep when dJdu is called during the direct solve.
      # VV = calcVV(mesh, sbp, adj, opts, t_nextstep)     # scalar only because our x-value of interest is unchanging
      VV = calcVV(mesh, sbp, adj, opts, t)     # scalar only because our x-value of interest is unchanging
      # NOTE 20170712: think I need to be doing t, not t_nextstep. this CN rev's i corresponds to eqn_fwd's i_fwd and t

      # println("       checking direct method: VV: ", VV)
      println("       checking direct method: size(VV): ", size(VV))
      check_directmethod = transpose(dJdu)*VV
      filename = string("check_directmethod-", i, ".dat")
      writedlm(filename, check_directmethod)

    end

    if neg_time == false
      ctx_residual = (physics_func, eqn, h, newton_data)
    else
      # i is the time step index in the reverse sweep.
      #   It moves forward from 2 to (t_steps+1) even though the adjoint is going backwards in time.
      #   The adjoint calculation requires data from the forward sweep to be loaded from disk.
      #   At a given time step in the adjoint solve, although the index i corresponds to this
      #     loop's time step, the index i does not correspond to the same i of the forward sweep.
      #   The adjustment is not just (t_steps - i) because the loop starts at 2 and ends at t_steps + 1.
      i_fwd = (t_steps + 1) - (i - 2)
      # TODO: eventual fix, issue #92, t_steps + 2. needs to be done in RK4 also
      println(" time step variables-  i: ", i, "  i_fwd: ", i_fwd, "  t_steps: ", t_steps)
      ctx_residual = (physics_func, adj, h, newton_data, i_fwd, dJdu)

      #-------------------------------------------------
      # Output times & indices
      println(" --- CN times and indices ---")
      println("  CN: i = $i")
      println("  CN: i_fwd = $i_fwd")
      println("  CN: t = $t")
      println("  CN: t_nextstep = $t_nextstep")
      #-------------------------------------------------

    end

    if neg_time == false
      # allow for user to select CN's internal Newton's method. Only supports dense FD Jacs, so only for debugging
      if opts["cleansheet_CN_newton"]
        # cnNewton: in cnNewton.jl
        cnNewton(mesh, sbp, opts, h, physics_func, eqn, eqn_nextstep, t)
      else
        @time newtonInner(newton_data, mesh, sbp, eqn_nextstep, opts, cnRhs, cnJac, jac, rhs_vec, ctx_residual, t)

        if opts["uadj_global"]
          ### dRdu check: fwd
          println(" GLOBAL: forming global dRdu, i = $i")
          # first, let's try forming the forward global dRdu without the state perturbed by WWW.
          blksz = mesh.numDof
          # blksz = 3 # testing 44
          row_ix_start = (i-2)*blksz + 1
          row_ix_end = (i-2)*blksz + blksz
          col_ix_start = (i-2)*blksz + 1
          col_ix_end = (i-2)*blksz + blksz
          println(" GLOBAL: diagonal block ix's: [",row_ix_start,":",row_ix_end,", ",col_ix_start,":", col_ix_end,"]")

          II = eye(blksz)

          # TODO: check eqn or eqn_nextstep
          newton_data_discard, jac_for_dRdu_global, rhs_vec_discard = setupNewton(mesh, mesh, sbp, eqn_nextstep, opts, physics_func)
          assert(opts["jac_method"] == 2)
          epsilon = opts["epsilon"]
          pert = complex(0, epsilon)
          calcJacobianComplex(newton_data_discard, mesh, sbp, eqn_nextstep, opts, physics_func, pert, jac_for_dRdu_global, t)
          # TODO check that this jacobian equals the one inside newtonInner
          println(" GLOBAL: norm of jac after newtonInner call: ", norm(jac_for_dRdu_global))

          # blk = II - 0.5*h*jac
          # blk = II - 0.5*h*jac[1:3,1:3]
          # blk = II - 0.5*h*jac_for_dRdu_global[1:3, 1:3]
          blk = II - 0.5*h*jac_for_dRdu_global
          # blk = ones(blksz, blksz)*44     # testing 44
          println(" GLOBAL: size of blk: ", size(blk))

          # this time step's actual portion of the global dRdu is the cnJac. so it has to be I-0.5*h*physicsJac
          jac_filename = string("jac_fwd_cJc_t-",t,".dat")
          writedlm(jac_filename, round(real(blk), 4))

          println(" GLOBAL: size(dRdu_global_fwd): ", size(dRdu_global_fwd))
          dRdu_global_fwd[row_ix_start:row_ix_end, col_ix_start:col_ix_end] = blk

          if i != t_steps + 1     # final time step has only one block
            row_ix_start = (i-2)*blksz + blksz + 1
            row_ix_end = (i-2)*blksz + 2*blksz
            col_ix_start = (i-2)*blksz + 1
            col_ix_end = (i-2)*blksz + blksz
            println(" GLOBAL: off-diagonal block ix's: [",row_ix_start,":",row_ix_end,", ",col_ix_start,":", col_ix_end,"]")

            # blk = -1.0*II - 0.5*h*jac
            # blk = -1.0*II - 0.5*h*jac[1:3,1:3]
            # blk = -1.0*II - 0.5*h*jac_for_dRdu_global[1:3, 1:3]
            blk = -1.0*II - 0.5*h*jac_for_dRdu_global
            # blk = ones(blksz, blksz)*55     # testing 44

            dRdu_global_fwd[row_ix_start:row_ix_end, col_ix_start:col_ix_end] = blk
          end   # end of "if != t_steps + 1"
        end # end of "if opts["uadj_global"]

      end   # end of else clause of "if opts["cleansheet_CN_newton"]"

      # NOTE: this has been commented out to remove Newton from the adjoint testing. Currently direct solving.
    # else      # call newtonInner using cnAdjJac and cnAdjRhs
      # @time newtonInner(newton_data, mesh, sbp, adj_nextstep, opts, cnAdjRhs, cnAdjJac, jac, rhs_vec, ctx_residual, t)

    else    # else clause of "if neg_time == false"

      # direct solve for psi_i
      # note: this needs to be called with t, not t_nextstep. 
      #   within cnAdjDirect, t_nextstep is calculated and used throughout.
      # (adj_nextstep.q_vec, jac) = cnAdjDirect(mesh, sbp, opts, adj, physics_func, jac, i_fwd, h, t)
      (adj_nextstep.q_vec, jac) = cnAdjDirect(mesh, sbp, opts, adj, physics_func, jac, i_fwd, h, t_steps, t, dRdu_global_rev)
      disassembleSolution(mesh, sbp, adj_nextstep, opts, adj_nextstep.q, adj_nextstep.q_vec)

    end

    # TODO TODO 20170711: dRdA timing: is this being contracted in the right time step?
    # advection adjoint check
    if neg_time == true
      # dRdA_CS = calcdRdA_CS(mesh, sbp, adj_nextstep, opts, t_nextstep)
      # dRdA_FD = calcdRdA_FD(mesh, sbp, adj_nextstep, opts, t_nextstep)
      dRdA_CS = calcdRdA_CS(mesh, sbp, eqn_fwd, opts, t)
      dRdA_FD = calcdRdA_FD(mesh, sbp, eqn_fwd, opts, t)
      dRdA = dRdA_CS
      # check_adjointmethod = transpose(adj_nextstep.q_vec)*dRdA
      # NOTE 20170712: think I need to be doing adj, not adj_nextstep. this CN rev's i corresponds to eqn_fwd's i_fwd and t
      # check_adjointmethod = transpose(adj.q_vec)*dRdA
      check_adjointmethod = transpose(adj.q_vec)*(-1.0*dRdA)
      filename = string("check_adjointmethod-", i, ".dat")
      writedlm(filename, check_adjointmethod)
    end

    # This allows the solution to be updated from _nextstep without a deepcopy.
    #   There are two memory locations used by eqn & eqn_nextstep, 
    #   and this flips the location of eqn & eqn_nextstep every time step
    if neg_time == false
      eqn_temp = eqn
      eqn = eqn_nextstep
      eqn_nextstep = eqn_temp
    else
      adj_temp = adj
      adj = adj_nextstep
      adj_nextstep = adj_temp
    end


    # adj.q_vec now contains the adjoint at time step i. 
    # Previously, adj_nextstep corresponded to time step i, and adj corresponded to time step i+1.

    #------- adjoint check
    # TODO this is wrong, what is this. seems like a first attempt at check_adjointmethod
    if neg_time == true
      # eqn_fwd = cnAdjLoadChkpt(mesh, sbp, opts, adj, physics_func, i_fwd, t)     # t has been updated, so no t_nextstep
      # jac = cnAdjCalcdRdu(mesh, sbp, opts, eqn_fwd, physics_func, t)
      omega = 1.0
      v_bc = sin(omega*t)
      dRdA = -1.0*jac*v_bc
      dJdA = transpose(adj.q_vec)*dRdA
      filename = string("dJdA_check-", i, ".dat")
      writedlm(filename, transpose(dJdA))
      # if i == 10
        # println("v_bc: ", v_bc)
        # println("jac: ", jac)
        # println("dRdA: ", dRdA)
        # println("adj.q_vec: ", adj.q_vec)
      # end
    end

    # Note: we now need to copy the updated q over for the initial newton guess
    if neg_time == false
      for dof_ix = 1:mesh.numDof
        eqn_nextstep.q_vec[dof_ix] = eqn.q_vec[dof_ix]
      end
      disassembleSolution(mesh, sbp, eqn_nextstep, opts, eqn_nextstep.q, eqn_nextstep.q_vec)
    else
      for dof_ix = 1:mesh.numDof
        adj_nextstep.q_vec[dof_ix] = adj.q_vec[dof_ix]
      end
      disassembleSolution(mesh, sbp, adj_nextstep, opts, adj_nextstep.q, adj_nextstep.q_vec)
    end


    t = t_nextstep        # update time step


  end   # end of t step loop

  i = i+1

  # Checkpoint final time step
  if neg_time == false
    # for adjoint_straight option: stores every time step's q to disk
    # Note: cannot store full eqn object without extending one of the julia write methods
    if store_u_to_disk == true
      filename = string("qvec_for_adj-", i, ".dat")
      writedlm(filename, eqn.q_vec)
      vis_filename = string("solution_storedtodisk_i-", i)
      saveSolutionToMesh(mesh, real(eqn.q_vec))
      writeVisFiles(mesh, vis_filename)

      time_filename = string("t-",i,".dat")
      writedlm(time_filename, t)              # make sure that this time value, whether t or t_nextstep,
                                              #   correctly corresponds to the q_vec stored just above
    end
  else
    # save every time step's adjoint to disk. This is for the final time step
    if opts["adjoint_saveall"]
      filename = string("adj-", i, ".dat")
      writedlm(filename, adj.q_vec)

      vis_filename = string("adj_i-", i)
      saveSolutionToMesh(mesh, real(adj.q_vec))
      writeVisFiles(mesh, vis_filename)
    end
  end

  if neg_time == false
    println(" -------------- eqn.q_vec of last time step. after CN loop. t = $t, i = $i --------------")
    # print_qvec_coords(mesh, sbp, eqn, opts)
    println(" -------------- End of CN --------------")
  else
    println(" -------------- adj.q_vec of last time step. after CN loop. t = $t, i = $i --------------")
    # print_qvec_coords(mesh, sbp, adj, opts)
    println(" -------------- End of CN --------------")
  end
  println(" ")


  #=
  # new: writing checkpoint data at start of time step
  if neg_time == false
    # for adjoint_straight option: stores every time step's q to disk
    # Note: cannot store full eqn object without extending one of the julia write methods
    if store_u_to_disk == true
      filename = string("qvec_for_adj-", i, ".dat")
      writedlm(filename, eqn.q_vec)
      vis_filename = string("solution_storedtodisk_i-", i)
      saveSolutionToMesh(mesh, real(eqn.q_vec))
      writeVisFiles(mesh, vis_filename)

      time_filename = string("t-",i,".dat")
      writedlm(time_filename, t)              # make sure that this time value, whether t or t_nextstep,
                                              #   correctly corresponds to the q_vec stored just above
    end
  else
    # save every time step's adjoint to disk
    if opts["adjoint_saveall"]
      filename = string("adj-", i, ".dat")
      writedlm(filename, adj.q_vec)
    end
  end
  =#

  if opts["uadj_global"]
    ### GLOBAL: saving global dRdu
    writedlm("global_dRdu_fwd.dat", dRdu_global_fwd)
    writedlm("global_dRdu_rev.dat", dRdu_global_rev)
  end

  # depending on how many timesteps we do, this may or may not be necessary
  #   usage: copy!(dest, src)   
  # this copy! is defined in ODLCommonTools
  if neg_time == false
    copy!(eqn, eqn_nextstep)      # copying eqn_nextstep to eqn
    writedlm("solution_final_inCN.dat", real(eqn.q_vec))
  else
    copy!(adj, adj_nextstep)      # copying adj_nextstep to eqn
    writedlm("adjoint_final_inCN.dat", real(adj.q_vec))
  end

  #=
  # For debugging: store integer set to mesh and save, just to visualize ordering of q_vec
  for i = 1:length(eqn.q_vec)
    eqn.q_vec[i] = convert(Float64, i)
  end
  filename = string("qvec_integers.dat")
  writedlm(filename, eqn.q_vec)
  vis_filename = string("qvec_integers")
  saveSolutionToMesh(mesh, real(eqn.q_vec))
  writeVisFiles(mesh, vis_filename)
  print_qvec_coords(mesh, sbp, eqn, opts)
  # End debugging integer ordering output
  =#

  if jac_type == 3      # if jac is a Petsc matrix, it needs to be freed when we're done using it
    # contents of ctx_newton: (jacp, x, b, ksp)
    NonlinearSolvers.destroyPetsc(jac, newton_data.ctx_newton...)
  end

  @debug1 println("============= end of CN: t = $t ===============")

  return t

end   # end of crank_nicolson function

function negativeZeroCheck(t)
  # prevent floating point errors from setting t to a negative number at t = 0.0
  if abs(t) < 1e-14 && t < 0.0
    println("Barely negative time detected (-1e-14 < t < 0.0) , setting t = 0.0.")
    t = 0.0
  end
  return t
end

function check_q_qvec_consistency(mesh, sbp, eqn, opts)

  if eqn.q != reshape(eqn.q_vec, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)
    error("q & q_vec consistency error")
  end

  return nothing

end
