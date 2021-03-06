# crank_nicolson.jl
# Crank-Nicolson implicit solver for PDEs

export crank_nicolson_uadj, negativeZeroCheck

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
# function crank_nicolson(physics_func::Function, h::AbstractFloat, t_max::AbstractFloat,
                        # mesh::AbstractMesh{Tmsh}, sbp::AbstractSBP, eqn::AbstractSolutionData{Tsol},
                        # opts, res_tol=-1.0; neg_time=false, obj_fn=obj_zero, store_u_to_disk=false) where {Tmsh, Tsol}
function crank_nicolson_uadj(physics_func::Function, h::AbstractFloat, t_max::AbstractFloat,
                        mesh::AbstractMesh{Tmsh}, sbp::AbstractSBP, eqn::AbstractSolutionData{Tsol},
                        opts, 
                        WWW, ZZZ, dRdu_global_fwd, dRdu_global_rev, dRdu_global_rev_PM,
                        res_tol=-1.0; neg_time=false, obj_fn=obj_zero, store_u_to_disk=false) where {Tmsh, Tsol}
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
  # time_of_final_step = (t_steps-1)*h
  time_of_final_step = (t_steps)*h
  # changed from (t_steps-1)*h, 20170908

  if neg_time == false    # negative time is for unsteady adjoint
    t = 0.0     # start time at 0.0
  else    
    t = time_of_final_step
    # opts["use_Minv_override_for_uadj"] = true
  end

  # NOTE: WWW, ZZZ, dRdu_global_fwd & rev formed in initialization.jl

  if neg_time == false
    # make a copy of the eqn object for storage of t_(n+1) information
    eqn_nextstep = eqn_deepcopy(mesh, sbp, eqn, opts)
    # TODO: copyForMultistage! does not give correct values.
    #     deepcopy works for now, but uses more memory than copyForMultistage!, if it worked
    # eqn_nextstep = copyForMultistage!(eqn)
    # eqn_nextstep.q = reshape(eqn_nextstep.q_vec, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)
    # eqn_nextstep.res = reshape(eqn_nextstep.res_vec, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)
  else
    # Initialize adjoint pde eqn objects
    adj = eqn_deepcopy(mesh, sbp, eqn, opts)
    # adj.q = reshape(adj.q_vec, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)
    # adj.res = reshape(adj.res_vec, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)

    adj_nextstep = eqn_deepcopy(mesh, sbp, eqn, opts)
    # adj_nextstep.q = reshape(adj_nextstep.q_vec, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)
    # adj_nextstep.res = reshape(adj_nextstep.res_vec, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)
  end

  println("========================= Starting CN =========================")
  println("   neg_time: ", neg_time)
  println("   t = ", round(t,3))
  println("   t_max = ", t_max)
  println("   t_steps = ", t_steps)
  println("   eqn.params.alpha_x: ",eqn.params.alpha_x)
  println("   eqn.params.alpha_y: ",eqn.params.alpha_y)


  #-------------------------------------------------------------------------------
  # allocate Jac outside of time-stepping loop
  # these jacs are for full CN or CN adj jac
  if neg_time == false
    newton_data, jac, rhs_vec = setupNewton(mesh, mesh, sbp, eqn, opts)
  else
    newton_data, jac, rhs_vec = setupNewton(mesh, mesh, sbp, adj, opts)
  end

  # Setting IC for reverse sweep
  if neg_time == true

    #----------------
    # this section:
    #   1) reads the checkpointed q_vec at the last time step of the forward sweep (n'th time step)
    #   2) uses calcJacobianComplex to calculate dRdu at time step n
    i_fwd = t_steps + 2  # index of the last time step's checkpoint. 
                         #  It was saved after the end of the time stepping loop during the forward sweep

    println("\n--- Adj IC: setting i_fwd. ---")
    println("     t_steps = ", t_steps)
    println("     i_fwd = ", i_fwd)
    println("     t = ", round(t,3))
    println("     t_max = ", t_max)
    println("     h = ", h)

    # load checkpoint to calculate dRdu at this time step
    println(" ")
    println("  Setting IC for reverse sweep, i_fwd (forward sweep time step index): ", i_fwd)
    eqn_fwd = cnAdjLoadChkpt(mesh, sbp, opts, adj, physics_func, i_fwd)
    check_q_qvec_consistency(mesh, sbp, eqn_fwd, opts)

    # 20170907: t being passed in is the wrong time. 
    t_ic = t_max

    jac = cnAdjCalcdRdu(mesh, sbp, opts, eqn_fwd, physics_func, i_fwd, t_ic)
    dRdu_n = jac      # TODO: check transpose
    # println(" size of dRdu: ", size(dRdu_n))
    #----------------

    dJdu_CS = calcdJdu_CS(mesh, sbp, eqn_fwd, opts, h, t_ic)  # obtain dJdu at time step n
    dJdu_FD = calcdJdu_FD(mesh, sbp, eqn_fwd, opts, h, t_ic)  # obtain dJdu at time step n
    dJdu = dJdu_CS
    dJdu_analytical = calcObjectiveFn(mesh, sbp, eqn_fwd, opts, h, t_ic, isDeriv=true)

    #---------------------------------------------------------------------------------------------
    # testing calcObjectiveFn
    eqn_poly = eqn_deepcopy(mesh, sbp, eqn_fwd, opts)
    for dof_ix = 1:mesh.numDof
      st = ind2sub((mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl), dof_ix)
      dofnum = st[1]
      nodenum = st[2]
      elnum = st[3]
      if getindex(eqn_poly.q, dofnum, nodenum, elnum) != eqn_poly.q_vec[dof_ix]
        error("problem with ind2sub")
      end
      coords_this_dof = getindex(mesh.coords, :, nodenum, elnum)
      x_coord = coords_this_dof[1]
      y_coord = coords_this_dof[2]

      # u = 2(x+y)
      eqn_poly.q_vec[dof_ix] = 2*(x_coord+y_coord)

      # u = 3
      # eqn_poly.q_vec[dof_ix] = 3.0
    end
    array3DTo1D(mesh, sbp, eqn, opts, eqn_poly.q, eqn_poly.q_vec)
    vis_filename = string("poly_soln")
    saveSolutionToMesh(mesh, real(eqn_poly.q_vec))
    writeVisFiles(mesh, vis_filename)
    J_poly = calcObjectiveFn(mesh, sbp, eqn_poly, opts, h, t_ic, isDeriv=false)
    println(";;;;;;;;;;;;;;;;;;;;;;;;;; J_poly = ", J_poly)
    #---------------------------------------------------------------------------------------------

    J_ic = calcObjectiveFn(mesh, sbp, eqn_fwd, opts, h, t_ic, isDeriv=false)
    print_qvec_coords(mesh, sbp, eqn, opts, to_file=true, filename="IC_qvec_coords.dat", other_field=eqn_fwd.q_vec)
    writedlm("IC_J.dat", J_ic)
    writedlm("IC_qvec.dat", eqn_fwd.q_vec)
    writedlm("IC_dJdu_CS.dat", dJdu_CS)
    writedlm("IC_dJdu_FD.dat", dJdu_FD)
    writedlm("IC_dJdu_analytical.dat", reshape(dJdu_analytical, (mesh.numDof, 1)))
    writedlm("IC_dRdu_n.dat", real(dRdu_n))

    # now that dRdu and dJdu at time step n has been obtained, we can now set the IC for the adjoint eqn
    # println("size of eqn_fwd.q_vec: ", size(eqn_fwd.q_vec))
    # println("size of dRdu_n: ", size(dRdu_n))
    # println("size of dJdu: ", size(dJdu))
    # println("size of adj.q_vec: ", size(adj.q_vec))

    # note: matrix calculus says that the derivative of a scalar by a vector is a row vector.
    #   therefore dJdu is a row vector. However, it is implemented in the code as a column vector.
    #   because of this, the final transpose of psi is not required.
    #
    # No! This is wrong. See 08/30/2017 notes for actual rigorous derivation.
    # I = eye(length(eqn_fwd.q_vec))
    # B = I - (0.5*h*dRdu_n)
    # psi = transpose(B\(-dJdu))
    # psi = B\(-dJdu)

    println(" ::::: in IC, norm(dRdu_n, 1): ", norm(dRdu_n,1))

    # 20170831 old, before xpose fix
    I = eye(length(eqn_fwd.q_vec))
    B = (I - (h/2) * (dRdu_n))
    psi = transpose(B)\(-dJdu)

    # println("size of psi: ", size(psi))
    # println("size of copy(psi): ", size(copy(psi)))
    adj.q_vec = copy(psi)
    writedlm("adj_ic.dat", adj.q_vec)

    array1DTo3D(mesh, sbp, adj, opts, adj.q_vec, adj.q)     # diarray3DTo1D: q_vec -> q

    println("--- Adj IC: end ---")

  end

  #--------------------------------------------------------------------------------------------------------------------------
  # Start of CN loop
  #--------------------------------------------------------------------------------------------------------------------------
  # TODO: Ideally should be t_steps + 2 for clarity, issue #92. Should be fixed for RK4 also.
  t_steps_end = t_steps + 1
  for i = 2:t_steps_end

    println(" ")
    println(" -------------- start of this CN iter. t = ", round(t,3),", i = $i, neg_time = ", neg_time, " --------------")
    @debug1 flush(eqn.params.f)

    # println("   eqn.params.alpha_x: ",eqn.params.alpha_x)
    # println("   eqn.params.alpha_y: ",eqn.params.alpha_y)

    # Write checkpoint data at start of time step
    if neg_time == false
      # for adjoint_straight option: stores every time step's q to disk
      # Note: cannot store full eqn object without extending one of the julia write methods
      if store_u_to_disk == true
        println(" >>>>> Writing solution to disk. i = ", i,", t = ", round(t,3))
        filename = string("qvec_for_adj-", i, ".dat")
        writedlm(filename, eqn.q_vec)
        vis_filename = string("solution_storedtodisk_i-", i)
        saveSolutionToMesh(mesh, real(eqn.q_vec))
        writeVisFiles(mesh, vis_filename)

        time_filename = string("t_during_fwd_sweep_i-",i,".dat")
        writedlm(time_filename, t)              # make sure that this time value, whether t or t_nextstep,
                                                #   correctly corresponds to the q_vec stored just above
      end
    else
      # save every time step's adjoint to disk
      if opts["adjoint_saveall"]
        filename = string("adj_irev-", i, ".dat")
        writedlm(filename, adj.q_vec)

        vis_filename = string("adj_i-", i)
        saveSolutionToMesh(mesh, real(adj.q_vec))
        writeVisFiles(mesh, vis_filename)
      end
    end

    if neg_time == false
      # println(" -------------- eqn.q_vec start of this CN iter. t = $t, i = $i --------------")
      # print_qvec_coords(mesh, sbp, eqn, opts)
      # println(" -------------- J for this eqn.q_vec --------------")
      J_arr = calcObjectiveFn(mesh, sbp, eqn, opts, h, t)
      # println(" -------------- eqn.q_bndry for this eqn.q_vec --------------")
      # print_qvec_coords(mesh, sbp, eqn, opts; bndry=true)
      J = J_arr[1]
      # println("  J: ", J)
      J_filename = string("J_during_fwd_sweep_i-", i,".dat")
      writedlm(J_filename, J)
      println("   Now using eqn.q_vec to compute eqn_nextstep.q_vec.")
    else
      # println(" -------------- adj.q_vec start of this CN iter. t = $t, i = $i --------------")
      # print_qvec_coords(mesh, sbp, adj, opts)
      println("   Now using adj.q_vec to compute adj_nextstep.q_vec.")
    end
    # println(" ")

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
    
      # NOTE: The below checks are being done for the 'current' time step. Not 'nextstep'. 
      #   So the time step needs to be the one before the current ??????????????????????????????????/

    #-------------
    # objective function section
    # 1. read option to indicate which obj fun
    # 2. call it, complex step it, and store it in dJdu
    if neg_time == true
      dJdu = zeros(Tsol, length(eqn.q_vec))
      dJdu_CS = calcdJdu_CS(mesh, sbp, eqn_fwd, opts, h, t)  # obtain dJdu at time step n
      dJdu_FD = calcdJdu_FD(mesh, sbp, eqn_fwd, opts, h, t)  # obtain dJdu at time step n
      dJdu = dJdu_CS
      dJdu_analytical = calcObjectiveFn(mesh, sbp, eqn_fwd, opts, h, t, isDeriv=true)

      filename = string("dJdu_irev-",i,"_CS.dat")
      writedlm(filename, dJdu_CS)
      filename = string("dJdu_irev-",i,"_FD.dat")
      writedlm(filename, dJdu_FD)
      filename = string("dJdu_irev-",i,"_analytical.dat")
      writedlm(filename, reshape(dJdu_analytical, (mesh.numDof, 1)))

      # println("       checking direct method: size(dJdu): ", size(dJdu))

      # VV is the algebraic v, which is dudA, calculated for the advection adjoint test

      # TODO: figure out if calcVV is to be called with t or t_nextstep. I think it's t_nextstep. 
      #       It is for sure t_nextstep when dJdu is called during the direct solve.
      # VV = calcVV(mesh, sbp, adj, opts, t_nextstep)     # scalar only because our x-value of interest is unchanging
      VV = calcVV(mesh, sbp, adj, opts, t)     # scalar only because our x-value of interest is unchanging
      # NOTE 20170712: think I need to be doing t, not t_nextstep. this CN rev's i corresponds to eqn_fwd's i_fwd and t

      # println("       checking direct method: VV: ", VV)
      # println("       checking direct method: size(VV): ", size(VV))
      check_directmethod = transpose(dJdu)*VV
      filename = string("check_directmethod_irev-", i, ".dat")
      writedlm(filename, check_directmethod)

    end

    # advection adjoint check
    if neg_time == true
      println("  calculating dRdA for adjoint check. irev = ", i, ", ifwd = ", i_fwd,", t = ", t)
      dRdA_CS = calcdRdA_CS(mesh, sbp, eqn_fwd, opts, i, t)
      dRdA_FD = calcdRdA_FD(mesh, sbp, eqn_fwd, opts, i, t)
      # dRdA = dRdA_CS
      dRdA = dRdA_FD
      filename = string("dRdA_CS_irev-",i,".dat")
      writedlm(filename, dRdA_CS)
      filename = string("dRdA_FD_irev-",i,".dat")
      writedlm(filename, dRdA_FD)
      println("  writing dRdA file: ", filename)
      # check_adjointmethod = transpose(adj_nextstep.q_vec)*dRdA
      # NOTE 20170712: think I need to be doing adj, not adj_nextstep. this CN rev's i corresponds to eqn_fwd's i_fwd and t
      # check_adjointmethod = transpose(adj.q_vec)*dRdA
      check_adjointmethod = transpose(adj.q_vec)*(-1.0*dRdA)
      filename = string("check_adjointmethod_irev-", i, ".dat")
      writedlm(filename, check_adjointmethod)
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
      println("\n   time step variables-  i: ", i, "  i_fwd: ", i_fwd, "  t_steps: ", t_steps)
      ctx_residual = (physics_func, adj, h, newton_data, i_fwd, dJdu)

      #-------------------------------------------------
      # Output times & indices
      println(" --- CN times and indices ---")
      println("    CN: i = ", i)
      println("    CN: i_fwd = ", i_fwd)
      println("    CN: t = ", round(t,3))
      println("    CN: t_nextstep = ", round(t_nextstep,3))
      #-------------------------------------------------

    end

    if neg_time == false
      # allow for user to select CN's internal Newton's method. Only supports dense FD Jacs, so only for debugging
      if opts["cleansheet_CN_newton"]
        # cnNewton: in cnNewton.jl
        cnNewton(mesh, sbp, opts, h, physics_func, eqn, eqn_nextstep, t)
      else
        @time newtonInner(newton_data, mesh, sbp, eqn_nextstep, opts, cnRhs_uadj, cnJac_uadj, jac, rhs_vec, ctx_residual, t)
        # filename = string("jac_ifwd-", i,".dat")
        # writedlm(filename, real(jac))

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
          newton_data_discard, jac_for_dRdu_global, rhs_vec_discard = setupNewton(mesh, mesh, sbp, eqn_nextstep, opts)
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

          # this time step's actual portion of the global dRdu is the cnJac_uadj. so it has to be I-0.5*h*physicsJac
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
      array1DTo3D(mesh, sbp, adj_nextstep, opts, adj_nextstep.q_vec, adj_nextstep.q)

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
      array1DTo3D(mesh, sbp, eqn_nextstep, opts, eqn_nextstep.q_vec,  eqn_nextstep.q)
    else
      for dof_ix = 1:mesh.numDof
        adj_nextstep.q_vec[dof_ix] = adj.q_vec[dof_ix]
      end
      array1DTo3D(mesh, sbp, adj_nextstep, opts, adj_nextstep.q_vec, adj_nextstep.q)
    end


    t = t_nextstep        # update time step


  end   # end of t step loop

  i = i+1

  # Checkpoint final time step
  if neg_time == false
    # for adjoint_straight option: stores every time step's q to disk
    # Note: cannot store full eqn object without extending one of the julia write methods
    if store_u_to_disk == true
      println(" >>>>> Writing solution to disk. i = ", i,", t = ", t)
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
      filename = string("adj_irev-", i, ".dat")
      writedlm(filename, adj.q_vec)

      vis_filename = string("adj_i-", i)
      saveSolutionToMesh(mesh, real(adj.q_vec))
      writeVisFiles(mesh, vis_filename)
    end
  end

  if neg_time == false
    println(" -------------- eqn.q_vec of last time step. after CN loop. t = $t, i = $i --------------")
    # print_qvec_coords(mesh, sbp, eqn, opts)
  else
    println(" -------------- adj.q_vec of last time step. after CN loop. t = $t, i = $i --------------")
    # print_qvec_coords(mesh, sbp, adj, opts)
  end
  println("\n========================= End of CN =========================")
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
  #   usage: copyForMultistage!(dest, src)   
  # this copyForMultistage! is defined in ODLCommonTools
  if neg_time == false
    copyForMultistage!(eqn, eqn_nextstep)      # copying eqn_nextstep to eqn
    writedlm("solution_final_inCN.dat", real(eqn.q_vec))
  else
    copyForMultistage!(adj, adj_nextstep)      # copying adj_nextstep to eqn
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

  cleanupNewton(newton_data, mesh, mesh, sbp, eqn, opts)

  # @debug1 println("============= end of CN: t = $t ===============")

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
