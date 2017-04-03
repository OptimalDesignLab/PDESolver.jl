# crank_nicolson.jl
# Crank-Nicolson implicit solver for PDEs

export crank_nicolson, cnResidual

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
function crank_nicolson{Tmsh, Tsol}(physics_func::Function, h::AbstractFloat, t_max::AbstractFloat,
                        mesh::AbstractMesh{Tmsh}, sbp::AbstractSBP, eqn::AbstractSolutionData{Tsol},
                        opts, res_tol=-1.0; neg_time=false, obj_fn=obj_zero, store_u_to_disk=false)
                        # NEWNEW: neg_time, obj_fn, store_u_to_disk
  #----------------------------------------------------------------------

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
 
  if neg_time == false    # negative time is for unsteady adjoint
    t = 0.0     # start time at 0.0
  else    
    t = t_max   # start time at t_max
  end

  # calculate t_steps, the number of time steps that CN will take
  t_steps = round(Int, t_max/h)

  if neg_time == false
    # make a copy of the eqn object for storage of t_(n+1) information
    eqn_nextstep = deepcopy(eqn)
    # TODO: copyForMultistage does not give correct values.
    #     deepcopy works for now, but uses more memory than copyForMultistage, if it worked
    # eqn_nextstep = copyForMultistage(eqn)
    eqn_nextstep.q = reshape(eqn_nextstep.q_vec, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)
    eqn_nextstep.res = reshape(eqn_nextstep.res_vec, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)
  else
    # Initialize adjoint pde eqn objects
    adj = deepcopy(eqn)
    adj.q = reshape(adj.q_vec, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)
    adj.res = reshape(adj.res_vec, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)

    adj_nextstep = deepcopy(eqn)
    adj_nextstep.q = reshape(adj_nextstep.q_vec, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)
    adj_nextstep.res = reshape(adj_nextstep.res_vec, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)
  end

  @debug1 println("============ In CN ============")

  #-------------------------------------------------------------------------------
  # allocate Jac outside of time-stepping loop
  # these jacs are for full CN or CN adj jac
  println("===== neg_time: ", neg_time, " =====")
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
    i_actual = t_steps + 1  # index during forward sweep of the n'th q_vec. +1 instead of +3-i because the loop adds 2

    # load checkpoint to calculate dRdu at this time step
    eqn_dummy = cnAdjLoadChkpt(mesh, sbp, opts, adj, physics_func, i_actual, t)

    jac = cnAdjCalcdRdu(mesh, sbp, opts, eqn_dummy, physics_func, t)
    dRdu_n = jac      # TODO: check transpose
    #----------------

    dJdu = calcdJdu_CS(mesh, sbp, eqn_dummy, opts)  # obtain dJdu at time step n
    # now that dRdu and dJdu at time step n has been obtained, we can now set the IC for the adjoint eqn
    I = eye(length(eqn_dummy.q_vec))
    B = (I - (h/2) * (dRdu_n))
    psi = transpose(B)\(-dJdu)
    adj.q_vec = psi

  end

  for i = 2:(t_steps + 1)

    @debug1 println(eqn.params.f, "====== CN: at the top of time-stepping loop, t = $t, i = $i")
    @debug1 flush(eqn.params.f)

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

    #-------------
    # objective function section
    # 1. read option to indicate which obj fun
    # 2. call it, complex step it, and store it in dJdu
    if neg_time == true
      dJdu = zeros(Tsol, length(eqn.q_vec))
      J = calcObjectiveFn(mesh, sbp, adj, opts; isDeriv=false)
      dJdu = calcdJdu_CS(mesh, sbp, adj, opts)
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
      i_actual = t_steps + 3 - i
      ctx_residual = (physics_func, adj, h, newton_data, i_actual, dJdu)
    end

    @debug1 println(fstdout, "in CN: before call to newtonInner")

    # time step update: h is passed in as argument to crank_nicolson
    if neg_time == false
      # need to add h in the forward time usage
      t_nextstep = t + h
    else
      # need to subtract h in the reverse time usage
      t_nextstep = t - h
    end

    if neg_time == false
      # allow for user to select CN's internal Newton's method. Only supports dense FD Jacs, so only for debugging
      if opts["cleansheet_CN_newton"]
        # cnNewton: in cnNewton.jl
        cnNewton(mesh, sbp, opts, h, physics_func, eqn, eqn_nextstep, t)
      else
        @time newtonInner(newton_data, mesh, sbp, eqn_nextstep, opts, cnRhs, cnJac, jac, rhs_vec, ctx_residual, t)
      end
    # else      # call newtonInner using cnAdjJac and cnAdjRhs
      # @time newtonInner(newton_data, mesh, sbp, adj_nextstep, opts, cnAdjRhs, cnAdjJac, jac, rhs_vec, ctx_residual, t)
    else    # direct solve for psi_i

      # TODO TODO TODO
      # adj_nextstep.q_vec = cnAdjDirect(mesh, sbp, opts, adj_nextstep, physics_func, jac, i_actual, h, t)
      adj_nextstep.q_vec = cnAdjDirect(mesh, sbp, opts, adj, physics_func, jac, i_actual, h, t)
      disassembleSolution(mesh, sbp, adj_nextstep, opts, adj_nextstep.q, adj_nextstep.q_vec)

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

    # for adjoint_straight option: stores every time step's q to disk
    # Note: cannot store full eqn object without extending one of the julia write methods
    if store_u_to_disk == true
      filename = string("qvec_for_adj-", i, ".dat")
      writedlm(filename, eqn.q_vec)
    end

    # Note: we now need to copy the updated q over for the initial newton guess
    if neg_time == false
      for i = 1:mesh.numDof
        eqn_nextstep.q_vec[i] = eqn.q_vec[i]
      end
      disassembleSolution(mesh, sbp, eqn_nextstep, opts, eqn_nextstep.q, eqn_nextstep.q_vec)
    else
      for i = 1:mesh.numDof
        adj_nextstep.q_vec[i] = adj.q_vec[i]
      end
      disassembleSolution(mesh, sbp, adj_nextstep, opts, adj_nextstep.q, adj_nextstep.q_vec)
    end

    t = t_nextstep        # update time step

  end   # end of t step loop

  # depending on how many timesteps we do, this may or may not be necessary
  #   usage: copy!(dest, src)
  if neg_time == false
    copy!(eqn, eqn_nextstep)      # copying eqn_nextstep to eqn
    writedlm("solution_final_inCN.dat", real(eqn.q_vec))
  else
    copy!(adj, adj_nextstep)      # copying adj_nextstep to eqn
    writedlm("adjoint_final_inCN.dat", real(adj.q_vec))
  end

  if jac_type == 3      # if jac is a Petsc matrix, it needs to be freed when we're done using it
    # contents of ctx_newton: (jacp, x, b, ksp)
    NonlinearSolvers.destroyPetsc(jac, newton_data.ctx_newton...)
  end

  @debug1 println("============= end of CN: t = $t ===============")
  return t

end   # end of crank_nicolson function
