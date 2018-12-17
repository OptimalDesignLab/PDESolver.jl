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

  # this loop is 2:(t_steps+1) when not restarting
  for i = istart:(t_steps + 1)

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

    # TODO: Allow for some kind of stage loop: ES-Dirk

    # TODO: output freq

    # NOTE: Must include a comma in the ctx tuple to indicate tuple
    # f is the physics function, like evalEuler

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
    end

    # do the callback using the current eqn object at time t
    eqn.majorIterationCallback(i, mesh, sbp, eqn, opts, BSTDOUT)

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

    # check stopping conditions
    if (res_norm < res_tol)  # ???
      if myrank == 0
        println(BSTDOUT, "breaking due to res_tol, res norm = $res_norm")
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


    # This allows the solution to be updated from _nextstep without a deepcopy.
    # There are two memory locations used by eqn & eqn_nextstep, 
    #   and this flips the location of eqn & eqn_nextstep every time step
    eqn_temp = eqn
    eqn = eqn_nextstep
    eqn_nextstep = eqn_temp

    # Note: we now need to copy the updated q over for the initial newton guess
    for i = 1:mesh.numDof
      eqn_nextstep.q_vec[i] = eqn.q_vec[i]
    end
    array1DTo3D(mesh, sbp, eqn_nextstep, opts, eqn_nextstep.q_vec, eqn_nextstep.q)

#    t = t_nextstep
    flush(BSTDOUT)

  end   # end of t step loop

  # final time update
  t += h

  # depending on how many timesteps we do, this may or may not be necessary
  #   usage: copyForMultistage!(dest, src)
  copyForMultistage!(eqn, eqn_nextstep)


  #TODO: return the NewtonData?
  free(newton_data)

  @debug1 println("============= end of CN: t = $t ===============")
  return t

  flush(BSTDOUT)

end   # end of crank_nicolson_ds function

