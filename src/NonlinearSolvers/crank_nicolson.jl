# crank_nicolson.jl
# Crank-Nicolson implicit solver for PDEs

export crank_nicolson, cnResidual

push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/solver/euler"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/NonlinearSolvers"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Debugging"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Utils"))

@doc """
crank_nicolson

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
"""->
function crank_nicolson(f::Function, h::AbstractFloat, t_max::AbstractFloat,
                        mesh::AbstractMesh, sbp::AbstractSBP, eqn::AbstractSolutionData,
                        opts, res_tol=-1.0, real_time=true)

  myrank = MPI.Comm_rank(MPI.COMM_WORLD)
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

  if jac_type == 4
    throw(ErrorException("CN not implemented for matrix-free ops. (jac_type cannot be 4)"))
  end

  @mpi_master if myrank == 0
    _f1 = open("convergence.dat", "a+")
    f1 = BufferedIO(_f1)
  end

  t = 0.0
  t_steps = round(Int, t_max/h)

  # eqn_nextstep = deepcopy(eqn)
  eqn_nextstep = eqn_deepcopy(mesh, sbp, eqn, opts)

  # TODO: copyForMultistage! does not give correct values.
  #     deepcopy works for now, but uses more memory than copyForMultistage!, if it worked
  # eqn_nextstep = copyForMultistage!(eqn)
  eqn_nextstep.q = reshape(eqn_nextstep.q_vec, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)
  eqn_nextstep.res = reshape(eqn_nextstep.res_vec, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)

  @debug1 println("============ In CN ============")

  #-------------------------------------------------------------------------------
  # allocate Jac outside of time-stepping loop
  newton_data, jac, rhs_vec = setupNewton(mesh, mesh, sbp, eqn, opts, f)

  for i = 2:(t_steps + 1)

    @mpi_master println(BSTDOUT, "\ni = ", i, ", t = ", t)
    @debug1 println(eqn.params.f, "====== CN: at the top of time-stepping loop, t = $t, i = $i")
    @debug1 flush(eqn.params.f)

    #----------------------------
    # zero out Jac
    #   this works for both PETSc and Julia matrices.
    #   when jac is a Julia matrix, this is effectively wrapping: fill!(jac, 0.0)
    PetscMatZeroEntries(jac)

    # TODO: Allow for some kind of stage loop: ES-Dirk

    # TODO: output freq

    # NOTE: Must include a comma in the ctx tuple to indicate tuple
    # f is the physics function, like evalEuler

    # NOTE: eqn_nextstep changed to eqn 20161013
    ctx_residual = (f, eqn, h, newton_data)

    t_nextstep = t + h

    # allow for user to select CN's internal Newton's method. Only supports dense FD Jacs, so only for debugging
    if opts["cleansheet_CN_newton"]
      cnNewton(mesh, sbp, opts, h, f, eqn, eqn_nextstep, t)
    else
      newtonInner(newton_data, mesh, sbp, eqn_nextstep, opts, cnRhs, cnJac, 
                  jac, rhs_vec, ctx_residual, t, itermax=30, step_tol=opts["step_tol"], 
                  res_abstol=opts["res_abstol"], res_reltol=opts["res_reltol"],                   res_reltol0=opts["res_reltol0"])
    end

    # do the callback using the current eqn object at time t
    eqn.majorIterationCallback(i, mesh, sbp, eqn, opts, STDOUT)

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
    if (res_norm < res_tol)
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
    disassembleSolution(mesh, sbp, eqn_nextstep, opts, eqn_nextstep.q, eqn_nextstep.q_vec)

    t = t_nextstep
    flush(BSTDOUT)

  end   # end of t step loop

  # depending on how many timesteps we do, this may or may not be necessary
  #   usage: copyForMultistage!(dest, src)
  copyForMultistage!(eqn, eqn_nextstep)


  if jac_type == 3
    # contents of ctx_newton: (jacp, x, b, ksp)
    NonlinearSolvers.destroyPetsc(jac, newton_data.ctx_newton...)
  end

  @debug1 println("============= end of CN: t = $t ===============")
  return t

  flush(BSTDOUT)

end   # end of crank_nicolson function

@doc """
###NonlinearSolvers.cnJac

  Jac of the CN calculation.
  Effectively a wrapper for physicsJac, because the CN Jac is:
    CN_Jac = I + dt/2 * physicsJac

  ctx:    
    physics func must be the first element, i.e. evalEuler
    eqn_nextstep must be the second element
    h must be the third element
    newton_data must be the fourth element
"""->
function cnJac(newton_data, mesh::AbstractMesh, sbp::AbstractSBP,
               eqn_nextstep::AbstractSolutionData, opts, jac, ctx, t)

  myrank = MPI.Comm_rank(MPI.COMM_WORLD)

  physics_func = ctx[1]
  # NOTE: eqn instead of eqn_nextstep, 20161013
  eqn = ctx[2]
  h = ctx[3]
  newton_data = ctx[4]

  jac_type = opts["jac_type"]

  t_nextstep = t + h

  # Forming the CN Jacobian:
  #   call physicsJac with eqn_nextstep & t_nextstep
  #   then form CN_Jac = I + dt/2 * physics_Jac

  NonlinearSolvers.physicsJac(newton_data, mesh, sbp, eqn_nextstep, opts, jac, ctx, t_nextstep)

  # need to flush assembly cache before performing the scale operation.
  #   These are defined for Julia matrices; they're just noops
  PetscMatAssemblyBegin(jac, PETSc.PETSC_MAT_FINAL_ASSEMBLY)
  PetscMatAssemblyEnd(jac, PETSc.PETSC_MAT_FINAL_ASSEMBLY)

  #--------------------------
  # applying dt/2 to jac
  # Jacobian is always 2D
  scale_factor = h*-0.5
  # make this into a petsc_scale_factor
  petsc_scale_factor = PetscScalar(scale_factor)

  # PetscMatScale is defined for all jac types, PETSc and Julia
  # when jac is julia array, this effectively does: scale!(jac, scale_factor)
  PetscMatScale(jac, petsc_scale_factor)

  # PetscMatAssembly___ not necessary here; PetscMatScale is provably local so it doesn't cache its stuff

  #--------------------------
  # adding identity
  ix_petsc_row = zeros(PetscInt, 1)
  ix_petsc_col = zeros(PetscInt, 1)
  value_to_add = zeros(PetscScalar, 1, 1)
  value_to_add[1,1] = 1.0
  flag = PETSc.PETSC_ADD_VALUES

  for i = 1:mesh.numDof
    ix_petsc_row[1,1] = i + mesh.dof_offset       # jac diag index + local to global offset
    ix_petsc_col[1,1] = i + mesh.dof_offset       # jac diag index + local to global offset

    # PETSc function: set_values1!(Jac, [2], [2], Jac[2,2] + 1)
    # args: array, row (as an array of len 1), col (as an array of len 1), new values (as a 2D array of len 1)
    #   set_values1! has different methods for both PETSc matrices and Julia matrices
    #   for Julia dense & sparse arrays, set_values1! effectively does this: jac[i,i] += 1
    #   in serial, mesh.dof_offset is set to 0 automatically
    set_values1!(jac, ix_petsc_row, ix_petsc_col, value_to_add, flag)
  end

  # set_values1! only caches the results; need to be assembled. This happens in petscSolve in petsc_funcs.jl
  #   (The assemble funcs are defined for Julia matrices; they're just noops)

  # jac is now I + dt/2 * physics_jac

  return nothing

end

@doc """
###NonlinearSolvers.cnRhs

  RHS of the CN calculation

  ctx:    
    physics func must be the first element, i.e. evalEuler
    eqn_nextstep must be the second element
    h must be the third element

"""->
function cnRhs(mesh::AbstractMesh, sbp::AbstractSBP, eqn_nextstep::AbstractSolutionData, opts, rhs_vec, ctx, t)

  # eqn comes in through ctx_residual, which is set up in CN before the newtonInner call

  physics_func = ctx[1]
  eqn = ctx[2]
  h = ctx[3]

  t_nextstep = t + h

  # these two flipped 20161219
  physics_func(mesh, sbp, eqn_nextstep, opts, t_nextstep)
  assembleSolution(mesh, sbp, eqn_nextstep, opts, eqn_nextstep.res, eqn_nextstep.res_vec)

  # TODO: can parallel_type != 2 and have CN work?
  if opts["parallel_type"] == 2 && mesh.npeers > 0
    startSolutionExchange(mesh, sbp, eqn, opts)
  end
  physics_func(mesh, sbp, eqn, opts, t)
  assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)

  #   what this is doing:
  #   u_(n+1) - 0.5*dt* (del dot G_(n+1)) - u_n - 0.5*dt* (del dot G_n)
  for i = 1:mesh.numDof

    temp1 = eqn_nextstep.q_vec[i] - 0.5*h*eqn_nextstep.res_vec[i]
    temp2 = eqn.q_vec[i] + 0.5*h*eqn.res_vec[i]

    rhs_vec[i] = temp1 - temp2 

    # NOTE: question: is there a sign problem here? should rhs_vec = -rhs_vec ?
    #     NO. this negative gets applied in newton.jl, where res_0[i] = -res_0[i]

  end

  return nothing

end     # end of function cnRhs

# the goal is to replace newton.jl.
# this will go into CN in the time-stepping loop
function cnNewton(mesh, sbp, opts, h, physics_func, eqn, eqn_nextstep, t)
  println("++++++++++++++++ clean sheet Newton being run ++++++++++")

  println("---- physics_func: ",physics_func)

  # Jac on eqn or eqn_nextstep?

  epsilon = 1e-8
  t_nextstep = t + h

  jac = zeros(mesh.numDof, mesh.numDof)

  # emulates physicsJac
  # so we need to step through the jacobian column wise.
  #   d res[1]/d q[1]   d res[1]/d q[2]   d res[1]/d q[3] ...
  #   d res[2]/d q[1]   d res[2]/d q[2]   d res[2]/d q[3] ...
  #   d res[3]/d q[1]   d res[3]/d q[2]   d res[3]/d q[3] ...
  #   ...               ...               ...

  newton_itermax = 2
  delta_q_vec = zeros(eqn_nextstep.q_vec)

  # newton_loop starting here?
  for newton_i = 1:newton_itermax

    #--------------------------
    # emulates physicsJac
    unperturbed_q_vec = copy(eqn_nextstep.q_vec)

    physics_func(mesh, sbp, eqn_nextstep, opts, t_nextstep)
    # needed b/c physics_func only updates eqn.res
    assembleSolution(mesh, sbp, eqn_nextstep, opts, eqn_nextstep.res, eqn_nextstep.res_vec)
    # Comment here about mass matrix inv multiplication TODO
    applyMassMatrixInv(mesh, eqn_nextstep, eqn_nextstep.res_vec)
    unperturbed_res_vec = copy(eqn_nextstep.res_vec)

    for i = 1:mesh.numDof
      eqn_nextstep.q_vec[i] = eqn_nextstep.q_vec[i] + epsilon

      physics_func(mesh, sbp, eqn_nextstep, opts, t_nextstep)
      assembleSolution(mesh, sbp, eqn_nextstep, opts, eqn_nextstep.res, eqn_nextstep.res_vec)
      applyMassMatrixInv(mesh, eqn, eqn_nextstep.res_vec)

      jac[:,i] = (eqn_nextstep.res_vec - unperturbed_res_vec)/epsilon

      eqn_nextstep.q_vec[i] = unperturbed_q_vec[i]

    end

    #--------------------------
    # emulates cnJac
    scale!(jac, -0.5*h)
    for i = 1:mesh.numDof
      jac[i,i] += 1
    end

    #--------------------------
    # emulates cnRhs
    #   what this is doing:
    #   u_(n+1) - 0.5*dt* (del dot G_(n+1)) - u_n - 0.5*dt* (del dot G_n)
    #=
    physics_func(mesh, sbp, eqn, opts, t)
    assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
    physics_func(mesh, sbp, eqn_nextstep, opts, t_nextstep)
    assembleSolution(mesh, sbp, eqn_nextstep, opts, eqn_nextstep.res, eqn_nextstep.res_vec)

    rhs_vec = zeros(eqn.q_vec)

    for i = 1:mesh.numDof
      rhs_vec[i] = eqn_nextstep.q_vec[i] - h*0.5*eqn_nextstep.res_vec[i] - eqn.q_vec[i] - h*0.5*eqn.res_vec[i]
    end
    =#
    physics_func(mesh, sbp, eqn, opts, t)
    assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
    applyMassMatrixInv(mesh, eqn, eqn.res_vec)
    current_t_step_contribution = zeros(eqn.q_vec)
    for i = 1:mesh.numDof
      current_t_step_contribution[i] = - eqn.q_vec[i] - h*0.5*eqn.res_vec[i]
    end

    # Test for 3D Minv results
    # this works!
#     res_vec_control = deepcopy(eqn.res_vec)
#     res_vec_test = deepcopy(eqn.res_vec)
#     res_control = deepcopy(eqn.res)
#     res_test = deepcopy(eqn.res)
#     applyMassMatrixInv3D(mesh, sbp, eqn, res_test)
#     assembleSolution(mesh, sbp, eqn, opts, res_test, res_vec_test)
#     println("=+=+=+ norm of diff btwn res_vec_test & res_vec_control: ", norm(res_vec_test - res_vec_control))

    physics_func(mesh, sbp, eqn_nextstep, opts, t_nextstep)
    assembleSolution(mesh, sbp, eqn_nextstep, opts, eqn_nextstep.res, eqn_nextstep.res_vec)
    applyMassMatrixInv(mesh, eqn_nextstep, eqn_nextstep.res_vec)
    next_t_step_contribution = zeros(eqn.q_vec)
    for i = 1:mesh.numDof
      next_t_step_contribution[i] = eqn_nextstep.q_vec[i] - h*0.5*eqn_nextstep.res_vec[i] 
    end

    rhs_vec = zeros(eqn.q_vec)

    for i = 1:mesh.numDof
      rhs_vec[i] = current_t_step_contribution[i] + next_t_step_contribution[i]
    end

    # TODO: check these args
    rhs_norm = calcNorm(eqn, rhs_vec, strongres=true)

    #--------------------------
    # start of actual Newton
    neg_rhs = scale(rhs_vec, -1.0)

    fill!(delta_q_vec, 0.0)
    delta_q_vec = jac\neg_rhs
    fill!(jac, 0.0)

    for i = 1:mesh.numDof
      eqn_nextstep.q_vec[i] += delta_q_vec[i]
    end

    rhs_norm_tol = 1e-6
    if rhs_norm < rhs_norm_tol
      println("=== cnNewton converged with rhs_norm under $rhs_norm_tol -- newton iters: $newton_i ===")
      return nothing
    end

  end   # end of newton iterations

  println("=== cnNewton did not converge ===")
  return nothing


end

# TODO: comment here
function applyMassMatrixInv(mesh, eqn, vec)

  for k = 1:mesh.numDof
    vec[k] = eqn.Minv[k] * vec[k]
  end

  return vec
end

# TODO: comment here
function applyMassMatrixInv3D(mesh, sbp, eqn, arr)

  for i = 1:mesh.numEl
    for j = 1:sbp.numnodes
      for k = 1:mesh.numDofPerNode
        arr[k, j, i] = eqn.Minv3D[k, j, i] * arr[k, j, i]
      end
    end
  end

  return arr
end


