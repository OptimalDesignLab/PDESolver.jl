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
"""->
function crank_nicolson(f::Function, h::AbstractFloat, t_max::AbstractFloat,
                        mesh::AbstractMesh, sbp::AbstractSBP, eqn::AbstractSolutionData,
                        opts, res_tol=-1.0, real_time=true; neg_time=false, obj_fn=obj_zero)
                        # NEWNEW: neg_time, obj_fn
  #----------------------------------------------------------------------
#   throw(ErrorException("Crank-Nicolson is in development. Exiting."))

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
    throw(ErrorException("CN not implemented for matrix-free ops. (jac_type cannot be 4)"))
  end

  if myrank == 0
    _f1 = open("convergence.dat", "a+")
    f1 = BufferedIO(_f1)
  end
 
  if neg_time == false
    # start time at 0.0
    t = 0.0
  else    # negative time for unsteady adjoint
    # start time at t_max
    t = t_max
  end

  # calculate t_steps, the number of time steps that CN will take
  t_steps = round(Int, t_max/h)

  # make a copy of the eqn object for storage of t_(n+1) information
  eqn_nextstep = deepcopy(eqn)
  # TODO: copyForMultistage does not give correct values.
  #     deepcopy works for now, but uses more memory than copyForMultistage, if it worked
#   eqn_nextstep = copyForMultistage(eqn)
  eqn_nextstep.q = reshape(eqn_nextstep.q_vec, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)
  eqn_nextstep.res = reshape(eqn_nextstep.res_vec, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)

  @debug1 println("============ In CN ============")

  #-------------------------------------------------------------------------------
   allocate Jac outside of time-stepping loop
  # NOTE 20161103: supplying eqn_nextstep does not work for x^2 + t^2 case, need to use eqn
  newton_data, jac, rhs_vec = setupNewton(mesh, mesh, sbp, eqn, opts, f)

  for i = 2:(t_steps + 1)

    @debug1 println(eqn.params.f, "====== CN: at the top of time-stepping loop, t = $t, i = $i")
    @debug1 flush(eqn.params.f)

    #----------------------------
    # zero out Jac
    #   this works for both PETSc and Julia matrices.
    #   when jac is a Julia matrix, this is effectively wrapping: fill!(jac, 0.0)
    PetscMatZeroEntries(jac)

    # TODO: Allow for some kind of stage loop: ES-Dirk

    # TODO: output freq

    # NOTE:
    # majorIterationCallback: called before every step of Newton's method
    # majorIterationCallback(itr, mesh::AbstractMesh, sbp::AbstractSBP, eqn::AbstractEulerData, opts)

    # NOTE: Must include a comma in the ctx tuple to indicate tuple
    # f is the physics function, like evalEuler

    # NOTE: eqn_nextstep changed to eqn 20161013
    ctx_residual = (f, eqn, h, newton_data)

    @debug1 println(fstdout, "in CN: before call to newtonInner")

    # time step update: h is passed in as argument to crank_nicolson
    if neg_time == false
      # need to add h in the forward time usage
      t_nextstep = t + h
    else
      # need to add h in the forward time usage
      t_nextstep = t - h
    end

    # allow for user to select CN's internal Newton's method. Only supports dense FD Jacs, so only for debugging
    if opts["cleansheet_CN_newton"]
      cnNewton(mesh, sbp, opts, h, f, eqn, eqn_nextstep, t)
    else
      @time newtonInner(newton_data, mesh, sbp, eqn_nextstep, opts, cnRhs, cnJac, jac, rhs_vec, ctx_residual, t)
    end

    # This allows the solution to be updated from _nextstep without a deepcopy.
    #   There are two memory locations used by eqn & eqn_nextstep, 
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

  end   # end of t step loop

  # depending on how many timesteps we do, this may or may not be necessary
  #   usage: copy!(dest, src)
  copy!(eqn, eqn_nextstep)
  writedlm("solution_final_inCN.dat", real(eqn.q_vec))


  if jac_type == 3
    # contents of ctx_newton: (jacp, x, b, ksp)
    NonlinearSolvers.destroyPetsc(jac, newton_data.ctx_newton...)
  end

  @debug1 println("============= end of CN: t = $t ===============")
  return t

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
  # NOTE: changed to eqn 20161013
  eqn = ctx[2]
  h = ctx[3]

  t_nextstep = t + h

  # these two flipped 20161219
  physics_func(mesh, sbp, eqn_nextstep, opts, t_nextstep)
  assembleSolution(mesh, sbp, eqn_nextstep, opts, eqn_nextstep.res, eqn_nextstep.res_vec)

  if opts["parallel_type"] == 2 && mesh.npeers > 0
    startDataExchange(mesh, opts, eqn.q, eqn.q_face_send, eqn.q_face_recv, eqn.params.f)
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

@doc """
### NonlinearSolvers.pde_pre_func

  The pre-function for solving partial differential equations with a physics
  module.  The only operation it performs is disassembling eqn.q_vec into
  eqn.q

  Inputs:
    mesh
    sbp
    eqn
    opts
"""->
function pde_pre_func(mesh, sbp, eqn, opts)
  
  eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
end


@doc """
### NonlinearSolvers.pde_post_func

  The post-function for solving partial differential equations with a physics
  module.  This function multiplies by A0inv, assembles eqn.res into
  eqn.res_vec, multiplies by the inverse mass matrix, and calculates
  the SBP approximation to the integral L2 norm

  Inputs:
    mesh
    sbp
    eqn
    opts

"""->
function pde_post_func(mesh, sbp, eqn, opts; calc_norm=true)
  eqn.multiplyA0inv(mesh, sbp, eqn, opts, eqn.res)
  eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
  for j=1:length(eqn.res_vec) eqn.res_vec[j] = eqn.Minv[j]*eqn.res_vec[j] end
  if calc_norm
    local_norm = calcNorm(eqn, eqn.res_vec)
    eqn.params.time.t_allreduce += @elapsed global_norm = MPI.Allreduce(local_norm*local_norm, MPI.SUM, mesh.comm)
    return sqrt(global_norm)
  end

   return nothing
end
