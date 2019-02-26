# newton.jl: function to do Newtons method, including calculating the Jacobian
# includes residual_evaluation.jl

export newton
import Utils.free

include("newton_setup.jl")

@doc """
### NonlinearSolvers.newton

  This function uses Newton's method to reduce the residual.  The Jacobian
  is calculated using one of several available methods.

  This two options that control how the jacobian is calculated and stored are 
  `jac_type` and `jac_method`.  `jac_type` = 1, 2, or 3 correspond to a 
  dense jacobian matrix, a SparseMatrixCSC matrix that uses direct solves

  jac_type | Meaning
     1        dense matrix
     2        SparseMatrixCSC 
     3        Petsc explictly stored matrix (using CSR internally)
     4        Petsc matrix shell (jacobian vector products)

  jac_method | Meaning
     1          use finite differences
     2          use complex step

  All combinations of jac_type and jac_method are supported except
  jac_type = 4 and jac_method = 2.
 
  The magnitude of perturbation for both finite differences and complex step is
  opts["epsilon"] 


  This function supports extensive debugging options. See the document 
  describing all options supported by PDESolver for details.

  The initial condition in eqn.q_vec should be in whatever variables
  the residual evaluation uses.

  **Arguments**
    * func  : function that evalutes the residual
    * mesh : mesh to use in evaluating the residual
    * sbp : sbp operator to be used to evaluate the residual
    * eqn : EulerData to use to evaluate the residual
    * opts : options dictonary
    * pmesh : mesh used for preconditioning, defaults to mesh

    **Optional Arguments**

    func must have the signature func(mesh, sbp, eqn, opts, t=0.0) 

"""->
function newton(func::Function, mesh::AbstractMesh, sbp::AbstractSBP,
                eqn::AbstractSolutionData, opts, pmesh=mesh, t=0.0)

  jac_type = opts["jac_type"]::Int  # jacobian sparse or dense

  # physicsRhs is defind in residual_evaluation.jl
  rhs_func = physicsRhs
  ctx_residual = (func, )

  pc, lo = getNewtonPCandLO(mesh, sbp, eqn, opts)
  ls = StandardLinearSolver(pc, lo, eqn.comm, opts)

  # allocate jac & rhs, and construct newton_data
  newton_data, rhs_vec = setupNewton(mesh, pmesh, sbp, eqn, opts, ls, alloc_rhs=false)

  newtonInner(newton_data, mesh, sbp, eqn, opts, rhs_func, ls,
              rhs_vec, ctx_residual, t)

  #TODO: return newton_data, don't free it
  free(newton_data)
end

"""
  This function contains the algorithm for Newton's method.  It is intended
  to be used by other functions rather than invoked as a solver directly.
  For example [`newton`](@ref) uses newtonInner to solve steady problems while
  [`crank_nicolson`](@ref) uses it to solve the nonlinear problem arising
  from the time discretization.

  For reference, Newton's method solves the problem R(q) = 0 using
    dR/dq delta q = -R(q)
  where R is the rsidual and q is the solution

  On entry, eqn.q_vec must contain the initial guess for q.  On exit, eqn.q_vec
  will contain the solution to f(q) = 0.  eqn.q will also be consistent with
  eqn.q_vec, as will the send and receive buffers in eqn.shared_data

  On exit, rhs_vec will have the residual corresponding to eqn.q_vec in it,
  with the imaginary part set to zero.

  The [`AbstractLO`](@ref) and [`AbstractPC`](@ref) supplied
  inside the [`LinearSolver`](@ref) object must match the `jac_type`.

  When doing inexact Newton-Krylov, `newonInner` modifies the tolerances
  of the linear solver.  Users calling `newtonInner` repeatedly with the
  same linear solver object should reset the initial tolerances as needed.

  **Inputs:**

   * newton_data: NewtonData object, typically obtained from setupNewton
   * mesh: a mesh object
   * sbp: an SBP operator
   * eqn: a solution data object
   * opts: options dictionary
   * ls: a LinearSolver with the preconditioner and linear operator fully
         initialized.  
   * rhs_vec: vector to store R(q) in
   * ctx_residual: extra data required by rhs_func


  The user must supply two functions, one to calculate the residual vector
  (referred to as rhs_vec), and another to compute the Jacobian.

  rhs_func should compute (eqn.q_vec) -> (rhs_vec) and have the signature

    rhs_func(mesh, sbp, eqn, opts, rhs_vec, ctx_residual, t=0.0)

  Note that the ctx_residual passed into newtonInner is passed directly
  to rhs_func.  The contents of ctx_residual can be whatever is needed by
  rhs_func to perform its computation.  See [`physicsRhs`](@ref) for
  the specific requirements on rhs_func and for an example implementation.

  The same ctx_residual passed into newtonInner is passed directly to
  [`calcLinearOperator`](@ref)..

  This function supports jacobian/preconditioner freezing using the
  prefix "newton".

  Aliasing restrictions: None.  In particular, rhs_vec *can* alias eqn.res_vec,
                         and this leads so some efficiency because it avoids
                         needlessly copying data.

"""
function newtonInner(newton_data::NewtonData, mesh::AbstractMesh, 
                     sbp::AbstractSBP, eqn::AbstractSolutionData,
                     opts, rhs_func, ls::LinearSolver,
                     rhs_vec, ctx_residual=(), 
                     t=0.0;)

  # println(BSTDOUT, " ---- DEBUG eqn.q_vec[1]: ", eqn.q_vec[1])   # TODO: DEBUG GUESS


  # verbose 4 = when newtonInner is used as an inner method for time marching 
  #             or something

  @debug1 println(eqn.params.f, "==== Entered newton")
  @debug1 flush(eqn.params.f)

  myrank = mesh.myrank
  reinitNewtonData(newton_data)
  verbose = newton_data.verbose

  # options
  jac_method = opts["jac_method"]::Int  # finite difference or complex step
  jac_type = opts["jac_type"]::Int  # jacobian sparse or dense
  epsilon = opts["epsilon"]::Float64
#  recalc_prec_freq = opts["recalc_prec_freq"]::Int

  @assert opts["parallel_type"] == 2

  @verbose5 @mpi_master println(BSTDOUT, "\nEntered Newtons Method")
  @verbose5 @mpi_master begin
    println("step_tol = ", newton_data.step_tol)
    println("res_abstol = ", newton_data.res_abstol)
    println("res_reltol = ", newton_data.res_reltol)
  end

  if jac_method == 1  # finite difference
    pert = epsilon
  elseif jac_method == 2  # complex step
    pert = complex(0, epsilon)
  end

  # unpack variables
  m = length(rhs_vec)
  Tsol = typeof(rhs_vec[1])
  res_0 = newton_data.res_0  # residual vector (real)
  res_0_norm = 0.0  # norm of res_0  #TODO: rename this to res_norm
  delta_q_vec = newton_data.delta_q_vec  # newton update


  # evaluating residual at initial condition
  @verbose5 @mpi_master println(BSTDOUT, "evaluating residual at initial condition"); flush(BSTDOUT)
  res_0_norm = rhs_func(mesh, sbp, eqn, opts, rhs_vec, ctx_residual, t)
  recordResNorm(newton_data, res_0_norm)

  # extract the real components to res_0
  for i=1:m
    res_0[i] = real(rhs_vec[i])         # not ok to remove call to real(). This is the reason for the CSR method
  end

  # if the user said to use the first residual for computing relative residuals
  if newton_data.set_rel_norm
    newton_data.res_norm_rel = res_0_norm
  end

  writeFiles(newton_data, mesh, sbp, eqn, opts)
  is_converged = checkConvergence(newton_data)

  # if initial residual satisfies tolerances, return
  if is_converged
    return nothing
  end

  # do Newton's method if not converged
  @verbose5 @mpi_master print(BSTDOUT, "\n")

  #----------------------------------------------------------------------------
  # Start of newton iteration loop
  eqn.params.time.t_newton += @elapsed for i=1:newton_data.itermax

    @verbose5 @mpi_master println(BSTDOUT, "===== newton iteration: ", i)
    newton_data.itr = i

    # recalculate PC and Jacobian if needed
    doRecalculation(newton_data.recalc_policy, i,
                    ls, mesh, sbp, eqn, opts, ctx_residual, t)

    #=             
    recalc_type = decideRecalculation(newton_data.recalc_policy, i)
    if recalc_type == RECALC_BOTH
      calcPCandLO(ls, mesh, sbp, eqn, opts, ctx_residual, t)
    elseif recalc_type == RECALC_PC
      calcPC(ls, mesh, sbp, eqn, opts, ctx_residual, t)
    elseif recalc_type == RECALC_LO
      calcLinearOperator(ls, mesh, sbp, eqn, opts, ctx_residual, t)
    end
    =# 
    #=
    if ((i % recalc_prec_freq)) == 0 || i == 1
      calcPCandLO(ls, mesh, sbp, eqn, opts, ctx_residual, t)
    else  # only recalculate the linear operator
      calcLinearOperator(ls, mesh, sbp, eqn, opts, ctx_residual, t)
    end
    =#

    # compute eigs, condition number, etc.
    doMatrixCalculations(newton_data, opts)

    # negate res
    for j=1:m
      res_0[j] = -res_0[j]
    end

    # println(BSTDOUT, " ---- DEBUG eqn.q_vec[1]: ", eqn.q_vec[1])   # TODO: DEBUG GUESS

    # calculate Newton step
    flush(BSTDOUT)
    tsolve = @elapsed linearSolve(ls, res_0, delta_q_vec, verbose)
    eqn.params.time.t_solve += tsolve
    step_norm = calcNorm(eqn, delta_q_vec)
    recordStepNorm(newton_data, step_norm)

    # perform Newton update

    for j=1:m
      eqn.q_vec[j] += newton_data.step_fac*delta_q_vec[j]
    end
  
    saveSolutionToMesh(mesh, real(eqn.q_vec))
    fname = string("newton_", i)
    writeVisFiles(mesh, fname)



    # calculate residual at updated location, used for next iteration rhs
    res_0_norm = rhs_func(mesh, sbp, eqn, opts, rhs_vec, ctx_residual, t)
    recordResNorm(newton_data, res_0_norm)
    
    # extract real component to res_0
    for j=1:m
      res_0[j] = real(rhs_vec[j])         # not ok to remove call to real(). This is the reason for the CSR method
    end

    writeFiles(newton_data, mesh, sbp, eqn, opts)
    is_converged = checkConvergence(newton_data)
    # @mpi_master println(BSTDOUT, " newton_data.itr: ", newton_data.itr,"  newton_data.step_norm_i: ", newton_data.step_norm_i)
    # DEBUG_CNTHES

    if is_converged
      # remove the imaginary part of rhs_vec before exiting
      for j=1:m
        rhs_vec[j] = real(rhs_vec[j])         # not ok to remove call to real(). This is the reason for the CSR method
      end
      flush(BSTDOUT)

      return nothing
    end

    updateKrylov(newton_data)
    @verbose5 @mpi_master print(BSTDOUT, "\n")
    
  end  # end loop over newton iterations

  # print final information before exiting (unconverged)
  @mpi_master begin
    println(BSTDOUT, "Warning: Newton iteration did not converge in ", newton_data.itermax, " iterations")
    println(BSTDOUT, "  Final step size: ", newton_data.step_norm_i)
    println(BSTDOUT, "  Final residual: ", res_0_norm)
    println(BSTDOUT, "  Final relative residual: ", res_0_norm/newton_data.res_norm_rel)
    @verbose5 close(newton_data.fconv); 
  end
  flush(BSTDOUT)

  # remove the imaginary part of the residual before exiting
  for j=1:m
    rhs_vec[j] = real(rhs_vec[j])         # not ok to remove call to real(). This is the reason for the CSR method
  end
  clearEulerConstants()

  return nothing
end               # end of function newton()

"""
  Do a bunch of file IO for Newtons method.  This function is called
  immediately after a residual evaluation is done.  It also calles
  the majorIterationCallback.

  On entry, res_0, eqn.q, eqn.q_vec, and eqn.res should have the values
  from the current iteration of the Newton loop.

  **Inputs**

   * newton_data: the NewtonDAta
   * mesh
   * sbp
   * eqn
   * opts

  **Options Keys**


"""
function writeFiles(newton_data::NewtonData, mesh, sbp, eqn, opts)

  verbose = newton_data.verbose
  myrank = newton_data.myrank
  itr = newton_data.itr

  res_norm = newton_data.res_norm_i
  res_norm_rel = newton_data.res_norm_rel
  step_norm = newton_data.step_norm_i

  @verbose4 @mpi_master begin
    println(BSTDOUT, "residual norm = ", res_norm)
    println(BSTDOUT, "relative residual ", res_norm/res_norm_rel)
  end

  # decide which files to write
  write_q = opts["writeq"]::Bool
  write_rhs = opts["write_rhs"]::Bool
  write_res = opts["write_res"]::Bool
  write_sol = opts["write_sol"]::Bool


  # write starting values for next iteration to file
  if write_sol
    fname = string("q_vec", itr, "_", myrank, ".dat")
    writedlm(fname, real(eqn.q_vec))
  end

  if write_q
    array1DTo3D(mesh, sbp, eqn, opts, eqn.q_vec)

    fname = string("q", itr, "_", myrank, ".dat")
    writedlm("q$i_$myrank.dat", eqn.q)
  end

  # write rhs to file
  if write_rhs
    fname = string("rhs", itr, "_", myrank, ".dat")
    writedlm(fname, newton_data.res_0)
  end

  if write_res
    fname = string("res", itr, "_", myrank, ".dat")
    writedlm(fname, eqn.res)
  end

  #TODO: have an option to control this
#    saveSolutionToMesh(mesh, eqn.q_vec)
#    writeVisFiles(mesh, "newton_$i")


  @verbose5 eqn.majorIterationCallback(itr, mesh, sbp, eqn, opts, BSTDOUT)

  # write to convergence file
  @verbose5 @mpi_master begin
    println(newton_data.fconv, itr, " ", res_norm, " ", step_norm)
    println(BSTDOUT, "printed to convergence.dat")
  end

  return nothing
end

"""
  Checks whether the Newton iteration has converged or not and prints
  appropriate messages.

  The residual and step norms for the current iteration must be saved to
  the NewtonData object before this function is called

  **Inputs**

   * newton_data: the NewtonData object

  **Outputs**

   * is_converged: true if converged, false otherwise
"""
function checkConvergence(newton_data::NewtonData)

  verbose = newton_data.verbose
  myrank = newton_data.myrank
  itr = newton_data.itr

  res_norm = newton_data.res_norm_i
  res_norm_rel = newton_data.res_norm_rel
  step_norm = newton_data.step_norm_i

  is_converged = false

  # absolute tolerance
  if res_norm < newton_data.res_abstol
    is_converged = true
    @mpi_master if itr == 0
      println(BSTDOUT, "Initial condition satisfies res_tol with residual norm ", res_norm, " < ", newton_data.res_abstol)
    else
      println(BSTDOUT, "Newton iteration converged with residual norm ", res_norm, " < ", newton_data.res_abstol)
    end
  end

  # relative tolerance
  rel_norm = res_norm/res_norm_rel
  if rel_norm < newton_data.res_reltol
    is_converged = true
    if itr == 0
      println(BSTDOUT, "Initial condition satisfied res_reltol with relative residual ", rel_norm, " < ", newton_data.res_reltol)
    else
      println(BSTDOUT, "Newton iteration converged with relative residual norm ", rel_norm, " < ", newton_data.res_reltol)
    end
  end

  # step tolerance
  if step_norm <= newton_data.step_tol && itr > 0
    println("step tolerance satisfied")
    is_converged = true
    @mpi_master println(BSTDOUT, "Newton iteration converged with step_norm = ", step_norm, " < ", newton_data.step_tol)
    @mpi_master println(BSTDOUT, "Final residual = ", res_norm)
  end

  if is_converged
    @mpi_master if itr == 0
      println(BSTDOUT, "Not entering Newton iteration loop")
    else
      println(BSTDOUT, "Iteration count: ", itr)
    end

    clearEulerConstants()
  end

  flush(BSTDOUT)
  return is_converged
end

"""
  Do calculations on the jacobian (eigenvalues, condition number, etc.) if
  requested by options

  **Inputs**

   * newton_data: the NewtonData object
"""
function doMatrixCalculations(newton_data::NewtonData, opts)

  # unpack options
  write_jac = opts["write_jac"]::Bool
  print_cond = opts["print_cond"]::Bool
  print_eigs = opts["print_eigs"]::Bool
  write_eigs = opts["write_eigs"]::Bool
  write_eigdecomp = opts["write_eigdecomp"]::Bool
  jac_type = opts["jac_type"]::Int  # jacobian sparse or dense

  ls = newton_data.ls
  itr = newton_data.itr
  myrank = newton_data.myrank

  # print as determined by options 
  if write_jac && jac_type != 4
    jac = getBaseLO(ls.lo).A
    if jac_type == 3
      MatAssemblyBegin(jac, MAT_FINAL_ASSEMBLY)
      MatAssemblyEnd(jac, MAT_FINAL_ASSEMBLY)
    end
    writedlm("jacobian$itr.dat", full(jac))
    @mpi_master println(BSTDOUT, "finished printing jacobian"); flush(BSTDOUT)
  end
  
  # calculate Jacobian condition number
  if print_cond && ( jac_type == 1 || jac_type == 2)
    jac = getBaseLO(ls.lo).A
    println(BSTDOUT, "calculating condition number of jacobian"); flush(BSTDOUT)
    cond_j = cond(full(jac))
    println(BSTDOUT, "Condition number of jacobian = ", cond_j); flush(BSTDOUT);
  end


  # if eigenvalues requested and we can calculate them
  if (( print_eigs || write_eigs) && (jac_type == 1 || jac_type == 2))
    jac = getBaseLO(ls.lo).A
    println(BSTDOUT, "calculating eigenvalues of jacobian")
    flush(BSTDOUT)
    if typeof(jac) <: Array
      eigs_i = reverse(eigvals(jac))
    else
      eigs_i = reverse(eigvals(full(jac)))
    end
    if print_eigs
      println(BSTDOUT, "eigenvalues =")
      for i=1:length(eigs_i)
        println(BSTDOUT, eigs_i[i])
      end
    end
    flush(BSTDOUT)

    if write_eigs
      writedlm("eigs_real$itr.dat", real(eigs_i))
      writedlm("eigs_imag$itr.dat", imag(eigs_i))
    end

  end

  # do the full eigenvalue decomposition
  # if requested and if julia owns the Jacobian matrix
  if write_eigdecomp && ( jac_type == 1 || jac_type == 2)
    jac = getBaseLO(ls.lo).A
    println(BSTDOUT, "doing eigen decomposition"); flush(BSTDOUT)
    # make a dense jacobian so we can get *all* the eigenvalues and vectors
    # the algorithm for sparse matrices can get all - 2 of them, and might
    # have problems for matrices with a wide range of eigenvalues
    jac_dense = full(jac)
    D, V = eig(jac_dense)
    writedlm("eigdecomp_real$itr.dat", real(D))
    writedlm("eigdecomp_imag$itr.dat", imag(D))
#      writedlm("eigdecomp_realvecs$i.dat", real(V))
#      writedlm("eigdecomp_imagvecs$i.dat", imag(V))

    max_val = typemin(Float64)
    min_val = typemax(Float64)
  elseif write_eigdecomp # && we can't calculate it
    println(BSTDERR, "Warning: not performing eigen decomposition for jacobian of type $jac_type")

  end

  return nothing
end
