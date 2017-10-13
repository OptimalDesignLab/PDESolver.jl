# newton.jl: function to do Newtons method, including calculating the Jacobian
# includes residual_evaluation.jl and petsc_funcs.jl

export newton
global const insert_freq = 1
@doc """
  This type holds all the data the might be needed for Newton's method,
  including globalization.  This simplifies the data handling and 
  the C callback used by Petsc
"""->
type NewtonData{Tsol, Tres}

  # MPI info
  myrank::Int
  commsize::Int

  # inexact Newton-Krylov parameters
  reltol::Float64
  abstol::Float64
  dtol::Float64
  itermax::Int
  krylov_gamma::Float64  # update parameter for krylov tolerance

  res_norm_i::Float64  # current step residual norm
  res_norm_i_1::Float64  # previous step residual norm
  # Pseudo-transient continuation Euler
  tau_l::Float64  # current pseudo-timestep
  tau_vec::Array{Float64, 1}  # array of solution at previous pseudo-timestep

  # temporary arrays used to for Petsc MatSetValues
  insert_idx::Int
  localsize::Int
  vals_tmp::Array{Float64, 2}
  idx_tmp::Array{PetscInt, 1}
  idy_tmp::Array{PetscInt, 1}

  # tuple of values needed to cleanup Petsc data structures
  ctx_newton

  # ctx needed by MatShell
  # TODO: get list of contents from jac-vec prod wrapper
  ctx_petsc

  #TODO; make this an outer constructor
  function NewtonData(mesh, sbp, eqn, opts)

    println("entered NewtonData constructor")
    println("typeof(eqn) = ", typeof(eqn))

    myrank = mesh.myrank
    commsize = mesh.commsize

    reltol = opts["krylov_reltol"]
    abstol = opts["krylov_abstol"]
    dtol = opts["krylov_dtol"]
    itermax = opts["krylov_itermax"]
    krylov_gamma = opts["krylov_gamma"]

    res_norm_i = 0.0
    res_norm_i_1 = 0.0
    if opts["newton_globalize_euler"]
      tau_l, tau_vec = initEuler(mesh, sbp, eqn, opts)
    else
      tau_l = opts["euler_tau"]
      tau_vec = []
    end

    local_size = mesh.numNodesPerElement*mesh.numDofPerNode*insert_freq
    vals_tmp = zeros(local_size, 1) # values
    idx_tmp = zeros(PetscInt, local_size)  # row index
    idy_tmp = zeros(PetscInt, 1)  # column indices
    localsize = mesh.numNodesPerElement*mesh.numDofPerNode

    # NOTE: we are leaving ctx_newton uninitialized because 
    #   createPetscData needs it, but ctx_newton depends on its returned values

    # initialize these to something, will be replaced in setupNewton
    ctx_newton = ()
    ctx_petsc = ()

    return new(myrank, commsize, reltol, abstol, dtol, 
                      itermax, krylov_gamma, 
                      res_norm_i, res_norm_i_1, tau_l, tau_vec, 1, localsize, vals_tmp, 
                      idx_tmp, idy_tmp, ctx_newton, ctx_petsc)
  end

end

include("residual_evaluation.jl")  # some functions for residual evaluation
include("jacobian.jl")
include("petsc_funcs.jl")  # Petsc related functions

@doc """
### NonlinearSolvers.setupNewton
  
  alloc_rhs: keyword arg to allocate a new object or not for rhs_vec
                true (default) allocates a new vector
                false will use eqn.res_vec

  rhs_func: only used for Petsc in matrix-free mode to do Jac-vec products
            should be the rhs_func passed into [`newtonInner`](@ref)
  ctx_residual: ctx_residual passed into [`newtonInner`](@ref)

  Allocates Jac & RHS

  See [`cleanupNewton`](@ref) to the cleanup function

"""->
function setupNewton{Tsol, Tres}(mesh, pmesh, sbp,
                     eqn::AbstractSolutionData{Tsol, Tres}, opts;
                     alloc_rhs=true)

  jac_type = opts["jac_type"]
  Tjac = typeof(real(eqn.res_vec[1]))  # type of jacobian, residual
  m = mesh.numDof

  # ctx_newton is not defined yet
  newton_data = NewtonData{Tsol, Tres}(mesh, sbp, eqn, opts)

  # Allocation of Jacobian, depending on type of matrix
  eqn.params.time.t_alloc += @elapsed if jac_type == 1  # dense
    jac = zeros(Tjac, m, m)  # storage of the jacobian matrix
    ctx_newton = ()
  elseif jac_type == 2  # sparse
    if typeof(mesh) <: AbstractCGMesh
      println("creating CG SparseMatrix")
      jac = SparseMatrixCSC(mesh.sparsity_bnds, Tjac)
    else
      println("Creating DG sparse matrix")
      jac = SparseMatrixCSC(mesh, Tjac)
    end
    ctx_newton = ()
  elseif jac_type == 3 || jac_type == 4 # petsc
    jac, jacp, x, b, ksp = createPetscData(mesh, pmesh, sbp, eqn, opts, newton_data)
    ctx_newton = (jac, jacp, x, b, ksp)
  end

  # now put ctx_newton into newton_data
  newton_data.ctx_newton = ctx_newton

  # For simple cases, especially for Newton's method as a steady solver,
  #   having rhs_vec and eqn.res_vec pointing to the same memory
  #   saves us from having to copy back and forth
  if alloc_rhs 
    rhs_vec = zeros(Tsol, size(eqn.res_vec))
  else
    rhs_vec = eqn.res_vec
  end

  # should be all zeros if alloc_rhs is true
#   writedlm("rhs_initial.dat", rhs_vec)

  return newton_data, jac, rhs_vec

end   # end of setupNewton

"""
  Cleans up after running Newton's method.  Petsc will leak memory if this
  function is not called

  **Inputs**

   * newton_data: the NewtonData object
   * mesh: the AbstractMesh
   * pmesh: the AbstractMesh used for preconditioning
   * sbp: the SBP operator
   * eqn: the AbstractSolutionData

"""
function cleanupNewton{Tsol, Tres}(newton_data::NewtonData, mesh, pmesh, sbp,
                                   eqn::AbstractSolutionData{Tsol, Tres}, opts)

  jac_type = opts["jac_type"]
  if jac_type == 3 || jac_type == 4
    NonlinearSolvers.destroyPetsc(newton_data.ctx_newton...)
  end

  return nothing
end

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

  Arguments:
    * func  : function that evalutes the residual
    * mesh : mesh to use in evaluating the residual
    * sbp : sbp operator to be used to evaluate the residual
    * eqn : EulerData to use to evaluate the residual
    * opts : options dictonary
    * pmesh : mesh used for preconditioning, defaults to mesh

    Optional Arguments
    * itermax : maximum number of Newton iterations
    * step_tol : step size stopping tolerance
    * res_tol : residual stopping tolerance
    * ctx_newton : 'context', i.e. the tuple of objects that is passed to func.
            The tuple is splatted before being passed to func.

    func must have the signature func(mesh, sbp, eqn, opts, t=0.0) 

"""->
function newton(func::Function, mesh::AbstractMesh, sbp::AbstractSBP, eqn::AbstractSolutionData, opts, pmesh=mesh, t=0.0; 
                itermax=10, step_tol=1e-6, res_abstol=1e-6, res_reltol=1e-6, res_reltol0=-1.0)

  jac_type = opts["jac_type"]::Int  # jacobian sparse or dense

  # physicsRhs and physicsJac are defined in newton.jl
  # They evaluate the basic Jac & RHS of only the physics function, such as evalEuler
  rhs_func = physicsRhs
  jac_func = physicsJac
  ctx_residual = (func, )

  # allocate jac & rhs, and construct newton_data
  newton_data, jac, rhs_vec = setupNewton(mesh, pmesh, sbp, eqn, opts, alloc_rhs=false)


  newtonInner(newton_data, mesh, sbp, eqn, opts, rhs_func, jac_func, jac,
              rhs_vec, ctx_residual, t, itermax=itermax, step_tol=step_tol,
              res_abstol=res_abstol, res_reltol=res_reltol, 
              res_reltol0=res_reltol0)

  #TODO; create a general destructor for newtonInner, the complement of 
  #      setupNewton

  cleanupNewton(newton_data, mesh, pmesh, sbp, eqn, opts)
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


  **Inputs:**

   * newton_data: NewtonData object, typically obtained from setupNewton
   * mesh: a mesh object
   * sbp: an SBP operator
   * eqn: a solution data object
   * opts: options dictionary
   * rhs_func: function to evaluate f(q)
   * jac_func: function to compute the Jacobian df/dq
   * rhs_vec: vector to store R(q) in
   * ctx_residual: extra data required by rhs_func


  The user must supply tow functions, one to calculate the residual vector
  (referred to as rhs_vec), and another to compute the Jacobian.

  rhs_func should compute (eqn.q_vec) -> (rhs_vec) and have the signature

    rhs_func(mesh, sbp, eqn, opts, rhs_vec, ctx_residual, t=0.0)

  Note that the ctx_residual passed into newtonInner is passed directly
  to rhs_func.  The contents of ctx_residual can be whatever is needed by
  rhs_func to perform its computation.  See [`physicsRhs`](@ref) for
  the specific requirements on rhs_func and for an example implementation.

  jac_func should compute the Jacobian and have the signature

    jac_func(newton_data::NewtonData, mesh, sbp, eqn, opts, jac, ctx_residual, t=0.0)

  The same ctx_residual passed into newtonInner is passed directly to jac_func.
  It is currently the same ctx_residual as is passed to rhs_func (TODO:
  stop that).  See [`physicsJac`](@ref) for details on the requirements
  of this function.


  Aliasing restrictions: None.  In particular, rhs_vec *can* alias eqn.res_vec,
                         and this leads so some efficiency because it avoids
                         needlessly copying data.

"""
function newtonInner(newton_data::NewtonData, mesh::AbstractMesh, sbp::AbstractSBP, eqn::AbstractSolutionData,
                     opts, rhs_func, jac_func, jac, rhs_vec, ctx_residual=(), t=0.0;
                     itermax=10, step_tol=1e-6, res_abstol=1e-6, res_reltol=1e-6, res_reltol0=-1.0 )

  # verbose 4 = when newtonInner is used as an inner method for time marching or something


  # this function drives the non-linear residual to some specified tolerance
  # using Newton's Method
  # the jacobian is formed using finite differences
  # the initial condition is stored in eqn.q_vec
  # itermax is the maximum number of iterations
  # this function is type unstable for certain variables, but that's ok
  # the dispatch to the backslash solver and possibly the jacobian calculation
  # function will be runtime dispatched

  @debug1 println(eqn.params.f, "==== Entered newton")
  @debug1 flush(eqn.params.f)

  myrank = mesh.myrank

  # options
  write_rhs = opts["write_rhs"]::Bool
  write_jac = opts["write_jac"]::Bool
  print_cond = opts["print_cond"]::Bool
  print_eigs = opts["print_eigs"]::Bool
  write_eigs = opts["write_eigs"]::Bool
  write_eigdecomp = opts["write_eigdecomp"]::Bool
  write_sol = opts["write_sol"]::Bool
  write_vis = opts["write_vis"]::Bool
  output_freq = opts["output_freq"]::Int
  write_qic = opts["write_qic"]::Bool
  write_res = opts["write_res"]::Bool
  write_q = opts["writeq"]::Bool
  jac_method = opts["jac_method"]::Int  # finite difference or complex step
  jac_type = opts["jac_type"]::Int  # jacobian sparse or dense
  epsilon = opts["epsilon"]::Float64
  globalize_euler = opts["newton_globalize_euler"]::Bool
  recalc_prec_freq = opts["recalc_prec_freq"]::Int
  use_jac_precond = opts["use_jac_precond"]::Bool
  verbose = opts["newton_verbosity"]::Int

  @assert opts["parallel_type"] == 2

  @verbose5 @mpi_master println(BSTDOUT, "\nEntered Newtons Method")
  @verbose5 @mpi_master begin
    println("step_tol = ", step_tol)
    println("res_abstol = ", res_abstol)
    println("res_reltol = ", res_reltol)
    println("res_reltol0 = ", res_reltol0)
  end

  if (t != 0.0) && (jac_type == 4)
    throw(ErrorException, "Matrix free Petsc cannot be used for solving unsteady problems, see TODO in calcJacVecProd_wrapper")
  end


  if jac_method == 1  # finite difference
    pert = epsilon
  elseif jac_method == 2  # complex step
    pert = complex(0, epsilon)
  end

  # set the ctx pointer for the matrix-free matrix
  if jac_type == 4
    ctx_petsc = (mesh, sbp, eqn, opts, newton_data, rhs_func, ctx_residual)
    newton_data.ctx_petsc = ctx_petsc
    ctx_ptr = pointer_from_objref(ctx_petsc)
    MatShellSetContext(jac, ctx_ptr)
  end

  if jac_type == 3 || jac_type == 4
    jacp = newton_data.ctx_newton[2]
  end

  Tjac = typeof(real(rhs_vec[1]))  # type of jacobian, residual
  m = length(rhs_vec)

  step_fac = 1.0 # step size limiter
#  jac_recal = 0  # number of iterations since jacobian was recalculated
  Tsol = typeof(rhs_vec[1])
  res_0 = zeros(Tjac, m)  # function evaluated at u0
  res_0_norm = 0.0  # norm of res_0
  delta_q_vec = zeros(Tjac, m)  # newton update
  step_norm = zero(Tjac)  # norm of newton update
  step_norm_1 = zero(Tjac) # norm of previous newton update


  ##### Write iteration 0 output #####
  # open file to write convergence data to
  # append to be on the safe side
  @verbose5 @mpi_master begin
    _fconv = open("convergence.dat", "a+")
    fconv = BufferedIO(_fconv)
  end
  # evaluating residual at initial condition
  @verbose5 @mpi_master println(BSTDOUT, "evaluating residual at initial condition"); flush(BSTDOUT)
  res_0_norm = newton_data.res_norm_i = rhs_func(mesh, sbp, eqn, opts, rhs_vec, ctx_residual, t)
  @verbose5 @mpi_master println(BSTDOUT, "res_0_norm = ", res_0_norm); flush(BSTDOUT)

#  saveSolutionToMesh(mesh, eqn.q_vec)
#  writeVisFiles(mesh, "newton_0")

  @verbose5 eqn.majorIterationCallback(0, mesh, sbp, eqn, opts, BSTDOUT)


  # extract the real components to res_0
  for i=1:m
    res_0[i] = real(rhs_vec[i])
  end
#   res_0 = deepcopy(real(rhs_vec))
  # TODO TODO TODO: changed to deepcopy 20161113
#   writedlm("rhs_before_newton_loop.dat", res_0)

  @verbose5 @mpi_master begin
    println(fconv, 0, " ", res_0_norm, " ", 0)
    flush(fconv)
  end

  # post-residual iteration 0 output
  if write_rhs
    writedlm("rhs0_$myrank.dat", res_0)
  end

  if write_qic
    writedlm("q0_$myrank.dat", eqn.q)
  end

  if write_res
    writedlm("res0_$myrank.dat", eqn.res)
  end

  # check if initial residual satisfied absolute or relative tolerances
  if res_0_norm < res_abstol || 
    (res_reltol0 > 0 && res_0_norm/res_reltol0 < res_reltol)

    # print which criteria was statisfied
    if res_0_norm/res_reltol0 < res_reltol
      @mpi_master println(BSTDOUT, "Initial condition satisfied res_reltol with relative residual ", res_0_norm/res_reltol0)
      @mpi_master println(BSTDOUT, "Residual ", res_0_norm)
    else
      @mpi_master println(BSTDOUT, "Initial condition satisfies res_tol with residual norm ", res_0_norm)
    end
    # no need to assemble q into q_vec because it never changed

    @verbose5 @mpi_master close(fconv)
    flush(BSTDOUT)

    @verbose5 @mpi_master println(BSTDOUT, "Not entering Newton iteration loop")
    flush(BSTDOUT)
    return nothing

  end  # end if tolerances satisfied

  # use the relative residual convergence criteria if specified (ie. if > 0)
  @verbose5 @mpi_master println(BSTDOUT, "res_reltol0 = ", res_reltol0)
  if res_reltol0 > 0  # use the supplied res_reltol0 value
    @verbose5 @mpi_master println(BSTDOUT, "using supplied value for relative residual")
    res_reltol_0 = res_reltol0
  else
    @verbose5 @mpi_master println(BSTDOUT, "using initial residual for relative residual")
    res_reltol_0 = res_0_norm
  end


  # do Newton's method if not converged
  @verbose5 @mpi_master print(BSTDOUT, "\n")

  #------------------------------------------------------------------------------------
  # Start of newton iteration loop
  eqn.params.time.t_newton += @elapsed for i=1:itermax

    @verbose5 @mpi_master println(BSTDOUT, "===== newton iteration: ", i)

    # Calculate Jacobian here
    jac_func(newton_data, mesh, sbp, eqn, opts, jac, ctx_residual, t)

    if use_jac_precond
      if ((i % recalc_prec_freq)) == 0 || i == 1
        jac_func(newton_data, pmesh, sbp, eqn, opts, jacp, ctx_residual, t)
      end
    end

    # apply globalization
    if globalize_euler

      if jac_type == 3 || jac_type == 4
        @verbose5 @mpi_master println(BSTDOUT, "applying Euler globalization to jacp")
        @verbose5 @mpi_master println(BSTDOUT, "tau = ", newton_data.tau_l)
        applyEuler(mesh, sbp, eqn, opts, newton_data, jacp)
      end

      if jac_type != 4
        @verbose5 @mpi_master println(BSTDOUT, "applying Euler gloablization to jac")
        applyEuler(mesh, sbp, eqn, opts, newton_data, jac)
      end
    end

#    checkJacVecProd(newton_data, mesh, sbp, eqn, opts, func, pert)

    # print as determined by options
    if write_jac && jac_type != 4
      if jac_type == 3
        PetscMatAssemblyBegin(jac, PETSC_MAT_FINAL_ASSEMBLY)
        PetscMatAssemblyEnd(jac, PETSC_MAT_FINAL_ASSEMBLY)
      end
      writedlm("jacobian$i.dat", full(jac))
      @mpi_master println(BSTDOUT, "finished printing jacobian"); flush(BSTDOUT)
    end
    
    # calculate Jacobian condition number
    if print_cond && ( jac_type == 1 || jac_type == 2)
      println(BSTDOUT, "calculating condition number of jacobian"); flush(BSTDOUT)
#      singular_vals = svdvals(full(jac))
#      writedlm("svd_vals.dat", singular_vals)
      cond_j = cond(full(jac))
      println(BSTDOUT, "Condition number of jacobian = ", cond_j); flush(BSTDOUT);
    end


    # if eigenvalues requested and we can calculate them
    if (( print_eigs || write_eigs) && (jac_type == 1 || jac_type == 2))
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
        writedlm("eigs_real$i.dat", real(eigs_i))
        writedlm("eigs_imag$i.dat", imag(eigs_i))
      end

    end

    # do the full eigenvalue decomposition
    # if requested and if julia owns the Jacobian matrix
    if write_eigdecomp && ( jac_type == 1 || jac_type == 2)
      println(BSTDOUT, "doing eigen decomposition"); flush(BSTDOUT)
      # make a dense jacobian so we can get *all* the eigenvalues and vectors
      # the algorithm for sparse matrices can get all - 2 of them, and might
      # have problems for matrices with a wide range of eigenvalues
      jac_dense = full(jac)
      D, V = eig(jac_dense)
      writedlm("eigdecomp_real$i.dat", real(D))
      writedlm("eigdecomp_imag$i.dat", imag(D))
#      writedlm("eigdecomp_realvecs$i.dat", real(V))
#      writedlm("eigdecomp_imagvecs$i.dat", imag(V))

      max_val = typemin(Float64)
      min_val = typemax(Float64)
    elseif write_eigdecomp # && we can't calculate it
      println(STDERR, "Warning: not performing eigen decomposition for jacobian of type $jac_type")

    end

    # negate res
    for j=1:m
      res_0[j] = -res_0[j]
    end

    # calculate Newton step
    flush(BSTDOUT)
    step_norm = matrixSolve(newton_data, eqn, mesh, opts, jac, delta_q_vec, res_0, BSTDOUT, verbose=verbose)
    #TODO: don't do this here, do it in the jac_func instead?
    if jac_type == 1 || jac_type == 2
      fill!(jac, 0.0)
    end
    
    # perform Newton update
    for j=1:m
      eqn.q_vec[j] += step_fac*delta_q_vec[j]
    end
   
#    saveSolutionToMesh(mesh, eqn.q_vec)
#    writeVisFiles(mesh, "newton_$i")

    @verbose5 eqn.majorIterationCallback(i, mesh, sbp, eqn, opts, BSTDOUT)


    # write starting values for next iteration to file
    if write_sol
      writedlm("q_vec$i_$myrank.dat", eqn.q_vec)
    end

    if write_q
      disassembleSolution(mesh, sbp, eqn, opts, eqn.q_vec)
      writedlm("q$i_$myrank.dat", eqn.q)
    end

    # calculate residual at updated location, used for next iteration rhs
    newton_data.res_norm_i_1 = newton_data.res_norm_i
    # extract real component to res_0
    res_0_norm = newton_data.res_norm_i = rhs_func(mesh, sbp, eqn, opts, rhs_vec, ctx_residual, t)
    for j=1:m
      res_0[j] = real(rhs_vec[j])
    end


    @verbose4 @mpi_master begin
      println(BSTDOUT, "residual norm = ", res_0_norm)
      println(BSTDOUT, "relative residual ", res_0_norm/res_reltol_0)
    end


    # write to convergence file
    @verbose5 @mpi_master begin
      flush(BSTDOUT)
      println(fconv, i, " ", res_0_norm, " ", step_norm)
      println(BSTDOUT, "printed to convergence.dat")
      flush(fconv)
    end
    flush(BSTDOUT)


#    tmp = i+1
    # write rhs to file
    if write_rhs
      fname = string("rhs", i, "_", myrank, ".dat")
      writedlm(fname, res_0)
    end

    if write_res
      fname = string("res", i, "_", myrank, ".dat")
      writedlm(fname, eqn.res)
    end


   if res_0_norm < res_abstol || res_0_norm/res_reltol_0 < res_reltol
     @mpi_master if res_0_norm < res_abstol 
       println(BSTDOUT, "Newton iteration converged with residual norm ", res_0_norm)
     end
     @mpi_master if res_0_norm/res_reltol_0 < res_reltol
      println(BSTDOUT, "Newton iteration converged with relative residual norm ", res_0_norm/res_reltol_0)
    end

     # put residual into rhs_vec
     for j=1:m
       rhs_vec[j] = res_0[j]
     end

     @mpi_master println("Iteration count: $i")
     @verbose5 @mpi_master close(fconv)

     flush(BSTDOUT)

     return nothing
    end  # end if tolerances satisfied

    if (step_norm < step_tol)
      @mpi_master println(BSTDOUT, "Newton iteration converged with step_norm = ", step_norm)
      @mpi_master println(BSTDOUT, "Final residual = ", res_0_norm)

      # put residual into rhs_vec
      for j=1:m
        rhs_vec[j] = res_0[j]
      end
      @verbose5 @mpi_master close(fconv)
      
      flush(BSTDOUT)

      return nothing
    end

#=
    # adjust step size limiter
    if (step_norm < step_norm_1)  # decreasing step size
      step_fac *= 1.2

      if step_fac > 1.0
        step_fac = 1.0
      end
    end

    if (step_norm > step_norm_1)
      step_fac /= 1.1
    end

    if step_norm < 0.001
      step_fac = 1.0
    end
=#
    # update globalization parameters
    if globalize_euler
      updateEuler(newton_data)
    end

    if jac_type == 3 || jac_type == 4
      updateKrylov(newton_data)
    end

    @verbose5 @mpi_master print(BSTDOUT, "\n")
    step_norm_1 = step_norm
    
  end  # end loop over newton iterations

  @mpi_master begin
    println(STDERR, "Warning: Newton iteration did not converge")

    println(BSTDOUT, "Warning: Newton iteration did not converge in ", itermax, " iterations")
    println(BSTDOUT, "  Final step size: ", step_norm)
    println(BSTDOUT, "  Final residual: ", res_0_norm)
    println(BSTDOUT, "  Final relative residual: ", res_0_norm/res_reltol_0)
    @verbose5 close(fconv); 
  end
  flush(BSTDOUT)

   # put residual into rhs_vec
   for j=1:m
     rhs_vec[j] = res_0[j]
   end
 
  return nothing

end               # end of function newton()
# TODO: more 'end of' comments inside newton()


