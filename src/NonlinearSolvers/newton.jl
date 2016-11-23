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

  # tuple to be passed to func
  ctx_newton

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
    vals_tmp = zeros(1, local_size) # values
    idx_tmp = zeros(PetscInt, local_size)  # row index
    idy_tmp = zeros(PetscInt, 1)  # column indices
    localsize = mesh.numNodesPerElement*mesh.numDofPerNode

    # NOTE: we are leaving ctx_newton uninitialized because 
    #   createPetscData needs it, but ctx_newton depends on its returned values

    return new(myrank, commsize, reltol, abstol, dtol, 
                      itermax, krylov_gamma, 
                      res_norm_i, res_norm_i_1, tau_l, tau_vec, 1, localsize, vals_tmp, 
                      idx_tmp, idy_tmp)
  end

end

include("residual_evaluation.jl")  # some functions for residual evaluation
include("petsc_funcs.jl")  # Petsc related functions

@doc """
### NonlinearSolvers.setupNewton
  
  alloc_rhs: keyword arg to allocate a new object or not for rhs_vec
                true (default) allocates a new vector
                false will use eqn.res_vec

  Allocates Jac & RHS

"""->
function setupNewton{Tsol, Tres}(mesh, pmesh, sbp, eqn::AbstractSolutionData{Tsol, Tres}, opts; alloc_rhs=true)

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
    ctx_newton = (jacp, x, b, ksp)
  end

  # now put ctx_newton into newton_data
  newton_data.ctx_newton = ctx_newton

  # For simple cases, especially for Newton's method as a steady solver,
  #   having rhs_vec and eqn.res_vec pointing to the same memory
  #   saves us from having to copy back and forth
  if alloc_rhs 
#     rhs_vec = deepcopy(eqn.res_vec)
    rhs_vec = zeros(eqn.res_vec)
  else
    rhs_vec = eqn.res_vec
  end

  # should be all zeros if alloc_rhs is true
  writedlm("rhs_initial.dat", rhs_vec)

  return newton_data, jac, rhs_vec

end   # end of setupNewton

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

  rhs_func = physicsRhs
  jac_func = physicsJac

  # allocate jac & rhs, and construct newton_data
  newton_data, jac, rhs_vec = setupNewton(mesh, pmesh, sbp, eqn, opts, alloc_rhs=false)

  ctx_residual = (func, )

  newtonInner(newton_data, mesh, sbp, eqn, opts, rhs_func, jac_func, jac, rhs_vec, ctx_residual, t,
                itermax=itermax, step_tol=step_tol, res_abstol=res_abstol, res_reltol=res_reltol, res_reltol0=res_reltol0)

  if jac_type == 3
    NonlinearSolvers.destroyPetsc(jac, newton_data.ctx_newton...)
  end

end

function newtonInner(newton_data::NewtonData, mesh::AbstractMesh, sbp::AbstractSBP, eqn::AbstractSolutionData,
                     opts, rhs_func, jac_func, jac, rhs_vec, ctx_residual=(), t=0.0;
                     itermax=10, step_tol=1e-6, res_abstol=1e-6, res_reltol=1e-6, res_reltol0=-1.0)

  # this function drives the non-linear residual to some specified tolerance
  # using Newton's Method
  # the jacobian is formed using finite differences
  # the initial condition is stored in eqn.q_vec
  # itermax is the maximum number of iterations
  # this function is type unstable for certain variables, but that's ok
  # the dispatch to the backslash solver and possibly the jacobian calculation
  # function will be runtime dispatched

  myrank = mesh.myrank
  fstdout = BufferedIO(STDOUT)

  @mpi_master println(fstdout, "\nEntered Newtons Method")
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
  @mpi_master begin
    println("step_tol = ", step_tol)
    println("res_abstol = ", res_abstol)
    println("res_reltol = ", res_reltol)
    println("res_reltol0 = ", res_reltol0)
  end

  if (t != 0.0) && (jac_type == 4)
    throw(ErrorException, "Petsc cannot be used for solving unsteady problems, see TODO in calcJacVecProd_wrapper")
  end


  if jac_method == 1  # finite difference
    pert = epsilon
  elseif jac_method == 2  # complex step
    pert = complex(0, epsilon)
  end


  if jac_type == 4
    # TODO: call Petsc MatShellSetContext her
    ctx_petsc = createPetscCtx(mesh, sbp, eqn, opts, newton_data, func)
  end

  if jac_type == 3 || jac_type == 4
    jacp = newton_data.ctx_newton[1]
  end

  Tjac = typeof(real(rhs_vec[1]))  # type of jacobian, residual
  m = length(rhs_vec)

  @mpi_master println(fstdout, "typeof(jac) = ", typeof(jac))

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
  @mpi_master begin
    _fconv = open("convergence.dat", "a+")
    fconv = BufferedIO(_fconv)
  end
  eqn.majorIterationCallback(0, mesh, sbp, eqn, opts, fstdout)

  # evaluating residual at initial condition
  @mpi_master println(fstdout, "evaluating residual at initial condition"); flush(fstdout)
  res_0_norm = newton_data.res_norm_i = calcResidual(mesh, sbp, eqn, opts, rhs_func, rhs_vec, ctx_residual, t)
  @mpi_master println(fstdout, "res_0_norm = ", res_0_norm); flush(fstdout)

  # extract the real components to res_0
  for i=1:m
    res_0[i] = real(rhs_vec[i])
  end
#   res_0 = deepcopy(real(rhs_vec))
  # TODO TODO TODO: changed to deepcopy 20161113
  writedlm("rhs_before_newton_loop.dat", res_0)

  @mpi_master begin
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
      @mpi_master println(fstdout, "Initial condition satisfied res_reltol with relative residual ", res_0_norm/res_reltol0)
      @mpi_master println(fstdout, "Residual ", res_0_norm)
    else
      @mpi_master println(fstdout, "Initial condition satisfies res_tol with residual norm ", res_0_norm)
    end
    # no need to assemble q into q_vec because it never changed

    @mpi_master close(fconv)
    flush(fstdout)

    @mpi_master println(fstdout, "Not entering Newton iteration loop")
    flush(fstdout)
    return nothing

  end  # end if tolerances satisfied

  @mpi_master println(fstdout, "res_reltol0 = ", res_reltol0)
  if res_reltol0 > 0  # use the supplied res_reltol0 value
    @mpi_master println(fstdout, "using supplied value for relative residual")
    res_reltol_0 = res_reltol0
    # Jared has a good explanation for why this is here
    # It has to do with choosing between user-supplied relative norm or the norm of the IC
  else
    @mpi_master println(fstdout, "using initial residual for relative residual")
    res_reltol_0 = res_0_norm
  end


  # do Newton's method if not converged
  @mpi_master print(fstdout, "\n")

  #------------------------------------------------------------------------------------
  # Start of newton iteration loop
  eqn.params.time.t_newton += @elapsed for i=1:itermax
    @mpi_master println(fstdout, "==================== Newton iteration: ", i)
    @mpi_master println(fstdout, "step_fac = ", step_fac)

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
        @mpi_master println(fstdout, "applying Euler globalization to jacp")
        @mpi_master println(fstdout, "tau = ", newton_data.tau_l)
        applyEuler(mesh, sbp, eqn, opts, newton_data, jacp)
      end

      if jac_type != 4
        @mpi_master println(fstdout, "applying Euler gloablization to jac")
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
      @mpi_master println(fstdout, "finished printing jacobian"); flush(fstdout)
    end
    
    # calculate Jacobian condition number
    if print_cond && ( jac_type == 1 || jac_type == 2)
      println(fstdout, "calculating condition number of jacobian"); flush(fstdout)
      cond_j = cond(full(jac))
      println(fstdout, "Condition number of jacobian = ", cond_j); flush(fstdout);
    end


    # if eigenvalues requested and we can calculate them
    if (( print_eigs || write_eigs) && (jac_type == 1 || jac_type == 2))
      eigs_i = reverse(eigvals(jac))
      if print_eigs
        println(fstdout, "eigenvalues =")
        for i=1:length(eigs_i)
          println(fstdout, eigs_i[i])
        end
      end
      flush(fstdout)

      if write_eigs
        writedlm("eigs$i.dat", eigs_i)
      end

    end

    # do the full eigenvalue decomposition
    # if requested and if julia owns the Jacobian matrix
    if write_eigdecomp && ( jac_type == 1 || jac_type == 2)
      println(fstdout, "doing eigen decomposition"); flush(fstdout)
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

#     println("========= in Newton, before delta_q_vec update. t = $t")
#     for dof_ix = 21461:21464
#       println("eqn.q_vec($dof_ix) = ", eqn.q_vec[dof_ix])
#     end
#     println(" ")
#     for dof_ix = 14997:15000
#       println("eqn.q_vec($dof_ix) = ", eqn.q_vec[dof_ix])
#     end

    println("----- in newtonInner, before linear solve -----")
    println("t: $t")
#     println("res_0[15]: ",res_0[15])
#     println("rhs_vec[15]: ",rhs_vec[15])
#     println("eqn.q_vec[15]: ",eqn.q_vec[15])
#     println("eqn.res_vec[15]: ",eqn.res_vec[15])
    println("-")
#     writedlm("jac_inside_newtonjl_iter_$i.dat", full(jac))
#     writedlm("rhs_inside_newtonjl_iter_$i.dat", full(res_0))

    # calculate Newton step
    flush(fstdout)
    if jac_type == 1 || jac_type == 2  # julia jacobian
      tmp, t_solve, t_gc, alloc = @time_all begin
        jac_f = factorize(jac)
        delta_q_vec[:] = jac_f\(res_0)  #  calculate Newton update
      end
      fill!(jac, 0.0)
#    @time solveMUMPS!(jac, res_0, delta_q_vec)
    elseif jac_type == 3 || jac_type == 4  # petsc jacobian
      # contents of ctx: (jacp, x, b, ksp)
      tmp, t_solve, t_gc, alloc = @time_all petscSolve(newton_data, jac, newton_data.ctx_newton..., opts, 
                                                       res_0, delta_q_vec, mesh.dof_offset)
    end

    println("=============+++++In newton+++++ i: ", i)
    println("=============+++++In newton+++++ t: ", t)
    println("=============+++++In newton+++++ norm(delta_q_vec): ", norm(delta_q_vec))
    println("=============+++++In newton+++++ norm(res_0): ", norm(res_0))

    eqn.params.time.t_solve += t_solve
    @mpi_master print(fstdout, "matrix solve: ")
    @mpi_master print_time_all(fstdout, t_solve, t_gc, alloc)
    step_norm = norm(delta_q_vec)
    #TODO: make this a regular reduce?
    step_norm = sqrt(MPI.Allreduce(step_norm*step_norm, MPI.SUM, mesh.comm))
    @mpi_master println(fstdout, "step_norm = ", step_norm)
    flush(fstdout)

    println("=============+++++In newton+++++ performing Newton update")
    # perform Newton update
    for j=1:m
      eqn.q_vec[j] += step_fac*delta_q_vec[j]
    end
    writedlm("q_vec_in_newtonjl_after_update-iter_$i.dat", eqn.q_vec)

#     println("========= in Newton, after delta_q_vec update. t = $t")
#     for dof_ix = 21461:21464
#       println("eqn.q_vec($dof_ix) = ", eqn.q_vec[dof_ix])
#     end
#     println(" ")
#     for dof_ix = 14997:15000
#       println("eqn.q_vec($dof_ix) = ", eqn.q_vec[dof_ix])
#     end
    
    eqn.majorIterationCallback(i, mesh, sbp, eqn, opts, fstdout)
 
    # write starting values for next iteration to file
    if write_sol
      writedlm("q_vec$i_$myrank.dat", eqn.q_vec)
    end

    if write_q
      disassembleSolution(mesh, sbp, eqn, opts, eqn.q_vec)
      writedlm("q$i_$myrank.dat", eqn.q)
    end

    println("=============+++++In newton+++++ calling calcResidual")
    # calculate residual at updated location, used for next iteration rhs
    newton_data.res_norm_i_1 = newton_data.res_norm_i
    # extract real component to res_0
    res_0_norm = newton_data.res_norm_i = calcResidual(mesh, sbp, eqn, opts, rhs_func, rhs_vec, ctx_residual, t)
    # TODO: should this be a res_vec instead of rhs_vec? 20161116
    for j=1:m
      res_0[j] = real(rhs_vec[j])
    end
#     res_0 = deepcopy(real(rhs_vec))
    # TODO TODO TODO: changed to deepcopy 20161113



    @mpi_master begin
      println(fstdout, "residual norm = ", res_0_norm)
      println(fstdout, "relative residual ", res_0_norm/res_reltol_0)
    end


    # write to convergence file
    @mpi_master begin
      flush(fstdout)
      println(fconv, i, " ", res_0_norm, " ", step_norm)
      println(fstdout, "printed to convergence.dat")
      flush(fconv)
    end
    flush(fstdout)


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
       println(fstdout, "Newton iteration converged with residual norm ", res_0_norm)
     end
     @mpi_master if res_0_norm/res_reltol_0 < res_reltol
      println(fstdout, "Newton iteration converged with relative residual norm ", res_0_norm/res_reltol_0)
    end

     # put residual into rhs_vec
     for j=1:m
       rhs_vec[j] = res_0[j]
     end

     @mpi_master close(fconv)

     flush(fstdout)

     return nothing
    end  # end if tolerances satisfied

    if (step_norm < step_tol)
      @mpi_master println(fstdout, "Newton iteration converged with step_norm = ", step_norm)
      @mpi_master println(fstdout, "Final residual = ", res_0_norm)

      # put residual into rhs_vec
      for j=1:m
        rhs_vec[j] = res_0[j]
      end
      @mpi_master close(fconv)
      
      flush(fstdout)

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

    @mpi_master print(fstdout, "\n")
    step_norm_1 = step_norm
    
  end  # end loop over newton iterations

  @mpi_master begin
    println(STDERR, "Warning: Newton iteration did not converge")

    println(fstdout, "Warning: Newton iteration did not converge in ", itermax, " iterations")
    println(fstdout, "  Final step size: ", step_norm)
    println(fstdout, "  Final residual: ", res_0_norm)
    println(fstdout, "  Final relative residual: ", res_0_norm/res_reltol_0)
    close(fconv); 
  end
  flush(fstdout)

   # put residual into rhs_vec
   for j=1:m
     rhs_vec[j] = res_0[j]
   end
 
  return nothing

end               # end of function newton()
# TODO: more 'end of' comments inside newton()

@doc """
###NonlinearSolver.physicsJac
  
  Jacobian (of the physics) calculation, separate from the Newton function

  ctx_residual: func must be the first element. func is the residual evaluation function, i.e. evalEuler
  t:            simulation time

  For solving steady problems, this function can be used directly by newtonInner 
    for the Jac calculation.

  For solving unsteady problems, this function can be used as a building block 
    for the calculation of the time-marching Jac.

"""->
function physicsJac(newton_data::NewtonData, mesh, sbp, eqn, opts, jac, ctx_residual, t=0.0; is_preconditioned::Bool=false)

#   DEBUG = false
  DEBUG = true

  fstdout = BufferedIO(STDOUT)

  println(fstdout, "in physicsJac")

  myrank = mesh.myrank
  fstdout = BufferedIO(STDOUT)

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

  if jac_method == 1  # finite difference
    pert = epsilon
  elseif jac_method == 2  # complex step
    pert = complex(0, epsilon)
  end

  # ctx_residual: func must be the first element
  func = ctx_residual[1]

  #----------------------------------------------------------------------
  # Calculate Jacobian using selected method 
  print_jacobian_timing = true
  eqn.params.time.t_jacobian += @elapsed if jac_method == 1
    @mpi_master println(fstdout, "calculating finite difference jacobian")

    if jac_type == 1  # dense jacobian
      @mpi_master println(fstdout, "calculating dense FD jacobian")
      # TODO: need to make q/q_vec and res/res_vec consistent
      arrToVecAssign(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
      # res_0 is the unperturbed res, and it needs to be passed in vector form
      assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
      res_copy_vec = copy(eqn.res_vec)
      #TODO: don't copy the giant vector!
      tmp, t_jac, t_gc, alloc = @time_all calcJacFD(newton_data, mesh, sbp, eqn, opts, func, res_copy_vec, pert, jac, t)

    elseif jac_type == 2  # Julia sparse jacobian
      @mpi_master println(fstdout, "calculating sparse FD jacobian")
      #TODO: don't copy the giant array!
      res_copy = copy(eqn.res)  # copy unperturbed residual

      tmp, t_jac, t_gc, alloc = @time_all calcJacobianSparse(newton_data, mesh, sbp, eqn, opts, func, res_copy, pert, jac, t)
    elseif jac_type == 3  # Petsc sparse jacobian
      @mpi_master println(fstdout, "calculating sparse FD jacobian")
      res_copy = copy(eqn.res)  # copy unperturbed residual
      tmp, t_jac, t_gc, alloc = @time_all calcJacobianSparse(newton_data, mesh, sbp, eqn, opts, func, res_copy, pert, jac, t)
    elseif jac_type == 4  # Petsc jacobian-vector product
      throw(ErrorException("No handling of jac_method = 1 and jac_type = 4: 
                           finite differencing isn't permitted for Petsc mat-free"))
    end

  elseif jac_method == 2
    @mpi_master println(fstdout, "calculating complex step jacobian")

    if jac_type == 1  # dense jacobian
      @mpi_master println(fstdout, "calculating dense complex step jacobian")
      tmp, t_jac, t_gc, alloc = @time_all calcJacobianComplex(newton_data, mesh, sbp, eqn, opts, func, pert, jac, t)
    elseif jac_type == 2  # Julia sparse jacobian 
      @mpi_master println(fstdout, "calculating sparse complex step jacobian")
      res_dummy = Array(Float64, 0, 0, 0)  # not used, so don't allocation memory
      tmp, t_jac, t_gc, alloc = @time_all calcJacobianSparse(newton_data, mesh, sbp, eqn, opts, func, res_dummy, pert, jac, t)
    elseif jac_type == 3 # Petsc sparse jacobian
      res_dummy = Array(Float64, 0, 0, 0)  # not used, so don't allocation memory
      @mpi_master println(fstdout, "calculating explicit Petsc jacobian")

      @mpi_master println(fstdout, "calculating main jacobain")
      tmp, t_jac, t_gc, alloc = @time_all calcJacobianSparse(newton_data, mesh, sbp, eqn, opts, func, res_dummy, pert, jac, t)

    elseif jac_type == 4 # Petsc jacobian-vector product
      # calculate preconditioner matrix only
      res_dummy = Array(Float64, 0, 0, 0)
      print_jacobian_timing = false

      # if jac_method == 2 (CS) and jac_type == 4 (Petsc mat-free), only calc the jac if it is a preconditioned jac
      if is_preconditioned
        print_jacobian_timing = true
        tmp, t_jac, t_gc, alloc = @time_all calcJacobianSparse(newton_data, mesh, sbp, eqn, opts, func, res_dummy, pert, jac, t)
      end

    end   # end of jac_type check

  end  # end of jac_method check

  # TODO: all printing should actually be handled outside of this function
  if print_jacobian_timing
    @mpi_master print(fstdout, "jacobian calculation: ")
    @mpi_master print_time_all(fstdout, t_jac, t_gc, alloc)
  end

  flush(fstdout)

  return nothing

end   # end of physicsJac function

#------------------------------------------------------------------------------
# jacobian vector product functions
#------------------------------------------------------------------------------
@doc """
###NonlinearSolver.calcJacVecProd

  This function calculates a Jacobian vector product Ax=b without explicitly 
  computing the Jacobian.

  The Jacobian refers to the Jacobian of the point stored in eqn.q_vec
  Inputs:
    newton_data:  NewtonData object
    mesh: AbstractMesh
    sbp:  SBP operator
    eqn:  AbstractEquation object
    opts: options dictionary
    pert: perturbation to use for the algorithmic differentiation.  Currently,
          only complex numbers are supported.
    func: residual evaluation function
    vec:  the x vector in Ax=b.  Can be AbstractVector type.
    b:    location to store the result (the b an Ax=b).  Can be any
          AbstractVector type

  Outputs:
    none


  Aliasing restrictions: vec, b, and eqn.q must not alias each other.

"""
function calcJacVecProd(newton_data::NewtonData, mesh, sbp, eqn, opts, pert, func::Function, vec::AbstractVector, b::AbstractVector, t=0.0)
# calculates the product of the jacobian with the vector vec using a directional
# derivative
# only intended to work with complex step
# might also work for finite differences
# vec is the vector the jacobian is multiplied by
# the result is stored in b
# the product uses the jacobian at the point stored in eqn.q_vec
# ie. J(eqn.q)*v = b = imag(J(u + pert*v))/pert
# pert is the perturbation, either real or complex for finite difference or 
# complex step, respectively
# func is the residual evaluation function
 
  itr = eqn.params.krylov_itr
  globalize_euler = opts["newton_globalize_euler"]::Bool

  epsilon = imag(pert)  # magnitude of perturbationa

  # apply perturbation
  for i=1:mesh.numDof
    eqn.q_vec[i] += pert*vec[i]
  end

  # scatter into eqn.q
  disassembleSolution(mesh, sbp, eqn, opts, eqn.q_vec) 
  func(mesh, sbp, eqn, opts, t)

  # gather into eqn.res_vec
  assembleResidual(mesh, sbp, eqn, opts, eqn.res_vec, assemble_edgeres=opts["use_edge_res"])
  
  # calculate derivatives, store into b
  calcJacCol(b, eqn.res_vec, epsilon)

  if globalize_euler
    applyEuler(mesh, sbp, eqn, opts, vec, newton_data, b)
  end

  # undo perturbation
  for i=1:mesh.numDof
    eqn.q_vec[i] -= pert*vec[i]
  end

  eqn.params.krylov_itr += 1

  return nothing
end


@doc """
### NonlinearSolvers.checkJacVecProd

  This function calculates a jacobian vector product, then computes the 
  Jacobian and multiplies it by the vector and compares the results.
  This is very expensive, only used for debugging.

  The Jacobian is calculated about the point stored in eqn.q_vec


  Inputs:
    newton_data:  NewtonData object
    mesh: AbstractMesh object
    sbp: SBP operator
    opts: options dictonary
    func: residual evaluation function
    pert: perturbation used to calculate the Jacobian.  Can be real (for 
          finite differences), or some AD datatype.

    Aliasing restrictions: none
"""->
function checkJacVecProd(newton_data::NewtonData, mesh, sbp, eqn, opts, func, pert, t=0.0)
  
  v = ones(mesh.numDof)
  result1 = zeros(mesh.numDof)
  calcJacVecProd(newton_data, mesh, sbp, eqn, opts, pert, func, v, result1, t)
  jac = SparseMatrixCSC(mesh.sparsity_bnds, Float64)
  res_dummy = []

  disassembleSolution(mesh, sbp, eqn,opts, eqn.q_vec)
  @time calcJacobianSparse(newton_data, mesh, sbp, eqn, opts, func, res_dummy, pert, jac, t)

  if opts["newton_globalize_euler"]
    applyEuler(mesh, sbp, eqn, opts, newton_data, jac)
  end
  result2 = jac*v

  cnt = 0
  for i=1:mesh.numDof
    if abs(result1[i] - result2[i]) > 1e-14
      cnt += 1
    end
  end

  if cnt != 0
    println(STDERR, "Warning: jacobian vector product check 1 failed")
    println("cnt = ", cnt)
    result_diff = result1 - result2
    diff_norm = calcNorm(eqn, result_diff)
    println("diff norm = ", diff_norm)
    println("result_diff = ", result_diff)
  else
    println("jacobian vector product check 1 passed")
  end

  # check some additional vectors
  for j=2:4
#  v2 = collect(1:mesh.numDof)
   v2 = linspace(j, j+1, mesh.numDof)
  result3 = jac*v2

  result4 = zeros(mesh.numDof)

  calcJacVecProd(newton_data, mesh, sbp, eqn, opts, pert, func, v2, result4, t)

  cnt = 0
  for i=1:mesh.numDof
    if abs(result3[i] - result4[i]) > 1e-14
      cnt += 1
    end
  end

  if cnt != 0
    println(STDERR, "Warning: jacobian vector product check $j failed")
    println("cnt = ", cnt)
    result_diff = result3 - result4
    diff_norm = calcNorm(eqn, result_diff)
    println("diff norm = ", diff_norm)
#    println("result_diff = ", result_diff)
  else
    println("jacobian vector product check $j passed")
  end

end


  return nothing
end


@doc """
### NonlinearSolvers.calcJacVecProd_wrapper

  This function is passed to Petsc so it can calculate Jacobian-vector products
  Ax=b.

  Inputs
    A:  PetscMat object 
    x:  PetscVec to multiply A with

  Inputs/Outputs:
    b:  PetscVec to store result in

  Aliasing restrictions: see Petsc documentation

"""->
function calcJacVecProd_wrapper(A::PetscMat, x::PetscVec, b::PetscVec)
# calculate Ax = b

  # TODO 20161102: this needs a 't' argument

#  println("entered calcJacVecProd wrapper")
  # get the context
  # the context is a pointer to a tuple of all objects needed
  # for a residual evaluation
  ctx = MatShellGetContext(A)
  tpl = unsafe_pointer_to_objref(ctx)
  # unpack the tuple (could use compact syntax)
  mesh = tpl[1]
  sbp = tpl[2]
  eqn = tpl[3]
  opts = tpl[4]
  newton_data = tpl[5]
  func = tpl[6]

  epsilon =  opts["epsilon"]
  jac_method = opts["jac_method"]

  # the perturbation had better match the type of the eqn object
  if jac_method == 1  # finite difference
    pert = epsilon
  elseif jac_method == 2  # complex step
    pert = complex(0, epsilon)
  end

  # get the arrays underlying x and b
  x_arr, xptr = PetscVecGetArrayRead(x)  # read only
  b_arr, bptr = PetscVecGetArray(b)  # writeable

  calcJacVecProd(newton_data, mesh, sbp, eqn, opts, pert, func, x_arr, b_arr)

#  println("finished calculating JacVecProd")
  PetscVecRestoreArrayRead(x, xptr)
  PetscVecRestoreArray(b, bptr)


  return PetscErrorCode(0)
end

#------------------------------------------------------------------------------
# Functions for calculating the Jacobian
#------------------------------------------------------------------------------
@doc """
### NonlinearSolvers.calcJacFD

  This function calculates the Jacobian using finite differences, perturbing
  one degree of freedom at a time.  This is slow and not very accurate.  
  The Jacobian is calculated about the point in eqn.q_vec.

  Inputs:
    newton_data:  NewtonData object
    mesh: AbstractMesh
    sbp:  SBP operator
    eqn:  AbstractEquation object
    opts: options dictionary
    pert: perturbation to use the finite differences.  Must be of type Tsol.
    func: residual evaluation function
    res_0: vector containing residual at the point the Jacobian is calculated

  Inputs/Outputs:
    jac:: Jacobian matrix to be populated.  Must be a dense matrix

  Aliasing restrictions: res_0 must not alias eqn.res_vec

  At the start, calcJacFD assumes:
    The Jacobian will be calculated at the state that is specified in eqn.q_vec .
    res_0 should have the residual at that state in it

  At exit, eqn.q_vec will have the same values as at the start.

  eqn.q and eqn.res will be overwritten in the course of this function.

"""->
function calcJacFD(newton_data::NewtonData, mesh, sbp, eqn, opts, func, res_0, pert, jac::DenseArray, t=0.0)
# calculate the jacobian using finite difference
  (m,n) = size(jac)
  entry_orig = zero(eltype(eqn.q_vec))
  epsilon = norm(pert)  # finite difference perturbation
  # calculate jacobian
  for j=1:m
    if j==1
      entry_orig = eqn.q_vec[j]
      eqn.q_vec[j] +=  epsilon
    else
      eqn.q_vec[j-1] = entry_orig # undo previous iteration pertubation
      entry_orig = eqn.q_vec[j]
      eqn.q_vec[j] += epsilon
    end


    disassembleSolution(mesh, sbp, eqn, opts, eqn.q_vec)
    # evaluate residual
    func(mesh, sbp, eqn, opts, t)

    assembleResidual(mesh, sbp, eqn, opts,  eqn.res_vec)
    calcJacCol(sview(jac, :, j), res_0, eqn.res_vec, epsilon)
    
  end

  # undo final perturbation
  eqn.q_vec[m] = entry_orig


  return nothing
end


@doc """
### NonlinearSolvers.calcJacComplex

  This function calculates the Jacobian (dense) using the complex step method, 
  perturbing one degree of freedom at a time.  This is very slow.  The jacobian
  is calculated about the point in eqn.q_vec.

  Inputs:
    newton_data:  NewtonData object
    mesh: AbstractMesh
    sbp:  SBP operator
    eqn:  AbstractEquation object
    opts: options dictionary
    pert: perturbation to use.  Must be of type Tsol.
    func: residual evaluation function

  Inputs/Outputs:
    jac:: Jacobian matrix to be populated.  Must be a dense matrix

  Aliasing restrictions: res_0 must not alias eqn.res_vec
"""->
function calcJacobianComplex(newton_data::NewtonData, mesh, sbp, eqn, opts, func, pert, jac, t=0.0)

  epsilon = norm(pert)  # complex step perturbation
  entry_orig = zero(eltype(eqn.q_vec))
  (m,n) = size(jac)
  # calculate jacobian
  for j=1:m
    if j==1
      entry_orig = eqn.q_vec[j]
      eqn.q_vec[j] +=  pert
    else
      eqn.q_vec[j-1] = entry_orig # undo previous iteration pertubation
      entry_orig = eqn.q_vec[j]
      eqn.q_vec[j] += pert
    end

    disassembleSolution(mesh, sbp, eqn, opts, eqn.q_vec)
    # evaluate residual
    func(mesh, sbp, eqn, opts, t)

    assembleResidual(mesh, sbp, eqn, opts, eqn.res_vec)
    calcJacCol(sview(jac, :, j), eqn.res_vec, epsilon)
    
  end  # end loop over rows of jacobian


  # undo final perturbation
  eqn.q_vec[m] = entry_orig
#

  return nothing
end



@doc """
### NonlinearSolvers.calcJacobianSparse

  This function calculate the Jacobian sparsely (only the entries 
    within the sparsity bounds), using either finite differences or algorithmic 
    differentiation.  The jacobian is calculated about the point stored in 
    eqn.q (not eqn.q_vec).

  Inputs:
    newton_data:  NewtonData object
    mesh: AbstractMesh
    sbp:  SBP operator
    eqn:  AbstractEquation object
    opts: options dictionary
    pert: perturbation to use for the algorithmic differentiation.  Currently,
          only complex numbers are supported.
    func: residual evaluation function
    res_0: element-based (3 dimensional) array containing the residual evaluated           at the point where the Jacobian is being calculated.
           This is only used for finite differences (can be a 0 x 0 x 0 array
           otherwise).

  Inputs/Outputs:
    jac:  Jacobian matrix.  Must be a sparse matrix type of some kind, 
          (including PetscMat).

  Aliasing restrictions: res_0 must not alias eqn.res

 
"""->
function calcJacobianSparse(newton_data::NewtonData, mesh, sbp, eqn, opts, func, res_0::Abstract3DArray, pert, jac::Union{SparseMatrixCSC, PetscMat}, t=0.0)
# res_0 is 3d array of unperturbed residual, only needed for finite difference
# pert is perturbation to apply
# this function is independent of perturbation type

#  filter_orig = eqn.params.use_filter  # record original filter state
#  eqn.params.use_filter = false  # don't repetatively filter

  epsilon = norm(pert)  # get magnitude of perturbation
  m = length(res_0)
  myrank = mesh.myrank
  f = eqn.params.f
  time = eqn.params.time
  time.t_color += @elapsed for color=1:mesh.maxColors  # loop over max colors, only do calculation for numColors
    for j=1:mesh.numNodesPerElement  # loop over nodes 
      for i=1:mesh.numDofPerNode  # loop over dofs on each node

	# apply perturbation to q
        if color <= mesh.numColors
          applyPerturbation(mesh, eqn.q, eqn.q_face_recv, color, pert, i, j,f)
  	  # evaluate residual
          time.t_func += @elapsed func(mesh, sbp, eqn, opts, t)
        end

        if !(color == 1 && j == 1 && i == 1) 
          PetscMatAssemblyEnd(jac, PETSC_MAT_FLUSH_ASSEMBLY)
        end

	# assemble res into jac
        if color <= mesh.numColors
          time.t_insert += @elapsed for k=1:mesh.numEl  # loop over elements in residual
            el_pert = mesh.pertNeighborEls[k, color] # get perturbed element
            #TODO: find a way to get rid of this if statement
            # Solution: make pertNeighbor Els only hold the perturbed elements
            if el_pert != 0   # if element was actually perturbed for this color

              col_idx = mesh.dofs[i, j, el_pert]  # = dof_pert
              #TODO: make an immutable type to hold the bookeeping info
              assembleElement(newton_data, mesh, eqn, eqn.res, res_0, k, el_pert, col_idx, epsilon, jac)
            end  # end if el_pert != 0
          end  # end loop over k

          # now do res_edge, if needed
          for edge = 1:size(eqn.res_edge, 4)
            res_edge = sview(eqn.res_edge, :, :, :, edge)
            for k=1:mesh.numEl  # loop over elements in residual
              el_pert = mesh.pertNeighborEls_edge[k, edge] # get perturbed element
              if el_pert != 0   # if element was actually perturbed for this color

                col_idx = mesh.dofs[i, j, el_pert] # = dof_pert
                #TODO: make an immutable type to hold the bookeeping info
                assembleElement(newton_data, mesh, eqn, res_edge, res_0, k, el_pert, col_idx, epsilon, jac)
              end  # end if el_pert != 0
            end  # end loop over k
          end  # end loop over local edges
        end

        PetscMatAssemblyBegin(jac, PETSC_MAT_FLUSH_ASSEMBLY)

        # undo perturbation
        if color <= mesh.numColors
          applyPerturbation(mesh, eqn.q, eqn.q_face_recv, color, -pert, i, j)
        end
      end  # end loop i
    end  # end loop j
  end  # end loop over colors

  PetscMatAssemblyEnd(jac, PETSC_MAT_FLUSH_ASSEMBLY)
#  flush(f)
  # now jac is complete
#  eqn.params.use_filter = filter_orig # reset filter
  return nothing

end  # end function


@doc """
### NonlinearSolvers.applyPerturbation

  This function applies a perturbation to a the specified degree of freedom
  on each element according to a mask.

  Inputs:
    mesh: an AbstractMesh
    color: the color to perturb
    pert: perturbation to apply.  Can be any datatype
    i: local degree of freedom number (in range 1:numDofPerNode) to perturb
    j: local node number (in range 1:numNodesPerElement) to perturb

  Inputs/Outputs:
    arr: element based (3D) array of values to perturb
    arr_shared: array of 3D arrays for ghost elements to be perturbed

  Aliasing restrictions: none
"""->
function applyPerturbation{T}(mesh::AbstractMesh, arr::Abstract3DArray, arr_shared::Array{Array{T, 3}, 1},  color::Integer, pert, i, j, f=STDOUT; perturb_shared=true)
  # applys perturbation pert to array arr according to a mask
  # color is the color currently being perturbed, used to select the mask
  # i, j specify the dof, node number within arr
  # the length of mask must equal the third dimension of arr
  # this function is independent of the type of pert

  @assert i <= size(arr, 1)
  @assert j <= size(arr, 2)

  (ndof, nnodes, numel) = size(arr)
  mask = mesh.color_masks[color]
  for k=1:numel
    arr[i, j, k] += pert*mask[k]
  end

  if perturb_shared
    for peer=1:mesh.npeers
      mask_i = mesh.shared_element_colormasks[peer][color]
      arr_i = arr_shared[peer]
      for k=1:length(mask_i)
        arr_i[i, j, k] += pert*mask_i[k]
      end
    end
  end

  return nothing
end




#------------------------------------------------------------------------------
# Helper functions: assembleElement, calcJacCol for finite difference
#------------------------------------------------------------------------------
@doc """
### NonlinearSolver.assembleElement
  
  This function adds the contribution of one element to the Petsc Jacobian, if
  finite differences were used to calculate the entries.

  Inputs:
    newton_data: a NewtonData object
    mesh:  AbstractMesh object
    eqn:  AbstractEquation
    res_arr: element-based (3D) array of perturbed residual values
    res_0:  element-based (3D) array of non-perturbed residual values
    el_res: element number of the element we are observing the change in
    el_pert: element number of the element that was perturbed
    dof_pert: the degree of freedom number of the perturbed dof
    epsilon: magnitude of perturbation

  Inputs/Outputs:
    jac: Petsc Jacobian object

  Aliasing restrictions: res_arr and res_0 must not alias each other.

"""->
# finite difference Petsc
function assembleElement{Tsol <: Real}(newton_data::NewtonData, mesh, eqn::AbstractSolutionData{Tsol}, res_arr, res_0, el_res::Integer, el_pert::Integer, dof_pert::Integer, epsilon, jac::PetscMat)
# assemble an element contribution into jacobian
# making this a separate function enables dispatch on type of jacobian
# el_res is the element in the residual to assemble
# el_pert is the element that was perturbed
# dof_pert is the dof number (global) of the dof that was perturbed
# nodenum_local is the local node number of the perturbed node
# dof_pert_local is the dofnumber local to the node of the perturbed dof
# typically either el_pert or dof_pert will be needed, not both


  # resize array
  # basically a no-op if array is already the right size
  local_size = PetscInt(mesh.numNodesPerElement*mesh.numDofPerNode)

  # get row number
  newton_data.idy_tmp[1] = dof_pert - 1 + mesh.dof_offset

  pos = 1
  for j_j = 1:mesh.numNodesPerElement
    for i_i = 1:mesh.numDofPerNode
      newton_data.idx_tmp[pos] = mesh.dofs[i_i, j_j, el_res] - 1 + mesh.dof_offset
  
      tmp = (res_arr[i_i,j_j, el_res] - res_0[i_i, j_j, el_res])/epsilon
      newton_data.vals_tmp[pos] = tmp
      pos += 1
    end
  end

  PetscMatSetValues(jac, newton_data.idx_tmp, newton_data.idy_tmp, newton_data.vals_tmp, PETSC_ADD_VALUES)

  return nothing

end


@doc """
### NonlinearSolver.assembleElement
  
  This function adds the contribution of one element to the Julia 
  SparseMatrixCSC Jacobian, if finite differences were used to calculate the 
  entries.

  Inputs:
    newton_data: a NewtonData object
    mesh:  AbstractMesh object
    eqn:  AbstractEquation
    res_arr: element-based (3D) array of perturbed residual values
    res_0:  element-based (3D) array of non-perturbed residual values
    el_res: element number of the element we are observing the change in
    el_pert: element number of the element that was perturbed
    dof_pert: the degree of freedom number of the perturbed dof
    epsilon: magnitude of perturbation

  Inputs/Outputs:
    jac: SparseMatrixCSC Jacobian object

  Aliasing restrictions: res_arr and res_0 must not alias each other.

"""->
#finite difference SparseMatrixCSC
function assembleElement{Tsol <: Real}(newton_data::NewtonData, mesh, eqn::AbstractSolutionData{Tsol}, res_arr, res_0, el_res::Integer, el_pert::Integer, dof_pert::Integer, epsilon, jac::SparseMatrixCSC)
# assemble an element contribution into jacobian
# making this a separate function enables dispatch on type of jacobian
# el_res is the element in the residual to assemble
# el_pert is the element that was perturbed
# dof_pert is the dof number (global) of the dof that was perturbed
# typically either el_pert or dof_pert will be needed, not both

  for j_j = 1:mesh.numNodesPerElement
    for i_i = 1:mesh.numDofPerNode
      row_idx = mesh.dofs[i_i, j_j, el_res] + mesh.dof_offset

      tmp = (res_arr[i_i,j_j, el_res] - res_0[i_i, j_j, el_res])/epsilon
      jac[row_idx, dof_pert + mesh.dof_offset] += tmp

    end
  end
  
  return nothing

end


@doc """
### NonlinearSolvers.calcJacCol

  This function extracts the entries for one column of the Jacobian from two residual evaluates that come from finite differences.

  Inputs:
    res_0: vector of unperturbed residual values
    res: vector of perturbed residual values
    epsilon: magnitude of perturbation

  Inputs/Outputs:
    jac_row = vector to be populated with the Jacobian entries

  Aliasing restrictions: res_0 and res cannot alias (obviously).

"""->
function calcJacCol{T <: Real}(jac_row, res_0, res::AbstractArray{T,1}, epsilon)
# calculate a row of the jacobian from res_0, the function evaluated 
# at the original point, and res, the function evaluated at a perturbed point

  m = length(res_0)
  for i=1:m
    jac_row[i] = (res[i] - res_0[i])/epsilon
  end

  return nothing

end



#------------------------------------------------------------------------------
# helper functions: assembleElement, calcJacCol for complex numbers
#------------------------------------------------------------------------------

@doc """
### NonlinearSolver.assembleElement
  
  This function adds the contribution of one element to the PetscMat Jacobian, 
  if algorthmic differentiation was used to calculate the entries.

  Inputs:
    newton_data: a NewtonData object
    mesh:  AbstractMesh object
    eqn:  AbstractEquation
    res_arr: element-based (3D) array of perturbed residual values
    res_0:  element-based (3D) array of non-perturbed residual values (not used)
    el_res: element number of the element we are observing the change in
    el_pert: element number of the element that was perturbed
    dof_pert: the degree of freedom number of the perturbed dof
    epsilon: magnitude of perturbation

  Inputs/Outputs:
    jac: PetscMat Jacobian object

  Aliasing restrictions: res_arr and res_0 must not alias each other.

"""->
# for complex numbers
function assembleElement{Tsol <: Complex}(newton_data::NewtonData, mesh, eqn::AbstractSolutionData{Tsol}, res_arr, res_0,  el_res::Integer, el_pert::Integer, dof_pert::Integer, epsilon, jac::PetscMat)
# assemble an element contribution into jacobian
# making this a separate function enables dispatch on type of jacobian
# el_res is the element in the residual to assemble
# el_pert is the element that was perturbed
# dof_pert is the dof number (global) of the dof that was perturbed
# typically either el_pert or dof_pert will be needed, not both

# get row number
newton_data.idy_tmp[1] = dof_pert - 1 + mesh.dof_offset
pos = 1
for j_j = 1:mesh.numNodesPerElement
  for i_i = 1:mesh.numDofPerNode
    newton_data.idx_tmp[pos] = mesh.dofs[i_i, j_j, el_res] - 1 + mesh.dof_offset
    row = newton_data.idx_tmp[pos]
    col = newton_data.idy_tmp[1]


    newton_data.vals_tmp[pos] = imag(res_arr[i_i,j_j, el_res])/epsilon
    val = newton_data.vals_tmp[pos]
    pos += 1
  end
end

PetscMatSetValues(jac, newton_data.idx_tmp, newton_data.idy_tmp, newton_data.vals_tmp, PETSC_ADD_VALUES)

return nothing

end


@doc """
### NonlinearSolver.assembleElement
  
  This function adds the contribution of one element to the Julia 
  SparseMatrixCSC Jacobian, if algorthmic differentiation was used to 
  calculate the entries.

  Inputs:
    newton_data: a NewtonData object
    mesh:  AbstractMesh object
    eqn:  AbstractEquation
    res_arr: element-based (3D) array of perturbed residual values
    res_0:  element-based (3D) array of non-perturbed residual values
    el_res: element number of the element we are observing the change in
    el_pert: element number of the element that was perturbed
    dof_pert: the degree of freedom number of the perturbed dof
    epsilon: magnitude of perturbation

  Inputs/Outputs:
    jac: SparseMatrixCSC Jacobian object

  Aliasing restrictions: res_arr and res_0 must not alias each other.

"""->
function assembleElement{Tsol <: Complex}(newton_data::NewtonData, mesh, eqn::AbstractSolutionData{Tsol}, res_arr, res_0,  el_res::Integer, el_pert::Integer, dof_pert::Integer, epsilon, jac::SparseMatrixCSC)
# assemble an element contribution into jacobian
# making this a separate function enables dispatch on type of jacobian
# el_res is the element in the residual to assemble
# el_pert is the element that was perturbed
# dof_pert is the dof number (global) of the dof that was perturbed
# typically either el_pert or dof_pert will be needed, not both

#println(" element $el_res res = ", sview(res_arr, :, :, el_res))

  for j_j = 1:mesh.numNodesPerElement
    for i_i = 1:mesh.numDofPerNode
      row_idx = mesh.dofs[i_i, j_j, el_res] + mesh.dof_offset
      jac[row_idx, dof_pert + mesh.dof_offset] += imag(res_arr[i_i,j_j, el_res])/epsilon
    end
  end

  return nothing

end


 @doc """
### NonlinearSolvers.calcJacCol

  This function extracts the entries for one column of the Jacobian from a 
  complex step residual evaluation

  Inputs:
    res: vector of perturbed residual values
    epsilon: magnitude of perturbation

  Inputs/Outputs:
    jac_row = vector to be populated with the Jacobian entries

  Aliasing restrictions: none

"""->
function calcJacCol{T <: Complex}(jac_row, res::AbstractArray{T, 1}, epsilon)
# calculate a row of the jacobian from res_0, the function evaluated 
# at the original point, and res, the function evaluated at a perturbed point

  m = length(res)

  for i=1:m
    jac_row[i] = imag(res[i])/epsilon
  end

  return nothing

end
