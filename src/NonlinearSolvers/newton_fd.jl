export newton, newton_check, newton_check_fd, initializeTempVariables, calcResidual

@doc """
  This type holds all the data the might be needed for Newton's method,
  including globalization.  This simplifies the data handling and 
  the C callback used by Petsc
"""->
type NewtonData{Tsol, Tres}
 
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



  # use default inner constructor
end

function NewtonData(mesh, sbp, eqn, opts)

  println("entered NewtonData constructor")
  println("typeof(eqn) = ", typeof(eqn))
  Tsol = eltype(eqn.q)
  Tres = eltype(eqn.res)

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

  println("creating NewtonData object, tau = ", tau_l)
  println("opts[euler_tau] = ", opts["euler_tau"])

  return NewtonData{Tsol, Tres}(reltol, abstol, dtol, itermax, krylov_gamma, res_norm_i, res_norm_i_1, tau_l, tau_vec)
end

@doc """
  This function uses Newton's method to reduce the residual.  The Jacobian
  is calculated using one of several available methods.

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

    func must have the signature func(mesh, sbp, eqn, opts, eqn.q_vec, eqn.res_vec) 

"""->
function newton(func::Function, mesh::AbstractMesh, sbp, eqn::AbstractSolutionData, opts, pmesh=mesh; itermax=200, step_tol=1e-6, res_abstol=1e-6,  res_reltol=1e-6, res_reltol0=-1.0)
  # this function drives the non-linear residual to some specified tolerance
  # using Newton's Method
  # the jacobian is formed using finite differences
  # the initial condition is stored in eqn.q_vec
  # itermax is the maximum number of iterations
  # this function is type unstable for certain variables, but thats ok
  # the dispatch to the backslash solver and possibly the jacobian calculation
  # function will be runtime dispatched

  println("\nEntered Newtons Method")
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

  println("write_rhs = ", write_rhs)
  println("write_res = ", write_res)
  println("step_tol = ", step_tol)
  println("res_abstol = ", res_abstol)
  println("res_reltol = ", res_reltol)
  println("res_reltol0 = ", res_reltol0)

  println("before printout")
  println("typeof(eqn) = ", typeof(eqn))
  println("after printout")
  newton_data = NewtonData(mesh, sbp, eqn, opts)

  if jac_method == 1  # finite difference
    pert = epsilon
  elseif jac_method == 2  # complex step
    pert = complex(0, epsilon)
  end


  Tjac = typeof(real(eqn.res_vec[1]))  # type of jacobian, residual
  m = length(eqn.res_vec)
#  jac = SparseMatrixCSC(mesh.sparsity_bnds, Tjac)


  if jac_type == 1  # dense
    jac = zeros(Tjac, m, m)  # storage of the jacobian matrix
  elseif jac_type == 2  # sparse
    jac = SparseMatrixCSC(mesh.sparsity_bnds, Tjac)
  elseif jac_type == 3 || jac_type == 4 # petsc
    jac, jacp, x, b, ksp, ctx = createPetscData(mesh, pmesh, sbp, eqn, opts, newton_data, func)
  end

  step_fac = 1.0 # step size limiter
#  jac_recal = 0  # number of iterations since jacobian was recalculated
  Tsol = typeof(eqn.res_vec[1])
  res_0 = zeros(Tjac, m)  # function evaluated at u0
  res_0_norm = 0.0  # norm of res_0
  delta_res_vec = zeros(Tjac, m)  # newton update
  step_norm = zero(Tjac)  # norm of newton update
  step_norm_1 = zero(Tjac) # norm of previous newton update


  ##### Write iteration 0 output #####
  # open file to write convergence data to
  # append to be on the safe side
  fconv = open("convergence.dat", "a+")
  eqn.majorIterationCallback(0, mesh, sbp, eqn, opts)

  # evaluating residual at initial condition
  println("evaluating residual at initial condition")
  res_0_norm = newton_data.res_norm_i = calcResidual(mesh, sbp, eqn, opts, func)
  println("res_0_norm = ", res_0_norm)

  # extract the real components to res_0
  for i=1:m
    res_0[i] = real(eqn.res_vec[i])
  end

  println(fconv, 0, " ", res_0_norm, " ", 0)
  flush(fconv)

  # post-residual iteration 0 output
  if write_rhs
    writedlm("rhs0.dat", res_0)
  end

  if write_qic
    writedlm("q0.dat", eqn.q)
  end

  if write_res
    writedlm("res0.dat", eqn.res)
  end

  # check if initial residual satisfied absolute or relative tolerances
  if res_0_norm < res_abstol || 
    (res_reltol0 > 0 && res_0_norm/res_reltol0 < res_reltol)

    # print which criteria was statisfied
    if res_0_norm/res_reltol0 < res_reltol
      println("Initial condition satisfied res_reltol with relative residual ", res_0_norm/res_reltol0)
      println("Residual ", res_0_norm)
    else
     println("Initial condition satisfies res_tol with residual norm ", res_0_norm)
    end
    # no need to assemble q into q_vec because it never changed

    close(fconv)

    if jac_type == 3
      destroyPetsc(jac, jacp, x, b, ksp)
    end


   return nothing
 end  # end if tolerances satisfied

 println("res_reltol0 = ", res_reltol0)
 if res_reltol0 > 0  # use the supplied res_reltol0 value
   println("using supplied value for relative residual")
   res_reltol_0 = res_reltol0
 else
   println("using initial residual for relative residual")
   res_reltol_0 = res_0_norm
 end


  # do Newton's method if not converged
  print("\n")


  for i=1:itermax
    println("Newton iteration: ", i)
    println("step_fac = ", step_fac)

    # calculate jacobian using selected method
    if jac_method == 1
#      println("calculating finite difference jacobian")

      if jac_type == 1  # dense jacobian
	println("calculating dense FD jacobian")
        @time calcJacFD(mesh, sbp, eqn, opts, func, res_0, pert, jac)

      elseif jac_type == 2  # Julia sparse jacobian
	println("calculating sparse FD jacobian")
        res_copy = copy(eqn.res)  # copy unperturbed residual
 
        @time calcJacobianSparse(mesh, sbp, eqn, opts, func, res_copy, pert, jac)
      elseif jac_type == 3  # Petsc sparse jacobian
	println("calculating sparse FD jacobian")
        res_copy = copy(eqn.res)  # copy unperturbed residual
        # use artificial dissipation for the preconditioner 
	use_dissipation_orig = eqn.params.use_dissipation
	use_edgestab_orig = eqn.params.use_edgestab
        eqn.params.use_dissipation = opts["use_dissipation_prec"]
	eqn.params.use_edgestab = opts["use_edgestab_prec"]

        @time calcJacobianSparse(pmesh, sbp, eqn, opts, func, res_copy, pert, jacp)
#        addDiagonal(mesh, sbp, eqn, jacp)        
	# use normal stabilization for the real jacobian
	eqn.params.use_dissipation = use_dissipation_orig
	eqn.params.use_edgestab = use_edgestab_orig
        @time calcJacobianSparse(mesh, sbp, eqn, opts, func, res_copy, pert, jac)
      end
      println("FD jacobian calculation @time printed above")

    elseif jac_method == 2
#      println("calculating complex step jacobian")

      if jac_type == 1  # dense jacobian
	      println("calculating dense complex step jacobian")
        @time calcJacobianComplex(mesh, sbp, eqn, opts, func, pert, jac)
      elseif jac_type == 2  # Julia sparse jacobian 
	      println("calculating sparse complex step jacobian")
        res_dummy = []  # not used, so don't allocation memory
        @time calcJacobianSparse(mesh, sbp, eqn, opts, func, res_dummy, pert, jac)
      elseif jac_type == 3 # Petsc sparse jacobian
        res_dummy = []  # not used, so don't allocation memory
	use_dissipation_orig = eqn.params.use_dissipation
	use_edgestab_orig = eqn.params.use_edgestab
        eqn.params.use_dissipation = opts["use_dissipation_prec"]
	eqn.params.use_edgestab = opts["use_edgestab_prec"]

        @time calcJacobianSparse(pmesh, sbp, eqn, opts, func, res_dummy, pert, jacp)
        
#        addDiagonal(mesh, sbp, eqn, jacp)        
	# use normal stabilization for the real jacobian
	eqn.params.use_dissipation = use_dissipation_orig
	eqn.params.use_edgestab = use_edgestab_orig
        @time calcJacobianSparse(mesh, sbp, eqn, opts, func, res_dummy, pert, jac)

      elseif jac_type == 4 # Petsc jacobian-vector product
	# calculate preconditioner matrix only
	res_dummy = []
	use_dissipation_orig = eqn.params.use_dissipation
	use_edgestab_orig = eqn.params.use_edgestab
        eqn.params.use_dissipation = opts["use_dissipation_prec"]
	eqn.params.use_edgestab = opts["use_edgestab_prec"]

	if ((i % recalc_prec_freq)) == 0 || i == 1

          @time calcJacobianSparse(pmesh, sbp, eqn, opts, func, res_dummy, pert, jacp)
	end
#        addDiagonal(mesh, sbp, eqn, jacp)        
	# use normal stabilization for the real jacobian
	eqn.params.use_dissipation = use_dissipation_orig
	eqn.params.use_edgestab = use_edgestab_orig
 
      end

      if ((i % recalc_prec_freq)) == 0 || i == 1
        println("complex step jacobian calculate @time printed above")
      end
    end

    # apply globalization
    if globalize_euler

      if jac_type == 3 || jac_type == 4
	println("applying Euler globalization to jacp")
	println("tau = ", newton_data.tau_l)
        applyEuler(mesh, sbp, eqn, opts, newton_data, jacp)
      end

      if jac_type != 4
	println("applying Euler gloablization to jac")
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
      println("finished printing jacobian")
    end
    #println(jac)

    
    # calculate Jacobian condition number
    if print_cond
      cond_j = cond(jac)
      println("Condition number of jacobian = ", cond_j)
    end


    # if eigenvalues requested and we can calculate them
    if (( print_eigs || write_eigs) && (jac_type == 1 || jac_type == 2))
      eigs_i = reverse(eigvals(jac))
      if print_eigs
	println("eigenvalues =")
	for i=1:length(eigs_i)
	  println(eigs_i[i])
	end
      end

      if write_eigs
	writedlm("eigs$i.dat", eigs_i)
      end
    end

    # do the full eigenvalue decomposition
    # if requested and if julia owns the Jacobian matrix
    if  write_eigdecomp && ( jac_type == 1 || jac_type == 2)
      println("doing eigen decomposition")
      # make a dense jacobian so we can get *all* the eigenvalues and vectors
      # the algorithm for sparse matrices can get all - 2 of them, and might
      # have problems for matrices with a wide range of eigenvalues
      jac_dense = full(jac)
      D, V = eig(jac_dense)
      writedlm("eigdecomp_real$i.dat", real(D))
      writedlm("eigdecomp_imag$i.dat", imag(D))
      writedlm("eigdecomp_realvecs$i.dat", real(V))
      writedlm("eigdecomp_imagvecs$i.dat", imag(V))
    elseif write_eigdecomp # && we can't calculate it
      println(STDERR, "Warning: not performing eigen decomposition for jacobian of type $jac_type")

    end

    # negate res
    for j=1:m
      res_0[j] = -res_0[j]
    end
   
    # calculate Newton step
    if jac_type == 1 || jac_type == 2  # julia jacobian
      @time delta_res_vec[:] = jac\(res_0)  #  calculate Newton update
      fill!(jac, 0.0)
#    @time solveMUMPS!(jac, res_0, delta_res_vec)
    elseif jac_type == 3 || jac_type == 4  # petsc jacobian
      @time petscSolve(newton_data, jac, jacp, x, b, ksp, opts, res_0, delta_res_vec)
    end
    
    println("matrix solve @time printed above")
    step_norm = norm(delta_res_vec)/m
    println("step_norm = ", step_norm)

    # perform Newton update
    for i=1:m
      eqn.q_vec[i] += step_fac*delta_res_vec[i]
    end
    
    eqn.majorIterationCallback(i, mesh, sbp, eqn, opts)
 
    # write starting values for next iteration to file
    if write_sol
      writedlm("q_vec$i.dat", eqn.q_vec)
    end

    if write_q
      disassembleSolution(mesh, sbp, eqn, opts, eqn.q_vec)
      writedlm("q$i.dat", eqn.q)
    end

    # write paraview files
    if write_vis && ((i % output_freq)) == 0 || i == 1
      vals = abs(real(eqn.q_vec))  # remove unneded imaginary part
      saveSolutionToMesh(mesh, vals)
      fname = string("solution_newton", i)
      writeVisFiles(mesh, fname)
    end
 

    # calculate residual at updated location, used for next iteration rhs
    newton_data.res_norm_i_1 = newton_data.res_norm_i
    res_0_norm = newton_data.res_norm_i = calcResidual(mesh, sbp, eqn, opts, func)
    # extract real component to res_0
    for i=1:m
      res_0[i] = real(eqn.res_vec[i])
    end

    println("residual norm = ", res_0_norm)
    println("relative residual ", res_0_norm/res_reltol_0)


    # write to convergence file
    println(fconv, i, " ", res_0_norm, " ", step_norm)
    println("printed to convergence.dat")
    flush(fconv)


#    tmp = i+1
    # write rhs to file
    if write_rhs
      writedlm("rhs$i.dat", res_0)
    end

    if write_res
      writedlm("res$i.dat", eqn.res)
    end


   if res_0_norm < res_abstol || res_0_norm/res_reltol_0 < res_reltol
     if res_0_norm < res_abstol 
       println("Newton iteration converged with residual norm ", res_0_norm)
     end
     if res_0_norm/res_reltol_0 < res_reltol
      println("Newton iteration converged with relative residual norm ", res_0_norm/res_reltol_0)
    end

     # put residual into eqn.res_vec
     for i=1:m
       eqn.res_vec[i] = res_0[i]
     end

     close(fconv)

     if jac_type == 3
	destroyPetsc(jac, jacp, x, b, ksp)
      end

    
     return nothing
   end  # end if tolerances satisfied

    if (step_norm < step_tol)
      println("Newton iteration converged with step_norm = ", step_norm)
      println("Final residual = ", res_0_norm)

      # put residual into eqn.res_vec
      for i=1:m
        eqn.res_vec[i] = res_0[i]
      end
      close(fconv)
      
      if jac_type == 3
	destroyPetsc(jac, jacp, x, b, ksp)
      end

      return nothing
    end


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

    # update globalization parameters
    if globalize_euler
      updateEuler(newton_data)
    end

    if jac_type == 3 || jac_type == 4
      updateKrylov(newton_data)
    end

    print("\n")
    step_norm_1 = step_norm
    
  end  # end loop over newton iterations

  println(STDERR, "Warning: Newton iteration did not converge")

  println("Warning: Newton iteration did not converge in ", itermax, " iterations")
  println("  Final step size: ", step_norm)
  println("  Final residual: ", res_0_norm)
  println("  Final relative residual: ", res_0_norm/res_reltol_0)
  close(fconv)


   # put residual into eqn.res_vec
   for i=1:m
     eqn.res_vec[i] = res_0[i]
   end
 

  if jac_type == 3
    destroyPetsc(jac, jacp,  x, b, ksp)
  end
  return nothing
end


function checkJacVecProd(newton_data, mesh, sbp, eqn, opts, func, pert)
  
  v = ones(mesh.numDof)
  result1 = zeros(mesh.numDof)
#  writedlm("qvec_before.dat", eqn.q_vec)
  calcJacVecProd(newton_data, mesh, sbp, eqn, opts, pert, func, v, result1)
#  writedlm("qvec_prod.dat", eqn.q_vec)
#  writedlm("q_prod.dat", eqn.q)
  jac = SparseMatrixCSC(mesh.sparsity_bnds, Float64)
  res_dummy = []

  disassembleSolution(mesh, sbp, eqn,opts, eqn.q_vec)
  @time calcJacobianSparse(mesh, sbp, eqn, opts, func, res_dummy, pert, jac)

  if opts["newton_globalize_euler"]
    applyEuler(mesh, sbp, eqn, opts, newton_data, jac)
  end
#  writedlm("qvec_explicit.dat", eqn.q_vec)
#  writedlm("q_explicit.dat", eqn.q)
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


  for j=2:4
#  v2 = collect(1:mesh.numDof)
   v2 = linspace(j, j+1, mesh.numDof)
  result3 = jac*v2

  result4 = zeros(mesh.numDof)

  calcJacVecProd(newton_data, mesh, sbp, eqn, opts, pert, func, v2, result4)

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

function assembleResidual{T}(mesh, sbp, eqn, opts, res_vec::AbstractArray{T, 1}; assemble_edgeres=true, zero_resvec=true)
# assembles all of the residuals into res_vec
# no aliaising concerns

  eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.res, res_vec)

  if assemble_edgeres

    for i=1:size(eqn.res_edge, 4)
      eqn.assembleSolution(mesh, sbp, eqn, opts, view(eqn.res_edge, :, :, :, i), res_vec, zero_resvec=zero_resvec)
    end
  end

  return nothing
end


function disassembleSolution{T}(mesh, sbp, eqn, opts, q_vec::AbstractArray{T, 1})
# scatters the q_vec to the 3d array eqn.q
# no aliasing concerns here
  eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, q_vec)

  return nothing
end


@doc """
### NonlinearSolvers.calcResidual

  This function takes the eqn object with the solution varaibles stored in
    q_vec, scatters them into q, evaluates the residual, and then gathers the 
    residual values into res_vec.

    Effectively, this is a wrapper around the physics module residual evaluation
      function (which performs eqn.q -> eqn.res) that performs eqn.q_vec ->
      eqn.res_vec.

    The norm of the residual (using the SBP norm) is calculated and returned

    Inputs:
      mesh:  an AbstractMesh object
      sbp:  an SBP operator
      eqn:  an AbstractSolutionData object
      opts:  options dictonary
      func: residual evaluation function
      res_0:  a temporary vector of length numDof.

    Outputs:
      res_0_norm:  norm of residual

    Aliasing restrictions: res_0 should not be eqn.res_vec if 
"""->
function calcResidual(mesh, sbp, eqn, opts, func)
# calculate the residual and its norm

  disassembleSolution(mesh, sbp, eqn, opts, eqn.q_vec)
  func(mesh, sbp, eqn, opts)
#  res_0[:] = real(eqn.res_vec)  # is there an unnecessary copy here?

  fill!(eqn.res_vec, 0.0)
  assembleResidual(mesh, sbp, eqn, opts, eqn.res_vec, assemble_edgeres=false)

#=
  for j=1:m
    res_0[j] = real(eqn.res_vec[j])
  end

  strongres = eqn.Minv.*res_0
=#
  res_0_norm = calcNorm(eqn, eqn.res_vec)
#  println("residual norm = ", res_0_norm)

 return res_0_norm
end


function calcJacVecProd_wrapper(A::PetscMat, x::PetscVec, b::PetscVec)
# calculate Ax = b

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




function calcJacVecProd(newton_data, mesh, sbp, eqn, opts, pert, func, vec::AbstractVector, b::AbstractVector)
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
#  writedlm("q_prod_after$itr.dat", eqn.q)
  func(mesh, sbp, eqn, opts)

#  writedlm("prod_res$itr.dat", eqn.res)
  fill!(eqn.res_vec, 0.0)
  # gather into eqn.res_vec
  assembleResidual(mesh, sbp, eqn, opts, eqn.res_vec, assemble_edgeres=opts["use_edge_res"])
  
#  writedlm("res_vec$itr.dat", eqn.res_vec)
  # calculate derivatives, store into b
  calcJacRow(b, eqn.res_vec, epsilon)

  if globalize_euler
    applyEuler(mesh, sbp, eqn, opts, vec, newton_data, b)
  end



#  writedlm("pprod$itr.dat", b)
#  println("b = ", b)
  # undo perturbation
  for i=1:mesh.numDof
    eqn.q_vec[i] -= pert*vec[i]
  end

#=
  cnt = 0
  for i=1:mesh.numDof
    if abs(b2[i] - b[i]) > 1e-14
      cnt += 1
#      println("  i = ", i)
    end
  end

  if cnt != 0
    println(STDERR, "Warning: jacobian vector product check failed")
    println("cnt = ", cnt)
    result_diff = b2 - b
    diff_norm = calcNorm(eqn, result_diff)
    println("diff norm = ", diff_norm)
#    println("result_diff = ", result_diff)
  else
    println("jacobian vector product check passed")
  end
=#
  eqn.params.krylov_itr += 1

#  print("\n")
  return nothing
end

  

function calcJacFD(mesh, sbp, eqn, opts, func, res_0, pert, jac)
# calculate the jacobian using finite difference
  #println(res_0)
  (m,n) = size(jac)
  entry_orig = zero(eltype(eqn.q_vec))
  epsilon = norm(pert)  # finite difference perturbation
  # calculate jacobian
  for j=1:m
#      println("  jacobian iteration ", j)
    if j==1
      entry_orig = eqn.q_vec[j]
      #println(eqn.res_vec)
      eqn.q_vec[j] +=  epsilon
    else
      eqn.q_vec[j-1] = entry_orig # undo previous iteration pertubation
      entry_orig = eqn.q_vec[j]
      eqn.q_vec[j] += epsilon
    end

    disassembleSolution(mesh, sbp, eqn, opts, eqn.q_vec)
    # evaluate residual
    func(mesh, sbp, eqn, opts)
#     println("column ", j, " of jacobian, res_vec = ", eqn.res_vec)


    fill!(eqn.res_vec, 0.0)
    assembleResidual(mesh, sbp, eqn, opts,  eqn.res_vec)
    #println(eqn.res_vec)
    calcJacRow(unsafe_view(jac, :, j), res_0, eqn.res_vec, epsilon)
#      println("res_vec norm = ", norm(res_vec)/m)
    
  end

  # undo final perturbation
  eqn.q_vec[m] = entry_orig


  return nothing
end


function addDiagonal(mesh, sbp, eqn, jac)
# add the mass matrix to the jacobian

  for i=1:mesh.numDof
     idx = PetscInt[i-1]
     idy = PetscInt[i-1]
     vals = [100*eqn.M[i]]

#     println("adding ", vals, " to jacobian entry ", i, ",", i)
     PetscMatSetValues(jac, idx, idy, vals, PETSC_ADD_VALUES)
   end

   return nothing

 end

function calcJacobianSparse(mesh, sbp, eqn, opts, func, res_0, pert, jac::Union(SparseMatrixCSC, PetscMat))
# res_0 is 3d array of unperturbed residual, only needed for finite difference
# pert is perturbation to apply
# this function is independent of perturbation type

  filter_orig = eqn.params.use_filter  # record original filter state
  eqn.params.use_filter = false  # don't repetatively filter

#  epsilon = 1e-6  # finite difference perturbation
  epsilon = norm(pert)  # get magnitude of perturbation
#  (m,n) = size(jac)
   m = length(res_0)
#  fill!(jac.nzval, 0.0)

  # for each color, store the perturbed element corresponding to each element
#  perturbed_els = zeros(eltype(mesh.neighbor_nums), mesh.numEl)

  # debugging: do only first color
  for color=1:mesh.numColors  # loop over colors
#    println("color = ", color)
#    getPertNeighbors(mesh, color, perturbed_els)
    for j=1:mesh.numNodesPerElement  # loop over nodes 
#      println("node ", j)
      for i=1:mesh.numDofPerNode  # loop over dofs on each node
#	println("dof ", i)
        # do perturbation for each residual here:


	# apply perturbation to q
#	println("  applying perturbation")
        applyPerturbation(eqn.q, mesh.color_masks[color], pert, i, j)
#	println("wrote imag(eqn.q)")
#	println("size(eqn.q) = ", size(eqn.q))

	# evaluate residual
#	println("  evaluating residual")
        func(mesh, sbp, eqn, opts)
#	
	# assemble res into jac
#        println("  assembling jacobian")
#	println("  local dof = ", i, ", local node = ", j)
#	println("  color_mask = ", mesh.color_masks[color])
#	println("  assembling standard residual")
	for k=1:mesh.numEl  # loop over elements in residual
#	  el_pert = perturbed_els[k] # get perturbed element

	  el_pert = mesh.pertNeighborEls[k, color] # get perturbed element
          if el_pert != 0   # if element was actually perturbed for this color

#	    println("el_pert = ", el_pert)
            col_idx = mesh.dofs[i, j, el_pert]
#	    println("dof_pert = ", col_idx)
	    #TODO: make an immutable type to hold the bookeeping info
	    assembleElement(mesh, eqn, eqn.res, res_0, k, el_pert, col_idx, epsilon, jac)
	 end  # end if el_pert != 0
       end  # end loop over k

       # now do res_edge, if needed
       #TODO: consider making eqn.res 4 dimensional, which is the 
       #      generalization of having multiple residuals 
        for edge = 1:size(eqn.res_edge, 4)
#	  println("assembling edge residual = ", edge)
	  res_edge = view(eqn.res_edge, :, :, :, edge)
	  for k=1:mesh.numEl  # loop over elements in residual
  #	  el_pert = perturbed_els[k] # get perturbed element

	    el_pert = mesh.pertNeighborEls_edge[k, edge] # get perturbed element
#	    println("el_pert = ", el_pert)
	    if el_pert != 0   # if element was actually perturbed for this color

	      col_idx = mesh.dofs[i, j, el_pert]
#	      println("dof_pert = ", col_idx)
	      #TODO: make an immutable type to hold the bookeeping info
	      assembleElement(mesh, eqn, res_edge, res_0, k, el_pert, col_idx, epsilon, jac)
	   end  # end if el_pert != 0
        end  # end loop over k
      end  # end loop over local edges



      # undo perturbation
      # is this the best way to undo the perturbation?
      # faster to just take the real part of every element?

      #      println("  undoing perturbation")
      applyPerturbation(eqn.q, mesh.color_masks[color], -pert, i, j)

      end  # end loop i
    end  # end loop j
  end  # end loop over colors

  # now jac is complete
  eqn.params.use_filter = filter_orig # reset filter
  return nothing

end

function initializeTempVariables(mesh::AbstractMesh)
# declare some static work variables
# this must be called before calling newton(...)
# this could  be a problem if we do multigrid or something
# where we solve different equations or orders of accuracy 
# in a single run

local_size = mesh.numNodesPerElement*mesh.numDofPerNode
global const vals_tmp = zeros(1, local_size) # values
global const idx_tmp = zeros(PetscInt, local_size)  # row index
global const idy_tmp = zeros(PetscInt, 1)  # column indices

return nothing

end



# finite difference
function assembleElement{Tsol <: Real}(mesh, eqn::AbstractSolutionData{Tsol}, res_arr, res_0, el_res::Integer, el_pert::Integer, dof_pert::Integer, epsilon, jac::PetscMat)
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
idy_tmp[1] = dof_pert - 1

pos = 1
for j_j = 1:mesh.numNodesPerElement
  for i_i = 1:mesh.numDofPerNode
    idx_tmp[pos] = mesh.dofs[i_i, j_j, el_res] - 1
#    col_idx = mesh.dofs[i, j, el_pert]

    tmp = (res_arr[i_i,j_j, el_res] - res_0[i_i, j_j, el_res])/epsilon
    vals_tmp[pos] = tmp
    pos += 1
  end
end

#println("assembling values ", vals_tmp, " into row ", dof_pert, " columns ", idy_tmp)

#println("idy_tmp = ", idy_tmp)
PetscMatSetValues(jac, idx_tmp, idy_tmp, vals_tmp, PETSC_ADD_VALUES)




return nothing

end

# finite difference Petsc
function assembleElement{Tsol <: Real}(mesh, eqn::AbstractSolutionData{Tsol}, res_arr, res_0, el_res::Integer, el_pert::Integer, dof_pert::Integer, epsilon, jac::SparseMatrixCSC)
# assemble an element contribution into jacobian
# making this a separate function enables dispatch on type of jacobian
# el_res is the element in the residual to assemble
# el_pert is the element that was perturbed
# dof_pert is the dof number (global) of the dof that was perturbed
# typically either el_pert or dof_pert will be needed, not both

for j_j = 1:mesh.numNodesPerElement
  for i_i = 1:mesh.numDofPerNode
    row_idx = mesh.dofs[i_i, j_j, el_res]
#    col_idx = mesh.dofs[i, j, el_pert]

    tmp = (res_arr[i_i,j_j, el_res] - res_0[i_i, j_j, el_res])/epsilon
    jac[row_idx, dof_pert] += tmp

  end
end

return nothing

end




function calcJacRow{T <: Real}(jac_row, res_0, res::AbstractArray{T,1}, epsilon)
# calculate a row of the jacobian from res_0, the function evaluated 
# at the original point, and res, the function evaluated at a perturbed point

m = length(res_0)
for i=1:m
  jac_row[i] = (res[i] - res_0[i])/epsilon
end
#println(jac_row)

return nothing

end





function calcJacobianComplex(mesh, sbp, eqn, opts, func, pert, jac)

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
    func(mesh, sbp, eqn, opts)

    fill!(eqn.res_vec, 0.0)
    assembleResidual(mesh, sbp, eqn, opts, eqn.res_vec)
#     println("column ", j, " of jacobian, res_vec = ", eqn.res_vec)
    calcJacRow(unsafe_view(jac, :, j), eqn.res_vec, epsilon)
#      println("res_vec norm = ", norm(res_vec)/m)
    
  end  # end loop over rows of jacobian


  # undo final perturbation
  eqn.q_vec[m] = entry_orig
#

  return nothing
end



# for complex numbers
function assembleElement{Tsol <: Complex}(mesh, eqn::AbstractSolutionData{Tsol}, res_arr, res_0,  el_res::Integer, el_pert::Integer, dof_pert::Integer, epsilon, jac::PetscMat)
# assemble an element contribution into jacobian
# making this a separate function enables dispatch on type of jacobian
# el_res is the element in the residual to assemble
# el_pert is the element that was perturbed
# dof_pert is the dof number (global) of the dof that was perturbed
# typically either el_pert or dof_pert will be needed, not both

# resize array
# basically a no-op if array is already the right size
local_size = PetscInt(mesh.numNodesPerElement*mesh.numDofPerNode)

# get row number
idy_tmp[1] = dof_pert - 1

pos = 1
for j_j = 1:mesh.numNodesPerElement
  for i_i = 1:mesh.numDofPerNode
    idx_tmp[pos] = mesh.dofs[i_i, j_j, el_res] - 1
#    col_idx = mesh.dofs[i, j, el_pert]

    vals_tmp[pos] = imag(res_arr[i_i,j_j, el_res])/epsilon

    pos += 1
  end
end

PetscMatSetValues(jac, idx_tmp, idy_tmp, vals_tmp, PETSC_ADD_VALUES)

return nothing

end




function assembleElement{Tsol <: Complex}(mesh, eqn::AbstractSolutionData{Tsol}, res_arr, res_0,  el_res::Integer, el_pert::Integer, dof_pert::Integer, epsilon, jac::SparseMatrixCSC)
# assemble an element contribution into jacobian
# making this a separate function enables dispatch on type of jacobian
# el_res is the element in the residual to assemble
# el_pert is the element that was perturbed
# dof_pert is the dof number (global) of the dof that was perturbed
# typically either el_pert or dof_pert will be needed, not both

#println(" element $el_res res = ", view(res_arr, :, :, el_res))

for j_j = 1:mesh.numNodesPerElement
  for i_i = 1:mesh.numDofPerNode
    row_idx = mesh.dofs[i_i, j_j, el_res]
#    col_idx = mesh.dofs[i, j, el_pert]
#=
     if row_idx == 1 && dof_pert == 5
       println("i_i = ", i_i)
       println("j_j = ", j_j)
       println("el_res = ", el_res)
       println("inserting value ", imag(res_arr[i_i,j_j, el_res])/epsilon, " into jacobian (1,5)")
     end
=#

    jac[row_idx, dof_pert] += imag(res_arr[i_i,j_j, el_res])/epsilon
  end
end

return nothing

end



function applyPerturbation(arr, mask, pert, i, j)
  # applys perturbation puert to array arr according to mask mask
  # i, j specify the dof, node number within arr
  # the length of mask must equal the third dimension of arr
  # this function is independent of the type of pert

  @assert size(arr,3) == length(mask)
  @assert i <= size(arr, 1)
  @assert j <= size(arr, 2)

  (ndof, nnodes, numel) = size(arr)

  for k=1:numel
    arr[i, j, k] += pert*mask[k]
  end

  return nothing
end

 
function calcJacRow{T <: Complex}(jac_row, res::AbstractArray{T, 1}, epsilon)
# calculate a row of the jacobian from res_0, the function evaluated 
# at the original point, and res, the function evaluated at a perturbed point

m = length(res)

for i=1:m
  jac_row[i] = imag(res[i])/epsilon
end

return nothing

end



function createPetscData(mesh::AbstractMesh, pmesh::AbstractMesh, sbp, eqn::AbstractSolutionData, opts, newton_data::NewtonData, func)
# initialize petsc and create Jacobian matrix A, and vectors x, b, and the
# ksp context
# serial only
# func residual evaluation function

jac_type = opts["jac_type"]::Int

# initialize Petsc
#PetscInitialize(["-malloc", "-malloc_debug", "-malloc_dump", "-ksp_monitor", "-pc_type", "ilu", "-pc_factor_levels", "4" ])

#PetscInitialize(["-malloc", "-malloc_debug", "-malloc_dump", "-sub_pc_factor_levels", "4", "ksp_gmres_modifiedgramschmidt", "-ksp_pc_side", "right", "-ksp_gmres_restart", "1000" ])
numDofPerNode = mesh.numDofPerNode
PetscInitialize(["-malloc", "-malloc_debug", "-malloc_dump", "-ksp_monitor", "-pc_type", "ilu", "-pc_factor_levels", "6", "ksp_gmres_modifiedgramschmidt", "-ksp_pc_side", "right", "-ksp_gmres_restart", "30" ])
comm = MPI.COMM_WORLD

obj_size = PetscInt(mesh.numDof)  # length of vectors, side length of matrices
println("creating b")
b = PetscVec(comm)
PetscVecSetType(b, VECSEQ)
PetscVecSetSizes(b, obj_size, obj_size)

println("creating x")
x = PetscVec(comm)
PetscVecSetType(x, VECSEQ)
PetscVecSetSizes(x, obj_size, obj_size)

# tuple of all the objects needed to evaluate the residual
# only used by Mat Shell
ctx = (mesh, sbp, eqn, opts, newton_data, func)  # tuple of all the objects needed to evalute the
                                    # residual
if jac_type == 3  # explicit sparse jacobian
  println("creating A")
  A = PetscMat(comm)
  PetscMatSetFromOptions(A)
  PetscMatSetType(A, PETSc.MATSEQBAIJ)
  PetscMatSetSizes(A, obj_size, obj_size, obj_size, obj_size)

elseif jac_type == 4  # jacobian-vector product
  # create matrix shell
  ctx_ptr = pointer_from_objref(ctx)  # make a pointer from the tuple
  A = MatCreateShell(comm, obj_size, obj_size, obj_size, obj_size, ctx_ptr)
  PetscMatSetFromOptions(A)  # necessary?
  PetscSetUp(A)

  # give PETSc the function pointer of Jacobian vector product function
  fptr = cfunction (calcJacVecProd_wrapper, PetscErrorCode, (PetscMat, PetscVec, PetscVec))
  MatShellSetOperation(A, PETSc.MATOP_MULT, fptr)

else
  println(STDERR, "Unsupported jacobian type requested")
  println(STDERR, "jac_type = ", jac_type)
end

println("type of A = ", MatGetType(A))

println("creating Ap")  # used for preconditioner
Ap = PetscMat(comm)
PetscMatSetFromOptions(Ap)
PetscMatSetType(Ap, PETSc.MATSEQBAIJ)
PetscMatSetSizes(Ap, PetscInt(mesh.numDof), PetscInt(mesh.numDof), PetscInt(mesh.numDof), PetscInt(mesh.numDof))

println("type of Ap = ", MatGetType(Ap))


println("preallocating Petsc Matrices")
# prellocate matrix
dnnz = zeros(PetscInt, mesh.numDof)  # diagonal non zeros per row
onnz = zeros(PetscInt, mesh.numDof)  # there is no off diagonal part for single proc case
dnnzu = zeros(PetscInt, 1)  # only needed for symmetric matrices
onnzu = zeros(PetscInt, 1)  # only needed for symmetric matrices
bs = PetscInt(mesh.numDofPerNode)  # block size

# calculate number of non zeros per row for A
for i=1:mesh.numNodes
  max_dof = mesh.sparsity_nodebnds[2, i]
  min_dof = mesh.sparsity_nodebnds[1, i]
  nnz_i = max_dof - min_dof + 1
  dnnz[i] = nnz_i
#  println("row ", i," has ", nnz_i, " non zero entries")
end

if jac_type == 3

  PetscMatXAIJSetPreallocation(A, bs, dnnz, onnz, dnnzu, onnzu)

  MatSetOption(A, PETSc.MAT_ROW_ORIENTED, PETSC_FALSE)
  PetscMatZeroEntries(A)
  matinfo = PetscMatGetInfo(A, Int32(1))
  println("A block size = ", matinfo.block_size)

  PetscMatAssemblyBegin(A, PETSC_MAT_FLUSH_ASSEMBLY)
  PetscMatAssemblyEnd(A, PETSC_MAT_FLUSH_ASSEMBLY)
end

# calculate number of nonzeros per row for A[
for i=1:mesh.numNodes
  max_dof = pmesh.sparsity_nodebnds[2, i]
  min_dof = pmesh.sparsity_nodebnds[1, i]
  nnz_i = max_dof - min_dof + 1
  dnnz[i] = nnz_i
#  println("row ", i," has ", nnz_i, " non zero entries")
end

PetscMatXAIJSetPreallocation(Ap, bs, dnnz, onnz, dnnzu, onnzu)

MatSetOption(Ap, PETSc.MAT_ROW_ORIENTED, PETSC_FALSE)
matinfo = PetscMatGetInfo(Ap, Int32(1))
println("Ap block size = ", matinfo.block_size)

# zero initialize the matrix just in case
println("zeroing Ap")
PetscMatZeroEntries(Ap)

PetscMatAssemblyBegin(Ap, PETSC_MAT_FLUSH_ASSEMBLY)
PetscMatAssemblyEnd(Ap, PETSC_MAT_FLUSH_ASSEMBLY)



# set some options
# MatSetValuesBlocked will interpret the array of values as being column
# major

# create KSP contex
ksp = KSP(comm)

KSPSetFromOptions(ksp)
KSPSetOperators(ksp, A, Ap)  # this was A, Ap

# set: rtol, abstol, dtol, maxits
#KSPSetTolerances(ksp, 1e-2, 1e-12, 1e5, PetscInt(1000))
#KSPSetUp(ksp)


pc = KSPGetPC(ksp)
pc_type = PCGetType(pc)
println("pc_type = ", pc_type)
#=
if pc_type == "bjacobi"
  n_local, first_local, ksp_arr = PCBJacobiGetSubKSP(pc)
  println("n_local = ", n_local, ", first_local = ", first_local)
  println("length(ksp_arr) = ", length(ksp_arr))

  sub_ksp = ksp_arr[first_local + 1]
  sub_pc = KSPGetPC(sub_ksp)
  pc_subtype = PCGetType(sub_pc)
  println("pc_subtype = ", pc_subtype)

  fill_level = PCFactorGetLevels(sub_pc)
  println("preconditioner using fill level = ", fill_level)
end

if pc_type == "ilu"
  fill_level = PCFactorGetLevels(sub_pc)
  println("preconditioner using fill level = ", fill_level)
end
=#


return A, Ap, x, b, ksp, ctx

end

function destroyPetsc(A::PetscMat, Ap::PetscMat, x::PetscVec,  b::PetscVec, ksp::KSP)
# destory Petsc data structures, finalize Petsc

PetscDestroy(A)
PetscDestroy(Ap)
PetscDestroy(x)
PetscDestroy(b)
PetscDestroy(ksp)
PetscFinalize()

return nothing

end


function petscSolve(newton_data::NewtonData, A::PetscMat, Ap::PetscMat, x::PetscVec, b::PetscVec, ksp::KSP, opts, res_0::AbstractVector, delta_res_vec::AbstractVector )

  # solve the system for the newton step, write it to delta_res_vec
  # writing it to delta_res_vec is an unecessary copy, becasue we could
  # write it directly to eqn.q, but for consistency we do it anyways

  jac_type = opts["jac_type"]::Int

  # copy res_0 into b
  # create the index array
  numDof = length(delta_res_vec)
  println("copying res_0 to b")
  idx = zeros(PetscInt, numDof)
  for i=1:numDof
    idx[i] = i - 1
  end

  # copy into Petsc and assemble
  # should do this by direct array access
  PetscVecSetValues(b, idx, res_0, PETSC_INSERT_VALUES)
  PetscVecAssemblyBegin(b)
  PetscVecAssemblyEnd(b)

  #=
  PetscVecSetValues(x, idx, res_0, PETSC_INSERT_VALUES)
  PetscVecAssemblyBegin(x)
  PetscVecAssemblyEnd(x)
  =#


  # assemble matrix
  # should this be overlapped with the vector assembly?
  if jac_type != 4
    println("assembling A")
    PetscMatAssemblyBegin(A, PETSC_MAT_FINAL_ASSEMBLY)
    PetscMatAssemblyEnd(A, PETSC_MAT_FINAL_ASSEMBLY)
  end

  PetscMatAssemblyBegin(Ap, PETSC_MAT_FINAL_ASSEMBLY)
  PetscMatAssemblyEnd(Ap, PETSC_MAT_FINAL_ASSEMBLY)

  if jac_type != 4
    matinfo = PetscMatGetInfo(A, MAT_LOCAL)
    println("number of mallocs for A = ", matinfo.mallocs)
    if matinfo.mallocs > 0.5  # if any mallocs
      println("Caution: non zero number of mallocs for A")
      println("  number of mallocs = ", matinfo.mallocs)
    end
  end

    matinfo = PetscMatGetInfo(Ap, MAT_LOCAL)

    if matinfo.mallocs > 0.5  # if any mallocs
      println("Caution: non zero number of mallocs for Ap")
      println("  number of mallocs = ", matinfo.mallocs)
    end
 

  # only call this first time?
  # what happens when A and Ap change?
  # this is not truely necessary for the common case, because KSPSolve
  # calls it if needed
  # it is necessary to call KSPSetUp before getting the preconditioner
  # context in some cases

  KSPSetTolerances(ksp, newton_data.reltol, newton_data.abstol, 
                   newton_data.dtol, PetscInt(newton_data.itermax))

  KSPSetUp(ksp)

  nx = PetscVecGetSize(x)
  nb = PetscVecGetSize(b)
  println("size of x = ", nx)
  println("size of y = ", nb)

  println("solving system")
  KSPSolve(ksp, b, x)
  reason = KSPGetConvergedReason(ksp)
  println("KSP converged reason = ", KSPConvergedReasonDict[reason])
  rnorm = KSPGetResidualNorm(ksp)
  println("Linear residual = ", rnorm)

  # copy solution from x into delta_res_vec
  arr, ptr_arr = PetscVecGetArray(x)
  for i=1:numDof
    delta_res_vec[i] = arr[i]
  end

  PetscVecRestoreArray(x, ptr_arr)

  # zero out the Jacobian for next use
  if opts["jac_type"] != 4
    PetscMatZeroEntries(A)
  end
  PetscMatZeroEntries(Ap)
  return nothing
end

@doc """
### newton_check

  Uses complex step to compare jacobian vector product to directional derivative.

"""->
function newton_check(func, mesh, sbp, eqn, opts)
  # this function drives the non-linear residual to some specified tolerance
  # using Newton's Method
  # the jacobian is formed using finite differences
  # the initial condition is stored in eqn.q_vec
  # itermax is the maximum number of iterations

  step_fac = 0.5  # step size limiter
  m = length(eqn.res_vec)
  Tsol = typeof(eqn.res_vec[1])
  Tjac = typeof(real(eqn.res_vec[1]))  # type of jacobian, residual
  jac = zeros(Tjac, m, m)  # storage of the jacobian matrix
  direction_der = zeros(mesh.numDof)
#  v = rand(mesh.numDof)
   v = readdlm("randvec.txt")

  epsilon = 1e-20  # complex step perturbation
  fill!(eqn.res_vec, 0.0)  # zero out res_vec
  # compute directional derivative
  for i=1:mesh.numDof
    eqn.q_vec[i] += complex(0, epsilon*v[i])  # apply perturbation
  end

  func(mesh, sbp, eqn, opts)
  println("evaluated directional derivative")

  # calculate derivative
  for i=1:mesh.numDof
    direction_der[i] = imag(eqn.res_vec[i])/epsilon
    eqn.q_vec[i] -= complex(0, epsilon*v[i])  # undo perturbation
  end



    println("Calculating Jacobian")

    # calculate jacobian
    for j=1:m
      println("\ncalculating column ", j, " of the jacobian")
      if j==1
	eqn.q_vec[j] +=  complex(0, epsilon)
      else
	eqn.q_vec[j-1] -= complex(0, epsilon) # undo previous iteration pertubation
	eqn.q_vec[j] += complex(0, epsilon)
      end

      # evaluate residual
      fill!(eqn.res_vec, 0.0)
      func(mesh, sbp, eqn, opts, eqn.q_vec, eqn.res_vec)
 #     println("column ", j, " of jacobian, res_vec = ", eqn.res_vec)
      calcJacRow(unsafe_view(jac, :, j), eqn.res_vec, epsilon)
#      println("res_vec norm = ", norm(res_vec)/m)
      
    end  # end loop over rows of jacobian

    # undo final perturbation
    eqn.q_vec[m] -= complex(0, epsilon)

    # now jac is complete

    fname = string("jacobian", ".dat")
    printMatrix(fname, jac)
    println("finished printing jacobian")

    cond_j = cond(jac)
    println("Condition number of jacobian = ", cond_j)
    svals = svdvals(jac)
    println("svdvals = \n", svals)

    jac_mult = jac*v

    # copy difference between directional derivative and
    # jacobian multiplication into res_vec for return

    for i=1:mesh.numDof
      eqn.res_vec[i] = direction_der[i] - jac_mult[i]
    end

    err_norm = norm(eqn.res_vec)/mesh.numDof
    println("step_norm = ", err_norm)
#    println("jac = ", jac)

    print("\n")

    println("finished newton_check")
  return nothing
end



@doc """
### newton_check

  This method calculates a single column of the jacobian with the complex step method.
"""->
function newton_check(func, mesh, sbp, eqn, opts, j)
# calculate a single column of hte jacobian
    
      jac_col = zeros(Float64, mesh.numDof)
      println("\ncalculating column ", j, " of the jacobian")

      epsilon = 1e-20

#      eqn.q_vec[j] += complex(0, epsilon)
      disassembleSolution(mesh, sbp, eqn, opts, eqn.q_vec)
      eqn.q[1, 2, 5] += complex(0, epsilon)
      writedlm("check_q.dat", imag(eqn.q))
#      eqn.q[1,1,1] += complex(0, epsilon)
      # evaluate residual
      func(mesh, sbp, eqn, opts)

      fill!(eqn.res_vec, 0.0)
      assembleResidual(mesh, sbp, eqn, opts, eqn.res_vec)
 #     println("column ", j, " of jacobian, res_vec = ", eqn.res_vec)
      calcJacRow(jac_col, eqn.res_vec, epsilon)
#      println("res_vec norm = ", norm(res_vec)/m)
      writedlm("check_res.dat", imag(eqn.res))

      return jac_col
end 

@doc """
### newton_check_fd

  This method calcualtes a single column of the jacobian with the finite difference method.

"""->
function newton_check_fd(func, mesh, sbp, eqn, opts, j)
# calculate a single column of hte jacobian
    
      jac_col = zeros(Float64, mesh.numDof)
      println("\ncalculating column ", j, " of the jacobian")

     disassembleSolution(mesh, sbp, eqn, opts, eqn.q_vec)
     func(mesh, sbp, eqn, opts, eqn.q_vec, eqn.res_vec)
     fill!(eqn.res_vec, 0.0)
     assembleResidual(mesh, sbp, eqn, opts, eqn.res_vec)
     res_0 = copy(eqn.res_vec)

      epsilon = 1e-6

      eqn.q_vec[j] += epsilon

      eqn.disassmbleSolution(mesh, sbp, eqn, opts, eqn.q_vec)

      # evaluate residual
      func(mesh, sbp, eqn, opts, eqn.q_vec, eqn.res_vec)

      fill!(eqn.res_vec, 0.0)
      assembleResidual(mesh, sbp, eqn, opts, eqn.res_vec)
 #     println("column ", j, " of jacobian, res_vec = ", eqn.res_vec)

      calcJacRow(jac_col, res_0, eqn.res_vec, epsilon)
#      println("res_vec norm = ", norm(res_vec)/m)

      return jac_col
end 
 
