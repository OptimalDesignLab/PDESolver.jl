export newton, newton_check, newton_check_fd, initializeTempVariables, calcResidual
@doc """
  This function uses Newton's method to reduce the residual.  The Jacobian
  is calculated using one of several available methods.

  Arguments:
    * func  : function that evalutes the residual
    * mesh : mesh to use in evaluating the residual
    * sbp : sbp operator to be used to evaluate the residual
    * eqn : EulerData to use to evaluate the residual
    * opts : options dictonary

    Optional Arguments
    * itermax : maximum number of Newton iterations
    * step_tol : step size stopping tolerance
    * res_tol : residual stopping tolerance

    func must have the signature func(mesh, sbp, eqn, opts, eqn.q_vec, eqn.res_vec) 

"""->
function newton(func, mesh, sbp, eqn, opts; itermax=200, step_tol=1e-6, res_abstol=1e-6,  res_reltol=1e-6, res_reltol0=-1.0)
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
  write_sol = opts["write_sol"]::Bool
  write_vis = opts["write_vis"]::Bool
  write_qic = opts["write_qic"]::Bool
  write_res = opts["write_res"]::Bool
  jac_method = opts["jac_method"]::Int  # finite difference or complex step
  jac_type = opts["jac_type"]::Int  # jacobian sparse or dense
  epsilon = opts["epsilon"]::Float64

  println("write_rhs = ", write_rhs)
  println("write_res = ", write_res)
  println("step_tol = ", step_tol)
  println("res_abstol = ", res_abstol)
  println("res_reltol = ", res_reltol)
  println("res_reltol0 = ", res_reltol0)

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
  elseif jac_type == 3  # petsc
    jac, jacp, x, b, ksp = createPetscData(mesh, eqn, opts)
  end



  step_fac = 1.0 # step size limiter
#  jac_recal = 0  # number of iterations since jacobian was recalculated
  Tsol = typeof(eqn.res_vec[1])
#  jac = zeros(Tjac, m, m)  # storage of the jacobian matrix
  res_0 = zeros(Tjac, m)  # function evaluated at u0
  res_0_norm = 0.0  # norm of res_0
  delta_res_vec = zeros(Tjac, m)  # newton update
  step_norm = zero(Tjac)  # norm of newton update
  step_norm_1 = zero(Tjac) # norm of previous newton update


  # open file to write convergence data to
  # append to be on the safe side
  fconv = open("convergence.dat", "a+")

  eqn.majorIterationCallback(0, mesh, sbp, eqn, opts)

  # evaluating residual at initial condition
  println("evaluating residual at initial condition")
  res_0_norm = calcResidual(mesh, sbp, eqn, opts, func, res_0)

  println(fconv, 0, " ", res_0_norm, " ", 0)
  flush(fconv)



  # write rhs to file
  if write_rhs
    writedlm("rhs1.dat", res_0)
  end

  if write_qic
    writedlm("qic.dat", eqn.q)
  end

  if write_res
    writedlm("res1.dat", eqn.res)
  end

  # check if initial residual satisfied absolute or relative tolerances
  if res_0_norm < res_abstol || (res_reltol0 > 0 && res_0_norm/res_reltol0 < res_reltol)

    if res_0_norm/res_reltol0 < res_reltol
      println("Initial condition satisfied res_reltol with relative residual ", res_0_norm/res_reltol0)
      println("Residual ", res_0_norm)
    else
     println("Initial condition satisfies res_tol with residual norm ", res_0_norm)
   end
#   println("writing to convergence.dat")
#   println(fconv, i, " ", res_0_norm, " ", 0.0)
   # no need to assemble q into q_vec because it never changed

   close(fconv)

    if jac_type == 3
      destroyPetsc(jac, jacp, x, b, ksp)
    end


   return nothing
 end

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
        eqn.params.use_dissipation = true
	eqn.params.use_edgestab = false

        @time calcJacobianSparse(mesh, sbp, eqn, opts, func, res_copy, pert, jacp)
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
      elseif jac_type == 3
        res_dummy = []  # not used, so don't allocation memory
	use_dissipation_orig = eqn.params.use_dissipation
	use_edgestab_orig = eqn.params.use_edgestab
        eqn.params.use_dissipation = true
	eqn.params.use_edgestab = false

        @time calcJacobianSparse(mesh, sbp, eqn, opts, func, res_dummy, pert, jacp)
        
#        addDiagonal(mesh, sbp, eqn, jacp)        
	# use normal stabilization for the real jacobian
	eqn.params.use_dissipation = use_dissipation_orig
	eqn.params.use_edgestab = use_edgestab_orig
        @time calcJacobianSparse(mesh, sbp, eqn, opts, func, res_dummy, pert, jac)
 
      end

      println("complex step jacobian calculate @time printed above")
    end

    # print as determined by options
    if write_jac
       # fname = string("jacobian", i, ".dat")
       # printMatrix(fname, jac)
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

    # negate res
    for j=1:m
      res_0[j] = -res_0[j]
    end
   
    # calculate Newton step
    if jac_type == 1 || jac_type == 2  # julia jacobian
      @time delta_res_vec[:] = jac\(res_0)  #  calculate Newton update
      fill!(jac, 0.0)
#    @time solveMUMPS!(jac, res_0, delta_res_vec)
    elseif jac_type == 3   # petsc
      @time petscSolve(jac, jacp, x, b, ksp, res_0, delta_res_vec)
    end
    
    println("matrix solve @time printed above")
    step_norm = norm(delta_res_vec)/m
    println("step_norm = ", step_norm)

    # perform Newton update
    eqn.q_vec[:] += step_fac*delta_res_vec  # update q_vec
    
    eqn.majorIterationCallback(i, mesh, sbp, eqn, opts)
 
    # write starting values for next iteration to file
    if write_sol
      writedlm("q_vec$i.dat", eqn.q_vec)
    end

    # write paraview files
    if write_vis
      vals = abs(real(eqn.q_vec))  # remove unneded imaginary part
      saveSolutionToMesh(mesh, vals)
      fname = string("solution_newton", i)
      writeVisFiles(mesh, fname)
    end
 

    # calculate residual at updated location, used for next iteration rhs
    res_0_norm = calcResidual(mesh, sbp, eqn, opts, func, res_0)
    println("relative residual ", res_0_norm/res_reltol_0)


    # write to convergence file
    println("i = ", i)
    println(fconv, i, " ", res_0_norm, " ", step_norm)
    println("printed to convergence.dat")
    flush(fconv)


    tmp = i+1
    # write rhs to file
    if write_rhs
      writedlm("rhs$tmp.dat", res_0)
    end

    if write_res
      writedlm("res$tmp.dat", eqn.res)
    end


   if res_0_norm < res_abstol || res_0_norm/res_reltol_0 < res_reltol
     if res_0_norm < res_abstol 
       println("Newton iteration converged with residual norm ", res_0_norm)
     end
     if res_0_norm/res_reltol_0 < res_reltol
      println("Newton iteration converged with relative residual norm ", res_0_norm/res_reltol_0)
    end

     # put solution into q_vec
#     fill!(eqn.q_vec, 0.0)
#     eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)

     # put residual into eqn.res_vec
     eqn.res_vec[:] = res_0
     close(fconv)

     if jac_type == 3
	destroyPetsc(jac, jacp, x, b, ksp)
      end

    
     return nothing
   end

    if (step_norm < step_tol)
      println("Newton iteration converged with step_norm = ", step_norm)
      println("Final residual = ", res_0_norm)

     # put solution into q_vec
#     fill!(eqn.q_vec, 0.0)
#     eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)

      # put residual into eqn.res_vec
      eqn.res_vec[:] = res_0
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

    print("\n")
    step_norm_1 = step_norm
    
  end  # end loop over newton iterations

  println(STDERR, "Warning: Newton iteration did not converge")

  println("Warning: Newton iteration did not converge in ", itermax, " iterations")
  println("  Final step size: ", step_norm)
  println("  Final residual: ", res_0_norm)
  println("  Final relative residual: ", res_0_norm/res_reltol_0)
  close(fconv)


   # put solution into q_vec
#   fill!(eqn.q_vec, 0.0)
#   eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)

   # put residual into eqn.res_vec
   eqn.res_vec[:] = res_0
 

  if jac_type == 3
    destroyPetsc(jac, jacp,  x, b, ksp)
  end
  return nothing
end


function calcResidual(mesh, sbp, eqn, opts, func, res_0)
# calculate the residual and its norm

  m = length(res_0)

  eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q_vec)
  func(mesh, sbp, eqn, opts)
#  res_0[:] = real(eqn.res_vec)  # is there an unnecessary copy here?

  fill!(eqn.res_vec, 0.0)
  eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)

  for j=1:m
    res_0[j] = real(eqn.res_vec[j])
  end

  strongres = eqn.Minv.*res_0
  res_0_norm = calcNorm(eqn, strongres)
  println("residual norm = ", res_0_norm)

 return res_0_norm
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

    eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q_vec)
    # evaluate residual
    func(mesh, sbp, eqn, opts)
#     println("column ", j, " of jacobian, res_vec = ", eqn.res_vec)


    fill!(eqn.res_vec, 0.0)
    eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.res,  eqn.res_vec)
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
	for k=1:mesh.numEl  # loop over elements in residual
#	  el_pert = perturbed_els[k] # get perturbed element

	  el_pert = mesh.pertNeighborEls[k, color] # get perturbed element
          if el_pert != 0   # if element was actually perturbed for this color

            col_idx = mesh.dofs[i, j, el_pert]
	    #TODO: make an immutable type to hold the bookeeping info
	    assembleElement(mesh, eqn, res_0, k, el_pert, col_idx, epsilon, jac)
	 end  # end if el_pert != 0
       end  # end loop over k

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
function assembleElement{Tsol <: Real}(mesh, eqn::AbstractSolutionData{Tsol}, res_0, el_res::Integer, el_pert::Integer, dof_pert::Integer, epsilon, jac::PetscMat)
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

    tmp = (eqn.res[i_i,j_j, el_res] - res_0[i_i, j_j, el_res])/epsilon
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
function assembleElement{Tsol <: Real}(mesh, eqn::AbstractSolutionData{Tsol}, res_0, el_res::Integer, el_pert::Integer, dof_pert::Integer, epsilon, jac::SparseMatrixCSC)
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

    tmp = (eqn.res[i_i,j_j, el_res] - res_0[i_i, j_j, el_res])/epsilon
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

    eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q_vec)
    # evaluate residual
    func(mesh, sbp, eqn, opts)

    fill!(eqn.res_vec, 0.0)
    eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
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
function assembleElement{Tsol <: Complex}(mesh, eqn::AbstractSolutionData{Tsol}, res_0,  el_res::Integer, el_pert::Integer, dof_pert::Integer, epsilon, jac::PetscMat)
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

    vals_tmp[pos] = imag(eqn.res[i_i,j_j, el_res])/epsilon

    pos += 1
  end
end

PetscMatSetValues(jac, idx_tmp, idy_tmp, vals_tmp, PETSC_ADD_VALUES)

return nothing

end




function assembleElement{Tsol <: Complex}(mesh, eqn::AbstractSolutionData{Tsol}, res_0,  el_res::Integer, el_pert::Integer, dof_pert::Integer, epsilon, jac::SparseMatrixCSC)
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

    jac[row_idx, dof_pert] += imag(eqn.res[i_i,j_j, el_res])/epsilon
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



function createPetscData(mesh::AbstractMesh, eqn::AbstractSolutionData, opts)
# initialize petsc and create Jacobian matrix A, and vectors x, b, and the
# ksp context
# serial only

# initialize Petsc
#PetscInitialize(["-malloc", "-malloc_debug", "-malloc_dump", "-ksp_monitor", "-pc_type", "ilu", "-pc_factor_levels", "4" ])

numDofPerNode = mesh.numDofPerNode
PetscInitialize(["-malloc", "-malloc_debug", "-malloc_dump", "-sub_pc_factor_levels", "4", "ksp_gmres_modifiedgramschmidt", "-ksp_pc_side", "right", "-ksp_gmres_restart", "1000" ])
comm = MPI.COMM_WORLD

println("creating b")
b = PetscVec(comm)
PetscVecSetType(b, VECMPI)
PetscVecSetSizes(b, PetscInt(mesh.numDof), PetscInt(mesh.numDof))

println("creating x")
x = PetscVec(comm)
PetscVecSetType(x, VECMPI)
PetscVecSetSizes(x, PetscInt(mesh.numDof), PetscInt(mesh.numDof))

println("creating A")
A = PetscMat(comm)
PetscMatSetFromOptions(A)
PetscMatSetType(A, MATMPIBAIJ)
PetscMatSetSizes(A, PetscInt(mesh.numDof), PetscInt(mesh.numDof), PetscInt(mesh.numDof), PetscInt(mesh.numDof))

println("creating Ap")  # used for preconditioner
Ap = PetscMat(comm)
PetscMatSetFromOptions(Ap)
PetscMatSetType(Ap, MATMPIBAIJ)
PetscMatSetSizes(Ap, PetscInt(mesh.numDof), PetscInt(mesh.numDof), PetscInt(mesh.numDof), PetscInt(mesh.numDof))



println("preallocating A")
# prellocate matrix
dnnz = zeros(PetscInt, mesh.numDof)  # diagonal non zeros per row
onnz = zeros(PetscInt, mesh.numDof)  # there is no off diagonal part for single proc case
dnnzu = zeros(PetscInt, 1)  # only needed for symmetric matrices
onnzu = zeros(PetscInt, 1)  # only needed for symmetric matrices
bs = PetscInt(mesh.numDofPerNode)  # block size

# calculate number of non zeros per row
for i=1:mesh.numNodes
  max_dof = mesh.sparsity_nodebnds[2, i]
  min_dof = mesh.sparsity_nodebnds[1, i]
  nnz_i = max_dof - min_dof + 1
  dnnz[i] = nnz_i
#  println("row ", i," has ", nnz_i, " non zero entries")
end

PetscMatXAIJSetPreallocation(A, bs, dnnz, onnz, dnnzu, onnzu)
PetscMatXAIJSetPreallocation(Ap, bs, dnnz, onnz, dnnzu, onnzu)
# zero initialize the matrix just in case
println("zeroing A")
PetscMatZeroEntries(A)
PetscMatZeroEntries(Ap)


# set some options
# MatSetValuesBlocked will interpret the array of values as being column
# major
MatSetOption(A, PETSc.MAT_ROW_ORIENTED, PETSC_FALSE)
MatSetOption(Ap, PETSc.MAT_ROW_ORIENTED, PETSC_FALSE)

matinfo = PetscMatGetInfo(A, Int32(1))
println("A block size = ", matinfo.block_size)

# create KSP contex
ksp = KSP(comm)
KSPSetOperators(ksp, A, Ap)  # this was A, Ap
KSPSetFromOptions(ksp)

# set: rtol, abstol, dtol, maxits
KSPSetTolerances(ksp, 1e-2, 1e-12, 1e5, PetscInt(1000))
KSPSetUp(ksp)


pc = KSPGetPC(ksp)
pc_type = PCGetType(pc)
println("pc_type = ", pc_type)

if pc_type != "none"
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



return A, Ap, x, b, ksp

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


function petscSolve(A::PetscMat, Ap::PetscMat, x::PetscVec, b::PetscVec, ksp::KSP, res_0::AbstractVector, delta_res_vec::AbstractVector)

  # solve the system for the newton step, write it to delta_res_vec
  # writing it to delta_res_vec is an unecessary copy, becasue we could
  # write it directly to eqn.q, but for consistency we do it anyways

  # copy res_0 into b
  # create the index array
  numDof = length(delta_res_vec)
  println("copying res_0 to b")
  idx = zeros(PetscInt, numDof)
  for i=1:numDof
    idx[i] = i - 1
  end

  # copy into Petsc and assemble
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
  println("assembling A")
  PetscMatAssemblyBegin(A, PETSC_MAT_FINAL_ASSEMBLY)
  PetscMatAssemblyEnd(A, PETSC_MAT_FINAL_ASSEMBLY)

  PetscMatAssemblyBegin(Ap, PETSC_MAT_FINAL_ASSEMBLY)
  PetscMatAssemblyEnd(Ap, PETSC_MAT_FINAL_ASSEMBLY)


  matinfo = PetscMatGetInfo(A, MAT_LOCAL)
  println("number of mallocs for A = ", matinfo.mallocs)
  if matinfo.mallocs > 0.5  # if any mallocs
    println("Caution: non zero number of mallocs for A")
    println("  number of mallocs = ", matinfo.mallocs)
  end
  matinfo = PetscMatGetInfo(Ap, MAT_LOCAL)

  if matinfo.mallocs > 0.5  # if any mallocs
    println("Caution: non zero number of mallocs for Ap")
    println("  number of mallocs = ", matinfo.mallocs)
  end


 


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
  PetscMatZeroEntries(A)
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
      eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q_vec)
      eqn.q[1, 2, 5] += complex(0, epsilon)
      writedlm("check_q.dat", imag(eqn.q))
#      eqn.q[1,1,1] += complex(0, epsilon)
      # evaluate residual
      func(mesh, sbp, eqn, opts)

      fill!(eqn.res_vec, 0.0)
      eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
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

     eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q_vec)
     func(mesh, sbp, eqn, opts, eqn.q_vec, eqn.res_vec)
     fill!(eqn.res_vec, 0.0)
     eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
     res_0 = copy(eqn.res_vec)

      epsilon = 1e-6

      eqn.q_vec[j] += epsilon

      eqn.disassmbleSolution(mesh, sbp, eqn, opts, eqn.q_vec)

      # evaluate residual
      func(mesh, sbp, eqn, opts, eqn.q_vec, eqn.res_vec)

      fill!(eqn.res_vec, 0.0)
      eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
 #     println("column ", j, " of jacobian, res_vec = ", eqn.res_vec)

      calcJacRow(jac_col, res_0, eqn.res_vec, epsilon)
#      println("res_vec norm = ", norm(res_vec)/m)

      return jac_col
end 
 
