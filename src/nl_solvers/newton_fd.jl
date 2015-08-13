export newton, newton_check, newton_check_fd, initializeTempVariables
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

    func must have the signature func(mesh, sbp, eqn, opts, eqn.SL0, eqn.SL) 

"""->
function newton(func, mesh, sbp, eqn, opts; itermax=200, step_tol=1e-6, res_tol=1e-6)
  # this function drives the non-linear residual to some specified tolerance
  # using Newton's Method
  # the jacobian is formed using finite differences
  # the initial condition is stored in eqn.SL0
  # itermax is the maximum number of iterations
  # this function is type unstable for certain variables, but thats ok
  # the dispatch to the backslash solver and possibly the jacobian calculation
  # function will be runtime dispatched

  # options
  write_rhs = opts["write_rhs"]::Bool
  write_jac = opts["write_jac"]::Bool
  print_cond = opts["print_cond"]::Bool
  write_sol = opts["write_sol"]::Bool
  write_vis = opts["write_vis"]::Bool
  jac_method = opts["jac_method"]::Int  # finite difference or complex step
  jac_type = opts["jac_type"]::Int  # jacobian sparse or dense
  epsilon = opts["epsilon"]::Float64

  if jac_method == 1  # finite difference
    pert = epsilon
  elseif jac_method == 2  # complex step
    pert = complex(0, epsilon)
  end


  Tjac = typeof(real(eqn.SL[1]))  # type of jacobian, residual
  m = length(eqn.SL)
#  jac = SparseMatrixCSC(mesh.sparsity_bnds, Tjac)


  if jac_type == 1  # dense
    jac = zeros(Tjac, m, m)  # storage of the jacobian matrix
  elseif jac_type == 2  # sparse
    jac = SparseMatrixCSC(mesh.sparsity_bnds, Tjac)
  elseif jac_type == 3  # petsc
    jac, x, b, ksp = createPetscData(mesh, eqn, opts)
  end



  step_fac = 1.0 # step size limiter
#  jac_recal = 0  # number of iterations since jacobian was recalculated
  Tsol = typeof(eqn.SL[1])
#  jac = zeros(Tjac, m, m)  # storage of the jacobian matrix
  res_0 = zeros(Tjac, m)  # function evaluated at u0
  res_0_norm = 0.0  # norm of res_0
  delta_SL = zeros(Tjac, m)  # newton update
  step_norm = zero(Tjac)  # norm of newton update
  step_norm_1 = zero(Tjac) # norm of previous newton update


  # open file to write convergence data to
  # append to be on the safe side
  fconv = open("convergence.dat", "a+")


  # evaluating residual at initial condition
  println("evaluating residual at initial condition")
  res_0_norm = calcResidual(mesh, sbp, eqn, opts, func, res_0)

  # write rhs to file
  if write_rhs
    writedlm("rhs1.dat", res_0)
  end

  if res_0_norm < res_tol
   println("Initial condition satisfies res_tol with residual norm ", res_0_norm)
   println("writing to convergence.dat")
   println(fconv, i, " ", res_0_norm, " ", 0.0)

   close(fconv)

    if jac_type == 3
      destroyPetsc(jac, x, b, ksp)
    end


   return nothing
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

      elseif jac_type == 2 || jac_type == 3  # sparse jacobian
	println("calculating sparse FD jacobian")
        res_copy = copy(eqn.res)  # copy unperturbed residual
        @time calcJacobianSparse(mesh, sbp, eqn, opts, func, res_copy, pert, jac)
      end
      println("FD jacobian calculation @time printed above")

    elseif jac_method == 2
#      println("calculating complex step jacobian")

      if jac_type == 1  # dense jacobian
	println("calculating dense complex step jacobian")
        @time calcJacobianComplex(mesh, sbp, eqn, opts, func, pert, jac)

      elseif jac_type == 2 || jac_type == 3  # sparse jacobian 
	println("calculating sparse complex step jacobian")
        res_dummy = []  # not used, so don't allocation memory
        @time calcJacobianSparse(mesh, sbp, eqn, opts, func, res_dummy, pert, jac)
      end

      println("complex set jacobian calculate @time printed above")
    end

    # print as determined by options
    if write_jac
#      fname = string("jacobian", i, ".dat")
#      printMatrix(fname, jac)

      if jac_type == 3
        PetscMatAssemblyBegin(jac, PETSC_MAT_FINAL_ASSEMBLY)
        PetscMatAssemblyEnd(jac, PETSC_MAT_FINAL_ASSEMBLY)
      end


      writedlm("jacobian$i.dat", full(jac))
      println("finished printing jacobian")
    end

    # calculate Jacobian condition number
    if print_cond
      cond_j = cond(jac)
      println("Condition number of jacobian = ", cond_j)
    end

    # negate res
    for j=1:m
      res_0[j] = -res_0[j]
    end
    
    # calculate Newton step
    if jac_type == 1 || jac_type == 2  # julia jacobian
    @time delta_SL[:] = jac\(res_0)  #  calculate Newton update
    fill!(jac.nzval, 0.0)
#    @time solveMUMPS!(jac, res_0, delta_SL)
    elseif jac_type == 3   # petsc
      @time petscSolve(jac, x, b, ksp, res_0, delta_SL)
    end
    println("matrix solve @time prined above")
    step_norm = norm(delta_SL)/m
    println("step_norm = ", step_norm)

    # perform Newton update
    eqn.SL0[:] += step_fac*delta_SL  # update SL0

    # write starting values for next iteration to file
    if write_sol
      writedlm("SL0$i.dat", eqn.SL0)
    end

    # write paraview files
    if write_vis
      vals = abs(real(eqn.SL0))  # remove unneded imaginary part
      saveSolutionToMesh(mesh, vals)
      fname = string("solution_newton", i)
      writeVisFiles(mesh, fname)
    end
 

    # write to convergence file
    println(fconv, i, " ", res_0_norm, " ", step_norm)
    flush(fconv)

    # calculate residual at updated location, used for next iteration rhs
    res_0_norm = calcResidual(mesh, sbp, eqn, opts, func, res_0)

    # write rhs to file
    if write_rhs
      tmp = i+1
      writedlm("rhs$tmp.dat", res_0)
    end



   if res_0_norm < res_tol
     println("Newton iteration converged with residual norm ", res_0_norm)
     eqn.SL[:] = res_0
     close(fconv)

     if jac_type == 3
	destroyPetsc(jac, x, b, ksp)
      end


     return nothing
   end

    if (step_norm < step_tol)
      println("Newton iteration converged with step_norm = ", step_norm)
      println("Final residual = ", res_0_norm)
      eqn.SL[:] = res_0
      close(fconv)
      
      if jac_type == 3
	destroyPetsc(jac, x, b, ksp)
      end

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
=#
#    if (step_norm > step_norm_1)
#      step_fac /= 1.1
#    end


    print("\n")
    step_norm_1 = step_norm
  end  # end loop over newton iterations

  println("Warning: Newton iteration did not converge in ", itermax, " iterations")
  println("  Final step size: ", step_norm)
  println("  Final residual: ", res_0_norm)
  close(fconv)

  if jac_type == 3
    destroyPetsc(jac, x, b, ksp)
  end
  return nothing
end


function calcResidual(mesh, sbp, eqn, opts, func, res_0)
# calculate the residual and its norm

  m = length(res_0)

  eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.SL0)
  func(mesh, sbp, eqn, opts)
#  res_0[:] = real(eqn.SL)  # is there an unnecessary copy here?

  fill!(eqn.SL, 0.0)
  eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.SL)

  for j=1:m
    res_0[j] = real(eqn.SL[j])
  end

  res_0_norm = norm(eqn.SL)/m
  println("residual norm = ", res_0_norm)

 return res_0_norm
end


function calcJacFD(mesh, sbp, eqn, opts, func, res_0, pert, jac)
# calculate the jacobian using finite difference

  (m,n) = size(jac)
  entry_orig = zero(eltype(eqn.SL0))
  epsilon = norm(pert)  # finite difference perturbation
  # calculate jacobian
  for j=1:m
#      println("  jacobian iteration ", j)
    if j==1
      entry_orig = eqn.SL0[j]
      eqn.SL0[j] +=  epsilon
    else
      eqn.SL0[j-1] = entry_orig # undo previous iteration pertubation
      entry_orig = eqn.SL0[j]
      eqn.SL0[j] += epsilon
    end

    eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.SL0)
    # evaluate residual
    func(mesh, sbp, eqn, opts)
#     println("column ", j, " of jacobian, SL = ", eqn.SL)


    fill!(eqn.SL, 0.0)
    eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.SL)
    calcJacRow(unsafe_view(jac, :, j), res_0, eqn.SL, epsilon)
#      println("SL norm = ", norm(SL)/m)
    
  end

  # undo final perturbation
  eqn.SL0[m] = entry_orig


  return nothing
end



function calcJacobianSparse(mesh, sbp, eqn, opts, func, res_0, pert, jac::Union(SparseMatrixCSC, PetscMat))
# res_0 is 3d array of unperturbed residual, only needed for finite difference
# pert is perturbation to apply
# this function is independent of perturbation type

#  epsilon = 1e-6  # finite difference perturbation
  epsilon = norm(pert)  # get magnitude of perturbation
#  (m,n) = size(jac)
   m = length(res_0)
#  fill!(jac.nzval, 0.0)

  # for each color, store the perturbed element corresponding to each element
  perturbed_els = zeros(eltype(mesh.neighbor_nums), mesh.numEl)

  # debugging: do only first color
  for color=1:mesh.numColors  # loop over colors
#    println("color = ", color)
    getPertNeighbors(mesh, color, perturbed_els)
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
	  el_pert = perturbed_els[k] # get perturbed element
          if el_pert != 0   # if element was actually perturbed for this color

            col_idx = mesh.dofs[i, j, el_pert]
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

return nothing

end





function calcJacobianComplex(mesh, sbp, eqn, opts, func, pert, jac)

  epsilon = norm(pert)  # complex step perturbation
  entry_orig = zero(eltype(eqn.SL0))
  (m,n) = size(jac)
  # calculate jacobian
  for j=1:m
    if j==1
      entry_orig = eqn.SL0[j]
      eqn.SL0[j] +=  pert
    else
      eqn.SL0[j-1] = entry_orig # undo previous iteration pertubation
      entry_orig = eqn.SL0[j]
      eqn.SL0[j] += pert
    end

    eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.SL0)
    # evaluate residual
    func(mesh, sbp, eqn, opts)

    fill!(eqn.SL, 0.0)
    eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.SL)
#     println("column ", j, " of jacobian, SL = ", eqn.SL)
    calcJacRow(unsafe_view(jac, :, j), eqn.SL, epsilon)
#      println("SL norm = ", norm(SL)/m)
    
  end  # end loop over rows of jacobian


  # undo final perturbation
  eqn.SL0[m] = entry_orig
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
idx_tmp[1] = dof_pert - 1

pos = 1
for j_j = 1:mesh.numNodesPerElement
  for i_i = 1:mesh.numDofPerNode
    idy_tmp[pos] = mesh.dofs[i_i, j_j, el_res] - 1
#    col_idx = mesh.dofs[i, j, el_pert]

    vals_tmp[pos] = imag(eqn.res[i_i,j_j, el_res])/epsilon

    pos += 1
  end
end

PetscMatSetValues(jac, PetscInt(1), idx_tmp, local_size , idy_tmp, vals_tmp, PETSC_ADD_VALUES)

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




function getPertNeighbors(mesh, color, arr)
# populate the array with the element that is perturbed for each element
# element number == 0 if no perturbation

#  println("getting neighbor list")

  num_neigh = size(mesh.neighbor_colors, 1)
#  fill!(arr, 0)
  for i=1:mesh.numEl
    # find out if current element or its neighbors have the current color
    pos = 0
    for j=1:num_neigh
      if color == mesh.neighbor_colors[j, i]
	pos = j
	break
      end
    end

    if pos != 0
      arr[i] = mesh.neighbor_nums[pos, i]
    else
       arr[i] = 0
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
PetscInitialize(["-malloc", "-malloc_debug", "-malloc_dump"])
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
PetscMatSetType(A, MATMPIAIJ)
PetscMatSetSizes(A, PetscInt(mesh.numDof), PetscInt(mesh.numDof), PetscInt(mesh.numDof), PetscInt(mesh.numDof))

println("preallocating A")
# prellocate matrix
dnnz = zeros(PetscInt, mesh.numDof)  # diagonal non zeros per row
onnz = zeros(PetscInt, mesh.numDof)  # there is no off diagonal part for single proc case
dnnzu = zeros(PetscInt, 1)  # only needed for symmetric matrices
onnzu = zeros(PetscInt, 1)  # only needed for symmetric matrices
bs = PetscInt(1)  # block size

# calculate number of non zeros per row
for i=1:mesh.numDof
  max_dof = mesh.sparsity_bnds[2, i]
  min_dof = mesh.sparsity_bnds[1, i]
  nnz_i = max_dof - min_dof + 1
  dnnz[i] = nnz_i
#  println("row ", i," has ", nnz_i, " non zero entries")
end

PetscMatXAIJSetPreallocation(A, bs, dnnz, onnz, dnnzu, onnzu)

# zero initialize the matrix just in case
println("zeroing A")
PetscMatZeroEntries(A)

# create KSP contex
ksp = KSP(comm)
KSPSetOperators(ksp, A, A)
KSPSetFromOptions(ksp)

# set: rtol, abstol, dtol, maxits
KSPSetTolerances(ksp, 1e-14, 1e-14, 1e5, PetscInt(1000))
KSPSetUp(ksp)


return A, x, b, ksp

end

function destroyPetsc(A::PetscMat, x::PetscVec,  b::PetscVec, ksp::KSP)
# destory Petsc data structures, finalize Petsc

PetscDestroy(A)
PetscDestroy(x)
PetscDestroy(b)
PetscDestroy(ksp)
PetscFinalize()

return nothing

end


function petscSolve(A::PetscMat, x::PetscVec, b::PetscVec, ksp::KSP, res_0::AbstractVector, delta_SL::AbstractVector)

  # solve the system for the newton step, write it to delta_SL
  # writing it to delta_SL is an unecessary copy, becasue we could
  # write it directly to eqn.q, but for consistency we do it anyways

  # copy res_0 into b
  # create the index array
  numDof = length(delta_SL)
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

  matinfo = PetscMatGetInfo(A, MAT_LOCAL)
  println("number of mallocs = ", matinfo.mallocs)


  println("solving system")
  KSPSolve(ksp, b, x)
  reason = KSPGetConvergedReason(ksp)
  println("KSP converged reason = ", KSPConvergedReasonDict[reason])
  rnorm = KSPGetResidualNorm(ksp)
  println("Linear residual = ", rnorm)

  # copy solution from x into delta_SL
  arr, ptr_arr = PetscVecGetArray(x)
  for i=1:numDof
    delta_SL[i] = arr[i]
  end

  PetscVecRestoreArray(x, ptr_arr)

  # zero out the Jacobian for next use
  PetscMatZeroEntries(A)
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
  # the initial condition is stored in eqn.SL0
  # itermax is the maximum number of iterations

  step_fac = 0.5  # step size limiter
  m = length(eqn.SL)
  Tsol = typeof(eqn.SL[1])
  Tjac = typeof(real(eqn.SL[1]))  # type of jacobian, residual
  jac = zeros(Tjac, m, m)  # storage of the jacobian matrix
  direction_der = zeros(mesh.numDof)
#  v = rand(mesh.numDof)
   v = readdlm("randvec.txt")

  epsilon = 1e-20  # complex step perturbation
  fill!(eqn.SL, 0.0)  # zero out SL
  # compute directional derivative
  for i=1:mesh.numDof
    eqn.SL0[i] += complex(0, epsilon*v[i])  # apply perturbation
  end

  func(mesh, sbp, eqn, opts)
  println("evaluated directional derivative")

  # calculate derivative
  for i=1:mesh.numDof
    direction_der[i] = imag(eqn.SL[i])/epsilon
    eqn.SL0[i] -= complex(0, epsilon*v[i])  # undo perturbation
  end



    println("Calculating Jacobian")

    # calculate jacobian
    for j=1:m
      println("\ncalculating column ", j, " of the jacobian")
      if j==1
	eqn.SL0[j] +=  complex(0, epsilon)
      else
	eqn.SL0[j-1] -= complex(0, epsilon) # undo previous iteration pertubation
	eqn.SL0[j] += complex(0, epsilon)
      end

      # evaluate residual
      fill!(eqn.SL, 0.0)
      func(mesh, sbp, eqn, opts, eqn.SL0, eqn.SL)
 #     println("column ", j, " of jacobian, SL = ", eqn.SL)
      calcJacRow(unsafe_view(jac, :, j), eqn.SL, epsilon)
#      println("SL norm = ", norm(SL)/m)
      
    end  # end loop over rows of jacobian

    # undo final perturbation
    eqn.SL0[m] -= complex(0, epsilon)

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
    # jacobian multiplication into SL for return

    for i=1:mesh.numDof
      eqn.SL[i] = direction_der[i] - jac_mult[i]
    end

    err_norm = norm(eqn.SL)/mesh.numDof
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

#      eqn.SL0[j] += complex(0, epsilon)
      eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.SL0)
      eqn.q[1, 2, 5] += complex(0, epsilon)
      writedlm("check_q.dat", imag(eqn.q))
#      eqn.q[1,1,1] += complex(0, epsilon)
      # evaluate residual
      func(mesh, sbp, eqn, opts)

      fill!(eqn.SL, 0.0)
      eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.SL)
 #     println("column ", j, " of jacobian, SL = ", eqn.SL)
      calcJacRow(jac_col, eqn.SL, epsilon)
#      println("SL norm = ", norm(SL)/m)
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

     eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.SL0)
     func(mesh, sbp, eqn, opts, eqn.SL0, eqn.SL)
     fill!(eqn.SL, 0.0)
     eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.SL)
     res_0 = copy(eqn.SL)

      epsilon = 1e-6

      eqn.SL0[j] += epsilon

      eqn.disassmbleSolution(mesh, sbp, eqn, opts, eqn.SL0)

      # evaluate residual
      func(mesh, sbp, eqn, opts, eqn.SL0, eqn.SL)

      fill!(eqn.SL, 0.0)
      eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.SL)
 #     println("column ", j, " of jacobian, SL = ", eqn.SL)

      calcJacRow(jac_col, res_0, eqn.SL, epsilon)
#      println("SL norm = ", norm(SL)/m)

      return jac_col
end 
 
