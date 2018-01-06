"""
  This function computes the Jacobian of an evalResidual-like function.
  Specifically, it computes \\partial (eqn.q_vec)/ \\partial (eqn.res_vec).

  Other users of [`newtonInner`](@ref) (for example, implicit time marching
  methods) will need to implemnt their own version of this function.  The
  documentation for this function describes the requirements for the other
  implementations.

  [`newtonInner](@ref) guarantees that eqn.q and eqn.q_vec will be consistent
  when this function is called.  When using finite differences, eqn.res and
  eqn.res_vec must contain the residual of the physics.

  
  **Inputs**:
   
   * mesh: an AbstractMesh
   * sbp: an SBP operator
   * eqn: an AbstractSolutionData, eqn.res and eqn.res_vec may be overwritten
   * opts: options dictonary
   * jac: the Jacobian, can be an Array, SparseMatrixCSC, or PetscMat
   * ctx_residual: a tuple of values.  ctx_residual[1] must be an 
                   evalResidual-like function (eqn.q -> eqn.res) function
                   with signature func(mesh, sbp, eqn, opts, t).
                   See [`physicsRhs`](@ref) for a more thorough description
                    of ctx_residual.
   * t: simulation time

  **Options Keys:**
 
  **Implementation Notes:**

  This function should not have to do any parallel communication. `newtonInner`
  ensures that the `rhs_func` is called before `jac_func`, and `rhs_func` 
  handles
  the parallel communication.

  Implementations of this function may perform either (`eqn.q -> jacobian`) or
  (`eqn.q_vec -> jacobian`).  The first may be more computationally efficient,
  but the second can be simpler for debugging.

  This function supportes several types of jacobians (dense arrays,
  SparseMatrixCSC, PetscMat), and several methods for calculating them
  (finite difference and complex step).  Any function calling this function
  should support them as well.
  
  It is strongly recommneded to
  use this function to compute the spatial jacobian and them modify the
  resulting matrix (this function zeros the Jacobian matrix)

  When using Petsc matrices, the function may do intermediate assemblies
  (`PETSC_FLUSH_ASSEMBLY`), but does not need to do the final assembly.

"""
function physicsJac(mesh, sbp, eqn, opts, jac::AbstractMatrix,
                    ctx_residual, t=0.0;)

  #TODO: get rid of is_preconditioned
  #TODO: add start_comm option
  verbose = opts["newton_verbosity"]::Int

  myrank = mesh.myrank

  #TODO: figure out which of these are actually needed
  jac_method = opts["jac_method"]::Int  # finite difference or complex step
  jac_type = opts["jac_type"]::Int  # jacobian sparse or dense
  epsilon = opts["epsilon"]::Float64
  calc_jac_explicit = opts["calc_jac_explicit"]::Bool

  if jac_method == 1  # finite difference
    pert = epsilon
  elseif jac_method == 2  # complex step
    pert = complex(0, epsilon)
    removeComplex(mesh, sbp, eqn, opts)
  end

  # ctx_residual: func must be the first element
  func = ctx_residual[1]

  #----------------------------------------------------------------------
  # Calculate Jacobian using selected method 
  fill_zero!(jac)
  print_jacobian_timing = true
  eqn.params.time.t_jacobian += @elapsed if jac_method == 1
    @verbose5 @mpi_master println(BSTDOUT, "calculating finite difference jacobian")

    if jac_type == 1  # dense jacobian
      @verbose5 @mpi_master println(BSTDOUT, "calculating dense FD jacobian")
      # TODO: need to make q/q_vec and res/res_vec consistent
      # TODO: q/q_vec are guaranteed to be consistent, and res/res_vec don't
      #       matter (recompute res_vec internally?
      arrToVecAssign(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
      # res_0 is the unperturbed res, and it needs to be passed in vector form
      assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
      res_copy_vec = copy(eqn.res_vec)
      #TODO: don't copy the giant vector!
      tmp, t_jac, t_gc, alloc = @time_all calcJacFD(mesh, sbp, eqn, opts, func, res_copy_vec, pert, jac, t)

    elseif jac_type == 2  # Julia sparse jacobian
      @verbose5 @mpi_master println(BSTDOUT, "calculating sparse FD jacobian")
      #TODO: don't copy the giant array!
      res_copy = copy(eqn.res)  # copy unperturbed residual

      tmp, t_jac, t_gc, alloc = @time_all calcJacobianSparse(mesh, sbp, eqn, opts, func, res_copy, pert, jac, t)
    elseif jac_type == 3  # Petsc sparse jacobian
      @verbose5 @mpi_master println(BSTDOUT, "calculating sparse FD jacobian")
      res_copy = copy(eqn.res)  # copy unperturbed residual
      tmp, t_jac, t_gc, alloc = @time_all calcJacobianSparse(mesh, sbp, eqn, opts, func, res_copy, pert, jac, t)
    elseif jac_type == 4  # Petsc jacobian-vector product
      throw(ErrorException("No handling of jac_method = 1 and jac_type = 4: 
                           finite differencing isn't permitted for Petsc mat-free"))
    end

  elseif jac_method == 2
    @verbose5 @mpi_master println(BSTDOUT, "calculating complex step jacobian")

    if jac_type == 1  # dense jacobian
      @verbose5 @mpi_master println(BSTDOUT, "calculating dense complex step jacobian")
      tmp, t_jac, t_gc, alloc = @time_all calcJacobianComplex(mesh, sbp, eqn, opts, func, pert, jac, t)
    elseif jac_type == 2  # Julia sparse jacobian 
      @verbose5 @mpi_master println(BSTDOUT, "calculating sparse complex step jacobian")

      if calc_jac_explicit
        assembler = _AssembleElementData(jac, mesh, sbp, eqn, opts)
        tmp, t_jac, t_gc, alloc = @time_all evalJacobian(mesh, sbp, eqn, opts, assembler, t)
      else
        res_dummy = Array(Float64, 0, 0, 0)  # not used, so don't allocation memory
        tmp, t_jac, t_gc, alloc = @time_all calcJacobianSparse(mesh, sbp, eqn, 
                                            opts, func, res_dummy, pert, jac, t)
      end
    elseif jac_type == 3 # Petsc sparse jacobian

       @verbose5 @mpi_master println(BSTDOUT, "calculating Petsc jacobian")
      if calc_jac_explicit
        assembler = _AssembleElementData(jac, mesh, sbp, eqn, opts)
        tmp, t_jac, t_gc, alloc = @time_all evalJacobian(mesh, sbp, eqn, opts, assembler, t)
      else
        res_dummy = Array(Float64, 0, 0, 0)  # not used, so don't allocation memory

        tmp, t_jac, t_gc, alloc = @time_all calcJacobianSparse(mesh, sbp, eqn,
                                            opts, func, res_dummy, pert, jac, t)
      end


    end   # end of jac_type check

  end  # end of jac_method check

  # TODO: all printing should actually be handled outside of this function
  @verbose5 if print_jacobian_timing
    @mpi_master print(BSTDOUT, "jacobian calculation: ")
    @mpi_master print_time_all(BSTDOUT, t_jac, t_gc, alloc)
  end

  flush(BSTDOUT)

  return nothing

end   # end of physicsJac function


#------------------------------------------------------------------------------
# Functions for calculating the Jacobian
#------------------------------------------------------------------------------
@doc """
### NonlinearSolvers.calcJacFD

  This function calculates the Jacobian using finite differences, perturbing
  one degree of freedom at a time.  This is slow and not very accurate.  
  The Jacobian is calculated about the point in eqn.q_vec.

  Inputs:
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
function calcJacFD(mesh, sbp, eqn, opts, func, res_0, pert, jac::DenseArray, t=0.0)
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
function calcJacobianComplex(mesh, sbp, eqn, opts, func, pert, jac, t=0.0)

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

global const insert_freq = 1

"""
  Helper type for [`calcJacobianSparse`](@ref)
"""
type AssembleData{T <: AbstractMatrix}
  A::T
  # temporary arrays used to for Petsc MatSetValues
  insert_idx::Int
  local_size::Int
  vals_tmp::Array{Float64, 2}
  idx_tmp::Array{PetscInt, 1}
  idy_tmp::Array{PetscInt, 1}
end

function AssembleData{T}(A::T, mesh, sbp, eqn, opts)

  insert_idx = 1
  local_size = mesh.numNodesPerElement*mesh.numDofPerNode*insert_freq
  vals_tmp = zeros(local_size, 1) # values
  idx_tmp = zeros(PetscInt, local_size)  # row index
  idy_tmp = zeros(PetscInt, 1)  # column indices


  return AssembleData{T}(A, insert_idx, local_size, vals_tmp, idx_tmp,
                                 idy_tmp)
end

@doc """
### NonlinearSolvers.calcJacobianSparse

  This function calculate the Jacobian sparsely (only the entries 
    within the sparsity bounds), using either finite differences or algorithmic 
    differentiation.  The jacobian is calculated about the point stored in 
    eqn.q (not eqn.q_vec).  A mesh coloring approach is used to compute the
    jacobian.  Both eqn.q and the MPI send and receive buffers are perturbed
    during this process.

  Inputs:
    mesh: AbstractMesh
    sbp:  SBP operator
    eqn:  AbstractEquation object
    opts: options dictionary
    pert: perturbation to use for the algorithmic differentiation.  Currently,
          only complex numbers are supported.
    func: residual evaluation function (eqn.q -> eqn.res)
    res_0: element-based (3 dimensional) array containing the residual evaluated           at the point where the Jacobian is being calculated.
           This is only used for finite differences (can be a 0 x 0 x 0 array
           otherwise).

  Inputs/Outputs:
    jac:  Jacobian matrix.  Must be a sparse matrix type of some kind, 
          (including PetscMat).

  Aliasing restrictions: res_0 must not alias eqn.res

 
"""->
function calcJacobianSparse(mesh, sbp, eqn, opts, func,
                            res_0::Abstract3DArray, pert, 
                            jac::Union{SparseMatrixCSC, PetscMat}, t=0.0)
# res_0 is 3d array of unperturbed residual, only needed for finite difference
# pert is perturbation to apply
# this function is independent of perturbation type

#  filter_orig = eqn.params.use_filter  # record original filter state
#  eqn.params.use_filter = false  # don't repetatively filter

  # hold misc. data needed for assemble functions
  helper = AssembleData(jac, mesh, sbp, eqn, opts)

  epsilon = norm(pert)  # get magnitude of perturbation
  m = length(res_0)
  myrank = mesh.myrank
  f = eqn.params.f
  time = eqn.params.time
  time.t_color += @elapsed for color=1:mesh.maxColors  # loop over max colors, 
                                                       # only do calculation for
                                                       # numColors
    for j=1:mesh.numNodesPerElement  # loop over nodes
      for i=1:mesh.numDofPerNode  # loop over dofs on each node

        # apply perturbation to q
        if color <= mesh.numColors
          applyPerturbation(mesh, eqn.q, eqn.shared_data, color, pert, i, j, f)
          # evaluate residual
          time.t_func += @elapsed func(mesh, sbp, eqn, opts, t)
        end

        if !(color == 1 && j == 1 && i == 1) 
          assembly_end(jac, MAT_FLUSH_ASSEMBLY)
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
              assembleElement(helper, mesh, eqn, eqn.res, res_0, k, el_pert,
                              col_idx, epsilon, jac)
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
                assembleElement(helper, mesh, eqn, res_edge, res_0, k, el_pert, col_idx, epsilon, jac)
              end  # end if el_pert != 0
            end  # end loop over k
          end  # end loop over local edges
        end

        assembly_begin(jac, MAT_FLUSH_ASSEMBLY)

        # undo perturbation
        if color <= mesh.numColors
          applyPerturbation(mesh, eqn.q, eqn.shared_data, color, -pert, i, j)
        end
      end  # end loop i
    end  # end loop j
  end  # end loop over colors

  assembly_end(jac, MAT_FLUSH_ASSEMBLY)
#  flush(f)
  # now jac is complete
#  eqn.params.use_filter = filter_orig # reset filter
  return nothing

end  # end function


@doc """
### NonlinearSolvers.applyPerturbation

  This function applies a perturbation to a the specified degree of freedom
  on each element according to a mask.

  Because this is element based perturbation, opts["parallel_data"] must
  be "element".

  Inputs:
    mesh: an AbstractMesh
    color: the color to perturb
    pert: perturbation to apply.  Can be any datatype
    i: local degree of freedom number (in range 1:numDofPerNode) to perturb
    j: local node number (in range 1:numNodesPerElement) to perturb

  Inputs/Outputs:
    arr: element based (3D) array of values to perturb
    shared_data: array of SharedFaceData for ghost elements to be perturbed
                 The receive buffers are perturbed according to the masks
                 in mesh.shared_element_colormasks, the send buffers are
                 perturbed consistently with arr.

  Aliasing restrictions: none
"""->
function applyPerturbation{T}(mesh::AbstractMesh, arr::Abstract3DArray,
                           shared_data::Array{SharedFaceData{T}, 1},  
                           color::Integer, pert, i, j, f=BSTDOUT; 
                           perturb_shared=true)
  # applys perturbation pert to array arr according to a mask
  # color is the color currently being perturbed, used to select the mask
  # i, j specify the dof, node number within arr
  # the length of mask must equal the third dimension of arr
  # this function is independent of the type of pert

  @assert i <= size(arr, 1)
  @assert j <= size(arr, 2)

  # check that element, not face, data is shared in parallel
  for peer=1:length(shared_data)
    @assert size(shared_data[peer].q_send, 2) == mesh.numNodesPerElement
  end

  (ndof, nnodes, numel) = size(arr)
  mask = mesh.color_masks[color]
  for k=1:numel
    arr[i, j, k] += pert*mask[k]
  end

  if perturb_shared
    for peer=1:mesh.npeers
      # perturb receive buffer
      mask_i = mesh.shared_element_colormasks[peer][color]
      recv_arr_i = shared_data[peer].q_recv
      for k=1:length(mask_i)
        recv_arr_i[i, j, k] += pert*mask_i[k]
      end

      # perturb the send buffer, using the mask for eqn.q
      send_arr_i = shared_data[peer].q_send
#      bndries_local = shared_data[peer].bndries_local
      elnums_local = mesh.local_element_lists[peer]

      for k=1:length(elnums_local)
        el_k = elnums_local[k]
        send_arr_i[i, j, k] += pert*mask[el_k]
      end
    end
  end

  return nothing
end




#------------------------------------------------------------------------------
# Helper functions: assembleElement, calcJacCol for finite difference
#------------------------------------------------------------------------------
"""
  This function assembles the jacobian contribution for a single element
  into the matrix, when the jacobian is computed using finite differences.
  This is used by coloring-based methods for computing the jacobian.

  **Inputs:**

   * helper: a AssembleData object
   * mesh:  AbstractMesh object
   * eqn:  AbstractEquation
   * res_arr: element-based (3D) array of perturbed residual values
   * res_0:  element-based (3D) array of non-perturbed residual values
   * el_res: element number of the element we are observing the change in
   * el_pert: element number of the element that was perturbed
   * dof_pert: the degree of freedom number of the perturbed dof
   * epsilon: magnitude of perturbation

  **Inputs/Outputs:**

   * jac: any kind of matrix (dense, SparseMatrixCSC, PetscMat)

  Aliasing restrictions: res_arr and res_0 must not alias each other.
"""
function assembleElement{Tsol <: Real}(helper::AssembleData, mesh,
                         eqn::AbstractSolutionData{Tsol}, res_arr, res_0,
                         el_res::Integer, el_pert::Integer, dof_pert::Integer,
                         epsilon, jac::AbstractMatrix)
#

  # resize array
  # basically a no-op if array is already the right size
  local_size = PetscInt(mesh.numNodesPerElement*mesh.numDofPerNode)

  # get row number
  helper.idy_tmp[1] = dof_pert + mesh.dof_offset

  pos = 1
  for j_j = 1:mesh.numNodesPerElement
    for i_i = 1:mesh.numDofPerNode
      helper.idx_tmp[pos] = mesh.dofs[i_i, j_j, el_res] + mesh.dof_offset
  
      tmp = (res_arr[i_i,j_j, el_res] - res_0[i_i, j_j, el_res])/epsilon
      helper.vals_tmp[pos] = tmp
      pos += 1
    end
  end

  set_values1!(jac, helper.idx_tmp, helper.idy_tmp, helper.vals_tmp, ADD_VALUES)
  
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
"""
  Same as other method, but for complex numbers.  See that method for
  details.  res_0 is not used in this case
"""
function assembleElement{Tsol <: Complex}(helper::AssembleData, mesh,
                         eqn::AbstractSolutionData{Tsol}, res_arr, res_0,
                         el_res::Integer, el_pert::Integer, dof_pert::Integer,
                         epsilon, jac::AbstractMatrix)
# assemble an element contribution into jacobian
# making this a separate function enables dispatch on type of jacobian
# el_res is the element in the residual to assemble
# el_pert is the element that was perturbed
# dof_pert is the dof number (global) of the dof that was perturbed
# typically either el_pert or dof_pert will be needed, not both

  # get row number
  helper.idy_tmp[1] = dof_pert + mesh.dof_offset
  pos = 1
  for j_j = 1:mesh.numNodesPerElement
    for i_i = 1:mesh.numDofPerNode
      helper.idx_tmp[pos] = mesh.dofs[i_i, j_j, el_res] + mesh.dof_offset
      row = helper.idx_tmp[pos]
      col = helper.idy_tmp[1]


      helper.vals_tmp[pos] = imag(res_arr[i_i,j_j, el_res])/epsilon
      val = helper.vals_tmp[pos]
      pos += 1
    end
  end

  set_values1!(jac, helper.idx_tmp, helper.idy_tmp, helper.vals_tmp, ADD_VALUES)
#  PetscMatSetValues(jac, helper.idx_tmp, helper.idy_tmp, helper.vals_tmp, ADD_VALUES)

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


#------------------------------------------------------------------------------
# explicitly computed jacobian

"""
  Helper object for assembling element and interface jacobians into the
  system matrix.

  **Fields**

   * A: the matrix (can be Array, SparseMatrixCSC, or PetscMat)
   * idx: temporary array for row indices, length numDofPerNode
   * idy: temporary array for column indices, length numDofPerNode
   * vals: temporary array for matrix entries, size numDofPerNode square
   * idx_i: temporary array for row indices when assembling interface 
            jacobians, length 2 x numDofPerNode
   * idy_i: like idx_i, but for column indices
   * vals_i: temporary array for storing matrix entries when assembling
             interface jacobians, size 2 x numDofPerNode square
"""
type _AssembleElementData{T <: AbstractMatrix} <: AssembleElementData
  A::T

  # temporary array for element jacobian assembly
  idx::Array{PetscInt, 1}
  idy::Array{PetscInt, 1}
  vals::Array{PetscScalar, 2}

  # temporary arrays for interface jacobian assembly
  idx_i::Array{PetscInt, 1}  # 2x normal length
  idy_i::Array{PetscInt, 1}  # 2x normal length
  vals_i::Array{PetscScalar, 2}  # 2x by 2x normal size
                                 # note: this must be sized according the
                                 #       sparsity of the interface jacobian

  # temporary arrays for shared face jacobian assembly
  vals_sf::Array{PetscScalar, 2}  # normal length x 2x normal length
end


"""
  Outer constructor for [`_AssembleElementData`](@ref)
"""
function _AssembleElementData(A::AbstractMatrix, mesh, sbp, eqn, opts)

  idx = zeros(PetscInt, mesh.numDofPerNode)
  idy = zeros(PetscInt, mesh.numDofPerNode)
  vals = zeros(PetscScalar, mesh.numDofPerNode, mesh.numDofPerNode)

  idx_i = zeros(PetscInt, 2*mesh.numDofPerNode)
  idy_i = zeros(PetscInt, 2*mesh.numDofPerNode)
  vals_i = zeros(PetscScalar, 2*mesh.numDofPerNode, 2*mesh.numDofPerNode)

  vals_sf = zeros(PetscScalar, mesh.numDofPerNode, 2*mesh.numDofPerNode)

  return _AssembleElementData{typeof(A)}(A, idx, idy, vals, idx_i, idy_i, vals_i,
                                        vals_sf)
end

function _AssembleElementData()
  A = Array(PetscScalar, 0, 0)
  idx = Array(PetscInt, 0)
  idy = Array(PetscInt, 0)
  vals = Array(PetscScalar, 0, 0)

  idx_i = Array(PetscInt, 0)
  idy_i = Array(PetscInt, 0)
  vals_i = Array(PetscScalar, 0, 0)

  vals_sf = Array(PetscScalar, 0, 0)

  return _AssembleElementData{typeof(A)}(A, idx, idy, vals, idx_i, idy_i, vals_i,
                                        vals_sf)
end

"""
  An empty [`_AssembleElementData`](@ref).  Useful for giving a default value
  to fields.
"""
const NullAssembleElementData = _AssembleElementData()

#TODO: specialize for sparsity of sbpface
#TODO: use Petsc block matrix
#TODO: do flush assembly occasionally?

"""
  jac contains the data for the jacobian of the volume terms for a given
  element.  jac[i, j, p, q] = \\partial R[i, p, elnum] / \\partial eqn.q[j, q, elnum].
  Its size is numDofPerNode x numDofPerNode x numNodesPerElement x numNodesPerElement.



"""
function assembleElement{T}(helper::_AssembleElementData, mesh::AbstractMesh,
                            elnum::Integer, jac::AbstractArray{T, 4})

  #TODO: make a specialized version of this for block Petsc matrices

  numNodesPerElement = size(jac, 4)
  numDofPerNode = size(jac, 1)

  for q=1:numNodesPerElement

    # get dofs for node q
    for j=1:numDofPerNode
      helper.idy[j] = mesh.dofs[j, q, elnum] + mesh.dof_offset
    end
    for p=1:numNodesPerElement

      # get dofs for node p
      for i=1:numDofPerNode
        helper.idx[i] = mesh.dofs[i, p, elnum] + mesh.dof_offset
      end

      # get values
      for j=1:numDofPerNode
        for i=1:numDofPerNode
          helper.vals[i, j] = real(jac[i, j, p, q])
        end
      end

      # assemble them into matrix
      set_values1!(helper.A, helper.idx, helper.idy, helper.vals, ADD_VALUES)

    end  # end loop p
  end  # end loop q

  return nothing
end

"""
  jacAB where A = L or R and B = L or R, is the jacobian of the residual of
  element A with respect to the solution of element B.

  jacAB has the same size/layout as `jac` in [`assembleElement`](@ref).
"""
function assembleInterface{T}(helper::_AssembleElementData, sbpface::DenseFace,
                              mesh::AbstractMesh, iface::Interface,
                              jacLL::AbstractArray{T, 4},
                              jacLR::AbstractArray{T, 4},
                              jacRL::AbstractArray{T, 4},
                              jacRR::AbstractArray{T, 4})


  numNodesPerElement = size(jacLL, 4)
  numDofPerNode = size(jacLL, 1)

  permL = sview(sbpface.perm, :, iface.faceL)
  permR = sview(sbpface.perm, :, iface.faceR)

  for q in permL  # =1:numNodesPerElement

    # get indices for q
    for j=1:numDofPerNode
      helper.idy_i[j] = mesh.dofs[j, q, iface.elementL] + mesh.dof_offset
      helper.idy_i[j + numDofPerNode] = mesh.dofs[j, q, iface.elementR] + mesh.dof_offset
    end

    for p in permR # =1:numNodesPerElement

      # get indices for p
      for i=1:numDofPerNode
        helper.idx_i[i] = mesh.dofs[i, p, iface.elementL] + mesh.dof_offset
        helper.idx_i[i + numDofPerNode] = mesh.dofs[i, p, iface.elementR] + mesh.dof_offset
      end

      # put values into 2 x 2 block matrix
      @simd for j=1:numDofPerNode
        @simd for i=1:numDofPerNode
          helper.vals_i[i,                 j]                 = real(jacLL[i, j, p, q])
          helper.vals_i[i + numDofPerNode, j]                 = real(jacRL[i, j, p, q])
          helper.vals_i[i,                 j + numDofPerNode] = real(jacLR[i, j, p, q])
          helper.vals_i[i + numDofPerNode, j + numDofPerNode] = real(jacRR[i, j, p, q])
        end
      end

      set_values1!(helper.A, helper.idx_i, helper.idy_i, helper.vals_i, ADD_VALUES)
    end  # end loop q
  end  # end loop p

  return nothing
end

"""
  Assemble one half of an interface, used by shared face integrals.
  See [`assembleInterface`](@ref).
"""
function assembleSharedFace{T}(helper::_AssembleElementData, sbpface::DenseFace,
                               mesh::AbstractMesh,
                               iface::Interface,
                               jacLL::AbstractArray{T, 4},
                               jacLR::AbstractArray{T, 4})

  numNodesPerElement = size(jacLL, 4)
  numDofPerNode = size(jacLL, 1)

  permL = sview(sbpface.perm, :, iface.faceL)
  permR = sview(sbpface.perm, :, iface.faceR)


  for q in permL  # =1:numNodesPerElement

    # get indices for q
    for j=1:numDofPerNode
      helper.idy_i[j] = mesh.dofs[j, q, iface.elementL] + mesh.dof_offset
      helper.idy_i[j + numDofPerNode] = mesh.dofs[j, q, iface.elementR] + mesh.dof_offset
    end

    for p in permR # =1:numNodesPerElement

      # get indices for p
      for i=1:numDofPerNode
        helper.idx[i] = mesh.dofs[i, p, iface.elementL] + mesh.dof_offset
      end

      # put values into 2 x 2 block matrix
      for j=1:numDofPerNode
        for i=1:numDofPerNode
          helper.vals_sf[i, j]                 = real(jacLL[i, j, p, q])
          helper.vals_sf[i, j + numDofPerNode] = real(jacLR[i, j, p, q])
        end
      end

      set_values1!(helper.A, helper.idx, helper.idy_i, helper.vals_sf, ADD_VALUES)
    end  # end loop q
  end  # end loop p

  return nothing
end

 function assembleBoundary{T}(helper::_AssembleElementData, sbpface::DenseFace,
                               mesh::AbstractMesh,
                               bndry::Boundary,
                               jac::AbstractArray{T, 4})

  elnum = bndry.element
  numNodesPerElement = size(jac, 4)
  numDofPerNode = size(jac, 1)

  permL = sview(sbpface.perm, :, bndry.face)


  for q in permL # =1:numNodesPerElement

    # get dofs for node q
    for j=1:numDofPerNode
      helper.idy[j] = mesh.dofs[j, q, elnum] + mesh.dof_offset
    end
    for p in permL  # =1:numNodesPerElement

      # get dofs for node p
      for i=1:numDofPerNode
        helper.idx[i] = mesh.dofs[i, p, elnum] + mesh.dof_offset
      end

      # get values
      for j=1:numDofPerNode
        for i=1:numDofPerNode
          helper.vals[i, j] = real(jac[i, j, p, q])
        end
      end

      # assemble them into matrix
      set_values1!(helper.A, helper.idx, helper.idy, helper.vals, ADD_VALUES)

    end  # end loop p
  end  # end loop q

  return nothing
end


