# functions for computing the jacobian


"""
  This function computes the Jacobian of an evalResidual-like function.
  Specifically, it computes \\partial (eqn.q_vec)/ \\partial (eqn.res_vec),
  overwriting the Jacobian by default

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

  **Keyword Arguments**

   * zero_jac: if true, overwrite the Jacobian matrix, otherwise add to it,
               default true.

  **Options Keys:**

   * jac_type
   * jac_method
   * epsilon
   * calc_jac_explicit
 
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
                    ctx_residual, t=0.0; zero_jac=true)

  #TODO: add start_comm option

  # this function should not be called when A is a matrix-free jacobian,
  # but it can be called for the preconditioner matrix of jac_type == 4
  if typeof(jac) <: PetscMat
    @assert MatGetType(jac) != PETSc2.MATSHELL
  end

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
    removeComplex(eqn)
  end

  # ctx_residual: func must be the first element
  func = ctx_residual[1]

  #----------------------------------------------------------------------
  # Calculate Jacobian using selected method
  if zero_jac
    fill_zero!(jac)
  end
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
      array3DTo1D(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
      res_copy_vec = copy(eqn.res_vec)
      #TODO: don't copy the giant vector!
      tmp, t_jac, t_gc, alloc = @time_all calcJacFD(mesh, sbp, eqn, opts, func, res_copy_vec, pert, jac, t)

    elseif jac_type == 2  # Julia sparse jacobian
      @verbose5 @mpi_master println(BSTDOUT, "calculating sparse FD jacobian")
      #TODO: don't copy the giant array!
      res_copy = copy(eqn.res)  # copy unperturbed residual

      tmp, t_jac, t_gc, alloc = @time_all calcJacobianSparse(mesh, sbp, eqn, opts, func, res_copy, pert, jac, t)
    elseif jac_type == 3 || jac_type == 4  # Petsc sparse jacobian
      @verbose5 @mpi_master println(BSTDOUT, "calculating sparse FD jacobian")
      res_copy = copy(eqn.res)  # copy unperturbed residual
      tmp, t_jac, t_gc, alloc = @time_all calcJacobianSparse(mesh, sbp, eqn, opts, func, res_copy, pert, jac, t)
#    elseif jac_type == 4  # Petsc jacobian-vector product
#      throw(ErrorException("No handling of jac_method = 1 and jac_type = 4: 
#                           finite differencing isn't permitted for Petsc mat-free"))
    end

  elseif jac_method == 2
    @verbose5 @mpi_master println(BSTDOUT, "calculating complex step jacobian")

    if jac_type == 1  # dense jacobian
      @verbose5 @mpi_master println(BSTDOUT, "calculating dense complex step jacobian")
      tmp, t_jac, t_gc, alloc = @time_all calcJacobianComplex(mesh, sbp, eqn, opts, func, pert, jac, t)
    elseif jac_type == 2  # Julia sparse jacobian 
      @verbose5 @mpi_master println(BSTDOUT, "calculating sparse complex step jacobian")

      if calc_jac_explicit
        println(BSTDOUT, "calculating jacobian explicitly")
        assembler = _AssembleElementData(jac, mesh, sbp, eqn, opts)
        tmp, t_jac, t_gc, alloc = @time_all evalJacobian(mesh, sbp, eqn, opts, assembler, t)
      else
        println(BSTDOUT, "calculating jacobian with coloring")
        res_dummy = Array{Float64}(0, 0, 0)  # not used, so don't allocation memory
        tmp, t_jac, t_gc, alloc = @time_all calcJacobianSparse(mesh, sbp, eqn, 
                                            opts, func, res_dummy, pert, jac, t)
      end
    elseif jac_type == 3 || jac_type == 4 # Petsc sparse jacobian

       @verbose5 @mpi_master println(BSTDOUT, "calculating Petsc jacobian")
      if calc_jac_explicit
        println("calculating Jacobian explicitly")
        assembler = _AssembleElementData(jac, mesh, sbp, eqn, opts)
        tmp, t_jac, t_gc, alloc = @time_all evalJacobian(mesh, sbp, eqn, opts, assembler, t)
      else
        res_dummy = Array{Float64}(0, 0, 0)  # not used, so don't allocation memory

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
### Jacobians.calcJacFD

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


    array1DTo3D(mesh, sbp, eqn, opts, eqn.q_vec)
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
### Jacobians.calcJacComplex

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

    array1DTo3D(mesh, sbp, eqn, opts, eqn.q_vec)
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
mutable struct AssembleData{T <: AbstractMatrix}
  A::T
  # temporary arrays used to for Petsc MatSetValues
  insert_idx::Int
  local_size::Int
  vals_tmp::Array{Float64, 2}
  idx_tmp::Array{PetscInt, 1}
  idy_tmp::Array{PetscInt, 1}
end

function AssembleData(A::T, mesh, sbp, eqn, opts) where T

  insert_idx = 1
  local_size = mesh.numNodesPerElement*mesh.numDofPerNode*insert_freq
  vals_tmp = zeros(local_size, 1) # values
  idx_tmp = zeros(PetscInt, local_size)  # row index
  idy_tmp = zeros(PetscInt, 1)  # column indices


  return AssembleData{T}(A, insert_idx, local_size, vals_tmp, idx_tmp,
                                 idy_tmp)
end

@doc """
### Jacobians.calcJacobianSparse

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
### Jacobians.applyPerturbation

  This function applies a perturbation to a the specified degree of freedom
  on each element according to a mask.

  Because this is element based perturbation, opts["parallel_data"] must
  be PARALLEL_DATA_ELEMENT.

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
function applyPerturbation(mesh::AbstractMesh, arr::Abstract3DArray,
                        shared_data::Array{SharedFaceData{T}, 1},  
                        color::Integer, pert, i, j, f=BSTDOUT; 
                        perturb_shared=true) where T
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
function assembleElement(helper::AssembleData, mesh,
           eqn::AbstractSolutionData{Tsol}, res_arr, res_0,
           el_res::Integer, el_pert::Integer, dof_pert::Integer,
           epsilon, jac::AbstractMatrix) where Tsol <: Real
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
### Jacobians.calcJacCol

  This function extracts the entries for one column of the Jacobian from two residual evaluates that come from finite differences.

  Inputs:
    res_0: vector of unperturbed residual values
    res: vector of perturbed residual values
    epsilon: magnitude of perturbation

  Inputs/Outputs:
    jac_row = vector to be populated with the Jacobian entries

  Aliasing restrictions: res_0 and res cannot alias (obviously).

"""->
function calcJacCol(jac_row, res_0, res::AbstractArray{T,1}, epsilon) where T <: Real
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
function assembleElement(helper::AssembleData, mesh,
        eqn::AbstractSolutionData{Tsol}, res_arr, res_0,
        el_res::Integer, el_pert::Integer, dof_pert::Integer,
        epsilon, jac::AbstractMatrix) where Tsol <: Complex
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
### Jacobians.calcJacCol

  This function extracts the entries for one column of the Jacobian from a 
  complex step residual evaluation

  Inputs:
    res: vector of perturbed residual values
    epsilon: magnitude of perturbation

  Inputs/Outputs:
    jac_row = vector to be populated with the Jacobian entries

  Aliasing restrictions: none

"""->
function calcJacCol(jac_row, res::AbstractArray{T, 1}, epsilon) where T <: Complex
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

global const assem_min_volume_nodes = 3  # minimum number of volume nodes
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
mutable struct _AssembleElementData{T <: AbstractMatrix} <: AssembleElementData
  A::T

  # temporary array for element jacobian assembly
  idx::Array{PetscInt, 1}
  idy::Array{PetscInt, 1}
  vals::Array{PetscScalar, 2}

  # temporary arrays for block element jacobian assembly
  idx_b::Array{PetscInt, 1}
  idy_b::Array{PetscInt, 1}
  vals_b::Array{PetscScalar, 2}

  # temporary arrays for interface jacobian assembly
  idx_i::Array{PetscInt, 1}  # 2x normal length
  idy_i::Array{PetscInt, 1}  # 2x normal length
  vals_i::Array{PetscScalar, 2}  # 2x by 2x normal size
                                 # note: this must be sized according the
                                 #       sparsity of the interface jacobian

  idx_ib::Array{PetscInt, 1}  # length 2, used for block assembly
  idy_ib::Array{PetscInt, 1}


  # temporary arrays for shared face jacobian assembly
  vals_sf::Array{PetscScalar, 2}  # normal length x 2x normal length

  # temporary arrays for block boundary face jacobian assembly
  idx_bb::Array{PetscInt, 1}
  idy_bb::Array{PetscInt, 1}

  nonstencil_nodes::Array{Int, 2}
end


"""
  Outer constructor for [`_AssembleElementData`](@ref)
"""
function _AssembleElementData(A::AbstractMatrix, mesh, sbp, eqn, opts)

  idx = zeros(PetscInt, mesh.numDofPerNode)
  idy = zeros(PetscInt, mesh.numDofPerNode)
  vals = zeros(PetscScalar, mesh.numDofPerNode, mesh.numDofPerNode)

  idx_b = zeros(PetscInt, assem_min_volume_nodes)
  idy_b = zeros(PetscInt, assem_min_volume_nodes)
  nentries = assem_min_volume_nodes*mesh.numDofPerNode
  vals_b = zeros(PetscInt, nentries, nentries)

  idx_i = zeros(PetscInt, 2*mesh.numDofPerNode)
  idy_i = zeros(PetscInt, 2*mesh.numDofPerNode)
  vals_i = zeros(PetscScalar, 2*mesh.numDofPerNode, 2*mesh.numDofPerNode)
  idx_ib = zeros(PetscInt, 2)
  idy_ib = zeros(PetscInt, 2)

  vals_sf = zeros(PetscScalar, mesh.numDofPerNode, 2*mesh.numDofPerNode)

  idx_bb = zeros(PetscInt, 1)
  idy_bb = zeros(PetscInt, 1)

  nonstencil_nodes = getNonStencilNodes(sbp, mesh.sbpface)

  return _AssembleElementData{typeof(A)}(A, idx, idy, vals, idx_b, idy_b, vals_b,
                                         idx_i, idy_i,
                                         vals_i, idx_ib, idy_ib, vals_sf,
                                         idx_bb, idy_bb,
                                         nonstencil_nodes)
end

function _AssembleElementData()
  A = Array{PetscScalar}(0, 0)
  idx = Array{PetscInt}(0)
  idy = Array{PetscInt}(0)
  vals = Array{PetscScalar}(0, 0)
  idx_b = Array{PetscInt}(0)
  idy_b = Array{PetscInt}(0)
  vals_b = Array{PetscScalar}(0, 0)

  idx_i = Array{PetscInt}(0)
  idy_i = Array{PetscInt}(0)
  vals_i = Array{PetscScalar}(0, 0)
  idx_ib = Array{PetscInt}(0)
  idy_ib = Array{PetscInt}(0)

  vals_sf = Array{PetscScalar}(0, 0)

  idx_bb = Array{PetscInt}(0)
  idy_bb = Array{PetscInt}(0)

  nonstencil_nodes = Array{Int}(0, 0)

  return _AssembleElementData{typeof(A)}(A, idx, idy, vals, idx_b, idy_b, vals_b,
                                         idx_i, idy_i,
                                         vals_i, idx_ib, idy_ib, vals_sf,
                                         idx_bb, idy_bb,
                                         nonstencil_nodes)
end

"""
  An empty [`_AssembleElementData`](@ref).  Useful for giving a default value
  to fields.
"""
const NullAssembleElementData = _AssembleElementData()

#TODO: do flush assembly occasionally?

"""
  **Inputs**

   * helper: an _AssembleElementData
   * mesh: a mesh
   * elnum: element number
   * jac: 4 dimensional array containing the jacobian of the element

  jac contains the data for the jacobian of the volume terms for a given
  element.  \$jac[i, j, p, q] = \\partial R[i, p, elnum] / \\partial eqn.q[j, q, elnum]\$.
  Its size is numDofPerNode x numDofPerNode x numNodesPerElement x numNodesPerElement.

"""
function assembleElement(helper::_AssembleElementData, mesh::AbstractMesh,
                         elnum::Integer, jac::AbstractArray{T, 4}) where T

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
      @simd for j=1:numDofPerNode
        @simd for i=1:numDofPerNode
          helper.vals[i, j] = real(jac[i, j, p, q])
        end
      end

      # assemble them into matrix
      set_values1!(helper.A, helper.idx, helper.idy, helper.vals, ADD_VALUES)

    end  # end loop p
  end  # end loop q

  return nothing
end

# non-tiled version using MatSetValuesBlocked
#=
function assembleElement(helper::_AssembleElementData{PetscMat}, mesh::AbstractMesh,
                            elnum::Integer, jac::AbstractArray{T, 4})  where {T}

  numNodesPerElement = size(jac, 4)
  numDofPerNode = size(jac, 1)

  for q=1:numNodesPerElement

    dof1 = mesh.dofs[1, q, elnum] + mesh.dof_offset
    helper.idy_bb[1] = div(dof1 - 1, numDofPerNode) + 1 - 1

    for p=1:numNodesPerElement
      dof1 = mesh.dofs[1, p, elnum] + mesh.dof_offset
      helper.idx_bb[1] = div(dof1 - 1, numDofPerNode) + 1 - 1


      # get values
      @simd for j=1:numDofPerNode
        @simd for i=1:numDofPerNode
          helper.vals[i, j] = real(jac[i, j, p, q])
        end
      end

      # assemble them into matrix

      MatSetValuesBlocked(helper.A, helper.idx_bb, helper.idy_bb, helper.vals, ADD_VALUES)
#      set_values1!(helper.A, helper.idx, helper.idy, helper.vals, ADD_VALUES)

    end  # end loop p
  end  # end loop q

  return nothing
end
=#

# This tiled version
function assembleElement(helper::_AssembleElementData{PetscMat}, mesh::AbstractMesh,
                         elnum::Integer, jac::AbstractArray{T, 4}) where T

  numNodesPerElement = size(jac, 4)
  numDofPerNode = size(jac, 1)

  # use a tiling approach to do fewer, larger calls to Petsc
  nblocks = div(numNodesPerElement, assem_min_volume_nodes)
  nrem = numNodesPerElement - nblocks*assem_min_volume_nodes

  for yblock=1:nblocks
    yoffset = (yblock-1)*assem_min_volume_nodes

    for j=1:assem_min_volume_nodes
      idx = yoffset + j
      dof1 = mesh.dofs[1, idx, elnum] + mesh.dof_offset
      helper.idy_b[j] = div(dof1 - 1, numDofPerNode) + 1 - 1

    end

    for xblock=1:nblocks
      xoffset = (xblock-1)*assem_min_volume_nodes

      for i=1:assem_min_volume_nodes
        idx = xoffset + i
        dof1 = mesh.dofs[1, idx, elnum] + mesh.dof_offset
        helper.idx_b[i] = div(dof1 - 1, numDofPerNode) + 1 - 1
      end

      # get values
      for q=1:assem_min_volume_nodes
        for p=1:assem_min_volume_nodes
          @simd for i=1:numDofPerNode
            i_idx = i + (q-1)*numDofPerNode
            @simd for j=1:numDofPerNode
              j_idx = j + (p-1)*numDofPerNode
              helper.vals_b[j_idx, i_idx] = real(jac[j, i, xoffset + p , yoffset + q])
            end
          end
        end
      end

      MatSetValuesBlocked(helper.A, helper.idx_b, helper.idy_b, helper.vals_b, ADD_VALUES)

    end
  end

  # cleanup loops
  for _q=1:nrem
    q = nblocks*assem_min_volume_nodes + _q

    # get dofs for node q
    dof1 = mesh.dofs[1, q, elnum] + mesh.dof_offset
    helper.idy_bb[1] = div(dof1 - 1, numDofPerNode) + 1 - 1


    for p=1:numNodesPerElement
      dof1 = mesh.dofs[1, p, elnum] + mesh.dof_offset
      helper.idx_bb[1] = div(dof1 - 1, numDofPerNode) + 1 - 1

      # get values
      @simd for j=1:numDofPerNode
        @simd for i=1:numDofPerNode
          helper.vals[i, j] = real(jac[i, j, p, q])
        end
      end


      MatSetValuesBlocked(helper.A, helper.idx_bb, helper.idy_bb, helper.vals, ADD_VALUES)
    end
  end

  for q=1:(numNodesPerElement-1)
    # get dofs for node q
    dof1 = mesh.dofs[1, q, elnum] + mesh.dof_offset
    helper.idy_bb[1] = div(dof1 - 1, numDofPerNode) + 1 - 1

    for _p=1:nrem
      p = nblocks*assem_min_volume_nodes + _p
      dof1 = mesh.dofs[1, p, elnum] + mesh.dof_offset
      helper.idx_bb[1] = div(dof1 - 1, numDofPerNode) + 1 - 1

      @simd for j=1:numDofPerNode
        @simd for i=1:numDofPerNode
          helper.vals[i, j] = real(jac[i, j, p, q])
        end
      end


      MatSetValuesBlocked(helper.A, helper.idx_bb, helper.idy_bb, helper.vals, ADD_VALUES)
    end
  end


  return nothing
end


"""
  Assembles the jacobian of an interface into the matrix.  Specialized
  versions take advantage of the sparsity of the `sbpface`.

  **Inputs**

   * helper: _AssembleElementData
   * sbpface: an SBP face object
   * mesh: a mesh
   * iface: an Interface object identify the interface to be assembled
   * jacLL: see below
   * jacLR:
   * jacRL
   * jacRR

  jacAB where A = L or R and B = L or R, is the jacobian of the residual of
  element A with respect to the solution of element B.

  jacAB has the same size/layout as `jac` in [`assembleElement`](@ref).
"""
function assembleInterface(helper::_AssembleElementData, 
                           sbpface::DenseFace,
                           mesh::AbstractMesh, iface::Interface,
                           jacLL::AbstractArray{T, 4},
                           jacLR::AbstractArray{T, 4},
                           jacRL::AbstractArray{T, 4},
                           jacRR::AbstractArray{T, 4}) where T


  numNodesPerElement = size(jacLL, 4)
  numDofPerNode = size(jacLL, 1)
  idx1 = 1:mesh.numDofPerNode
  idx2 = (mesh.numDofPerNode+1):(2*mesh.numDofPerNode)

  permL = sview(sbpface.perm, :, iface.faceL)
  permR = sview(sbpface.perm, :, iface.faceR)

#  for q in permL  # =1:numNodesPerElement
  for q=1:sbpface.stencilsize
    qL = sbpface.perm[q, iface.faceL]
    qR = sbpface.perm[q, iface.faceR]

    # get indices for q
    for j=1:numDofPerNode
      helper.idy_i[j]                 = mesh.dofs[j, qL, iface.elementL] + mesh.dof_offset
      helper.idy_i[j + numDofPerNode] = mesh.dofs[j, qR, iface.elementR] + mesh.dof_offset
    end

#    for p in permR # =1:numNodesPerElement
    for p=1:sbpface.stencilsize
      pL = sbpface.perm[p, iface.faceL]
      pR = sbpface.perm[p, iface.faceR]


      # get indices for p
      for i=1:numDofPerNode
        helper.idx_i[i]                 = mesh.dofs[i, pL, iface.elementL] + mesh.dof_offset
        helper.idx_i[i + numDofPerNode] = mesh.dofs[i, pR, iface.elementR] + mesh.dof_offset
      end

      # put values into 2 x 2 block matrix
      @simd for j=1:numDofPerNode
        @simd for i=1:numDofPerNode
          helper.vals_i[i,                 j]                 = real(jacLL[i, j, pL, qL])
          helper.vals_i[i + numDofPerNode, j]                 = real(jacRL[i, j, pR, qL])
          helper.vals_i[i,                 j + numDofPerNode] = real(jacLR[i, j, pL, qR])
          helper.vals_i[i + numDofPerNode, j + numDofPerNode] = real(jacRR[i, j, pR, qR])
        end
      end

      set_values1!(helper.A, helper.idx_i, helper.idy_i, helper.vals_i, ADD_VALUES)
    end  # end loop q
  end  # end loop p

  return nothing
end


function assembleInterface(helper::_AssembleElementData, 
                           sbpface::SparseFace,
                           mesh::AbstractMesh, iface::Interface,
                           jacLL::AbstractArray{T, 4},
                           jacLR::AbstractArray{T, 4},
                           jacRL::AbstractArray{T, 4},
                           jacRR::AbstractArray{T, 4}) where T


  numNodesPerElement = size(jacLL, 4)
  numDofPerNode = size(jacLL, 1)

#  vals_before = copy(helper.vals_i)
#  vals_after = copy(helper.vals_i)
  for p=1:sbpface.numnodes
#    println("face node ", p)
    # row indices
    iR = sbpface.nbrperm[p, iface.orient]
    pL = sbpface.perm[p, iface.faceL]
    pR = sbpface.perm[iR, iface.faceR]

    # column indices
    # this matches the way SBP puts the data into the jac arrays,
    # but its a little weird
    qL = pL
    qR = pR

    # get indices for p
    for i=1:numDofPerNode
      helper.idx_i[i]                 = mesh.dofs[i, pL, iface.elementL] + mesh.dof_offset
      helper.idx_i[i + numDofPerNode] = mesh.dofs[i, pR, iface.elementR] + mesh.dof_offset

      helper.idy_i[i]                 = mesh.dofs[i, qL, iface.elementL] + mesh.dof_offset
      helper.idy_i[i + numDofPerNode] = mesh.dofs[i, qR, iface.elementR] + mesh.dof_offset
    end

    # put values into 2 x 2 block matrix
    @simd for j=1:numDofPerNode
      @simd for i=1:numDofPerNode
        helper.vals_i[i,                 j]                 = real(jacLL[i, j, pL, qL])
        helper.vals_i[i + numDofPerNode, j]                 = real(jacRL[i, j, pR, qL])
        helper.vals_i[i,                 j + numDofPerNode] = real(jacLR[i, j, pL, qR])
        helper.vals_i[i + numDofPerNode, j + numDofPerNode] = real(jacRR[i, j, pR, qR])
      end
    end

#    get_values1!(helper.A, helper.idx_i, helper.idy_i, vals_before)
    set_values1!(helper.A, helper.idx_i, helper.idy_i, helper.vals_i, ADD_VALUES)
#    get_values1!(helper.A, helper.idx_i, helper.idy_i, vals_after)
#    println("vals_before = \n", vals_before)
#    println("new values = \n", helper.vals_i)
#    println("vals after = \n", vals_after)
  end  # end loop p

  return nothing
end

"""
  Assemble an interface for the viscous terms when T4 = 0.  Same signature as
  `assembleInterface`.

  Note that this does not special case diagonal E operators because that may not
  be correct for all viscoscity modes (ie. ShockDiffusion).
"""
function assembleInterfaceVisc(helper::_AssembleElementData, 
                           sbpface::AbstractFace,
                           mesh::AbstractMesh, iface::Interface,
                           jacLL::AbstractArray{T, 4},
                           jacLR::AbstractArray{T, 4},
                           jacRL::AbstractArray{T, 4},
                           jacRR::AbstractArray{T, 4}) where T

  elL = iface.elementL; elR = iface.elementR
  stencilL = sview(sbpface.perm, :, iface.faceL)
  stencilR = sview(sbpface.perm, :, iface.faceR)
  nonstencilL = sview(helper.nonstencil_nodes, :, iface.faceL)
  nonstencilR = sview(helper.nonstencil_nodes, :, iface.faceR)

  idx = helper.idx_i
  idy = helper.idy_i
  vals = helper.vals_i

  # do the square blocks
  for _q=1:length(stencilL)
    qL = stencilL[_q]
    qR = stencilR[_q]

    for j=1:mesh.numDofPerNode
      idy[j                     ] = mesh.dofs[j, qL, elL]
      idy[j + mesh.numDofPerNode] = mesh.dofs[j, qR, elR]
    end


    for _p =1:length(stencilL)  # assume stencilL and R are the same length
      pL = stencilL[_p]
      pR = stencilR[_p]
      for i=1:mesh.numDofPerNode
        idx[i                     ] = mesh.dofs[i, pL, elL]
        idx[i + mesh.numDofPerNode] = mesh.dofs[i, pR, elR]
      end

      for j=1:mesh.numDofPerNode
        for i=1:mesh.numDofPerNode
          vals[i,                      j                     ] = real(jacLL[i, j, pL, qL])
          vals[i + mesh.numDofPerNode, j                     ] = real(jacRL[i, j, pR, qL])
          vals[i,                      j + mesh.numDofPerNode] = real(jacLR[i, j, pL, qR])
          vals[i + mesh.numDofPerNode, j + mesh.numDofPerNode] = real(jacRR[i, j, pR, qR])
        end
      end  # end j

      set_values1!(helper.A, idx, idy, vals, ADD_VALUES)
    end  # end p
  end  # end q

  # do the remaining rectangles

  # columns
  for _q = 1:length(stencilL)
    qL = stencilL[_q]
    qR = stencilR[_q]

    for j=1:mesh.numDofPerNode
      idy[j                     ] = mesh.dofs[j, qL, elL]
      idy[j + mesh.numDofPerNode] = mesh.dofs[j, qR, elR]
    end

    for _p=1:length(nonstencilL)  # assume nonstencilL and R are the same length
      pL = nonstencilL[_p]
      pR = nonstencilR[_p]

      for i=1:mesh.numDofPerNode
        idx[i                     ] = mesh.dofs[i, pL, elL]
        idx[i + mesh.numDofPerNode] = mesh.dofs[i, pR, elR]
      end

      for j=1:mesh.numDofPerNode
        for i=1:mesh.numDofPerNode
          vals[i,                      j                     ] = real(jacLL[i, j, pL, qL])
          vals[i + mesh.numDofPerNode, j                     ] = real(jacRL[i, j, pR, qL])
          vals[i,                      j + mesh.numDofPerNode] = real(jacLR[i, j, pL, qR])
          vals[i + mesh.numDofPerNode, j + mesh.numDofPerNode] = real(jacRR[i, j, pR, qR])
        end
      end  # end j

      set_values1!(helper.A, idx, idy, vals, ADD_VALUES)
    end  # end p
  end  # end q

  # rows
  for _p = 1:length(stencilL)
    pL = stencilL[_p]
    pR = stencilR[_p]

    for i=1:mesh.numDofPerNode
      idx[i                     ] = mesh.dofs[i, pL, elL]
      idx[i + mesh.numDofPerNode] = mesh.dofs[i, pR, elR]
    end

    for _q=1:length(nonstencilL)
      qL = nonstencilL[_q]
      qR = nonstencilR[_q]

      for j=1:mesh.numDofPerNode
        idy[j                     ] = mesh.dofs[j, qL, elL]
        idy[j + mesh.numDofPerNode] = mesh.dofs[j, qR, elR]

        for i=1:mesh.numDofPerNode
          vals[i,                      j                     ] = real(jacLL[i, j, pL, qL])
          vals[i + mesh.numDofPerNode, j                     ] = real(jacRL[i, j, pR, qL])
          vals[i,                      j + mesh.numDofPerNode] = real(jacLR[i, j, pL, qR])
          vals[i + mesh.numDofPerNode, j + mesh.numDofPerNode] = real(jacRR[i, j, pR, qR])
        end
      end  # end j

      set_values1!(helper.A, idx, idy, vals, ADD_VALUES)
    end  # end q
  end  # end p

  # remaining on diagonal entries
  # This is only needed for shock capturing (elements where dLambda/dq != 0)
  idx = helper.idx
  idy = helper.idy
  vals = helper.vals

  for _q=1:length(nonstencilL)
    qL = nonstencilL[_q]
    qR = nonstencilR[_q]

    # block indices (zero-based)
    for j=1:mesh.numDofPerNode
      idy[j] = mesh.dofs[j, qL, elL]
    end

    for _p=1:length(nonstencilL)
      pL = nonstencilL[_p]

      # block indices (zero-based)
      for i=1:mesh.numDofPerNode
        idx[i] = mesh.dofs[i, pL, elL]
      end

      @simd for j=1:mesh.numDofPerNode
        @simd for i=1:mesh.numDofPerNode
          vals[i, j] = real(jacLL[i, j, pL, qL])
        end
      end

      set_values1!(helper.A, idx, idy, vals, ADD_VALUES)
    end
    
    for j=1:mesh.numDofPerNode
      idy[j] = mesh.dofs[j, qR, elR]
    end

    for _p=1:length(nonstencilR)
      pR = nonstencilR[_p]

      for i=1:mesh.numDofPerNode
        idx[i] = mesh.dofs[i, pR, elR]
      end

      @simd for j=1:mesh.numDofPerNode
        @simd for i=1:mesh.numDofPerNode
          vals[i, j] = real(jacRR[i, j, pR, qR])
        end
      end

      set_values1!(helper.A, idx, idy, vals, ADD_VALUES)
    end  # end p
    
  end  # end q



  return nothing
end





function assembleInterface(helper::_AssembleElementData{PetscMat}, 
                           sbpface::DenseFace,
                           mesh::AbstractMesh, iface::Interface,
                           jacLL::AbstractArray{T, 4},
                           jacLR::AbstractArray{T, 4},
                           jacRL::AbstractArray{T, 4},
                           jacRR::AbstractArray{T, 4}) where T

  numNodesPerElement = size(jacLL, 4)
  numDofPerNode = size(jacLL, 1)

  permL = sview(sbpface.perm, :, iface.faceL)
  permR = sview(sbpface.perm, :, iface.faceR)

#  for q in permL  # =1:numNodesPerElement
  for q=1:sbpface.stencilsize
    qL = sbpface.perm[q, iface.faceL]
    qR = sbpface.perm[q, iface.faceR]

    # compute block indices for q (zero-based)
    dof1 = mesh.dofs[1, qL, iface.elementL] + mesh.dof_offset
    helper.idy_ib[1] = div(dof1 - 1, numDofPerNode) + 1 - 1
    dof1 = mesh.dofs[1, qR, iface.elementR] + mesh.dof_offset
    helper.idy_ib[2] = div(dof1 - 1, numDofPerNode) + 1 - 1


#    for p in permR # =1:numNodesPerElement
    for p=1:sbpface.stencilsize
      pL = sbpface.perm[p, iface.faceL]
      pR = sbpface.perm[p, iface.faceR]

      dof1 = mesh.dofs[1, pL, iface.elementL] + mesh.dof_offset
      helper.idx_ib[1] = div(dof1 - 1, numDofPerNode) + 1 - 1
      dof1 = mesh.dofs[1, pR, iface.elementR] + mesh.dof_offset
      helper.idx_ib[2] = div(dof1 - 1, numDofPerNode) + 1 - 1

      # put values into 2 x 2 block matrix
      @simd for j=1:numDofPerNode
        @simd for i=1:numDofPerNode
          helper.vals_i[i,                 j]                 = real(jacLL[i, j, pL, qL])
          helper.vals_i[i + numDofPerNode, j]                 = real(jacRL[i, j, pR, qL])
          helper.vals_i[i,                 j + numDofPerNode] = real(jacLR[i, j, pL, qR])
          helper.vals_i[i + numDofPerNode, j + numDofPerNode] = real(jacRR[i, j, pR, qR])
        end
      end


      MatSetValuesBlocked(helper.A, helper.idx_ib, helper.idy_ib, helper.vals_i, ADD_VALUES)
#      set_values1!(helper.A, helper.idx_i, helper.idy_i, helper.vals_i, ADD_VALUES)
    end  # end loop q
  end  # end loop p

  return nothing
end


function assembleInterface(helper::_AssembleElementData{PetscMat}, 
                           sbpface::SparseFace,
                           mesh::AbstractMesh, iface::Interface,
                           jacLL::AbstractArray{T, 4},
                           jacLR::AbstractArray{T, 4},
                           jacRL::AbstractArray{T, 4},
                           jacRR::AbstractArray{T, 4}) where T

  numNodesPerElement = size(jacLL, 4)
  numDofPerNode = size(jacLL, 1)

  permL = sview(sbpface.perm, :, iface.faceL)
  permR = sview(sbpface.perm, :, iface.faceR)

  for p=1:sbpface.numnodes
    iR = sbpface.nbrperm[p, iface.orient]
    pL = sbpface.perm[p, iface.faceL]
    pR = sbpface.perm[iR, iface.faceR]

    qL = pL
    qR = pR

    dof1 = mesh.dofs[1, pL, iface.elementL] + mesh.dof_offset
    helper.idx_ib[1] = div(dof1 - 1, numDofPerNode) + 1 - 1
    dof1 = mesh.dofs[1, pR, iface.elementR] + mesh.dof_offset
    helper.idx_ib[2] = div(dof1 - 1, numDofPerNode) + 1 - 1

    # compute block indices for q (zero-based)
    dof1 = mesh.dofs[1, qL, iface.elementL] + mesh.dof_offset
    helper.idy_ib[1] = div(dof1 - 1, numDofPerNode) + 1 - 1
    dof1 = mesh.dofs[1, qR, iface.elementR] + mesh.dof_offset
    helper.idy_ib[2] = div(dof1 - 1, numDofPerNode) + 1 - 1


    # put values into 2 x 2 block matrix
    @simd for j=1:numDofPerNode
      @simd for i=1:numDofPerNode
        helper.vals_i[i,                 j]                 = real(jacLL[i, j, pL, qL])
        helper.vals_i[i + numDofPerNode, j]                 = real(jacRL[i, j, pR, qL])
        helper.vals_i[i,                 j + numDofPerNode] = real(jacLR[i, j, pL, qR])
        helper.vals_i[i + numDofPerNode, j + numDofPerNode] = real(jacRR[i, j, pR, qR])
      end
    end


    MatSetValuesBlocked(helper.A, helper.idx_ib, helper.idy_ib, helper.vals_i, ADD_VALUES)
#      set_values1!(helper.A, helper.idx_i, helper.idy_i, helper.vals_i, ADD_VALUES)
  end  # end loop p

  return nothing
end

function assembleInterfaceVisc(helper::_AssembleElementData{PetscMat}, 
                           sbpface::AbstractFace,
                           mesh::AbstractMesh, iface::Interface,
                           jacLL::AbstractArray{T, 4},
                           jacLR::AbstractArray{T, 4},
                           jacRL::AbstractArray{T, 4},
                           jacRR::AbstractArray{T, 4}) where T

  elL = iface.elementL; elR = iface.elementR
  stencilL = sview(sbpface.perm, :, iface.faceL)
  stencilR = sview(sbpface.perm, :, iface.faceR)
  nonstencilL = sview(helper.nonstencil_nodes, :, iface.faceL)
  nonstencilR = sview(helper.nonstencil_nodes, :, iface.faceR)
  numDofPerNode = mesh.numDofPerNode

  idx = helper.idx_ib
  idy = helper.idy_ib
  vals = helper.vals_i

  # do the square blocks

  # LL and RL
  for _q=1:length(stencilL)
    qL = stencilL[_q]
    qR = stencilR[_q]

    # block indices (zero-based)
    dof1 = mesh.dofs[1, qL, elL]
    idy[1] = div(dof1 - 1, numDofPerNode) + 1 - 1
    dof1 = mesh.dofs[1, qR, elR]
    idy[2] = div(dof1 - 1, numDofPerNode) + 1 - 1

    for _p =1:length(stencilL)  # assume stencilL and R are the same length
      pL = stencilL[_p]
      pR = stencilR[_p]

      # block indices (zero-based)
      dof1 = mesh.dofs[1, pL, elL]
      idx[1] = div(dof1 - 1, numDofPerNode) + 1 - 1
      dof1 = mesh.dofs[1, pR, elR]
      idx[2] = div(dof1 - 1, numDofPerNode) + 1 - 1

      @simd for j=1:mesh.numDofPerNode
        @simd for i=1:mesh.numDofPerNode
          vals[i,                      j                     ] = real(jacLL[i, j, pL, qL])
          vals[i + mesh.numDofPerNode, j                     ] = real(jacRL[i, j, pR, qL])
          vals[i,                      j + mesh.numDofPerNode] = real(jacLR[i, j, pL, qR])
          vals[i + mesh.numDofPerNode, j + mesh.numDofPerNode] = real(jacRR[i, j, pR, qR])
        end
      end  # end j

      MatSetValuesBlocked(helper.A, idx, idy, vals, ADD_VALUES)
    end  # end p
  end  # end q

  # do the remaining rectangles

  # columns
  for _q = 1:length(stencilL)
    qL = stencilL[_q]
    qR = stencilR[_q]

    # block indices (zero-based)
    dof1 = mesh.dofs[1, qL, elL]
    idy[1] = div(dof1 - 1, numDofPerNode) + 1 - 1
    dof1 = mesh.dofs[1, qR, elR]
    idy[2] = div(dof1 - 1, numDofPerNode) + 1 - 1


    for _p=1:length(nonstencilL)  # assume nonstencilL and R are the same length
      pL = nonstencilL[_p]
      pR = nonstencilR[_p]

      # block indices (zero-based)
      dof1 = mesh.dofs[1, pL, elL]
      idx[1] = div(dof1 - 1, numDofPerNode) + 1 - 1
      dof1 = mesh.dofs[1, pR, elR]
      idx[2] = div(dof1 - 1, numDofPerNode) + 1 - 1

      @simd for j=1:mesh.numDofPerNode
        @simd for i=1:mesh.numDofPerNode
          vals[i,                      j                     ] = real(jacLL[i, j, pL, qL])
          vals[i + mesh.numDofPerNode, j                     ] = real(jacRL[i, j, pR, qL])
          vals[i,                      j + mesh.numDofPerNode] = real(jacLR[i, j, pL, qR])
          vals[i + mesh.numDofPerNode, j + mesh.numDofPerNode] = real(jacRR[i, j, pR, qR])
        end
      end  # end j

      MatSetValuesBlocked(helper.A, idx, idy, vals, ADD_VALUES)
    end  # end p
  end  # end q

  # rows
  for _p = 1:length(stencilL)
    pL = stencilL[_p]
    pR = stencilR[_p]

    # block indices (zero-based)
    dof1 = mesh.dofs[1, pL, elL]
    idx[1] = div(dof1 - 1, numDofPerNode) + 1 - 1
    dof1 = mesh.dofs[1, pR, elR]
    idx[2] = div(dof1 - 1, numDofPerNode) + 1 - 1


    for _q=1:length(nonstencilL)
      qL = nonstencilL[_q]
      qR = nonstencilR[_q]

      # block indices (zero-based)
      dof1 = mesh.dofs[1, qL, elL]
      idy[1] = div(dof1 - 1, numDofPerNode) + 1 - 1
      dof1 = mesh.dofs[1, qR, elR]
      idy[2] = div(dof1 - 1, numDofPerNode) + 1 - 1

      @simd for j=1:mesh.numDofPerNode
        @simd for i=1:mesh.numDofPerNode
          vals[i,                      j                     ] = real(jacLL[i, j, pL, qL])
          vals[i + mesh.numDofPerNode, j                     ] = real(jacRL[i, j, pR, qL])
          vals[i,                      j + mesh.numDofPerNode] = real(jacLR[i, j, pL, qR])
          vals[i + mesh.numDofPerNode, j + mesh.numDofPerNode] = real(jacRR[i, j, pR, qR])
        end
      end  # end j

      MatSetValuesBlocked(helper.A, idx, idy, vals, ADD_VALUES)
    end  # end q
  end  # end p

  # remaining on diagonal entries
  # This is only needed for shock capturing (elements where dLambda/dq != 0)
  idx = helper.idx_bb
  idy = helper.idy_bb
  vals = helper.vals

  for _q=1:length(nonstencilL)
    qL = nonstencilL[_q]
    qR = nonstencilR[_q]

    # block indices (zero-based)
    dof1 = mesh.dofs[1, qL, elL]
    idy[1] = div(dof1 - 1, numDofPerNode) + 1 - 1

    for _p=1:length(nonstencilL)
      pL = nonstencilL[_p]

      # block indices (zero-based)
      dof1 = mesh.dofs[1, pL, elL]
      idx[1] = div(dof1 - 1, numDofPerNode) + 1 - 1

      @simd for j=1:mesh.numDofPerNode
        @simd for i=1:mesh.numDofPerNode
          vals[i, j] = real(jacLL[i, j, pL, qL])
        end
      end

      MatSetValuesBlocked(helper.A, idx, idy, vals, ADD_VALUES)
    end
    
    dof1 = mesh.dofs[1, qR, elR]
    idy[1] = div(dof1 - 1, numDofPerNode) + 1 - 1

    for _p=1:length(nonstencilR)
      pR = nonstencilR[_p]
      dof1 = mesh.dofs[1, pR, elR]
      idx[1] = div(dof1 - 1, numDofPerNode) + 1 - 1

      @simd for j=1:mesh.numDofPerNode
        @simd for i=1:mesh.numDofPerNode
          vals[i, j] = real(jacRR[i, j, pR, qR])
        end
      end

      MatSetValuesBlocked(helper.A, idx, idy, vals, ADD_VALUES)
    end  # end p
    
  end  # end q

  return nothing
end


"""
  Assemble one half of an interface, used by shared face integrals.
  See [`assembleInterface`](@ref).

  **Inputs**

   * helper: _AssembleElementData
   * sbpface: an SBP face object
   * mesh: a mesh
   * jacLL
   * jacLR
"""
function assembleSharedFace(helper::_AssembleElementData{PetscMat}, sbpface::DenseFace,
                            mesh::AbstractMesh,
                            iface::Interface,
                            jacLL::AbstractArray{T, 4},
                            jacLR::AbstractArray{T, 4}) where T

  numNodesPerElement = size(jacLL, 4)
  numDofPerNode = size(jacLL, 1)

#  for q in permL  # =1:numNodesPerElement
  for q=1:sbpface.stencilsize
    qL = sbpface.perm[q, iface.faceL]
    qR = sbpface.perm[q, iface.faceR]

    # compute block indices for q (zero-based)
    dof1 = mesh.dofs[1, qL, iface.elementL] + mesh.dof_offset
    helper.idy_ib[1] = div(dof1 - 1, numDofPerNode) + 1 - 1
    dof1 = mesh.dofs[1, qR, iface.elementR] + mesh.dof_offset
    helper.idy_ib[2] = div(dof1 - 1, numDofPerNode) + 1 - 1


#    for p in permR # =1:numNodesPerElement
    for p=1:sbpface.stencilsize
      pL = sbpface.perm[p, iface.faceL]

      dof1 = mesh.dofs[1, pL, iface.elementL] + mesh.dof_offset
      helper.idx_bb[1] = div(dof1 - 1, numDofPerNode) + 1 - 1

      # put values into 2 x 2 block matrix
      for j=1:numDofPerNode
        for i=1:numDofPerNode
          helper.vals_sf[i, j]                 = real(jacLL[i, j, pL, qL])
          helper.vals_sf[i, j + numDofPerNode] = real(jacLR[i, j, pL, qR])
        end
      end

      MatSetValuesBlocked(helper.A, helper.idx_bb, helper.idy_ib, helper.vals_sf, ADD_VALUES)
#      set_values1!(helper.A, helper.idx, helper.idy_i, helper.vals_sf, ADD_VALUES)
    end  # end loop q
  end  # end loop p

  return nothing
end


function assembleSharedFace(helper::_AssembleElementData{PetscMat}, sbpface::SparseFace,
                            mesh::AbstractMesh,
                            iface::Interface,
                            jacLL::AbstractArray{T, 4},
                            jacLR::AbstractArray{T, 4}) where T

  numNodesPerElement = size(jacLL, 4)
  numDofPerNode = size(jacLL, 1)

  for p=1:sbpface.numnodes
    iR = sbpface.nbrperm[p, iface.orient]
    qL = sbpface.perm[p, iface.faceL]
    qR = sbpface.perm[iR, iface.faceR]

    pL = qL

    # compute block indices (zero-based)
    dof1 = mesh.dofs[1, pL, iface.elementL] + mesh.dof_offset
    helper.idx_bb[1] = div(dof1 - 1, numDofPerNode) + 1 - 1

    dof1 = mesh.dofs[1, qL, iface.elementL] + mesh.dof_offset
    helper.idy_ib[1] = div(dof1 - 1, numDofPerNode) + 1 - 1
    dof1 = mesh.dofs[1, qR, iface.elementR] + mesh.dof_offset
    helper.idy_ib[2] = div(dof1 - 1, numDofPerNode) + 1 - 1


    # put values into 2 x 2 block matrix
    for j=1:numDofPerNode
      for i=1:numDofPerNode
        helper.vals_sf[i, j]                 = real(jacLL[i, j, pL, qL])
        helper.vals_sf[i, j + numDofPerNode] = real(jacLR[i, j, pL, qR])
      end
    end

    MatSetValuesBlocked(helper.A, helper.idx_bb, helper.idy_ib, helper.vals_sf, ADD_VALUES)
#      set_values1!(helper.A, helper.idx, helper.idy_i, helper.vals_sf, ADD_VALUES)
  end  # end loop p

  return nothing
end



"""
  Assembles the jacobian of a boundary integral into the matrix.
  Specialized versions take advantage of the sparsity of the `sbpface`, ie.
  it only assembles nodes that are in the stencil of R.  See
  [`AssembleBoundaryFull`](@ref) for assembling the entire block.

  **Inputs**

   * helper: _AssembleElementData
   * sbpface: an SBP face object
   * mesh: a mesh
   * jac: a jac, same layout as [`assembleElement`](@ref)
"""
function assembleBoundary(helper::_AssembleElementData, sbpface::DenseFace,
                            mesh::AbstractMesh,
                            bndry::Boundary,
                            jac::AbstractArray{T, 4}) where T
# MatSetValuesBlocked not implemented for boundaries because the performance
# gain isn't that great

  elnum = bndry.element
  numNodesPerElement = size(jac, 4)
  numDofPerNode = size(jac, 1)

#  for q=1:numNodesPerElement
#    qL = q
  for q=1:sbpface.stencilsize
    qL = sbpface.perm[q, bndry.face]

    # get dofs for node q
    for j=1:numDofPerNode
      helper.idy[j] = mesh.dofs[j, qL, elnum] + mesh.dof_offset
    end

#    for p=1:numNodesPerElement
#      pL = p
    for p=1:sbpface.stencilsize
      pL = sbpface.perm[p, bndry.face]

      # get dofs for node p
      for i=1:numDofPerNode
        helper.idx[i] = mesh.dofs[i, pL, elnum] + mesh.dof_offset
      end

      # get values
      for j=1:numDofPerNode
        for i=1:numDofPerNode
          helper.vals[i, j] = real(jac[i, j, pL, qL])
        end
      end

      # assemble them into matrix
      set_values1!(helper.A, helper.idx, helper.idy, helper.vals, ADD_VALUES)
    end  # end loop p
  end  # end loop q

  return nothing
end

function assembleBoundary(helper::_AssembleElementData, sbpface::SparseFace,
                            mesh::AbstractMesh,
                            bndry::Boundary,
                            jac::AbstractArray{T, 4}) where T

  elnum = bndry.element
  numNodesPerElement = size(jac, 4)
  numDofPerNode = size(jac, 1)

  for p=1:sbpface.numnodes
    pL = sbpface.perm[p, bndry.face]
    qL = pL

    # get dofs for node p
    for i=1:numDofPerNode
      helper.idx[i] = mesh.dofs[i, pL, elnum] + mesh.dof_offset
      helper.idy[i] = mesh.dofs[i, qL, elnum] + mesh.dof_offset
    end

    # get values
    for j=1:numDofPerNode
      for i=1:numDofPerNode
        helper.vals[i, j] = real(jac[i, j, pL, qL])
      end
    end

    # assemble them into matrix
    set_values1!(helper.A, helper.idx, helper.idy, helper.vals, ADD_VALUES)

  end  # end loop p

  return nothing
end


#------------------------------------------------------------------------------
# Generic methods (internally they call one of the AssembleElementData
# implementation-specific functions

"""
  This function assembles the entire element Jacobian into the sparse matrix
  for a given boundary.

  **Inputs**

   * helper: [`AssembleElementData`](@ref)
   * mesh
   * bndry: a `Boundary` object
   * jac: the element jacobian
"""
function assembleBoundaryFull(helper::AssembleElementData,
                              mesh::AbstractMesh, bndry::Boundary,
                              jac::AbstractArray{T, 4}) where T

  assembleElement(helper, mesh, bndry.element, jac)

  return nothing
end


"""
  This function assembles all 4 element Jacobians into the sparse matrix for
  a given interface (rather than only those in the stencil of R).

  **Inputs**

   * helper: an [`AssembleElementData`](@ref)
   * mesh
   * iface: an `Interface` object
   * jacLL
   * jacLR
   * jacRL
   * jacRR
"""
function assembleInterfaceFull(helper::AssembleElementData, 
                           mesh::AbstractMesh, iface::Interface,
                           jacLL::AbstractArray{T, 4},
                           jacLR::AbstractArray{T, 4},
                           jacRL::AbstractArray{T, 4},
                           jacRR::AbstractArray{T, 4}) where T

  fullface = FullFace(mesh.numNodesPerElement, mesh.dim)
  assembleInterface(helper, fullface, mesh, iface, jacLL, jacLR,
                                                          jacRL, jacRR)

  return nothing
end

"""
  Similar to [`assembleInterfaceFull`](@ref), but for shared faces.
"""
function assembleSharedFaceFull(helper::AssembleElementData,
                            mesh::AbstractMesh,
                            iface::Interface,
                            jacLL::AbstractArray{T, 4},
                            jacLR::AbstractArray{T, 4}) where T

  fullface = FullFace(mesh.numNodesPerElement, mesh.dim)
  assembleSharedFace(helper, fullface, mesh, iface, jacLL, jacLR)

  return nothing
end

#------------------------------------------------------------------------------
# Internal functions

"""
  Get the list of nodes that are not in the stencil of the interpolation
  operator R for a given face.

  **Inputs**

   * sbpface: `AbstractFace`
   * face: the local face number
   * numNodesPerElement: the total number of nodes in the element

  **Inputs/Outputs**

   * nodes: an array to be populated with the non-stencil nodes.  This array
            may be longer than required.  The first `n` indices will be
            populated.

  **Outputs**

   * n: the number of non-stencil nodes
"""
function getNonStencilNodes(sbpface::AbstractFace, face::Integer,
                            numNodesPerElement::Integer, nodes::AbstractVector)

  stencil_nodes = sview(sbpface.perm, :, face)
  idx = 1
  for i=1:numNodesPerElement
    if !(i in stencil_nodes)
      nodes[idx] = i
      idx += 1
    end
  end

  return idx - 1
end


"""
  Construct a table of the non-stencil nodes for all faces of the element.

  **Inputs**

   * sbp: `AbstractOperator` or the number of nodes per element
   * sbpface: `AbstractFace`

  **Outputs**

   * table: num_nonstencil x number of faces per element.
"""
function getNonStencilNodes(numNodesPerElement::Integer, sbpface::AbstractFace)

  num_nonstencil = numNodesPerElement - size(sbpface.perm, 1)
  table = Array{Int}(num_nonstencil, size(sbpface.perm, 2))
  for i=1:size(sbpface.perm, 2)
    nodes_i = sview(table, :, i)
    n = getNonStencilNodes(sbpface, i, numNodesPerElement, nodes_i)
    @assert n == num_nonstencil
  end

  return table
end


function getNonStencilNodes(sbp::AbstractOperator, sbpface::AbstractFace)

  return getNonStencilNodes(sbp.numnodes, sbpface)
end
