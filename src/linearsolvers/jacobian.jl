"""
  This function computes the Jacobian of an evalResidual-like function.
  Specifically, it computes \\partial (eqn.q_vec)/ \\partial (eqn.res_vec).

  Other users of [`newtonInner`](@ref) (for example, implicit time marching
  methods) will need to implemnt their own version of this function.  The
  documentation for this function describes the requirements for the other
  implementations.

  [`newtonInner](@ref) guarantees that eqn.q and eqn.q_vec will be consistent
  when this function is called.  No guarantees are made about the contents
  of eqn.res or eqn.res_vec.

  
  **Inputs**:
   
   * lo: an [`AbstractLinearOperator`](@ref) (a matrix-explicit one)
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

  This function should not have to do any parallel communication. newtonInner
  ensures that the rhs_func is called before jac_func, and rhs_func handles
  the parallel communication.

  Implementations of this function may perform either (eqn.q -> jacobian) or
  (eqn.q_vec -> jacobian).  The first may be more computationally efficient,
  but the second can be simpler for some time-marching methods.

  This function supportes several types of jacobians (dense arrays,
  SparseMatrixCSC, PetscMat), and several methods for calculating them
  (finite difference and complex step).  All implementations of this function
  should support them as well.  For this reason, it is strongly recommneded to
  use this function to compute the spatial jacobian and them modify the
  resulting matrix.

  When using Petsc matrices, the function may do intermediate assemblies
  (PETSC_FLUSH_ASSEMBLY), but does not need to do the final assembly.

"""
function physicsJac(newton_data::NewtonData, mesh, sbp, eqn, opts, jac, ctx_residual, t=0.0; is_preconditioned::Bool=false)

  verbose = opts["newton_verbosity"]::Int

  myrank = mesh.myrank

  #TODO: figure out which of these are actually needed
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
      tmp, t_jac, t_gc, alloc = @time_all calcJacFD(newton_data, mesh, sbp, eqn, opts, func, res_copy_vec, pert, jac, t)

    elseif jac_type == 2  # Julia sparse jacobian
      @verbose5 @mpi_master println(BSTDOUT, "calculating sparse FD jacobian")
      #TODO: don't copy the giant array!
      res_copy = copy(eqn.res)  # copy unperturbed residual

      tmp, t_jac, t_gc, alloc = @time_all calcJacobianSparse(newton_data, mesh, sbp, eqn, opts, func, res_copy, pert, jac, t)
    elseif jac_type == 3  # Petsc sparse jacobian
      @verbose5 @mpi_master println(BSTDOUT, "calculating sparse FD jacobian")
      res_copy = copy(eqn.res)  # copy unperturbed residual
      tmp, t_jac, t_gc, alloc = @time_all calcJacobianSparse(newton_data, mesh, sbp, eqn, opts, func, res_copy, pert, jac, t)
    elseif jac_type == 4  # Petsc jacobian-vector product
      throw(ErrorException("No handling of jac_method = 1 and jac_type = 4: 
                           finite differencing isn't permitted for Petsc mat-free"))
    end

  elseif jac_method == 2
    @verbose5 @mpi_master println(BSTDOUT, "calculating complex step jacobian")

    if jac_type == 1  # dense jacobian
      @verbose5 @mpi_master println(BSTDOUT, "calculating dense complex step jacobian")
      tmp, t_jac, t_gc, alloc = @time_all calcJacobianComplex(newton_data, mesh, sbp, eqn, opts, func, pert, jac, t)
    elseif jac_type == 2  # Julia sparse jacobian 
      @verbose5 @mpi_master println(BSTDOUT, "calculating sparse complex step jacobian")
      res_dummy = Array(Float64, 0, 0, 0)  # not used, so don't allocation memory
      tmp, t_jac, t_gc, alloc = @time_all calcJacobianSparse(newton_data, mesh, sbp, eqn, opts, func, res_dummy, pert, jac, t)
    elseif jac_type == 3 # Petsc sparse jacobian
      res_dummy = Array(Float64, 0, 0, 0)  # not used, so don't allocation memory
      @verbose5 @mpi_master println(BSTDOUT, "calculating explicit Petsc jacobian")

      @verbose5 @mpi_master println(BSTDOUT, "calculating main jacobain")
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
  @verbose5 if print_jacobian_timing
    @mpi_master print(BSTDOUT, "jacobian calculation: ")
    @mpi_master print_time_all(BSTDOUT, t_jac, t_gc, alloc)
  end

  flush(BSTDOUT)

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

  The implict Euler globalization is supported in matrix-free mode.

  **Inputs**:

   * newton_data:  NewtonData object
   * mesh: AbstractMesh
   * sbp:  SBP operator
   * eqn:  AbstractEquation object
   * opts: options dictionary
   * pert: perturbation to use for the algorithmic differentiation.  Currently,
   *       only complex numbers are supported.
   * rhs_func: rhs_func from [`newtonInner`](@ref)
   * ctx_residual: ctx_residual from [`newtonInner`](@ref)
   * vec:  the x vector in Ax=b.  Can be AbstractVector type.
   * b:    location to store the result (the b an Ax=b).  Can be any
           AbstractVector type

  Outputs:
    none


  Aliasing restrictions: vec, b, and eqn.q must not alias each other.

"""
function calcJacVecProd(newton_data::NewtonData, mesh, sbp, eqn, opts, pert,
                        rhs_func::Function, ctx_residual, vec::AbstractVector,
                        b::AbstractVector, t=0.0)
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
 
#  itr = newton_data.krylov_itr
  globalize_euler = opts["newton_globalize_euler"]::Bool

  epsilon = imag(pert)  # magnitude of perturbationa

  # apply perturbation
  for i=1:mesh.numDof
    eqn.q_vec[i] += pert*vec[i]
  end

  rhs_func(mesh, sbp, eqn, opts, eqn.res_vec, ctx_residual, t)
#=
  # scatter into eqn.q
  disassembleSolution(mesh, sbp, eqn, opts, eqn.q_vec) 
  func(mesh, sbp, eqn, opts, t)

  # gather into eqn.res_vec
  assembleResidual(mesh, sbp, eqn, opts, eqn.res_vec, assemble_edgeres=opts["use_edge_res"])
=#  
  # calculate derivatives, store into b
  calcJacCol(b, eqn.res_vec, epsilon)

  if globalize_euler
    applyEuler(mesh, sbp, eqn, opts, vec, newton_data, b)
  end

  # undo perturbation
  for i=1:mesh.numDof
    eqn.q_vec[i] -= pert*vec[i]
  end

#  eqn.params.krylov_itr += 1

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
  calcJacobianSparse(newton_data, mesh, sbp, eqn, opts, func, res_dummy, pert, jac, t)

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

  This function does not pass t correctly.

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
  ctx_ptr = MatShellGetContext(A)
  ctx_petsc = unsafe_pointer_to_objref(ctx_ptr)
  @assert length(ctx_petsc) == 8
  # unpack the tuple (could use compact syntax)
  mesh = ctx_petsc[1]
  sbp = ctx_petsc[2]
  eqn = ctx_petsc[3]
  opts = ctx_petsc[4]
  newton_data = ctx_petsc[5]
  func = ctx_petsc[6]  # rhs_func from newtonInner
  ctx_residual = ctx_petsc[7]
  t = ctx_petsc[8]

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

  calcJacVecProd(newton_data, mesh, sbp, eqn, opts, pert, func, ctx_residual, x_arr, b_arr)

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
    eqn.q (not eqn.q_vec).  A mesh coloring approach is used to compute the
    jacobian.  Both eqn.q and the MPI send and receive buffers are perturbed
    during this process.

  Inputs:
    newton_data:  NewtonData object
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
function calcJacobianSparse(newton_data::NewtonData, mesh, sbp, eqn, opts, func,
                            res_0::Abstract3DArray, pert, 
                            jac::Union{SparseMatrixCSC, PetscMat}, t=0.0)
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
  time.t_color += @elapsed for color=1:mesh.maxColors  # loop over max colors, 
                                                       # only do calculation for numColors
    for j=1:mesh.numNodesPerElement  # loop over nodes 
      for i=1:mesh.numDofPerNode  # loop over dofs on each node

        # apply perturbation to q
        if color <= mesh.numColors
          applyPerturbation(mesh, eqn.q, eqn.shared_data, color, pert, i, j, f)
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
              assembleElement(newton_data, mesh, eqn, eqn.res, res_0, k, el_pert,
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
                assembleElement(newton_data, mesh, eqn, res_edge, res_0, k, el_pert, col_idx, epsilon, jac)
              end  # end if el_pert != 0
            end  # end loop over k
          end  # end loop over local edges
        end

        PetscMatAssemblyBegin(jac, PETSC_MAT_FLUSH_ASSEMBLY)

        # undo perturbation
        if color <= mesh.numColors
          applyPerturbation(mesh, eqn.q, eqn.shared_data, color, -pert, i, j)
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
                           color::Integer, pert, i, j, f=STDOUT; 
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

   * newton_data: a NewtonData object
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
function assembleElement{Tsol <: Real}(newton_data::NewtonData, mesh,
                         eqn::AbstractSolutionData{Tsol}, res_arr, res_0,
                         el_res::Integer, el_pert::Integer, dof_pert::Integer,
                         epsilon, jac::AbstractMatrix)
#

  # resize array
  # basically a no-op if array is already the right size
  local_size = PetscInt(mesh.numNodesPerElement*mesh.numDofPerNode)

  # get row number
  newton_data.idy_tmp[1] = dof_pert + mesh.dof_offset

  pos = 1
  for j_j = 1:mesh.numNodesPerElement
    for i_i = 1:mesh.numDofPerNode
      newton_data.idx_tmp[pos] = mesh.dofs[i_i, j_j, el_res] + mesh.dof_offset
  
      tmp = (res_arr[i_i,j_j, el_res] - res_0[i_i, j_j, el_res])/epsilon
      newton_data.vals_tmp[pos] = tmp
      pos += 1
    end
  end

  set_values1!(jac, newton_data.idx_tmp, newton_data.idy_tmp, newton_data.vals_tmp, PETSC_ADD_VALUES)
  
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
function assembleElement{Tsol <: Complex}(newton_data::NewtonData, mesh,
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
  newton_data.idy_tmp[1] = dof_pert + mesh.dof_offset
  pos = 1
  for j_j = 1:mesh.numNodesPerElement
    for i_i = 1:mesh.numDofPerNode
      newton_data.idx_tmp[pos] = mesh.dofs[i_i, j_j, el_res] + mesh.dof_offset
      row = newton_data.idx_tmp[pos]
      col = newton_data.idy_tmp[1]


      newton_data.vals_tmp[pos] = imag(res_arr[i_i,j_j, el_res])/epsilon
      val = newton_data.vals_tmp[pos]
      pos += 1
    end
  end

  set_values1!(jac, newton_data.idx_tmp, newton_data.idy_tmp, newton_data.vals_tmp, PETSC_ADD_VALUES)
#  PetscMatSetValues(jac, newton_data.idx_tmp, newton_data.idy_tmp, newton_data.vals_tmp, PETSC_ADD_VALUES)

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


"""
  This function performs a matrix solve x = inv(A)*b.  A can be a dense
  matrix, a SparseMatrixCSC, or a PetscMatrix.

  This function should be used anywhere a matrix solve is needed that
  should work with any type of matrix

  Inputs:
    newton_data: a NewtonData object
    eqn: an AbstractSolutionData object (only eqn.params is really needed)
    jac: should be the matrix corresponding to the newton_data object
    b: the right hand side vector
    
  Inputs/Outputs:
    x: the output vector

  Aliasing restrictions: x and b cannot alias
"""
function matrixSolve(newton_data::NewtonData, eqn::AbstractSolutionData, 
                     mesh::AbstractDGMesh, opts,
                     jac::AbstractMatrix, x::AbstractVector, b::AbstractVector,                      fstdout::IO;
                     verbose=5)

  jac_type = typeof(jac)
  myrank = mesh.myrank
  tmp, t_solve, t_gc, alloc = @time_all if jac_type <: Array || jac_type <: SparseMatrixCSC
      jac_f = factorize(jac)
      tmp2 = jac_f\b

      copy!(x, tmp2)

  elseif  jac_type <: PetscMat
    petscSolve(newton_data, newton_data.ctx_newton..., opts, b, x, 
               mesh.dof_offset)

  else
    throw(ErrorException("Unsupported jac_type $jac_type"))
  end

  eqn.params.time.t_solve += t_solve

  @verbose5 @mpi_master print(fstdout, "matrix solve: ")
  @verbose5 @mpi_master print_time_all(fstdout, t_solve, t_gc, alloc)
  step_norm = norm(x)
  step_norm = sqrt(MPI.Allreduce(step_norm*step_norm, MPI.SUM, mesh.comm))
  @verbose5 @mpi_master println(fstdout, "step_norm = ", step_norm)
  flush(fstdout)

  return step_norm
end


