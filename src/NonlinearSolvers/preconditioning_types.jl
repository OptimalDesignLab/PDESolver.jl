# types used for preconditioning
# The type definitions and the functions that use them need to be in separate
# files because they need to be included in a particular order

#TODO: it should be possible to move this into preconditioning.jl again


#------------------------------------------------------------------------------
# NewtonBDiagPC

"""
  This type holds all the data needed to calculate a preconditioner based only
  on the element-block diagonal of the Jacobian (the
  `mesh.numDofPerNode*mesh.numNodesPerElement` square block).
  This preconditioner is classified as matrix-free, although it does store the
  Jacobian diagonal.  This preconditioner only works for physics that have
  explicit Jacobian calculation.


  This preconditioner is not suitable for use as an inner preconditioner for
  other methods in NonlinearSolvers, but
  the functions [`calcBDiagPC`](@ref), [`factorBDiagPC`](@ref), and
  [`applyBDiagPC`](@ref) can be easily used to make a new PC.

  TODO: this doesn't work when implicit Euler globalization is used, because
        it calls `evalJacobian` not `physicsJac`.  `physicsJac` needs to
        be more like `evalJacobian`.

  **Fields**

   * pc_inner: a [`PetscMatFreePC`](@ref) 
   * assem: an [`AssembleDiagJacData`](@ref)
   * diag_jac: a [`DiagJac`]
   * ipiv: permutation information, `bs` x `numEl`
   * bs: the block size, `mesh.numDofPerNode*mesh.numNodesPerElement`
   * is_factored: if true, volume_jac has been factored (LU with partial
                  pivoting), false otherwise
   * evalJacobian: function to evaluate the Jacobian, must have same
                   signature as `PDESolver.evalJacobian`.
"""
mutable struct NewtonBDiagPC <: NewtonMatFreePC
  pc_inner::PetscMatFreePC
  assem::AssembleDiagJacData{Float64}
  diag_jac::DiagJac{Float64}  # same DiagJac used in assem
  ipiv::Array{BlasInt, 2}
  bs::Int  # block size
  is_factored::Bool
  evalJacobian::Function
end  # end type definition

function needParallelData(pc::NewtonBDiagPC)
  return true
end

"""
  Regular constructor

  **Inputs**

   * mesh: the mesh
   * sbp: the SBP operator
   * eqn: AbstractSolutionData
   * opts: options dictionary
"""
function NewtonBDiagPC(mesh::AbstractMesh, sbp, eqn::AbstractSolutionData, opts)


  pc_inner = PetscMatFreePC(mesh, sbp, eqn, opts)
  bs = mesh.numDofPerNode*mesh.numNodesPerElement
  numEl = mesh.numEl

  diag_jac = DiagJac(Float64, bs, numEl)
  assem = AssembleDiagJacData(mesh, sbp, eqn, opts, diag_jac)
  ipiv = Array{BlasInt}(bs, numEl)
  is_factored = false
  eval_jacobian = evalJacobian

  return NewtonBDiagPC(pc_inner, assem, diag_jac, ipiv, bs, is_factored,
                       eval_jacobian)
end

"""
  Default constructor (think synthesized default constructor in C++).
  Returns an object where all arrays have size zero

  **Inputs**

    none
"""
function NewtonBDiagPC()

  pc_inner = PetscMatFreePC(mesh, sbp, eqn, opts)
  bs = 0
  diag_jac = DiagJac(Float64, bs, numEl)
  assem = AssembleDiagJacData(mesh, sbp, eqn, opts, diag_jac)
  ipiv = Array{BlasInt}(0, 0)
  is_factored = false

  return NewtonBDiagPC(pc_inner, bs, diag_jac, assem, ipiv, is_factored)
end

"""
  Function to set a new `evalJacobian` function.

  **Inputs**

   * pc: the `NewtonBDiagPC
   * func: the new function
"""
function setEvalJacobian(pc::NewtonBDiagPC, func::Function)

  pc.evalJacobian = func
end


#------------------------------------------------------------------------------
# NewtonBJacobiPC

"""
  A block-jacobi smoother.  Uses the element diagonal
  (`mesh.numDofPerNode*mesh.numNodesPerElement`) rather than the strict
  diagonal.

  This preconditioner requires the solution and residual to be complex
  numbers.

  This is a stationary smoother, it does not update eqn.q_vec, and can be
  used to smooth any quantity.

  TODO: this function doesn't work when implicit euler globalization is
        used, similar reasons as `NewtonBDIagPC`.
  **Fields**

   * pc_inner: `PetscMatFreePC`
   * diag_pc: a [`NewtonBDiagPC`](@ref)
   * itermax: maximum number of iterations, default 10
   * res_tol: tolerance for Euclidian norm of linear residual, default 1e-4
   * verbose: Bool determining if output is done, default false
   * evalJacVecProduct: function to evaluate Jacobian-vector products,
                        must have same signature as
                        `PDESolver.evaldRdqProduct`
   * evalJacTVecProduct: function to evaluate trasposed Jacobian-vector
                         product.  Must have same signature as
                         `PDESolver.evalResidual_revq`
"""
mutable struct NewtonBJacobiPC <: NewtonMatFreePC
  pc_inner::PetscMatFreePC
  diag_pc::NewtonBDiagPC
  itermax::Int
  res_tol::Float64
  verbose::Bool
  evalJacVecProduct::Function
  evalJacTVecProduct::Function
end


function NewtonBJacobiPC(mesh::AbstractMesh, sbp, eqn::AbstractSolutionData,
                         opts; itermax=10, res_tol=1e-4, verbose=false)

  pc_inner = PetscMatFreePC(mesh, sbp, eqn, opts)
  diag_pc = NewtonBDiagPC(mesh, sbp, eqn, opts)
  evalJacVecProduct = evaldRdqProduct
  evalJacTVecProduct = evalResidual_revq
  return NewtonBJacobiPC(pc_inner, diag_pc, itermax, res_tol, verbose,
                         evalJacVecProduct, evalJacTVecProduct)
end


function free(obj::NewtonBJacobiPC)

  # remove any references to the diagonal jacobian so the GC can free it
  etype = eltype(obj.diag_pc.diag_jac)
  obj.diag_pc.diag_jac = NullDiagJac
  obj.diag_pc.assem = NullAssembleDiagJacData(etype)
end


"""
  Sets the functions that compute the Jacobian and matrix-free products.

  **Inputs**

   * pc: [`NewtonBJacobiPC`](@ref)
   * eval_jac: function that evaluates the Jacobian (see `NewtonBDiagPC`)
   * eval_jacT_vec: evaluates matrix-free Jacobian-vector products
   * eval_jacT_vec: evaluates transposed Jacobian-vector products,
                    this argument is only required if the transposed
                    preconditioner will be applied.
"""
function setEvalJacobian(pc::NewtonBJacobiPC, eval_jac::Function,
                         eval_jac_vec::Function,
                         eval_jacT_vec::Function=evalResidual_revq)

  setEvalJacobian(pc.diag_pc, eval_jac)
  pc.evalJacVecProduct = eval_jac_vec
  pc.evalJacTVecProduct = eval_jacT_vec

  return nothing
end
