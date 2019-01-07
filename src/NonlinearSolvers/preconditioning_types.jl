# types used for preconditioning
# The type definitions and the functions that use them need to be in separate
# files because they need to be included in a particular order

#TODO: it should be possible to move this into preconditioning.jl again

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
   
  **Fields**

   * pc_inner: a [`PetscMatFreePC`](@ref) 
   * assem: an [`AssembleDiagJacData`](@ref)
   * diag_jac: a [`DiagJac`]
   * ipiv: permutation information, `bs` x `numEl`
   * bs: the block size, `mesh.numDofPerNode*mesh.numNodesPerElement`
   * is_factored: if true, volume_jac has been factored (LU with partial
                  pivoting), false otherwise
"""
mutable struct NewtonBDiagPC <: AbstractPetscMatFreePC
  pc_inner::PetscMatFreePC
  assem::AssembleDiagJacData{Float64}
  diag_jac::DiagJac{Float64}  # same DiagJac used in assem
  ipiv::Array{BlasInt, 2}
  bs::Int  # block size
  is_factored::Bool
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

  return NewtonBDiagPC(pc_inner, assem, diag_jac, ipiv, bs, is_factored)
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


