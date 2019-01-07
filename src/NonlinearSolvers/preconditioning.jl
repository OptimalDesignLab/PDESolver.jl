# functions for calculating various preconditioners
using Base.LinAlg.BLAS  # trmv!

#------------------------------------------------------------------------------
# AbstractPC Interface
function calcPC(pc::NewtonBDiagPC, mesh::AbstractMesh, sbp::AbstractOperator,
                eqn::AbstractSolutionData, opts::Dict, ctx_residual, t)

  setPCCtx(pc, mesh, sbp, eqn, opts, ctx_residual, t)  # unnecessary?

  evalJacobian(mesh, sbp, eqn, opts, pc.assem)
  factorBDiagPC(pc, mesh, sbp, eqn, opts)
end

function applyPC(pc::NewtonBDiagPC, mesh::AbstractMesh, sbp::AbstractOperator,
                 eqn::AbstractSolutionData, opts::Dict, t, b::AbstractVector, 
                 x::AbstractVector)

  applyBDiagPC(pc, mesh, sbp, eqn, opts, b, x)

  return nothing
end

# we could define applyPCTranspose, but this PC is ineffective so don't bother

#------------------------------------------------------------------------------
# Volume preconditioner implementation


"""
  This function factors the volume preconditioner.  It uses Lapack to do an
  LU factorization with partial pivoting.  The user is allowed to modify
  `pc.diag_jac`_jac after calling [`calcBDiagPC`](@ref)
  but before calling this routine.

  **Inputs**
  
   * mesh
   * sbp
   * eqn
   * opts

  **Inputs/Outputs**

   * newton_data: vol_prec is updated. `pc.diag_jac` is factored in-place and
                 `pc.ipiv` is overwritten with the permutation

"""
function factorBDiagPC(pc::NewtonBDiagPC, mesh, sbp, eqn, opts)



  for i=1:mesh.numEl
    jac_i = sview(pc.diag_jac.A, :, :, i)
    ipiv_i = sview(pc.ipiv, :, i)

    # call Lapack here
    info = getrf!(jac_i, ipiv_i)
    @assert info == 0
  end

  pc.is_factored = true

  return nothing
end

"""
  Apply the preconditioner, ie. do inv(A)*x = b, where A is calculated by
  [`calcBDiagPC`](@ref).

  This is not likely to be a sensible preconditioner for Contiuous Galerkin
  discretizations.

  **Inputs**

   * pc: a [`BDiagPC `](@ref) object, already factored by [`factorBDiagPC`](@ref).
   * mesh
   * sbp
   * eqn
   * opts
   * x: the vector to multiply against

  **Inputs/Outputs**

   * b: the output vector (overwritten)

"""
function applyBDiagPC(pc::NewtonBDiagPC, mesh, sbp, eqn, opts, x::AbstractVector, b::AbstractVector)

  # we need to do inv(A)*x = b --> solve A*b = x using the factorization A.

  @assert pc.is_factored

  bs = pc.bs
  workvec = zeros(Float64, bs)  # hold values passed into LAPACK

  for i=1:mesh.numEl
    jacf_i = sview(pc.diag_jac.A, :, :, i)
    ipiv_i = sview(pc.ipiv, :, i)

    getValues(mesh, x, i, workvec)

    # call Lapack GETRS to solve for b (in-place)
    getrs2!('N', jacf_i, ipiv_i, workvec)

    setValues(mesh, workvec, i, b)

  end  # end loop i

  return nothing
end





"""
  Computes the action of the inverse of the preconditioner on a vector.
  This is not required by the AbstractPC interface, but is useful for
  constructing smoothers.

  Note that the preconditioner is defined as the inverse of some
  some approximation to the Jacobian. so the inverse of the preconditioner
  is the approximate Jacobian itself.

  **Inputs**

   * pc: [`NewtonBDiagPC`](@ref).  Can be factored or not
   * mesh
   * sbp
   * eqn
   * opts
   * x: vector to multiply against

  **Inputs/Outputs**

   * b: vector to overwrite with the result.
"""
function applyBDiagPCInv(pc::NewtonBDiagPC, mesh, sbp, eqn, opts, x::AbstractVector, b::AbstractVector)


  if !pc.is_factored
    diagMatVec(pc.diag_mat, mesh, x, b)
  else

    bs = pc.bs
    workvec = zeros(Float64, bs)  # hold values passed into LAPACK

    for i=1:mesh.numEl
      jacf_i = sview(pc.diag_jac.A, :, :, i)
      ipiv_i = sview(pc.ipiv, :, i)

      getValues(mesh, x, i, workvec)

      # do b = P*L*U*x
      trmv!('U', 'N', 'N', jacf_i, workvec)
      trmv!('L', 'N', 'U', jacf_i, workvec)
      applyIpiv!(ipiv_i, workvec)

      setValues(mesh, workvec, i, b)

    end  # end loop i

  end  # end if

  return nothing
end



