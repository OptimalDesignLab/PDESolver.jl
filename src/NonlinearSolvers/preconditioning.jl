# functions for calculating various preconditioners
using Base.LinAlg.BLAS  # trmv!

### NewtonBDiagPC
#------------------------------------------------------------------------------
# AbstractPC Interface
function calcPC(pc::NewtonBDiagPC, mesh::AbstractMesh, sbp::AbstractOperator,
                eqn::AbstractSolutionData, opts::Dict, ctx_residual, t)

  setPCCtx(pc, mesh, sbp, eqn, opts, ctx_residual, t)  # unnecessary?

  pc.evalJacobian(mesh, sbp, eqn, opts, pc.assem)
  factorBDiagPC(pc, mesh, sbp, eqn, opts)
end

function applyPC(pc::NewtonBDiagPC, mesh::AbstractMesh, sbp::AbstractOperator,
                 eqn::AbstractSolutionData, opts::Dict, t, b::AbstractVector, 
                 x::AbstractVector)

  applyBDiagPC(pc, mesh, sbp, eqn, opts, b, x)

  return nothing
end


function applyPCTranspose(pc::NewtonBDiagPC, mesh::AbstractMesh, sbp::AbstractOperator,
                 eqn::AbstractSolutionData, opts::Dict, t, b::AbstractVector, 
                 x::AbstractVector)

  applyBDiagPC(pc, mesh, sbp, eqn, opts, b, x, trans=true)

  return nothing
end



# we could define applyPCTranspose, but this PC is ineffective so don't bother

#------------------------------------------------------------------------------
# block diagonal preconditioner implementation


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

  **Keyword Arguments**

   * trans: if true, apply the transposed operation, default false

"""
function applyBDiagPC(pc::NewtonBDiagPC, mesh, sbp, eqn, opts, x::AbstractVector, b::AbstractVector; trans::Bool=false)

  # we need to do inv(A)*x = b --> solve A*b = x using the factorization A.

  @assert pc.is_factored
  if trans
    tchar = 'T'
  else
    tchar = 'N'
  end


  bs = pc.bs
  workvec = zeros(Float64, bs)  # hold values passed into LAPACK

  fill!(b, 0)

  for i=1:mesh.numEl
    jacf_i = sview(pc.diag_jac.A, :, :, i)
    ipiv_i = sview(pc.ipiv, :, i)

    getValues(mesh, x, i, workvec)

    # call Lapack GETRS to solve for b (in-place)
    getrs2!(tchar, jacf_i, ipiv_i, workvec)

    setValues(mesh, workvec, i, b)

  end  # end loop i

  return nothing
end

# these names are too generic to export
import Jacobian: getValues, setValues


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

  **Keyword Arguments**

   * trans: if true, apply the transpose, default false

   Aliasing Restrictions: `x` and `b` cannot alias
"""
function applyBDiagPCInv(pc::NewtonBDiagPC, mesh, sbp, eqn, opts,
                         x::AbstractVector, b::AbstractVector;
                         zero_output=true, trans::Bool=false)


  if !pc.is_factored
    diagMatVec(pc.diag_mat, mesh, x, b, zero_output=zero_output, trans=trans)
  else

    if zero_output
      fill!(b, 0)
    end

    _trans::Bool = trans

    bs = pc.bs
    workvec = zeros(Float64, bs)  # hold values passed into LAPACK

    for i=1:mesh.numEl
      jacf_i = sview(pc.diag_jac.A, :, :, i)
      ipiv_i = sview(pc.ipiv, :, i)

      getValues(mesh, x, i, workvec)

      if trans
        laswp!(workvec, 1, bs, ipiv_i)
        trmv!('L', 'T', 'U', jacf_i, workvec)
        trmv!('U', 'T', 'N', jacf_i, workvec)
      else
        # do b = P*L*U*x
        trmv!('U', 'N', 'N', jacf_i, workvec)
        trmv!('L', 'N', 'U', jacf_i, workvec)
        applyIpiv!(ipiv_i, workvec)
      end

      setValues(mesh, workvec, i, b)
    end  # end loop i

  end  # end if

  return nothing
end


### NewtonBJacobiPC
#------------------------------------------------------------------------------
# AbstractPC Interface
function calcPC(pc::NewtonBJacobiPC, mesh::AbstractMesh, sbp::AbstractOperator,
                eqn::AbstractSolutionData, opts::Dict, ctx_residual, t)

  calcPC(pc.diag_pc, mesh, sbp, eqn, opts, ctx_residual, t)
  setPCCtx(pc, mesh, sbp, eqn, opts, ctx_residual, t)  # unnecessary?
end


"""
  Applies the `NewtonBJacobiPC`.  Note that the initial value of `x` is
  used as the initial guess.
"""
function applyPC(pc::NewtonBJacobiPC, mesh::AbstractMesh, sbp::AbstractOperator,
                 eqn::AbstractSolutionData, opts::Dict, t, b::AbstractVector, 
                 x::AbstractVector{T}; trans::Bool=false) where {T}

  @assert length(x) == length(b)
  @assert length(x) == mesh.numDof

  t1 = zeros(T, length(x))
  t2 = zeros(T, length(x))
  _trans::Bool = trans

  exit_status = 1  # 1 = itermax, 2 = res_tol
  for i=1:pc.itermax

    # t1 = R*x
    fill!(t1, 0)
    computeRProduct(pc, mesh, sbp, eqn, opts, x, t1, trans=_trans)

    # t1 = b - R*x
    @simd for j=1:length(t1)
      t1[j] = b[j] - t1[j]
    end
    
    if pc.res_tol > 0 || pc.verbose  # save the expense of the mat-vec if residual will
                       # not be checked
      # compute residual: b - A*x = t1 - D*x
      applyBDiagPCInv(pc.diag_pc, mesh, sbp, eqn, opts, x, t2; trans=_trans)
      @simd for j=1:length(t2)
        t2[j] = t1[j] - t2[j]
      end

      res_norm = calcEuclidianNorm(eqn.comm, t2)
      if pc.verbose
        println(BSTDOUT, "  iteration ", i, " linear residual norm: ", res_norm)
        println(BSTDOUT,  "    max res = ", maximum(abs.(t2)))
      end

      if res_norm < pc.res_tol
        exit_status = 2
        break
      end
    end  # end if res_tol > 0

    # x = D^-1 * (b - R*x) = D^-1 * t1
    if _trans
      applyPCTranspose(pc.diag_pc, mesh, sbp, eqn, opts, 0.0, t1, x)
    else
      applyPC(pc.diag_pc, mesh, sbp, eqn, opts, 0.0, t1, x)
    end
  end

  if pc.verbose
    if exit_status == 1
      println(BSTDOUT, "exited Block Jacobi smoother due to itermax")
      # res_norm is one iteration old, so don't print it here
    else
      println(BSTDOUT, "exited Block Jacobi smoother due to residual norm")

    end
  end

  flush(BSTDOUT)

  return nothing
end


function applyPCTranspose(pc::NewtonBJacobiPC, mesh::AbstractMesh,
                 sbp::AbstractOperator,
                 eqn::AbstractSolutionData, opts::Dict, t, b::AbstractVector, 
                 x::AbstractVector{T}) where {T}

  applyPC(pc, mesh, sbp, eqn, opts, t, b, x, trans=true)
end


"""
  Multiply R (= A - D), where A is the Jacobian and D is the block diagonal), by
  a vector.

  **Inputs**

   * pc
   * mesh
   * sbp
   * eqn
   * opts
   * x: vector to multiply against

  **Inputs/Outputs**

   * b: vector to overwrite with result

  **Inputs/Outputs**

   * trans: if true, compute R.'x, default false
"""
function computeRProduct(pc::NewtonBJacobiPC, mesh, sbp, eqn, opts, x, b; trans::Bool=false)

  # compute R*x as A*x - D*x
  if trans
    pc.evalJacTVecProduct(mesh, sbp, eqn, opts, x, b)
  else
    pc.evalJacVecProduct(mesh, sbp, eqn, opts, x, b)
  end
  scale!(b, -1)
  applyBDiagPCInv(pc.diag_pc, mesh, sbp, eqn, opts, x, b, zero_output=false, trans=trans)
  scale!(b, -1)
  return nothing
end
