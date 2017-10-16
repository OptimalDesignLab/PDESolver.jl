# functions for calculating various preconditioners
"""
  This function computes the jacobian of the volume integrals for use as a
  preconditioner.  Boundary integrals are also included because it is easy
  to do so.

  This function uses eqn.q and eqn.res as temporary arrays.  On exit, 
  eqn.q will have the same value as on entry, and the real part of eqn.res will
  be consistent with the real part of eqn.q (but the complex part is undefined).
  On entry, the imaginary part of eqn.q must be zero.

  Complex step method only!

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * pert: the perturbation to use
   * func: evalResidual-like function (eqn.q -> eqn.res), same signature
    * t: current time

  **Inputs/Outputs**

   * newton_data: volume_prec is modified

  **Implementation Notes:**

  This function is actually doing a distance-0 coloring.  The number of residual
  evaluations is mesh.numDofPerNode*mesh.numNodesPerElement
"""
function calcVolumePreconditioner(newton_data::NewtonData, mesh, sbp, eqn, opts, pert, func::Function, t=0.0)

  # get original values of keys
#  addVolumeIntegrals = opts["add_volume_integrals"]
#  addBoundaryIntegrals = opts["add_boundaryIntegrals"]  #TODO: include this?
  addFaceIntegrals = opts["addFaceIntegrals"]
  # leave addStabilization alone

  opts["addFaceIntegrals"] = false
  newton_data.volume_prec.is_factored = false

  h = imag(pert)

  volume_jac = newton_data.volume_prec.volume_jac

  col = 1  # column of each element jacobian
  for i=1:mesh.numNodesPerElement
    for j=1:mesh.numDofPerNode

      # apply perturbation
      # volume integrals are local, so no need ot perturb parallel buffers
      for k=1:mesh.numEl
        eqn.q[i, j, k] += pert
      end

      func(mesh, sbp, eqn, opts, t)

      # extract jacobian of each element
      # also undo the perturbation
      for k=1:mesh.numEl
        pos = 1
        eqn.q[i, j, k] -= pert

        for p=1:mesh.numNodesPerElement
          for m=1:mesh.numDofPerNode
            volume_jac[pos, col, k] = imag(eqn.res[m, p, k])/h
            pos += 1
          end
        end
      end  # end loop k

      col += 1  # advance to next column
    end  # end loop j
  end  # end loop i

  # re-enable the face integrals

  opts["addFaceIntegrals"] = addFaceIntegrals

  return nothing
end

"""
  This function factors the volume preconditioner.  It uses Lapack to do an
  LU factorization with partial pivoting.

  **Inputs**
  
   * mesh
   * sbp
   * eqn
   * opts

  **Inputs/Outputs**

   * newton_data: volume_prec is updated.  volume_jac is factored in-place and
                  ipiv is overwritten with the permutation

"""
function factorVolumePreconditioner(newton_data::newtonData, mesh, sbp, eqn, opts)


  volume_prec = newton_data.volume_prec

  for i=1:mesh.numEl
    jac_i = sview(volume_prec.volume_jac, :, :, i)
    ipiv_i = sview(volume_prec.volume_jac, :, i)

    # call Lapack here

    # GETRF
  end

  volume_prec.is_factored = true

  return nothing
end


"""
  Solves Ax = b where A is the preconditioner matrix.
  This is not likely to be a sensible preconditioner for Contiuous Galerkin
  discretizations.

  **Inputs**

   * newton_data: NewtonData object containing a [`VolumePreconditioner `](@ref)
                  object, already factored
   * mesh
   * sbp
   * eqn
   * opts
   * b: the right hand side vector

  **Inputs/Outputs**

   * x: the output vector
"""
function applyVolumePreconditioner(newton_data::newtonData, mesh, sbp, eqn, opts, x::AbstractVector, b::AbstractVector)

  volume_prec = newton_data.volume_prec
  @assert volume_prec.is_factored

  for i=1:mesh.numEl 
    jacf_i = sview(volume_prec.volume_jac, :, :, i)

    # use mesh.dofs to get the right entries from b
    # the order has to be consistent with the order in which the dofs were
    # perturbed when calculating hte jacobian.

    # call Lapack GETRS to solve for x

    # put entries back into x using mesh.dofs

  end





