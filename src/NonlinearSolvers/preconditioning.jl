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

   * newton_data: vol_prec is modified

  **Implementation Notes:**

  This function is actually doing a distance-0 coloring.  The number of residual
  evaluations is mesh.numDofPerNode*mesh.numNodesPerElement
"""
function calcVolumePreconditioner(newton_data::NewtonData, mesh, sbp, eqn, opts, pert, func::Function, t=0.0)

  println("computing volume PC")
  if opts["jac_method"] != 2
    throw(ErrorException("jac method must be 2 for calcVolumePreconditioner"))
  end

  # get original values of keys
#  addVolumeIntegrals = opts["add_volume_integrals"]
#  addBoundaryIntegrals = opts["add_boundaryIntegrals"]  #TODO: include this?
  addFaceIntegrals = opts["addFaceIntegrals"]
  # leave addStabilization alone

  opts["addFaceIntegrals"] = false
  newton_data.vol_prec.is_factored = false

  h = imag(pert)

  volume_jac = newton_data.vol_prec.volume_jac

  col = 1  # column of each element jacobian
  for i=1:mesh.numNodesPerElement
    for j=1:mesh.numDofPerNode

      # apply perturbation
      # volume integrals are local, so no need ot perturb parallel buffers
      for k=1:mesh.numEl
        eqn.q[j, i, k] += pert
      end

      func(mesh, sbp, eqn, opts, t)

      # extract jacobian of each element
      # also undo the perturbation
      for k=1:mesh.numEl
        pos = 1
        eqn.q[j, i, k] -= pert

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
  LU factorization with partial pivoting.  The user is allowed to modify
  newton_data.vol_prec.volume_jac after calling [`calcVolumePreconditioner`](@ref)
  but before calling this routine.

  **Inputs**
  
   * mesh
   * sbp
   * eqn
   * opts

  **Inputs/Outputs**

   * newton_data: vol_prec is updated.  volume_jac is factored in-place and
                  ipiv is overwritten with the permutation

"""
function factorVolumePreconditioner(newton_data::NewtonData, mesh, sbp, eqn, opts)


  volume_prec = newton_data.vol_prec

  for i=1:mesh.numEl
    jac_i = sview(volume_prec.volume_jac, :, :, i)
    ipiv_i = sview(volume_prec.ipiv, :, i)

    # call Lapack here
    info = getrf!(jac_i, ipiv_i)
    @assert info == 0
  end

  volume_prec.is_factored = true

  return nothing
end

"""
  This function is a wrapper around [`applyVolumePreconditioner`](@ref)
  that is passed into Petsc as a callback

  The PC ctx must be: (mesh, sbp, eqn, opts, newton_data, func, ctx_residual,
  t), where ctx_residual is the ctx used by all the CN functions and
  func is the right hand side function passed into [`newtonInner`](@ref)

"""
function applyVolumePC_wrapper(pc::PC, x::PetscVec, b::PetscVec)

  println("applying volume PC")
  ctx_petsc_pc_ptr = PCShellGetContext(pc)
  ctx_petsc_pc = unsafe_pointer_to_objref(ctx_petsc_pc_ptr)

  # the ctx_petsc_pc is the same as for a shell matrix for convenience
  mesh = ctx_petsc_pc[1]
  sbp = ctx_petsc_pc[2]
  eqn = ctx_petsc_pc[3]
  opts = ctx_petsc_pc[4]
  newton_data = ctx_petsc_pc[5]
  func = ctx_petsc_pc[6]  # rhs_func from newtonInner
  ctx_residual = ctx_petsc_pc[7]
  t = ctx_petsc_pc[8]

  # get the arrays underlying x and b
  x_arr, xptr = PetscVecGetArrayRead(x)  # read only
  b_arr, bptr = PetscVecGetArray(b)  # writeable

  applyVolumePreconditioner(newton_data, mesh, sbp, eqn, opts, x_arr, b_arr)

  println("diffnorm = ", norm(x_arr - b_arr))
  PetscVecRestoreArrayRead(x, xptr)
  PetscVecRestoreArray(b, bptr)

  return PetscErrorCode(0)
end




"""
  Apply the preconditioner, ie. do inv(A)*x = b, where A is calculated by
  [`calcVolumePreconditioner`](@ref).

  This is not likely to be a sensible preconditioner for Contiuous Galerkin
  discretizations.

  **Inputs**

   * newton_data: NewtonData object containing a [`VolumePreconditioner `](@ref)
                  object, already factored by [`factorVolumePreconditioner`](@ref).
   * mesh
   * sbp
   * eqn
   * opts
   * x: the vector to multiply against

  **Inputs/Outputs**

   * b: the output vector

"""
function applyVolumePreconditioner(newton_data::NewtonData, mesh, sbp, eqn, opts, x::AbstractVector, b::AbstractVector)

  # we need to do inv(A)*x = b --> solve A*b = x using the factorization A.

  volume_prec = newton_data.vol_prec
  @assert volume_prec.is_factored

  jac_size = size(volume_prec.volume_jac, 1)
  workvec = zeros(Float64, jac_size)  # hold values passed into LAPACK

  for i=1:mesh.numEl
    jacf_i = sview(volume_prec.volume_jac, :, :, i)
    ipiv_i = sview(volume_prec.ipiv, :, i)
#    println("jacf = \n", jacf_i)

    # get the x values
    # this could be faster if we assume the dofs on each element are numbered
    # sequentially
    pos = 1
    for j=1:mesh.numNodesPerElement
      for k=1:mesh.numDofPerNode
        dof_i = mesh.dofs[k, j, i]
        workvec[pos] = x[dof_i]
        pos += 1
      end
    end

    # use mesh.dofs to get the right entries from b
    # the order has to be consistent with the order in which the dofs were
    # perturbed when calculating the jacobian.

    # call Lapack GETRS to solve for b (in-place)
    getrs2!('N', jacf_i, ipiv_i, workvec)

    # put entries back into b using mesh.dofs
    pos = 1
    for j=1:mesh.numNodesPerElement
      for k=1:mesh.numDofPerNode
        dof_i = mesh.dofs[k, j, i]
        b[dof_i] = workvec[pos]
        pos += 1
      end
    end



  end  # end loop i

  return nothing
end
