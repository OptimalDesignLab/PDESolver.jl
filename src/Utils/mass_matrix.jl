# function for computing and applying the mass matrix and its inverse

@doc """
  This function calculates the inverse mass matrix and returns it.
  Because we use SBP operators, the mass matrix is diagonal, so it is stored
  in a vector.  mesh.dofs is used to put the components of the inverse
  mass matrix in the same place as the corresponding values in eqn.res_vec

  Arguments:
    mesh: AbstractMesh
    sbp: SBP operator
    eqn: an implementation of AbstractSolutionData. Does not have to be fully initialized.

  Outputs:
    Minv: vector containing inverse mass matrix

"""->
# used by AbstractSolutionData Constructor
# mid level functions
function calcMassMatrixInverse{Tmsh,  Tsol, Tres}(mesh::AbstractMesh{Tmsh}, 
                                                  sbp::AbstractSBP, 
                                                  eqn::AbstractSolutionData{Tsol, Tres})
# calculate the inverse mass matrix so it can be applied to the entire solution vector
# mass matrix is diagonal, stores in vector eqn.Minv

  Minv = zeros(Tmsh, mesh.numDof)

  for i=1:mesh.numEl
    for j=1:sbp.numnodes
      for k=1:mesh.numDofPerNode
        dofnum_k = mesh.dofs[k,j,i]
        # multiplication is faster than division, so do the divisions here
        # and then multiply solution vector times Minv
        Minv[dofnum_k] += (sbp.w[j]/mesh.jac[j,i])
      end
    end
  end

  for i=1:mesh.numDof
    Minv[i] = 1/Minv[i]
  end

  return Minv

end     # end of calcMassMatrixInverse function

@doc """
  This function calculate the mass matrix and returns it.
  Beause w are using SBP operators, the mass matrix is diagonal, so it is
  stored in a vector.

  Arguments:
    mesh: AbstractMesh
    sbp: SBP operator
    eqn: an implementation of AbstractSolutionData. Does not have to be fully initialized.

  Outputs:
    M: vector containing mass matrix

"""->
function calcMassMatrix{Tmsh,  Tsol, Tres}(mesh::AbstractMesh{Tmsh}, 
                                           sbp::AbstractSBP, 
                                           eqn::AbstractSolutionData{Tsol, Tres})
# calculate the (diagonal) mass matrix as a vector
# return the vector M

  M = zeros(Tmsh, mesh.numDof)

  for i=1:mesh.numEl
    for j=1:sbp.numnodes
      for k=1:mesh.numDofPerNode
        dofnum_k = mesh.dofs[k,j,i]
        # multiplication is faster than division, so do the divions here
        # and then multiply solution vector times M
        M[dofnum_k] += (sbp.w[j]/mesh.jac[j,i])
      end
    end
  end

  return M

end     # end of calcMassMatrix function

@doc """

  This function calculates the mass matrix and returns it, in a 3D array format suitable for
  application to eqn.res.

  Arguments:
    mesh: AbstractMesh
    sbp: SBP operator
    eqn: an implementation of AbstractSolutionData. Does not have to be fully initialized.

  Outputs:
    M: vector containing mass matrix

"""->
# calcMassMatrixInverse3D: 
#   calculates the inverse mass matrix, returning it as a 3D array suitable for application to eqn.res
function calcMassMatrixInverse3D{Tmsh,  Tsol, Tres}(mesh::AbstractMesh{Tmsh}, 
                                                  sbp::AbstractSBP, 
                                                  eqn::AbstractSolutionData{Tsol, Tres})

  Minv3D = zeros(Tmsh, mesh.numDofPerNode, sbp.numnodes, mesh.numEl)

  for i=1:mesh.numEl
    for j=1:sbp.numnodes
      for k=1:mesh.numDofPerNode
        dofnum_k = mesh.dofs[k,j,i]
        # multiplication is faster than division, so do the divisions here
        # and then multiply solution vector times Minv
        Minv3D[k, j, i] += (sbp.w[j]/mesh.jac[j,i])
      end
    end
  end

  for i=1:mesh.numEl
    for j=1:sbp.numnodes
      for k=1:mesh.numDofPerNode
        Minv3D[k, j, i] = 1/Minv3D[k, j, i]
      end
    end
  end

  return Minv3D

end 

@doc """
  This function multiplies eqn.res_vec (the residual in vector form  by eqn.Minv,
  the diagonal mass matrix.  This is a very fast function because all values
  are precomputed and stored linearly in memory.

  This is a mid level function, and does the correct thing regardless of the
  dimension of the equation.

  Aliasing restrictions: none
"""->
# mid level function (although it doesn't really need to Tdim)
function applyMassMatrixInverse{Tsol, Tres}(eqn::AbstractSolutionData{Tsol, Tres}, 
                                            res_vec::AbstractVector{Tsol})
  # apply the inverse mass matrix stored eqn to res_vec

  ndof = length(res_vec)
  for i=1:ndof
    res_vec[i] *= eqn.Minv[i]
  end

  return nothing
end



