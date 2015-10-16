# getMass
#=
@doc"""
### getMass

Calculates the masss matrix and given the mesh and SBP operator

*  operator: SBP operator
*  mesh: Mesh object.

"""-> =#

function getMass(sbp::SBPOperator, mesh::PumiMesh2)
  # assemble mesh
  numnodes = mesh.numNodes  # number of dofs
  numdof = numnodes*mesh.numDofPerNode
  mass_matrix = zeros(numdof, numdof)

  for i=1:mesh.numEl
    dofnums_i = getGlobalNodeNumbers(mesh, i)
    nnodes = size(dofnums_i)[2]  # number of nodes
    for j=1:nnodes
      for k=1:mesh.numDofPerNode
        dofnum_k = dofnums_i[k,j]
        mass_matrix[dofnum_k, dofnum_k] += sbp.w[j]
      end
    end
  end

  return mass_matrix
end
