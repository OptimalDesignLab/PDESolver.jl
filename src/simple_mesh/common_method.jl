# the common method for developing code

function getElementVertCoords(mesh::simpleMesh)
  # Get IEN matrix
  # Getting Vertex coordinates
  coords = zeros(Float64, 3, 3, mesh.nel)
  for i = 1:mesh.nel
    coords[1:2,:,i] = mesh.vtx_loc[:,m.ien[1:3,i],i]
  end
  return mesh.ien, coords
end

function getShapeFunctionOrder(mesh::simpleMesh)
  order::Int = mesh.nnpe - 1
  return order
end

function getGlobalNodeNumber(mesh::simpleMesh, el_num::Integer, local_node_num::Integer)
  # globalNumber::Int = mesh.ien(local_node_num,el_num)
  # return globalNumber
  return mesh.ien[local_node_num,el_num]
end

function getNumEl(mesh::simpleMesh)
  # nel::Int = mesh.nel
  return mesh.nel
end

function numEdges(mesh::simpleMesh)
  # numEdges::Int = mesh.nedges
  # return numEdges
  return mesh.nedges
end

function getNumVerts(mesh::simpleMesh)
  # nVertices = mesh.nvertex
  return mesh.nvertex
end

function getNumNodes(mesh::simpleMesh)
  return mesh.nnp
end

function getBoundaryEdge(mesh::simpleMesh)
  return mesh.HBedges, mesh.VBedges, 
end