module SimpleMesh

include("simple_mesh_function.jl")

export simpleMesh

# abstract abstractMesh{T <: FloatingPoint}

@doc """
### SimpleMesh

simpleMesh is a mesh type which generates a rectangular mesh in 2D with triangular elements. It has the following fields

**Fields**

*  `lengthx` : Length of domain on the X-axis
*  `lengthy` : Length of domain on the Y-axis
*  `nedx`    : Number of element edges along X-axis
*  `nedy`    : Number of element edges along Y-axis
*  `nnpe`    : Number of nodes per element edge (Same for all edges)
*  `numNodes`  : Number of nodal points
*  `numEl`     : Number of elements mesh
*  `ien`     : Element information matrix
*  `vtx_loc` : Location coordinates of the vertices
*  `elem_edge_lengthx : Length of element edge along the X-axis
*  `elem_edge_lengthy : Length of element edge along the Y-axis
*  `NodeEdgex` : A 2*xnodes matrix that stores nodes on the boundaries along the X-axis. First row stores the nodes on edge 1 and the scond row stores nodes on edge 3 of the entire geometry. The nodes are stored from left to right.
*  `NodeEdgey` : A 2*ynodes matrix that stores nodes on the boundaries along the Y-axis.First row stores nodes on edge 2 and second row stores nodes on edge 4 of the entire geometry. The nodes are stored from bottom to top.
*  `HBedges` : A 2*nedx array which stores the edge numbers on the boundary along the X-axis. The first row stores the edges along the base of the rectangle. The second row stores edges on the top of the rectangle.
*  `VBedges` : A 2*nedy array which stores the edge numbers on the boundary along the Y-axis. The first row stores the vertical edges along the right side of the rectangle. The second row stores the edges on the left side of the rectangle.
*  `numEdge`  : Total number of mesh edges in the geometry.

"""->

# type simpleMesh{T} <: abstractMesh{T}
type simpleMesh{T <: FloatingPoint}

	lengthx::FloatingPoint
	lengthy::FloatingPoint
	nedx::Int
	nedy::Int
	nnpe::Int
	order::Int
	numNodes::Int
	numVert::Int
	numEdge::Int
	numEl::Int
	numDof::Int
	numDofPerNode::Int
	numBoundaryEdges::Int
	ien::Array{Int}
	vtx_loc::Array{T}
	elem_edge_lengthx::FloatingPoint
	elem_edge_lengthy::FloatingPoint
	NodeEdgex::Array{Int}
	NodeEdgey::Array{Int}
	HBedges::Array{Int}
	VBedges::Array{Int}

	# dofs::Array{Int,3} # Store dof numbers of solution array for speedy assembly

	# numBC::Int  # number of boundary conditions
	

	function simpleMesh(lengthx,lengthy,nedx,nedy,nnpe,numDofPerNode)
		order = nnpe-1 # Calculate the mesh order

		# calculating the members of the type using functions
		(ien,numNodes,numEl) = createMesh(lengthx,lengthy,nedx,nedy,nnpe)

		numVert = (nedx+1)*(nedy+1)
		numDof = numNodes * numDofPerNode
		numBoundaryEdges = 2*(nedx+nedy)

		(vtx_loc,elem_edge_lengthx,elem_edge_lengthy) = nodeLocation(lengthx,lengthy,nedx,nedy,nnpe)

		(NodeEdgex, NodeEdgey) = boundaryNodeInfo(nedx,nedy,nnpe)

		(HBedges, VBedges, numEdge) = boundaryEdgeInfo(nedx,nedy)


		new(lengthx,lengthy,nedx,nedy,nnpe,order,numNodes,numVert,numEdge,numEl,
			numDof,numDofPerNode,numBoundaryEdges,ien,vtx_loc,elem_edge_lengthx,
			elem_edge_lengthy,NodeEdgex,NodeEdgey,HBedges,VBedges)

	end # Ends the function simpleMesh
end # Ends the type simpleMesh



#=
function getElementVertCoords(mesh::simpleMesh)
  vtx_coord = zeros(Float64, 2, 3, mesh.numEl)
  for i = 1:mesh.numEl
    # coords[1:2,:,i] = mesh.vtx_loc[mesh.ien[1:3,i],i]
    # coords[1,1,i] = mesh.vtx_loc[1,]
    vtx_coord[:,:,i] = mesh.vtx_loc[:,mesh.ien[1:3,i]]
  end
  return vtx_coord
end

function getShapeFunctionOrder(mesh::simpleMesh)
  order::Int = mesh.nnpe - 1
  return order
end

function getGlobalNodeNumber(mesh::simpleMesh, el_num::Integer, local_node_num::Integer)
  return mesh.ien[local_node_num,el_num]
end

function getNumEl(mesh::simpleMesh)
  return mesh.numEl
end

function numEdges(mesh::simpleMesh)
  return mesh.numEdge
end

function getNumVerts(mesh::simpleMesh)
  return mesh.numVert
end

function getNumNodes(mesh::simpleMesh)
  return mesh.nnp
end

function getBoundaryEdge(mesh::simpleMesh)
  return mesh.HBedges, mesh.VBedges 
end

# end # Ends the abstract type abstractMesh
=#
end # Ends the module