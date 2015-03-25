module SimpleMesh

include("simple_mesh_function.jl")

export simpleMesh

abstract abstractMesh{T <: FloatingPoint}

@doc """
### SimpleMesh

simpleMesh is a mesh type which generates a rectangular mesh in 2D with triangular elements. It has the following fields

**Fields**

*  `lengthx` : Length on the X-axis
*  `lengthy` : Length on the Y-axis
*  `nedx`    : Number of element edges along X-axis
*  `nedy`    : Number of element edges along Y-axis
*  `nnpe`    : Number of nodes per element edge (Same for all edges)
*  `nnp`     : Number of nodal points
*  `nel`     : Number of elements mesh
*  `ien`     : Element information matrix
*  `vtx_loc` : Location coordinates of the vertices
*  `elem_edge_lengthx : Length of element edge along the X-axis
*  `elem_edge_lengthy : Length of element edge along the Y-axis
*  `NodeEdgex` : A 2*xnodes matrix that stores nodes on the boundaries along the X-axis. First row stores the nodes on edge 1 and the scond row stores nodes on edge 3 of the entire geometry. The nodes are stored from left to right.
*  `NodeEdgey` : A 2*ynodes matrix that stores nodes on the boundaries along the Y-axis.First row stores nodes on edge 2 and second row stores nodes on edge 4 of the entire geometry. The nodes are stored from bottom to top.
*  `HBedges` : A 2*nedx array which stores the edge numbers on the boundary along the X-axis. The first row stores the edges along the base of the rectangle. The second row stores edges on the top of the rectangle.
*  `VBedges` : A 2*nedy array which stores the edge numbers on the boundary along the Y-axis. The first row stores the vertical edges along the right side of the rectangle. The second row stores the edges on the left side of the rectangle.
*  `nedges`  : Total number of mesh edges in the geometry.

"""->

type simpleMesh{T} <: abstractMesh{T}

	lengthx::FloatingPoint
	lengthy::FloatingPoint
	nedx::Int
	nedy::Int
	nnpe::Int
	nnp::Int
	nel::Int
	ien::Array{Int}
	vtx_loc::Array{T}
	elem_edge_lengthx::FloatingPoint
	elem_edge_lengthy::FloatingPoint
	NodeEdgex::Array{Int}
	NodeEdgey::Array{Int}
	HBedges::Array{Int}
	VBedges::Array{Int}
	nedges::Int

	function simpleMesh(lengthx,lengthy,nedx,nedy,nnpe)

		# calculating the members of the type using functions
		(ien,nnp,nel) = createMesh(lengthx,lengthy,nedx,nedy,nnpe)
		(vtx_loc,elem_edge_lengthx,elem_edge_lengthy) = nodeLocation(lengthx,lengthy,nedx,nedy,nnpe)
		(NodeEdgex, NodeEdgey) = boundaryNodeInfo(nedx,nedy,nnpe)
		(HBedges, VBedges, nedges) = boundaryEdgeInfo(nedx,nedy)
		# new(ien,nnp,nel,vtx_loc,elem_edge_lengthx,elem_edge_lengthy,NodeEdgex, NodeEdgey,HBedges, VBedges)
		new(lengthx,lengthy,nedx,nedy,nnpe,nnp,nel,ien,vtx_loc,elem_edge_lengthx,
			elem_edge_lengthy,NodeEdgex,NodeEdgey,HBedges,VBedges,nedges)

	end # Ends the function simpleMesh
end # Ends the type simpleMesh

end # Ends the abstract type abstractMesh
# end # Ends the module

