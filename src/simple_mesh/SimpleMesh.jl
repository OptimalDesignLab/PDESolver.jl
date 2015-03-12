module SimpleMesh

include("simple_mesh_function.jl")

export simpleMesh

abstract abstractMesh{T <: FloatingPoint}

type simpleMesh{T} <: abstractMesh{T}

	lengthx::FloatingPoint
	lengthy::FloatingPoint
	nedx::Int
	nedy::Int
	nnpe::Int
	nnp::Int
	nel::Int
	IEN::Array{Int}
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
		(IEN,nnp,nel) = createMesh(lengthx,lengthy,nedx,nedy,nnpe)
		(vtx_loc,elem_edge_lengthx,elem_edge_lengthy) = nodeLocation(lengthx,lengthy,nedx,nedy,nnpe)
		(NodeEdgex, NodeEdgey) = boundaryNodeInfo(nedx,nedy,nnpe)
		(HBedges, VBedges, nedges) = boundaryEdgeInfo(nedx,nedy)
		# new(IEN,nnp,nel,vtx_loc,elem_edge_lengthx,elem_edge_lengthy,NodeEdgex, NodeEdgey,HBedges, VBedges)
		new(lengthx,lengthy,nedx,nedy,nnpe,nnp,nel,IEN,vtx_loc,elem_edge_lengthx,
			elem_edge_lengthy,NodeEdgex,NodeEdgey,HBedges,VBedges,nedges)

	end # Ends the function simpleMesh
end # Ends the type simpleMesh

end # Ends the abstract type abstractMesh
# end # Ends the module

