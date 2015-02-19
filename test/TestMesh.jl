# module TestMesh

using PDESolver
include("../src/simple_mesh/SimpleMeshType.jl")


using .SimpleMeshType

# Make sure that the elemnt size is a factor of the length and width
# [ ] = MeshType.createMesh(Length,Width,Element Size)

(xnodes,ynodes,nnp,nel,IEN) = SimpleMeshType.createMesh(4,4,4,4)
(NodeEdgex,NodeEdgey,ElemEdgex,ElemEdgey) = SimpleMeshType.boundaryInfo(xnodes,ynodes)
println(xnodes)
println(ynodes)





# end  # module
