# module TestMesh

using PDESolver
include("../src/simple_mesh/SimpleMeshType.jl")


using .SimpleMeshType

# Make sure that the elemnt size is a factor of the length and width
# [ ] = MeshType.createMesh(Length,Width,Element Size)

lengthx = 3;
lengthy = 2;
nedx = 3;
nedy = 2;
nnpe = 5;

(xnodes,ynodes,nnp,nel,IEN) = SimpleMeshType.createMesh(lengthx,lengthy,nedx,nedy,nnpe)
println(xnodes)
println(ynodes)
(NodeEdgex,NodeEdgey,HBedges,VBedges) = SimpleMeshType.boundaryInfo(xnodes,ynodes,nedx,nedy,nnpe)
println(NodeEdgex)
println(NodeEdgey)
println(HBedges)
println(VBedges)





# end  # module
