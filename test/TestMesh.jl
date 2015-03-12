# module TestMesh

push!(LOAD_PATH,"../src/simple_mesh/")

using PDESolver
using SimpleMesh
# include("../src/simple_mesh/SimpleMeshType.jl")
# include("../src/simple_mesh/SimpleMesh.jl")

# using .SimpleMeshType


lengthx = 3.0
lengthy = 2.0
nedx = 3;
nedy = 2;
nnpe = 4;
m = simpleMesh{Float64}(lengthx,lengthy,nedx,nedy,nnpe)
println(m.lengthx)
println(m.lengthy)
println(m.nedx)
println(m.nedy)
println(m.nnpe)
println(m.IEN)
println(m)


# The following part has been commented out
#=
(xnodes,ynodes,nnp,nel,IEN) = SimpleMeshType.createMesh(lengthx,lengthy,nedx,nedy,nnpe)
println(xnodes)
println(ynodes)
(NodeEdgex,NodeEdgey,HBedges,VBedges) = SimpleMeshType.boundaryInfo(xnodes,ynodes,nedx,nedy,nnpe)
println(NodeEdgex)
println(NodeEdgey)
println(HBedges)
println(VBedges)
=#
#=
m::simpleMesh
m.lengthx = 3
m.lengthy = 2
m.nedx = 3
m.nedy = 2
m.nnpe = 4
=#





# end  # module
