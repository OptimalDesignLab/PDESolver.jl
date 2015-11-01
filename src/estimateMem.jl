vals = readdlm(ARGS[1], Int)

numV = vals[1]
numEdge = vals[2]
numEl = vals[3]
numBoundaryEdges = vals[4]
numInterfaces = numEdge - numBoundaryEdges
numVertNodes = vals[5]
numEdgeNodes = vals[6]
numFaceNodes = vals[7]
numFaceNodes = 2*numVertNodes + numEdgeNodes
numNodes = 3*numVertNodes + 3*numEdgeNodes + numFaceNodes # per element
numDofPerNode  = vals[8]
numNodes_total = numVertNodes*numV + numEdgeNodes*numEdge + numFaceNodes*numEl

numDof = numNodes_total*numDofPerNode
jac_mem = vals[9]

println("numInterfaces = ", numInterfaces)
println("numFaceNodes = ", numFaceNodes)
println("numNodes per element = ", numNodes)
println("total number of nodes = ", numNodes_total)
println("numDof = ", numDof)


# add up the memory usage of all large arrays

verts = numV*64
edges = numEdge*64
elements = numEl*64
bndryfaces = 40*numBoundaryEdges
interfaces = 96*numInterfaces
interface_normals = 2*2*numFaceNodes*numInterfaces*64
coords = 2*numNodes*numEl*64
dxidx = 2*2*numNodes*numEl*64
jac = numNodes*numEl*64
dofs = numDofPerNode*numNodes*numEl*32
sparsity_bnds = 2*numDof*32
color_masks = 10*numEl*1
neighbor_colors = 4*numEl*8
neighbor_nums = 4*numEl*32

mesh_sum = verts + edges + elements + bndryfaces + interfaces + interface_normals + coords + dxidx + jac + dofs + sparsity_bnds + color_masks + neighbor_colors + neighbor_nums

println("mesh memory usage estimate: ", mesh_sum, " bits")
println("                          = ", mesh_sum/8, " bytes")
println("                          = ", mesh_sum/(8*1024*1024), " megabytes")

q = numDofPerNode*numNodes*numEl*64
aux_vars = 1*numEl*64
flux_parametric = numDofPerNode*numNodes*numEl*2*64
res = numDofPerNode*numNodes*numEl*64
res_vec = numDof*64
q_vec = numDof*64
edgestab_alpha = 2*2*numNodes*numEl*64
bndryflux = numDofPerNode*numFaceNodes*numBoundaryEdges*64
stabscale = numNodes*numInterfaces*64
Minv = numDof*64

eqn_sum = q + aux_vars + flux_parametric + res + res_vec + q_vec + edgestab_alpha + bndryflux + stabscale + Minv

println("solution data memory usage estimate: ", eqn_sum, " bits")
println("                                   = ", eqn_sum/8, " bytes")
println("                                   = ", eqn_sum/(8*1024*1024), " megabytes")

total_sum = eqn_sum + mesh_sum

println("total memory usage estimate: ", total_sum, " bits")
println("                           = ", total_sum/8, " bytes")
println("                           = ", total_sum/(8*1024*1024), " megabytes")

println("Jacobian memory usage estimate: ", jac_mem, " bits")
println("                           = ", jac_mem/8, " bytes")
println("                           = ", jac_mem/(8*1024*1024), " megabytes")




