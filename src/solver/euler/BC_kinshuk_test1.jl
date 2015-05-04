# Function Called for evaluating boundary conditions

push!(LOAD_PATH, "../../../../SummationByParts/src/")
# push!(LOAD_PATH,"../src/simple_mesh/")
using SummationByParts
# using SimpleMesh

function bndryflux(u, dxidx, nrm)
  return u*sum(nrm.'*dxidx)
end

p = 1
sbp = TriSBP{Float64}(degree = p);
x = zeros(Float64,(2,sbp.numnodes,2));
vtx = [0. 0.; 1. 0.; 0. 1.]
x[:,:,1] = calcnodes(sbp, vtx)
vtx = [1. 0.; 1. 1.; 0. 1.]
x[:,:,2] = calcnodes(sbp, vtx)
dxidx = zeros(Float64, (2,2,sbp.numnodes,2))
jac = zeros(Float64, (sbp.numnodes,2))
mappingjacobian!(sbp, x, dxidx, jac)
bndryfaces = Array(Boundary, 4)
bndryfaces[1] = Boundary(1,1) # element 1 edge 1
bndryfaces[2] = Boundary(1,3) # element 1 edge 3
bndryfaces[3] = Boundary(2,1) # element 2 edge 1
bndryfaces[4] = Boundary(2,2) # element 2 edge 2

u = zeros(Float64, (sbp.numnodes,2)) # first dmension = numnodes, dimension 2 = elemNum
println(size(u))
println(u,'\n')
for d = 0:p
  for j = 0:d
    i = d-j
    # k = (x[1,:,:].^i).*(x[2,:,:].^j)
    # println(round(k,2), '\n', size(k))
    u = squeeze((x[1,:,:].^i).*(x[2,:,:].^j), 1)
    println(size(u))
    # println("Type of u = ", typeof(u), '\n')
    println(round(u,2))
    res = zeros(u)
    boundaryintegrate!(sbp, bndryfaces, u, dxidx, bndryflux, res)
    println("Result = ", '\n', round(res,2), '\n')
    println(typeof(res))
  end
end