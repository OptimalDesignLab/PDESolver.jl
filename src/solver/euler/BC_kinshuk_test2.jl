# Test Boundary Conditions with Boundary flux

push!(LOAD_PATH, "../../../../SummationByParts/src/")
push!(LOAD_PATH, "../../../../PUMI")
# push!(LOAD_PATH, "../../../src/simple_mesh/")
using SummationByParts
# using SimpleMesh
using PumiInterface # pumi interface
using PdePumiInterface  # common mesh interface - pumi

include("../../equation/Equation.jl")  # equation types
include("../../rk4/rk4.jl")  # timestepping
include("./euler.jl")  # solver functions
include("./ic.jl")  # initial conditions functions
include("./BC_Roe_Kinshuk.jl")


sbp = TriSBP{Float64}()  # create linear sbp operator
eqn = EulerEquation(sbp)

# timestepping parameters
delta_t = 0.5
t_max = 1.00

# create mesh
dmg_name = ".null"
smb_name = "tri8l.smb"
mesh = PumiMesh2(dmg_name, smb_name, 1; dofpernode=4)  #create linear mesh with 4 dof per node

# Nodal Coordinates
x = zeros(Float64,(2,sbp.numnodes,getNumEl(mesh))); # nodal Coordinates of the marix
for i = 1:getNumEl(mesh)
  vtxcoord = getElementVertCoords(mesh, [i]);
  vtxcoord = squeeze(vtxcoord,3);
  vtxcoord = vtxcoord[1:2,:]
  vtxcoord = vtxcoord'
  # println(vtxcoord)
  x[:,:,i] = calcnodes(sbp, vtxcoord);
end

dxidx = zeros(Float64, (2,2,sbp.numnodes,getNumEl(mesh))); # Jacobian Matrix
dofJacobian = zeros(Float64, (2,2,4*sbp.numnodes,getNumEl(mesh))) # Jacobian Matrix needed for boundary integrate
jac = zeros(Float64, (sbp.numnodes,getNumEl(mesh))); # Determinant of the Jacobian Matrix
mappingjacobian!(sbp, x, dxidx, jac) # Get the Jocabian for transformation between actual and iso-parametric space

bndryfaces = Array(Boundary, mesh.numBoundaryEdges)
bndryfaces = getBoundaryArray(mesh)

# Calculate the intital condition
u0 = ones(Float64, 4, sbp.numnodes, mesh.numBoundaryEdges)
# u0 = ones(Float64, 4, sbp.numnodes, 2)

# res = zeros(Float64, 4,sbp.numnodes,mesh.numBoundaryEdges)
res = zeros(u0)

println("type of u0: ", typeof(u0))
println("size of dxidx: ", size(dxidx,4))
println("size of u0: ", size(u0,3))
println("size of res: ", size(res,3))
println("size of x: ", size(x,3))

boundaryintegrate!(sbp, bndryfaces, u0, x, dxidx, isentropicVortexBC, res)
