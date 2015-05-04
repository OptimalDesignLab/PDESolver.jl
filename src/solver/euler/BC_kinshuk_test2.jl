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

# Initialize solutions
u0 = ones(mesh.numDof) # Previous time step solution
u = zeros(mesh.numDof) # Current time step solution


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

println('\n', "Face Normal \n ", sbp.facenormal)

calcIsentropicVortex(x[:,1,1], eqn, ug)
println('\n', "ug = ", ug)

dxidx = zeros(Float64, (2,2,sbp.numnodes,getNumEl(mesh))); # Jacobian Matrix
dofJacobian = zeros(Float64, (2,2,4*sbp.numnodes,getNumEl(mesh))) # Jacobian Matrix needed for boundary integrate
jac = zeros(Float64, (sbp.numnodes,getNumEl(mesh))); # Determinant of the Jacobian Matrix
mappingjacobian!(sbp, x, dxidx, jac) # Get the Jocabian for transformation between actual and iso-parametric space

bndryfaces = Array(Boundary, mesh.numBoundaryEdges)
bndryfaces = getBoundaryArray(mesh)


res = zeros(Float64, 4,sbp.numnodes,)



#=
ug1 = zeros(Float64, 4)
ug2 = zeros(Float64, 4)
ug3 = zeros(Float64, 4)
ug = zeros(Float64, 4*sbp.numnodes) # the actual solution to the problem
u = zeros(Float64, 4*sbp.numnodes) 
for i = 1:getNumEl(mesh)
  f1g = zeros(Float64, 4*sbp.numnodes) # linear finite element with 3 nodes(at vertices) and 4 dof pernode
  f2g = zeros(Float64, 4*sbp.numnodes)
  calcIsentropicVortex(x[:,1,i], eqn, ug1) # gets the ghost variable
  calcIsentropicVortex(x[:,2,i], eqn, ug2)
  calcIsentropicVortex(x[:,3,i], eqn, ug3)
  ug = append!(ug1, ug2, ug3)
  getF1(mesh, sbp, eqn, ug, i, f1g)
  getF2(mesh, sbp, eqn, ug, i, f2g)
  m = 1
  for j = 1:sbp.numnodes
    normal = zeros(Float64,2)
    normal = dxidx[:,:,j,i]*sbp.facenormal[:,j]
    eulerRoeSAT(sgn::Int, tau::double, normal, u, u1g, sat)
  end
end


# calculate boundary integrate individually for all dof
q =#