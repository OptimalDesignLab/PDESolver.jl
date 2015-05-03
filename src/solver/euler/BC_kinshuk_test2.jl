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
include("./ic_kinshuk.jl")  # initial conditions functions

function bndryflux(u, dxidx, nrm)
  return u*sum(nrm.'*dxidx)
end

sbp = TriSBP{Float64}()  # create linear sbp operator

# timestepping parameters
delta_t = 0.5
t_max = 1.00

# create mesh
dmg_name = ".null"
smb_name = "tri8l.smb"
mesh = PumiMesh2(dmg_name, smb_name, 1; dofpernode=4)  #create linear mesh with 4 dof per node

eqn = EulerEquation(sbp)

# Initialize solutions
u0 = ones(mesh.numDof) # Previous time step solution
u = zeros(mesh.numDof) # Current time step solution

# for i = 1:16
#   u0[i] = i
# end

ICZero(mesh, sbp, eqn, u0) # populate u0 with initial condition

# Boundary integrals
# evalBoundaryIntegrals(mesh::AbstractMesh, operator::SBPOperator, eqn::EulerEquation, u::AbstractVector, u0::AbstractVector)

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
# for i = 1:size(bnd_array,1)
  # println(bnd_array[i].face)
# end

# bndryJacobian = zeros(Float64, (2,2,))

# bndryfaces = Array(Boundary, 4)
# bndryfaces[1] = Boundary(1,1) # first element edge 1
# bndryfaces[2] = Boundary(1,2) # element 1 edge 3
# bndryfaces[3] = Boundary(2,2) # element 2 edge 1
# bndryfaces[4] = Boundary(2,3) # element 2 edge 2

F1 = zeros(Float64, 2,sbp.numnodes,getNumEl(mesh))
F2 = zeros(Float64, 2,sbp.numnodes,getNumEl(mesh))
F3 = zeros(Float64, 2,sbp.numnodes,getNumEl(mesh))
F4 = zeros(Float64, 2,sbp.numnodes,getNumEl(mesh))

for i = 1:getNumEl(mesh)
  f1 = zeros(Float64, 4*sbp.numnodes) # linear finite element with 3 nodes(at vertices) and 4 dof pernode
  f2 = zeros(Float64, 4*sbp.numnodes)
  println(typeof(f1))
  println(typeof(f2))
  getF1(mesh, sbp, eqn, u0, i, f1)
  getF2(mesh, sbp, eqn, u0, i, f2)
  m::Int = 1
  for j = 1:sbp.numnodes
    F1[1,j,i] = f1[m]
    F1[2,j,i] = f2[m]
    F2[1,j,i] = f1[m+1]
    F2[2,j,i] = f2[m+1]
    F3[1,j,i] = f1[m+2]
    F3[2,j,i] = f2[m+2]
    F4[1,j,i] = f1[m+3]
    F4[2,j,i] = f2[m+3]
    m = m+4
  end
end

println(F1, '\n')
println(F2, '\n')
println(F3, '\n')
println(F4, '\n')

# println('\n', F)
res1 = zeros(F1)
res2 = zeros(F2)
res3 = zeros(F3)
res4 = zeros(F4)


# calculate boundary integrate individually for all dof

boundaryintegrate!(sbp, bndryfaces, F1, dxidx, bndryflux, res1)
boundaryintegrate!(sbp, bndryfaces, F2, dxidx, bndryflux, res2)
boundaryintegrate!(sbp, bndryfaces, F3, dxidx, bndryflux, res3)
boundaryintegrate!(sbp, bndryfaces, F4, dxidx, bndryflux, res4)

println('\n', "res1 = ", '\n', res1)