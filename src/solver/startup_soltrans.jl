# startup script for solving an equation

push!(LOAD_PATH, "/users/creanj/julialib_fork/PUMI.jl")
using PumiInterface # pumi interface
using PdePumiInterface  # common mesh interface - pumi
using SummationByParts  # SBP operators
include("../equation/Equation.jl")  # equation types
include("../rk4/rk4.jl")  # timestepping
include("./euler/euler.jl")  # solver functions
include("./euler/ic.jl")  # initial conditions functions
include("./euler/adaptfuncs.jl")  # mesh adapt function
#=
# timestepping parameters
delta_t = 0.5
t_max = 1.00
=#

# create operator
sbp = TriSBP{Float64}()  # create linear sbp operator
println("sbp.numnodes = ", sbp.numnodes)


function checkPtr(sbp)
  println("in main script")
  nnodes_before = sbp.numnodes
  sbp_ptr = pointer_from_objref(sbp)
  println("typeof(sbp_ptr) = ", sbp_ptr)
  println("sizeof(sbp_ptr) = ", sizeof(sbp_ptr))
  println("sbp_ptr = ", sbp_ptr)
  sbp2 = unsafe_pointer_to_objref(sbp_ptr)
  nnodes_after= sbp2.numnodes
  println("nnodes_before = ", nnodes_before, " nnodes_after = ", nnodes_after)
end



#=
nnodes_before = sbp.numnodes
sbp_ptr = pointer_from_objref(sbp)
println("typeof(sbp_ptr) = ", sbp_ptr)
sbp2 = unsafe_pointer_to_objref(sbp_ptr)
nnodes_after= sbp2.numnodes

println("nnodes_before = ", nnodes_before, " nnodes_after = ", nnodes_after)
=#


# create mesh
dmg_name = ".null"
#smb_name = "tri2l.smb"
#smb_name = "tri18l.smb"
smb_name = "adapt_big.smb"
mesh = PumiMesh2(dmg_name, smb_name, 1; dofpernode=4)  #create linear mesh with 1 dof per node

# create euler equation
eqn = EulerEquation(sbp)
println("eqn.bigQT_xi = \n", eqn.bigQT_xi)
println("eqn.bigQT_eta = \n", eqn.bigQT_eta)
println("sbp.QT_xi' = \n", sbp.Q[:,:,1].')
println("sbp.QT_eta' = \n", sbp.Q[:,:,2].')


# create vectors to hold solution at current, previous timestep
u0 = zeros(mesh.numDof)  # solution at previous timestep
u = zeros(mesh.numDof) # solution at current timestep



# populate u0 with initial condition
#ICLinear(mesh, sbp, eqn, u0) # change this
ICsmoothHeaviside(mesh, sbp, eqn, u0)
println("u0 = \n", u0)
saveSolutionToMesh(mesh, u0)

cfunc3 = cfunction(shockRefine2, Void, (Ptr{Void}, Ptr{Float64}, Ptr{Float64}, Ptr{Void}, Ptr{Void}, Ptr{TriSBP{Float64}}))

createAnisoFunc(mesh.m_ptr, cfunc3, mesh.f_ptr, sbp)
writeVtkFiles("output_pre", mesh.m_ptr)
checkPtr(sbp)

  runAnisoAdapt(mesh.m_ptr)

writeVtkFiles("output_post", mesh.m_ptr)
reinitPumiMesh2(mesh)
u = zeros(mesh.numDof)
retrieveSolutionFromMesh(mesh, u);
println("u = \n", u)

# test code, please ignore
#=
# test getF1
for i=1:16
  u0[i] = i
end

f1 = zeros(12)
getF1(mesh, sbp, eqn, u0, 1, f1)
f2 = zeros(12)
getF2(mesh, sbp, eqn, u0, 1, f2)


# test assembleU
#assembleU(f1, 1, u)

vec = [f1[1], f1[5], f1[9]]
assembleU(vec, 1, 1, u)


edgenum_local = getBoundaryEdgeLocalNum(mesh, 1)
println("edgenum_local = ", edgenum_local)
=#

#=
function evalEuler(t, x)
# this function is called by time stepping algorithm
# t is the current time
# x is the solution value at the previous timestep
# u = output, the function value at the current timestep
# u is declared outside this function to avoid reallocating memory

u[:] = 0.0  # zero out u before starting
evalVolumeIntegrals(mesh, sbp, eqn, u, x)
evalBoundaryIntegrals(mesh, sbp, eqn, u, x)
applyMassMatrixInverse(mesh, sbp, eqn, u, x)

return u

end  # end evalEuler


# call timestepper
rk4(evalEuler, delta_t, u0, t_max)
=#

function checkPtr(sbp)
  nnodes_before = sbp.numnodes
  sbp_ptr = pointer_from_objref(sbp)
  println("typeof(sbp_ptr) = ", sbp_ptr)
  sbp2 = unsafe_pointer_to_objref(sbp_ptr)
  nnodes_after= sbp2.numnodes
  println("nnodes_before = ", nnodes_before, " nnodes_after = ", nnodes_after)
end


