# startup script for solving an equation

 push!(LOAD_PATH, "/users/creanj/julialib_fork/PUMI.jl")
push!(LOAD_PATH, "../../../../PUMI")
using PumiInterface # pumi interface
using PdePumiInterface  # common mesh interface - pumi
using SummationByParts  # SBP operators
include("../../equation/Equation.jl")  # equation types
include("../../rk4/rk4.jl")  # timestepping
include("./euler.jl")  # solver functions
include("./ic.jl")  # initial conditions functions
include("./output.jl")  # printing results to files
# include("./euler/addEdgeStabilize.jl")  # printing results to files

# timestepping parameters
delta_t = 0.01
t_max = 20.0

# create operator
sbp = TriSBP{Float64}()  # create linear sbp operator

# create mesh
dmg_name = ".null"
#smb_name = "../../mesh_files/quarter_vortex3l.smb"
# smb_name = "../../mesh_files/quarter_vortex8l.smb"
smb_name = "../../mesh_files/tri2l.smb"
mesh = PumiMesh2(dmg_name, smb_name, 1; dofpernode=4)  #create linear mesh with 1 dof per node

# create euler equation
eqn = EulerEquation(sbp)
# println("eqn.bigQT_xi = \n", eqn.bigQT_xi)
# println("eqn.bigQT_eta = \n", eqn.bigQT_eta)
println("sbp.QT_xi' = \n", sbp.Q[:,:,1].')
println("sbp.QT_eta' = \n", sbp.Q[:,:,2].')


# create vectors to hold solution at current, previous timestep
SL0 = zeros(mesh.numDof)  # solution at previous timestep
SL = zeros(mesh.numDof) # solution at current timestep


# populate u0 with initial condition
# ICZero(mesh, sbp, eqn, SL0)
# ICLinear(mesh, sbp, eqn, SL0)
# ICIsentropicVortex(mesh, sbp, eqn, SL0)
ICRho1E2(mesh, sbp, eqn, SL0)
#ICIsentropicVortex(mesh, sbp, eqn, SL0)

SL_exact = deepcopy(SL0)

# ICIsentropicVortexWithNoise(mesh, sbp, eqn, SL0)

# more test code
#=
coords = [2, 1.5]
sol = zeros(4)
calcIsentropicVortex(coords, eqn, sol)
println("at (2, 1.5), sol = ", sol)
=#
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


function evalEuler(t, SL0)
  println("\n\n")
# this function is called by time stepping algorithm
# t is the current time
# x is the solution value at the previous timestep
# u = output, the function value at the current timestep
# u is declared outside this function to avoid reallocating memory

# SL[:] = 0.0  # zero out u before starting
SL = zeros(SL0)
println("SL0 = ", SL0)
# u, x, dxidx, jac, res, interface = dataPrep(mesh, sbp, eqn, SL, SL0)
# println("u = ", u)
# println("x = ", x)
# println("dxidx = ", dxidx)
# println("jac = ", jac)
# println("res = ", res)
# println("interface = ", interface)

evalVolumeIntegrals(mesh, sbp, eqn, SL, SL0)
println("VOLVOLVOL SL = ", SL)
evalBoundaryIntegrals(mesh, sbp, eqn, SL, SL0)
println("BCBCBCBC SL = ", SL)
SL_sum = sum(SL)
println("BCBCBCBC SL_sum: ",SL_sum)



addEdgeStabilize(mesh, sbp, eqn, SL, SL0)
# println("STABSTABSTAB SL = ", SL)
applyMassMatrixInverse(mesh, sbp, eqn, SL, SL0)
println("MASSMASSMASS SL = ", SL)

#=
# These two calls are for TESTING ONLY, delete in production code
println("SL0 = ", SL0,"\n\n")
println("SL = ", SL)
print("\n")
println("Running evalVolumeIntegrals")
evalVolumeIntegrals(mesh, sbp, eqn, SL, SL0)
print("\n")
println("SL0 = ", SL0,"\n\n")
println("SL = ", SL)
print("\n")
println("Running evalBoundaryIntegrals")
evalBoundaryIntegrals(mesh, sbp, eqn, SL, SL0)
print("\n")
println("SL0 = ", SL0,"\n\n")
println("SL = ", SL)
print("\n")
println("Running addEdgeStabilize")
addEdgeStabilize(mesh, sbp, eqn, SL, SL0)
print("\n")
println("at end: SL0 = ", SL0,"\n\n")
println("at end: SL = ", SL)
=#

println("+++++++++ SL +++++++++:\n",SL)

return SL

end  # end evalEuler

#=
# These two calls are for TESTING ONLY, delete in production code
println("SL0 = ", SL0,"\n\n")
println("SL = ", SL)
print("\n")
println("Running evalVolumeIntegrals")
evalVolumeIntegrals(mesh, sbp, eqn, SL, SL0)
print("\n")
println("SL0 = ", SL0,"\n\n")
println("SL = ", SL)
print("\n")
println("Running evalBoundaryIntegrals")
evalBoundaryIntegrals(mesh, sbp, eqn, SL, SL0)
print("\n")
println("SL0 = ", SL0,"\n\n")
println("SL = ", SL)
print("\n")
println("Running addEdgeStabilize")
addEdgeStabilize(mesh, sbp, eqn, SL, SL0)
print("\n")
println("at end: SL0 = ", SL0,"\n\n")
println("at end: SL = ", SL)
=#

# call timestepper
SL, SL_hist = rk4(evalEuler, delta_t, SL0, t_max)

SL_diff = SL - SL_exact
SL_norm = norm(SL_diff)/mesh.numDof
SL_side_by_side = [SL_exact  SL]

println("\n\n\n")
println("SL_diff: \n")
for i=1:size(SL_diff)[1]
  println(SL_diff[i,:])
end
println("SL_side_by_side: \n")
for i=1:size(SL_side_by_side)[1]
  println(SL_side_by_side[i,:])
end
println("SL_norm: \n",SL_norm,"\n")


saveSolutionToMesh(mesh, SL)
printSolution(mesh, SL)
printCoordinates(mesh)
writeVtkFiles("solution_done",mesh.m_ptr)

