# startup script for solving an equation

 push!(LOAD_PATH, "/users/creanj/julialib_fork/PUMI.jl")
 push!(LOAD_PATH, "../../equation")
 push!(LOAD_PATH, "../../common")
 push!(LOAD_PATH, "/users/creanj/.julia/v0.4/PDESolver/src/solver/euler")
#push!(LOAD_PATH, "../../../../PUMI")
using PDESolverCommon
using PumiInterface # pumi interface
using PdePumiInterface  # common mesh interface - pumi
using SummationByParts  # SBP operators
using EulerEquationMod
include("../../rk4/rk4.jl")  # timestepping
include("./output.jl")  # printing results to files


function runtest()
# timestepping parameters
delta_t = 0.005
#t_max = 0.025
t_max = 5.00
#t_max = 1.0

# types of the mesh, SBP, Equation objects
Tmsh= Float64
Tsbp = Float64
Tsol = Float64

# create operator
sbp = TriSBP{Tsbp}()  # create linear sbp operator

# create mesh
dmg_name = ".null"
smb_name = "../../mesh_files/quarter_vortex3l.smb"
#smb_name = "../../mesh_files/quarter_vortex8l.smb"
#smb_name = "../../mesh_files/tri30l.smb"
mesh = PumiMesh2{Tmsh}(dmg_name, smb_name, 1, sbp; dofpernode=4)  #create linear mesh with 4 dof per node


# create euler equation
eqn = EulerEquation{Tsol}(mesh, sbp)


SL0 = zeros(Tsol, mesh.numDof)  # solution at previous timestep
SL = zeros(Tsol, mesh.numDof) # solution at current timestep


# populate u0 with initial condition
# ICZero(mesh, sbp, eqn, SL0)
# ICLinear(mesh, sbp, eqn, SL0)
 ICIsentropicVortex(mesh, sbp, eqn, SL0)
#ICRho1E2(mesh, sbp, eqn, SL0)
#ICRho1E2U3(mesh, sbp, eqn, SL0)

#ICVortex(mesh, sbp, eqn, SL0)
#ICIsentropicVortex(mesh, sbp, eqn, SL0)

SL_exact = deepcopy(SL0)


saveSolutionToMesh(mesh, SL0)

writeVtkFiles("solution_ic",mesh.m_ptr)



# call timestepper

extra_args = (mesh, sbp, eqn)

SL = rk4(evalEuler, delta_t,SL,  SL0, t_max, extra_args)
#println("rk4 @time printed above")

SL_diff = SL - SL_exact
SL_norm = norm(SL_diff)/mesh.numDof
#SL_side_by_side = [SL_exact  SL]

#=
println("\n\n\n")
println("SL_diff: \n")
for i=1:size(SL_diff)[1]
  println(SL_diff[i,:])
end
println("SL_side_by_side: \n")
for i=1:size(SL_side_by_side)[1]
  println(i, " ", SL_side_by_side[i,:])
end
=#
println("SL_norm: \n",SL_norm,"\n")


saveSolutionToMesh(mesh, SL)
printSolution(mesh, SL)
printCoordinates(mesh)
writeVtkFiles("solution_done",mesh.m_ptr)

end


runtest()
