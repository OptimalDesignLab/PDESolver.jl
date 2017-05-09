# Test for computing the derivative using the design variables. This includes
# using FFD and mesh movement. The disgn variables are angle of attack and the
# coordinates of the control points


resize!(ARGS, 1)
ARGS[1] = "../input_vals_vortex.jl"
# ARGS[1] = "../input_vals_airfoil.jl"

push!(LOAD_PATH, joinpath(Pkg.dir("FFD.jl"), "src"))

using FreeFormDeformation
using MeshMovement
using PdePumiInterface
import ODLCommonTools.sview

# Get the adjoint vector
include("../startup.jl")
objective = EulerEquationMod.createObjectiveFunctionalData(mesh, sbp, eqn, opts)
EulerEquationMod.evalFunctional(mesh, sbp, eqn, opts, objective)
adjoint_vec = zeros(Complex128, mesh.numDof)
EulerEquationMod.calcAdjoint(mesh, sbp, eqn, opts, objective, adjoint_vec)

# Write VTK files
PdePumiInterface.saveSolutionToMesh(mesh, real(adjoint_vec))
PdePumiInterface.writeVisFiles(mesh, "adjoint_field")

# Initialize FFD and MeshWarping
geom_faces = opts["BC4"]
# geom_faces = opts["BC2"]
include("initMeshMotion.jl")

dpsiTRdXs = zeros(Float64, 2*3*sum(nwall_faces))
EulerEquationMod.evalLagrangianDerivative(mesh, sbp, eqn, opts, objective,
                                          adjoint_vec, dpsiTRdXs, ffd_map)

#=
# Verify against finite difference
eqn.params.aoa = opts["aoa"]
fill!(eqn.res, 0.0)
fill!(eqn.res_vec, 0.0)
res_norm = NonlinearSolvers.calcResidual(mesh, sbp, eqn, opts, evalEuler)
orig_res_vec = copy(eqn.res_vec)
eqn.params.aoa += 1e-6
fill!(eqn.res, 0.0)
fill!(eqn.res_vec, 0.0)
res_norm = NonlinearSolvers.calcResidual(mesh, sbp, eqn, opts, evalEuler)
dRdAlpha_FD = (eqn.res_vec - orig_res_vec)/1e-6

f = open("dRdAlpha_error.dat", "w")
for i = 1:length(eqn.res_vec)
  println(f, dRdAlpha_FD[i] - dRdAlpha[i])
end
close(f)
println("dRdAlpha error norm = ", norm(dRdAlpha_FD - dRdAlpha, 2))
=#









# Warp The mesh
warpMesh(param, volNodes, wallCoords)

# Redistribute it to mesh.coords
for i = 1:mesh.numEl
  for j = 1:mesh.numNodesPerElement
    local_vertnum = mesh.element_vertnums[j,i]
    mesh.vert_coords[2,j,i] = volNodes[2,local_vertnum]
    mesh.vert_coords[1,j,i] = volNodes[1,local_vertnum]
  end
end
for i = 1:mesh.numEl
  update_coords(mesh, i, mesh.vert_coords[:,:,i])
end  # End



PumiInterface.writeVtkFiles("warped_mesh", mesh.m_ptr)
gc()







#=
dLdx_adjoint = dJdAlpha + dot(adjoint_vec, dRdAlpha)
#=
#----- Finite Differencing -----#
# Get the design variable array
x_dv = reshape(ffd_map.cp_xyz, length(ffd_map.cp_xyz))
x_dv = append!([eqn.params.aoa], x_dv)
=#
# Finite difference derivative of Lagrangian wrt aoa, which is x_dv[1]

pert = 1e-6 # FD perturbation
eqn.params.aoa = opts["aoa"] + pert
fill!(eqn.q, 0.0)
fill!(eqn.q_vec, 0.0)
fill!(eqn.res, 0.0)
fill!(eqn.res_vec, 0.0)

# Rerun with the perturbed value
ICfunc_name = opts["IC_name"]
ICfunc = EulerEquationMod.ICDict[ICfunc_name]
ICfunc(mesh, sbp, eqn, opts, eqn.q_vec)
pmesh = mesh
EulerEquationMod.init(mesh, sbp, eqn, opts, pmesh)
call_nlsolver(mesh, sbp, eqn, opts, pmesh)

EulerEquationMod.evalFunctional(mesh, sbp, eqn, opts, objective)
dLdx = (real(objective.lift_val) - orig_val)/pert
println("orig_val = $orig_val, new_val = $(real(objective.lift_val)), dLdx = $dLdx")
println("dLdx_adjoint = $dLdx_adjoint")

errfd_norm = norm(dLdx - dLdx_adjoint,2)
println("err dLdx = $errfd_norm")
println("dJdAlpha = $dJdAlpha, dJdAlpha_comp = $dJdAlpha_comp")
println("dJdAlpha error = ", norm(dJdAlpha - dJdAlpha_comp,2))
=#
