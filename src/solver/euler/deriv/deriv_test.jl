# Test for computing the derivative using the design variables. This includes
# using FFD and mesh movement. The disgn variables are angle of attack and the
# coordinates of the control points


resize!(ARGS, 1)
# ARGS[1] = "../input_vals_vortex.jl"
# ARGS[1] = "../input_vals_airfoil.jl"
ARGS[1] = "../input_vals_airfoil_parallel.jl"


push!(LOAD_PATH, joinpath(Pkg.dir("FFD.jl"), "src"))

using FreeFormDeformation
using MeshMovement
using PdePumiInterface
import ODLCommonTools.sview

# Get the adjoint vector
include("../startup.jl")
objective = EulerEquationMod.createObjectiveFunctionalData(mesh, sbp, eqn, opts)
EulerEquationMod.evalFunctional(mesh, sbp, eqn, opts, objective)

lift = objective.lift_val

adjoint_vec = zeros(Complex128, mesh.numDof)
EulerEquationMod.calcAdjoint(mesh, sbp, eqn, opts, objective, adjoint_vec)


dJdaoa = EulerEquationMod.eval_dJdaoa(mesh, sbp, eqn, opts, objective, "lift", adjoint_vec)

# Check complete derivatives w.r.t alpha using finite difference
pert = 1e-6
eqn.params.aoa += pert
EulerEquationMod.solve_euler(mesh, sbp, eqn, opts, mesh)
EulerEquationMod.evalFunctional(mesh, sbp, eqn, opts, objective)
lift_pert = objective.lift_val
drag_pert = objective.drag_val
eqn.params.aoa -= pert

dLiftdaoa_fd = (lift_pert - lift)/pert

println("dJdaoa = $dJdaoa, dLiftdaoa_fd = $dLiftdaoa_fd")
deriv_err = norm(dLiftdaoa_fd - dJdaoa, 2)
println("Error = $deriv_err")
