# testing derivatives w.r.t states
resize!(ARGS, 1)
ARGS[1] = "input_vals_airfoil.jl"


#----  Initialize EulerEquationMod for all the global variables necessary  ----#
include("startup.jl")

# Create the objective function data object
objective = EulerEquationMod.OptimizationData{Tsol}(mesh, sbp, opts)

# Copy all the original values
disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
orig_q_vec = deepcopy(eqn.q_vec)
original_res_vec = deepcopy(eqn.res_vec)


EulerEquationMod.evalFunctional(mesh, sbp, eqn, opts, objective)
orig_Ju = deepcopy(objective.val) # Copy the original objective value

#---- Test \patial J/\partial u ----#
# 3D array into which func_deriv_arr gets interpolated
func_deriv_arr = zeros(eqn.q)
func_deriv = zeros(eqn.q_vec)
functional_edges = opts["geom_faces_objective"]
functional_name = EulerEquationMod.FunctionalDict[opts["objective_function"]]

boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, eqn.q, eqn.q_bndry)
EulerEquationMod.calcFunctionalDeriv(mesh, sbp, eqn, opts, functional_name, functional_edges,
                      objective, func_deriv_arr)  # populate df_dq_bndry
assembleSolution(mesh, sbp, eqn, opts, func_deriv_arr, func_deriv)

rand_vec = rand(length(eqn.q_vec))
contract_val = dot(rand_vec,func_deriv)

# Check with finite difference
eqn.q_vec += 1e-6*rand_vec
disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
println("orig_Ju = $orig_Ju")
println("objective.val before = ", objective.val)
EulerEquationMod.evalFunctional(mesh, sbp, eqn, opts, objective)
println("objective.val after = ", objective.val)
dJdu_fd = (objective.val-orig_Ju)/1e-6
eqn.q_vec -= 1e-6*rand_vec


#---- Test \partial R/\partial u ----#
rand_vec = rand(length(eqn.q_vec))
fill!(eqn.res, 0.0)
fill!(eqn.res_vec, 0.0)
res_norm = NonlinearSolvers.calcResidual(mesh, sbp, eqn, opts, evalEuler)
boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, eqn.q, eqn.q_bndry)
res_jac, jacData = EulerEquationMod.calcResidualJacobian(mesh, sbp, eqn, opts)

contract_vec = res_jac*rand_vec

# Check with FD
copy!(eqn.q_vec, orig_q_vec)
for i = 1:length(q_vec)
  eqn.q_vec[i] += 1e-6*rand_vec[i]
end
disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, eqn.q, eqn.q_bndry)
fill!(eqn.res, 0.0)
fill!(eqn.res_vec, 0.0)
res_norm = NonlinearSolvers.calcResidual(mesh, sbp, eqn, opts, evalEuler)
assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
partialRpartialu = (eqn.res_vec - original_res_vec)/1e-6

# Check error
f = open("directionaldRdu.dat", "w")
for i = 1:length(partialRpartialu)
  println(f,real(contract_vec[i] - partialRpartialu[i]))
end
close(f)
println("error partialRpartialu = ", norm(contract_vec - eqn.res_vec,2))

err_vals = readdlm("directionaldRdu.dat")
println("max val = ", maximum(err_vals), " min vals = ", minimum(err_vals))
println("dJdu error = ", norm(dJdu_fd - contract_val, 2))


# Finalize MPI if it hasn't been finalized
if MPI.Initialized()
  MPI.Finalize()
end