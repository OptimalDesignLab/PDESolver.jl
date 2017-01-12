# testing derivatives w.r.t states
resize!(ARGS, 1)
ARGS[1] = "input_vals_airfoil.jl"


#----  Initialize EulerEquationMod for all the global variables necessary  ----#
include("startup.jl")

# Create the objective function data object
objective = EulerEquationMod.OptimizationData{Tsol}(mesh, sbp, opts)

# Copy all the original values
disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
orig_q = deepcopy(eqn.q)
orig_q_vec = deepcopy(eqn.q_vec)
original_res_vec = copy(eqn.res_vec)

#---- Test \partial R/\partial u ----#
rand_vec = rand(length(eqn.q_vec))
fill!(eqn.res, 0.0)
fill!(eqn.res_vec, 0.0)
res_norm = NonlinearSolvers.calcResidual(mesh, sbp, eqn, opts, evalEuler)
boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, eqn.q, eqn.q_bndry)
res_jac, jacData = EulerEquationMod.calcResidualJacobian(mesh, sbp, eqn, opts)

contract_vec = res_jac*rand_vec

# Check with FD
eqn.q_vec = copy(orig_q_vec)
# copy!(eqn.q_vec, orig_q_vec)
eqn.q_vec += 1e-6*rand_vec
#=
for i = 1:length(q_vec)
  eqn.q_vec[i] += 1e-6*rand_vec[i]
end=#
# disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
# boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, eqn.q, eqn.q_bndry)
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

# Finalize MPI if it hasn't been finalized
if MPI.Initialized()
  MPI.Finalize()
end