# test parts
# include("pressure.jl")

resize!(ARGS, 1)
ARGS[1] = "input_vals_airfoil.jl"


#----  Initialize EulerEquationMod for all the global variables necessary  ----#
include("startup.jl")



#----  Write the surface pressure coefficient for the concerned edge  ----#
#=
g_edges = [3]
# Create data structure for storing coefficient of pressure
nface_arr = zeros(Int,length(g_edges))
for i = 1:length(g_edges)
  g_edge_number = g_edges[i] # Extract geometric edge number
  nface_arr[i] = EulerEquationMod.getnFaces(mesh, g_edge_number)
end
pressCoeffArrWrite = Array(Array{Float64,2},length(g_edges))
for i = 1:length(g_edges)
  pressCoeffArrWrite[i] = zeros(Float64, mesh.sbpface.numnodes, nface_arr[i])
end
EulerEquationMod.writeSurfacePressureCoeff(mesh, sbp, eqn, opts, g_edges, pressCoeffArrWrite)
=#


# Create the objective function data object
objective = EulerEquationMod.OptimizationData{Tsol}(mesh, sbp, opts)

#=
#----  Read in the target surface pressure coefficients  ----#
EulerEquationMod.readSurfacePressureCoeff(mesh, sbp, eqn, opts, g_edges,
objective.pressCoeff_obj.targetCp_arr)

# Check if the pressure coefficients are being read in correctly

ctr = 0
for i = 1:length(g_edges)
  for j = 1:size(pressCoeffArrWrite[i],2)
  	for k = 1:size(pressCoeffArrWrite[i],1)
      err = objective.pressCoeff_obj.targetCp_arr[i][k,j] - pressCoeffArrWrite[i][k,j]
      if abs(err) > 1e-12
        ctr += 1
        println("something is wrong")
      end
  	end
  end
end

if ctr == 0
  println("all good")
end
=#

#=
#---- Get the objective function value  ----#
EulerEquationMod.evalFunctional(mesh, sbp, eqn, opts, objective)
println("Objective function name = ", opts["objective_function"])
println("objective.val = $(objective.val)")
=#
#----------------------------- Test derivatives -------------------------------#

# Copy all the original values
# orig_Ju = deepcopy(objective.val) # Copy the original objective value
disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
orig_q = deepcopy(eqn.q)
orig_q_vec = deepcopy(eqn.q_vec)
original_res_vec = copy(eqn.res_vec)

#---- Get the objective function value  ----#
boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, eqn.q, eqn.q_bndry)
EulerEquationMod.evalFunctional(mesh, sbp, eqn, opts, objective)
orig_Ju = deepcopy(objective.val) # Copy the original objective value

#---- Test \patial J/\partial u ----#

#=
# 3D array into which func_deriv_arr gets interpolated
func_deriv_arr = zeros(eqn.q)
func_deriv = zeros(eqn.q_vec)
functional_edges = opts["geom_faces_objective"]
functional_name = EulerEquationMod.FunctionalDict[opts["objective_function"]]

boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, eqn.q, eqn.q_bndry)
EulerEquationMod.calcFunctionalDeriv(mesh, sbp, eqn, opts, functional_name, functional_edges,
                      objective, func_deriv_arr)  # populate df_dq_bndry
assembleArray(mesh, sbp, eqn, opts, func_deriv_arr, func_deriv) # Assemble into func_deriv

rand_vec = rand(length(func_deriv))
contract_val = dot(rand_vec,func_deriv)

# Check with finite difference
eqn.q_vec += 1e-6*rand_vec
disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
println("orig_Ju = $orig_Ju")
println("objective.val before = ", objective.val)
EulerEquationMod.evalFunctional(mesh, sbp, eqn, opts, objective)
println("objective.val after = ", objective.val)
dJdu_fd = (objective.val-orig_Ju)/1e-6

println("dJdu error = ", norm(dJdu_fd - contract_val, 2))
=#

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
eqn.q_vec += 1e-6*rand_vec
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
  println(f,contract_vec[i] - partialRpartialu[i])
end
close(f)
println("error partialRpartialu = ", norm(contract_vec - eqn.res_vec,2))


#=
# Calculate the adjoint vector
adjoint_vec = zeros(Tsol, mesh.numDof)
calcAdjoint(mesh, sbp, eqn, opts, objective, adjoint_vec)

saveSolutionToMesh(mesh, real(adjoint_vec))
writeVisFiles(mesh, "adjoint_field")
=# 
if MPI.Initialized()
  MPI.Finalize()
end
