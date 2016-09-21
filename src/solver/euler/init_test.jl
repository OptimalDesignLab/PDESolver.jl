# test parts
include("input_vals_airfoil.jl")
# include("pressure.jl")

resize!(ARGS, 1)
ARGS[1] = "input_vals_vortex_parallel.jl"


#----  Initialize EulerEquationMod for all the global variables necessary  ----#
include("startup.jl")



#----  Write the surface pressure coefficient for the concerned edge  ----#
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


#----  Read in the target surface pressure coefficients  ----#
# Create the objective function data object
objective = EulerEquationMod.OptimizationData{Tsol}(mesh, sbp, opts)
EulerEquationMod.readSurfacePressureCoeff(mesh, sbp, eqn, opts, g_edges,
objective.pressCoeff_obj.targetCp_arr)

# Check if the pressure coefficients are being read in correctly
#=
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

#---- Get the objective function value  ----#
EulerEquationMod.evalFunctional(mesh, sbp, eqn, opts, objective, is_objective_fn=true)
println("objective.val = $(objective.val)")

if MPI.Initialized()
  MPI.Finalize()
end
