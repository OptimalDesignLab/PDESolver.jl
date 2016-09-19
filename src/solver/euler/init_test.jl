# test parts
include("input_vals_airfoil.jl")
# include("pressure.jl")

resize!(ARGS, 1)
ARGS[1] = "input_vals_vortex.jl"

include("startup.jl")
# using PressureMod

g_edges = [3]
nface_arr = zeros(Int,length(g_edges))

for i = 1:length(g_edges)
  g_edge_number = g_edges[i] # Extract geometric edge number
  nface_arr[i] = getnFaces(mesh, g_edge_number)
end

pressCoeffArrWrite = Array(Array{Float64,2},length(g_edges))
for i = 1:length(g_edges)
  pressCoeffArrWrite[i] = zeros(Float64, mesh.sbpface.numnodes, nface_arr[i])
end

pressureInfo = EulerEquationMod.PressureData{Float64}(mesh, g_edges, nface_arr)

EulerEquationMod.writeSurfacePressureCoeff(mesh, sbp, eqn, opts, g_edges, pressCoeffArrWrite)

# Create data structure for storing coefficient of pressure
# Step 1: get nfaces on every geometric edge

EulerEquationMod.readSurfacePressureCoeff(mesh, sbp, eqn, opts, g_edges, pressInfo.targetCp_arr)

#=
ctr = 0
for i = 1:length(g_edges)
  for j = 1:size(pressCoeffArrRead[i],2)
  	for k = 1:size(pressCoeffArrRead[i],1)
      err = pressCoeffArrRead[i][k,j] - pressCoeffArrWrite[i][k,j]
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

EulerEquationMod.evalFunctional(mesh, sbp, eqn, opts)

if MPI.Initialized()
  MPI.Finalize()
end
