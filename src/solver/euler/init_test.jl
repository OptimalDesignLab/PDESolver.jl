# test parts
include("input_vals_airfoil.jl")
# include("pressure.jl")

resize!(ARGS, 1)
ARGS[1] = "input_vals_vortex_parallel.jl"

include("startup.jl")

g_edges = [3]
nface_arr = zeros(Int,length(g_edges))

for i = 1:length(g_edges)
  g_edge_number = g_edges[i] # Extract geometric edge number
  # get the boundary array associated with the geometric edge
  itr2 = 0
  for itr2 = 1:mesh.numBC
    if findfirst(mesh.bndry_geo_nums[itr2],g_edge_number) > 0
      break
    end
  end

  start_index = mesh.bndry_offsets[itr2]
  end_index = mesh.bndry_offsets[itr2+1]
  idx_range = start_index:(end_index-1)
  bndry_facenums = sview(mesh.bndryfaces, idx_range) # faces on geometric edge i
  nface_arr[i] = length(bndry_facenums)
end

pressCoeffArrWrite = Array(Array{Float64,2},length(g_edges))
pressCoeffArrRead = Array(Array{Float64,2},length(g_edges))
for i = 1:length(g_edges)
  pressCoeffArrWrite[i] = zeros(Float64, mesh.sbpface.numnodes, nface_arr[i])
  pressCoeffArrRead[i] = zeros(Float64, mesh.sbpface.numnodes, nface_arr[i])
end

EulerEquationMod.writeSurfacePressureCoeff(mesh, sbp, eqn, opts, g_edges, pressCoeffArrWrite)
# EulerEquationMod.evalFunctional(mesh, sbp, eqn, opts)

# Create data structure for storing coefficient of pressure
# Step 1: get nfaces on every geometric edge



EulerEquationMod.readSurfacePressureCoeff(mesh, sbp, eqn, opts, g_edges, pressCoeffArrRead)

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

if MPI.Initialized()
  MPI.Finalize()
end