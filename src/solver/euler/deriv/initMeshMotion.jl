# Mesh Warp and FFD

# Step 1: FFD Initialize
  #   Free Form deformation parameters
  ndim = 2
  order = [2,2,2]  # Order of B-splines in the 3 directions
  nControlPts = [2,2,2]
  Tmsh = Float64 # Complex128

  # Create Mapping object
  ffd_map = PumiMapping{Tmsh}(ndim, order, nControlPts, mesh, full_geom=false,
                              geom_faces=geom_faces)

  # Create knot vector
  calcKnot(ffd_map)
  println("\nlength controlPoint = $(length(ffd_map.cp_xyz))")

  # Create Bounding box
  offset = [0., 0., 0.5] # No offset in the X & Y direction
  ffd_box = PumiBoundingBox{Tmsh}(ffd_map, mesh, sbp, offset)

  # Control points
  controlPoint(ffd_map, ffd_box)

  # Populate map.xi
  calcParametricMappingNonlinear(ffd_map, ffd_box, mesh, geom_faces)

MPI.Barrier(mesh.comm)

# Step 2: Prep MeshWarping
volNodes, wallCoords, nwall_faces, param, mpiVar = initPumiWarp(mesh, geom_faces)

# Populate entries of param
param.aExp = 2.0
param.bExp = 2.0
param.LdefFact = 1.0
param.alpha = 0.2
param.symmTol = 1e-4
param.errTol = 1e-4
param.cornerAngle = 30.0
param.zeroCornerRotations = false
param.useRotations = true
param.evalMode = 5
param.bucket_size = convert(Int32,8)

symmetryPlanes = zeros(Float64,3,2)
flatWarpSurfPts = reshape(wallCoords, 3*size(wallCoords,2))
faceSizes = 4*ones(Int32,sum(nwall_faces)) # An extruded 2D edge for a linear
                                           # discretization will result in a
                                           # "element face" with 4 vertiecs.

initializeParameters(param, mpiVar, symmetryPlanes)
initializeWarping(param, volNodes, wallCoords, faceSizes, mesh.comm)
