# Mesh Warp and FFD

# Step 1: FFD Initialize
  #   Free Form deformation parameters
  ndim = 2
  order = [4,4,2]  # Order of B-splines in the 3 directions
  nControlPts = [6,6,2]
  Tmsh = Float64 # Complex128

  # Create Mapping object
  ffd_map = PumiMapping{Tmsh}(ndim, order, nControlPts, mesh, full_geom=false,
                              geom_faces=geom_faces)

  # Create knot vector
  calcKnot(ffd_map)

  # Create Bounding box
  offset = [0., 0., 0.5] # No offset in the X & Y direction
  ffd_box = PumiBoundingBox{Tmsh}(ffd_map, mesh, sbp, offset)

  # Control points
  controlPoint(ffd_map, ffd_box)

  # Populate map.xi
  calcParametricMappingNonlinear(ffd_map, ffd_box, mesh, geom_faces)

MPI.Barrier(mesh.comm)

# Step 2: Prep MeshWarping
volNodes, surfaceVtx, wallCoords, nwall_faces, param, mpiVar = initPumiWarp(mesh, geom_faces)
# Populate entries of param
param.aExp = 2.
param.bExp = 2.
param.LdefFact = 1
param.alpha = 0.2
param.symmTol = 1e-4
param.errTol = 1e-4
param.cornerAngle = 30.
param.zeroCornerRotations = false
param.useRotations = false
param.evalMode = 5

symmetryPlanes = zeros(Float64,3,2)
flatWarpSurfPts = reshape(wallCoords, 3*size(wallCoords,2))
faceSizes = mesh.dim*ones(Int32,sum(nwall_faces))

initializeParameters(param, mpiVar, symmetryPlanes)
initializeWarping(param, volNodes, wallCoords, faceSizes, mesh.comm)
