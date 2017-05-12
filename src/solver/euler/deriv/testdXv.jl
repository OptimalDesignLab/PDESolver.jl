# Check dxv test for functionals

resize!(ARGS, 1)
# ARGS[1] = "../input_vals_test_reversemode.jl"
ARGS[1] = "../input_vals_vortex.jl"

push!(LOAD_PATH, joinpath(Pkg.dir("FFD.jl"), "src"))

using FreeFormDeformation
using MeshMovement
using PdePumiInterface
import ODLCommonTools.sview

# Get the adjoint vector
include("../startup.jl")
# EulerEquationMod.dataPrep(mesh, sbp, eqn, opts)

function dJdXv_rev(mesh, sbp, eqn, opts, objective)

  EulerEquationMod.evalFunctional_revm(mesh, sbp, eqn, opts, objective, "lift")

  # Checking interpolation and entities from PumiInterface
  interpolateMapping_rev(mesh)
  getVertCoords_rev(mesh, sbp)
  vert_coords_bar2 = zeros(mesh.dim, mesh.numVert)
  PdePumiInterface.accumulateAtVerts(mesh, mesh.vert_coords_bar, vert_coords_bar2)

  return vert_coords_bar2
end

function bypassMeshMovement{Tmsh}(mesh::AbstractMesh{Tmsh}, sbp, eqn, opts, geom_faces, vert_coords_bar2)

  vtx_per_face = mesh.dim # Only for simplex elements
  nwall_faces = MeshMovement.getnWallFaces(mesh, geom_faces)
  nWallCoords = sum(nwall_faces)*vtx_per_face
  # wallCoords = zeros(3, nWallCoords)
  Xs_bar = zeros(Tmsh, 3, nWallCoords) # This has to be 3D always since it gets used by FFD
  ctr = 1
  for itr = 1:length(geom_faces)
    geom_face_number = geom_faces[itr]
    # get the boundary array associated with the geometric edge
    itr2 = 0
    for itr2 = 1:mesh.numBC
      if findfirst(mesh.bndry_geo_nums[itr2],geom_face_number) > 0
        break
      end
    end
    start_index = mesh.bndry_offsets[itr2]
    end_index = mesh.bndry_offsets[itr2+1]
    idx_range = start_index:(end_index-1)
    bndry_facenums = sview(mesh.bndryfaces, idx_range) # faces on geometric edge i
    nfaces = length(bndry_facenums)
    for i = 1:nfaces
      bndry_i = bndry_facenums[i]
      vtx_arr = mesh.topo.face_verts[:,bndry_i.face]
      for j = 1:length(vtx_arr)
        local_vertnum = mesh.element_vertnums[vtx_arr[j],bndry_i.element]
        for k = 1:mesh.dim
          Xs_bar[k, ctr] = vert_coords_bar2[k, local_vertnum]
          # println("radius = ", norm(mesh.vert_coords[:,vtx_arr[j],bndry_i.element],2))
        end
        ctr += 1
      end
    end # End for i = 1:nfaces
  end # End for itr = 1:length(geom_faces)

  return vec(Xs_bar) # Return 1D form of Xs_bar
end # End function bypassMeshMovement


Tmsh, Tsol, Tres = EulerEquationMod.getTypeParameters(mesh, eqn)
objective = EulerEquationMod.createObjectiveFunctionalData(mesh, sbp, eqn, opts)
geom_faces = objective.geom_faces_functional
EulerEquationMod.evalFunctional(mesh, sbp, eqn, opts, objective)
orig_lift_val = objective.lift_val

# FFD initialization
ndim = 2
order = [2,2,2]  # Order of B-splines in the 3 directions
nControlPts = [2,2,2]
bounding_box_offset = [0., 0., 0.5] # No offset in the X & Y direction
ffd_map = initializeFFD(mesh, sbp, order, nControlPts, Tmsh, bounding_box_offset,
                       false, geom_faces)


# Check reverse mode
fill!(mesh.dxidx_bndry_bar, 0.0)
fill!(mesh.dxidx_bar, 0.0)
fill!(mesh.vert_coords_bar, 0.0)
vert_coords_bar2 = dJdXv_rev(mesh, sbp, eqn, opts, objective)
Xs_bar = bypassMeshMovement(mesh, sbp, eqn, opts, geom_faces, vert_coords_bar2)
evaldXdControlPointProduct(ffd_map, mesh, Xs_bar)

dJdx = zeros(Tmsh, length(ffd_map.cp_xyz))
dJdx_sum = zero(Tmsh)
itr = 1
for i = 1:size(ffd_map.work, 4)
  for j = 1:size(ffd_map.work, 3)
    for k = 1:size(ffd_map.work, 2)
      for l = 1:3 # mesh.dim
        dJdx[itr] = ffd_map.work[l,k,j,i]
        itr += 1
      end
      dJdx_sum += ffd_map.work[2,k,j,i]
    end
  end
end

# Check against finite difference
pert = 1e-6

for i = 1:size(ffd_map.work, 4)
  for j = 1:size(ffd_map.work, 3)
    for k = 1:size(ffd_map.work, 2)
        ffd_map.cp_xyz[2,k,j,i] += pert
    end
  end
end
evalSurface(ffd_map, mesh)
for j = 1:mesh.numEl
  update_coords(mesh, j, mesh.vert_coords[:,:,j])
end
commit_coords(mesh, sbp)
EulerEquationMod.evalFunctional(mesh, sbp, eqn, opts, objective)
dJdx_fd= (objective.lift_val - orig_lift_val)/pert
error = norm(dJdx_fd - dJdx_sum, 2)
ctr = 0
if error > 1e-5
  println("dJdx_fd = $(dJdx_fd), dJdx_sum = $(dJdx_sum), error = $error")
  ctr += 1
end
println("ctr = $ctr")
#=
for i = 1:length(ffd_map.cp_xyz)
  ffd_map.cp_xyz[i] += pert
  evalSurface(ffd_map, mesh)
  for j = 1:mesh.numEl
    update_coords(mesh, j, mesh.vert_coords[:,:,j])
  end
  commit_coords(mesh, sbp)
  EulerEquationMod.evalFunctional(mesh, sbp, eqn, opts, objective)
  dJdx_fd_i = (objective.lift_val - orig_lift_val)/pert
  error = norm(dJdx[i] - dJdx_fd_i, 2)
  if error > 1e-5
    println("dJdx[$i] = $(dJdx[i]), dJdx_fd_i = $(real(dJdx_fd_i)), error = $error")
    ctr += 1
  end
  ffd_map.cp_xyz[i] -= pert
end
println("ctr = $ctr, length(ffd_map.cp_xyz) = $(length(ffd_map.cp_xyz))")
=#
#= # This test Passes!!
# Check against finite difference
volNodes = zeros(Tmsh, mesh.dim, mesh.numVert)
for i = 1:mesh.numEl
  for j = 1:size(mesh.vert_coords,2)
    # Get the vertex numbering on the portion of mesh owned by the processor
    local_vertnum = mesh.element_vertnums[j,i]
    volNodes[:, local_vertnum] = mesh.vert_coords[:,j,i] # mesh.element_vertnums
  end
end
pert = 1e-6
ctr = 0
error = zeros(volNodes)
dXv_fd = zeros(volNodes)
for i = 1:length(volNodes)
  volNodes[i] += pert
  sub_tuple = ind2sub(mesh.vert_coords,i)
  # for j = 1:mesh.numEl
  #   mesh.coords[:,:,j] = PdePumiInterface.vertToVolumeCoords(mesh, sbp, mesh.vert_coords[:,:,j])
  # end
  # PdePumiInterface.getMetrics(mesh, sbp)
  # PdePumiInterface.interpolateCoordinatesAndMetrics(mesh)
  # PdePumiInterface.getFaceNormals(mesh, sbp)

  for j = 1:mesh.numEl
    for k = 1:size(mesh.vert_coords,2)
      local_vertnum = mesh.element_vertnums[k,j]
      mesh.vert_coords[:,k,j] = volNodes[:, local_vertnum]
    end
    update_coords(mesh, j, mesh.vert_coords[:,:,j])
  end
  commit_coords(mesh, sbp)
  EulerEquationMod.evalFunctional(mesh, sbp, eqn, opts, objective)
  # println("lift_val = $(objective.lift_val), orig_lift_val = $orig_lift_val")
  dXv_fd[i] = (objective.lift_val - orig_lift_val)/pert
  error[i] = norm(dXv_fd[i] - vert_coords_bar2[i],2)

  if error[i] > 1e-5
    println("dXv_fd[$i] = $(dXv_fd[i]), vert_coords_bar2[$i] =$(mesh.vert_coords_bar[i]), objective.lift_val = $(objective.lift_val)")
    ctr += 1
  end

  # mesh.vert_coords[i] -= pert
  volNodes[i] -= pert

end
println("ctr = $ctr")
=#
#=
# Write error vtk
vtk_error_arr = zeros(4, size(error,2), size(error,3))
for i = 1:size(error,3)
  for j = 1:size(error, 2)
    for k = 1:size(error, 1)
      vtk_error_arr[k,j,i] = error[k,j,i]
    end
  end
end

PdePumiInterface.saveNodalSolution(mesh, vtk_error_arr)
writeVisFiles(mesh, "error_plot")
=#

#=
# THIS TEST PASSES!!!
coords_bar = zeros(mesh.coords)
for i = 1:length(coords_bar)
  coords_bar[i] = randn()
end
PdePumiInterface.volumeToVertCoords_rev(mesh, sbp, mesh.vert_coords_bar, coords_bar)

# Check against finite difference
coords_orig = zeros(mesh.coords)
for i = 1:mesh.numEl
  coords_orig[:,:,i] = PdePumiInterface.vertToVolumeCoords(mesh, sbp, mesh.vert_coords[:,:,i])
end
pert = 1e-6
ctr = 0
coords_pert = zeros(mesh.coords)
# coords_bar = reshape(coords_bar, length(coords_bar))
for i = 1:length(mesh.vert_coords)
  mesh.vert_coords[i] += pert
  for j = 1:mesh.numEl
    coords_pert[:,:,j] = PdePumiInterface.vertToVolumeCoords(mesh, sbp, mesh.vert_coords[:,:,j])
  end
  dXcoord_dVertcoords_i = (coords_pert - coords_orig)/pert

  prod = dot(vec(dXcoord_dVertcoords_i), vec(coords_bar))
  error = norm(mesh.vert_coords_bar[i] - prod, 2)
  if error > 1e-5
    # println("i = $i, error = $error")
    println("i = $i, prod = $prod, mesh.vert_coords_bar[$i] = $(mesh.vert_coords_bar[i])")
    ctr += 1
  end
  mesh.vert_coords[i] -= pert
end
println("ctr = $ctr, length(mesh.vert_coords) = $(length(mesh.vert_coords))")
=#
