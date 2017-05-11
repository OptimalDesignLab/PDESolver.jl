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

function dmv_dXv_rev(mesh, sbp, eqn, opts)

  getVertCoords_rev(mesh, sbp)

  return nothing
end
#=
for i = 1:length(mesh.dxidx_bar)
  mesh.dxidx_bar[i] = randn()
end

dmv_dXv_rev(mesh, sbp, eqn, opts)
# Check against finite difference
pert = 1e-6
ctr = 0
orig_dmv_dx = deepcopy(mesh.dxidx)
jacobian = zeros(length(mesh.vert_coords), length(mesh.dxidx))
for i = 1:length

for i = 1# :length(mesh.vert_coords)
  sub_tuple = ind2sub(mesh.vert_coords,i)

  mesh.vert_coords[i] += pert
  for j = 1:mesh.numEl
    update_coords(mesh, j, mesh.vert_coords[:,:,j])
  end
  commit_coords(mesh, sbp)

  dmv_dx_pert_i = deepcopy(mesh.dxidx)
  println("dmv_dx_pert_i = $dmv_dx_pert_i")
  deriv = (dmv_dx_pert_i - orig_dmv_dx)/pert
  println("deriv = $deriv")
  println("orig_dmv_dx = $orig_dmv_dx")
  val = dot(vec(mesh.dxidx_bar), vec(deriv))

  error = norm(val - mesh.vert_coords_bar[i], 2)

  if error > 1e-5
    println("val = $(val), mesh.vert_coords_bar[$i] =$(mesh.vert_coords_bar[i])")
    println("vert_coords = $(mesh.vert_coords[:, sub_tuple[2], sub_tuple[3]])")
  end

  mesh.vert_coords[i] -= pert
end
=#

Tmsh, Tsol, Tres = EulerEquationMod.getTypeParameters(mesh, eqn)
objective = EulerEquationMod.createObjectiveFunctionalData(mesh, sbp, eqn, opts)
EulerEquationMod.evalFunctional(mesh, sbp, eqn, opts, objective)
orig_lift_val = objective.lift_val


# Check reverse mode
fill!(mesh.dxidx_bndry_bar, 0.0)
fill!(mesh.dxidx_bar, 0.0)
fill!(mesh.vert_coords_bar, 0.0)
vert_coords_bar2 = dJdXv_rev(mesh, sbp, eqn, opts, objective)

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
    # println("error at $i dof, = $error")
    println("dXv_fd[$i] = $(dXv_fd[i]), vert_coords_bar2[$i] =$(mesh.vert_coords_bar[i]), objective.lift_val = $(objective.lift_val)")
    # println("vert_coords = $(mesh.vert_coords[:, sub_tuple[2], sub_tuple[3]])")
    # println("vert radius = $(norm(mesh.vert_coords[:, sub_tuple[2], sub_tuple[3]]))")
    ctr += 1
  end

  # mesh.vert_coords[i] -= pert
  volNodes[i] -= pert

end
println("ctr = $ctr")
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
