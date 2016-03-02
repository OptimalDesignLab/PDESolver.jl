# Calculate boundary "forces" in advection

function calc_numerical_bndryforces()

  # Specify the boundary conditions for the edge on which the force needs to be computed 
  # separately in the input dictionary. Use that boundary number to access the boundary 
  # offset array. Then proceed the same as bndryflux to get the forces using 
  # boundaryintegrate!

  g_edge_number = 1 # Geometric boundary edge on which the force needs to be computed
  start_index = mesh.bndry_offsets[g_edge_number]
  end_index = mesh.bndry_offsets[g_edge_number+1]
  bndry_facenums = view(mesh.bndryfaces, start_index:(end_index - 1)) # faces on geometric edge i
  

  nfaces = length(bndry_facenums)
  boundary_press = zeros(Tsol, Tdim, sbp.numfacenodes, nfaces)
  boundary_force = zeros(Tsol, Tdim, sbp.numfacenodes, nfaces)
  analytical_force = zeros(Tsol, sbp.numfacenodes, nfaces)

  for i = 1:nfaces
  	bndry_i = bndry_facenums[i]
  	for j = 1:sbp.numfacenodes
      k = sbp.facenodes[j, bndry_i.face]
      q = eqn.q[1,k,bndryi.element]
      x = view(mesh.coords, :, k, bndry_i.element)
      dxidx = view(mesh.dxidx, :, :, k, bndry_i.element)
      nrm = view(sbp.facenormal, :, bndry_i.face)
      nx = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
      ny = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]
  	end
  end

  return nothing
end