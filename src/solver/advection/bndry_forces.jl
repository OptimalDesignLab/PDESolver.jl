# Calculate boundary "forces" in advection

#=
function calcForceErrorNorm(mesh, sbp, eqn, g_edge_number)

  

  return norm_force_error
end
=#
function calcBndryforces{Tmsh, Tsol}(mesh::AbstractMesh{Tmsh},sbp::AbstractSBP,
                         eqn::AdvectionData{Tsol}, opts, g_edge_number)

  # Specify the boundary conditions for the edge on which the force needs to be computed 
  # separately in the input dictionary. Use that boundary number to access the boundary 
  # offset array. Then proceed the same as bndryflux to get the forces using 
  # boundaryintegrate!

  # println("mesh.bndry_offsets = \n", mesh.bndry_offsets)
  start_index = mesh.bndry_offsets[g_edge_number]
  # println("start_index = ", start_index)
  end_index = mesh.bndry_offsets[g_edge_number+1]
  # println("end_index = ", end_index)
  bndry_facenums = view(mesh.bndryfaces, start_index:(end_index - 1)) # faces on geometric edge i
  # println("bndry_facenums = ", bndry_facenums)

  nfaces = length(bndry_facenums)
  boundary_press = zeros(Tsol, 1, sbp.numfacenodes, nfaces)
  boundary_force = zeros(Tsol, 1, sbp.numnodes, mesh.numEl)
  analytical_force = zeros(Tsol, sbp.numfacenodes, nfaces)

  for i = 1:nfaces
  	bndry_i = bndry_facenums[i]
  	for j = 1:sbp.numfacenodes
      k = sbp.facenodes[j, bndry_i.face]
      q = eqn.q[1,k,bndry_i.element]
      #println("\nq = ", q)
      x = view(mesh.coords, :, k, bndry_i.element)
      #println("coordinates = ", x)
      dxidx = view(mesh.dxidx, :, :, k, bndry_i.element)
      #println("dxidx = \n", round(dxidx,3))
      nrm = view(sbp.facenormal, :, bndry_i.face)
      #println("sbp.facenormal = ", nrm)
      alpha_x = eqn.alpha_x[1, k, bndry_i.element]
      alpha_y = eqn.alpha_y[1, k, bndry_i.element]
      nx = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
      ny = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]
      #println("normal = ", round([nx, ny],3))
      # analytical_force[j,i] = calcAnalyticalForce(alpha_x, alpha_y, [nx, ny], x)
      #println("bndry_i.element = ", bndry_i.element)
      boundary_press[1,j,i] = (alpha_x*nx + alpha_y*ny)*q # Boundary Flux
      # boundary_press[1,j,i] = 1
      #println("boundary_press[1,j,i] = ", boundary_press[1,j,i], '\n')
  	end
  end

  boundaryintegrate!(sbp, mesh.bndryfaces[start_index:(end_index - 1)], boundary_press, boundary_force)

  #=
  for (bindex, bndry) in enumerate(mesh.bndryfaces[start_index:(end_index - 1)])
    for i = 1:sbp.numnodes
      @printf("boundary_force[1,%d,%d] = %f\n", i, bndry.element, boundary_force[1,i,bndry.element])
    end
    println("\nboundary_press[1,:,bndry.element] = ", boundary_press[1,:,bindex])
  end # end enumerate
  =#
  # Add all boundary_force nodal values along the edge to get the nodal force value
  numerical_force = 0.0::Tsol
  for (bindex, bndry) in enumerate(mesh.bndryfaces[start_index:(end_index - 1)])
    for i = 1:sbp.numfacenodes
      k = sbp.facenodes[i, bndry.face]
      numerical_force += boundary_force[1,k,bndry.element]
    end
  end  # end enumerate

  println("numerical_force = ", numerical_force)
  # Compare the residual of the solution
  # error_force = zeros(boundary_force)


  return numerical_force
end

function calcAnalyticalForce(alpha_x, alpha_y, nrm, coords)

  # specialized function to do integration over edge 1 (Integration along dx)
  # for exp(x+y)

  x = coords[1]
  y = coords[2]
  # int_value = (alpha_x*nrm[1] + alpha_y*nrm[2])*exp(x+y)
  int_value = (alpha_x*nrm[1] + alpha_y*nrm[2])*((x^6)/6 + x*(y^5))

  return int_value
end