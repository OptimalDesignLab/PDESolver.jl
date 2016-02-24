# Calculate the analytical force on the inner boundary of the isentropic vortex

function calc_analytical_forces{Tmsh}(mesh::AbstractMesh{Tmsh}, params::ParamType{2},
	                                  coords::AbstractArray{Tmsh})

  q = zeros(mesh.numDofPerNode)
  calcIsentropicVortex(coords, params, q)  # Get analytical q ath the coordinates
  p = calcPressure(params, q) # Get the analytical pressure
  r = sqrt(coords[1]*coords[1] + coords[2]*coords[2]) # get the curcular radius
  force = 0.5*pi*r*p

  return force
end

function calc_numerical_forces{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh}, 
                               sbp::SBPOperator, 
                               eqn::EulerData{Tsol, Tres, Tdim}, opts)

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
      q = view(eqn.q, :, k, bndry_i.element)
      aux_vars = view(eqn.aux_vars, :, k, bndry_i.element)
      x = view(mesh.coords, :, k, bndry_i.element)
      dxidx = view(mesh.dxidx, :, :, k, bndry_i.element)
      nrm = view(sbp.facenormal, :, bndry_i.face)

      analytical_force[k,bndry_i.element] = calc_analytical_forces(mesh, eqn.params, x)
      nx = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
      ny = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]

      # Calculate euler flux for the current iteration
      euler_flux = zeros(Tsol, mesh.numDofPerNode)
      calcEulerFlux(eqn.params, q, aux_vars, [nx, ny], euler_flux)
      
      # Boundary pressure in "ndimensions" direcion
      boundary_press[:,j,i] =  euler_flux[2:3]
    end # end for j = 1:sbp.numfacenodes
  end   # end for i = 1:nfaces
  boundaryintegrate!(sbp, bndry_facenums, boundary_press, boundary_force)


  return nothing
end

function calcForceError()

  #--- Getting nodal errors
  numfacenodes, nfaces, ndimensions = size(boundary_force)
  force_error_mag = zeros(analytical_force)

  for i = 1:nfaces
    for j = 1:numfacenodes
      boundary_force_mag = sqrt(boundary_force[j,i,1]*boundary_force[j,i,1] +
                                boundary_force[j,i,2]*boundary_force[j,i,2])
      force_error_mag[j,i] = analytical_force[j,i] - boundary_force_mag
    end
  end


  #--- Calculate the integral norm of the force error

  # Get a 1D force error array


  M = Array(Tmsh, 1)
  # Get the corresponding mass matrix
  for i = 1:nfaces
    for j = 1:numfacenodes
      for k = 1:mesh.numDofPerNode
        dofnum_k = mesh.dofs[k,j,i]
        M[dofnum_k] += (sbp.w[j]/mesh.jac[j,i])
      end
    end
  end

  return nothing
end

function calcPhysicalEulerFlux{Tsol}(params::ParamType{2}, q::AbstractArray{Tsol,1}, 
                               F::AbstractArray{Tsol, 2})

  u = q[2]/q[1]
  v = q[3]/q[1]
  p = calcPressure(params, q)

  # Calculate Euler Flux in X-direction
  F[1,1] = q[2]
  F[2,1] = q[2]*u + p
  F[3,1] = q[2]*v
  F[4,1] = u*(q[4] + p)

  # Calculate Euler Flux in Y-direction

  F[1,2] = q[3]
  F[2,2] = q[3]*u
  F[3,2] = q[3]*v + p
  F[4,2] = v*(q[4] + p)

  return nothing
end