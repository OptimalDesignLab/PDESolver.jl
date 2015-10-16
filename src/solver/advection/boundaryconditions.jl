# boundaryconditions.jl

function flux1(u_sbp_, dxi_dx_, nrm, mesh::PumiMesh2)
    # u_sbp_ is the entry from u_sbp for this node
    # dxi_dx is the jacobian for this node
    # nrm is the normal vector

  u_xi = mesh.dxidx[1,1]*alpha_x + mesh.dxidx[1,2]*alpha_y
  u_eta = mesh.dxidx[2,1]*alpha_x + mesh.dxidx[2,2]*alpha_y
  mag_x = u_xi*nrm[1]  # x flow magnitude
  mag_y = u_eta*nrm[2]  # y flow magnitude
  mag = mag_x + mag_y

  # because mag is scalar,not vector, can combine these if statements
  if mag < 0  # inflow condition, added factor of nrm[2]
    flx_xi = u_xi*u_bc
  else
    flx_xi = u_xi*u_sbp_
  end

  if mag < 0  # inflow condition, added factor of nrm[1]
    flx_eta = u_eta*u_bc
  else
    flx_eta = u_eta*u_sbp_
  end

  net_flux = flx_xi*nrm[1] + flx_eta*nrm[2]
  cntr += 1
  return net_flux
  
end # end function flux1