# bc_solvers.jl for advection

@doc """
flux1

Calculates the boundary flux for the advection equation. It works at the nodal
level.

**Inputs**

*  `u_sbp_`: The entry from u_sbp for this node
*  `dxidx` : The jacobian for this node
*  `nrm`   : nrm is the normal vector
*  `net_flux`:
*  `params`: the equation ParamType

**Outputs**

*  None
"""->

function flux1(u_sbp_, dxidx, nrm, net_flux, params::ParamType2)
  # This function works at the nodal level  
  # u_sbp_ is the entry from u_sbp for this node
  # dxi_dx is the jacobian for this node
  # nrm is the normal vector

  u_xi = dxidx[1,1]*alpha_x + dxidx[1,2]*alpha_y
  u_eta = dxidx[2,1]*alpha_x + dxidx[2,2]*alpha_y
  mag_x = u_xi*nrm[1]  # x flow magnitude
  mag_y = u_eta*nrm[2]  # y flow magnitude
  mag = mag_x + mag_y
  
  u_bc = sin(-1)  # boundary condition (for all sides)
  # println("u_bc = ", u_bc)

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
  
  return nothing
end # end function flux1

@doc """
### AdvectionEquationMod.RoeSolver

Roe solver for the advection equations. It determines the boundary flux on 
each boundary. It is called at the nodal level

**Inputs**

*  `u`    : Solution of advection equation at a particular node
*  `u_bc` : Prescribed solution value at the boundary
*  `params`: the equation ParamType object
*  `nrm`  : Summation-by-parts face normal vector
*  `dxidx`: Mapping jacobian at a particular node

**Outputs**

*  `bndryflux` : Boundary flux at the particular node

"""->
function RoeSolver{Tsol, Tmsh}(u::Tsol, u_bc, params::ParamType2, nrm, 
                               dxidx::AbstractArray{Tmsh,2})
    alpha_x = params.alpha_x
    alpha_y = params.alpha_y
    alpha_xi = dxidx[1,1]*alpha_x + dxidx[1,2]*alpha_y
    alpha_eta = dxidx[2,1]*alpha_x + dxidx[2,2]*alpha_y
    alpha_n  = alpha_xi*nrm[1] + alpha_eta*nrm[2]
    bndryflux = 0.5*alpha_n.*(u_bc + u) - 0.5*absvalue(alpha_n).*(u_bc - u)

  return bndryflux
end # end function RoeSolver

function RoeSolver{Tsol, Tmsh}(u::Tsol, u_bc, params::ParamType3, nrm, 
                               dxidx::AbstractArray{Tmsh,2})
    alpha_x = params.alpha_x
    alpha_y = params.alpha_y
    alpha_z = params.alpha_z
    alpha_xi = dxidx[1,1]*alpha_x + dxidx[1,2]*alpha_y + dxidx[1,3]*alpha_z
    alpha_eta = dxidx[2,1]*alpha_x + dxidx[2,2]*alpha_y + dxidx[2,3]*alpha_z
    alpha_psi = dxidx[3,1]*alpha_x + dxidx[3,2]*alpha_y + dxidx[3,3]*alpha_z

    alpha_n  = alpha_xi*nrm[1] + alpha_eta*nrm[2] + alpha_psi*nrm[3]
    bndryflux = 0.5*alpha_n.*(u_bc + u) - 0.5*absvalue(alpha_n).*(u_bc - u)

  return bndryflux
end # end function RoeSolver
