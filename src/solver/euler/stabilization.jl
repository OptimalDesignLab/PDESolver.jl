# this file contains all functions related to performing edge stabilization
# edge stabilization is is executed from euler.jl

# this function is going to be deprecated soon
function stabscale{T}(u::AbstractArray{T,1}, dxidx::AbstractArray{T,2}, nrm::AbstractArray{T,1}, mesh::AbstractMesh, eqn::EulerEquation)

#     println("==== entering stabscale ====")

    # grabbing conserved variables
    rho = u[1]
    vel_x = u[2]/rho
    vel_y = u[3]/rho
    Energy = u[4]

    # from JC's code below, eqn should still be in scope
    pressure = calcPressure(u, eqn)

    # solved eqn for e: E = rho*e + (1/2)*rho*u^2
    vel_squared = vel_x^2 + vel_y^2
    energy = Energy/rho - (1/2)*vel_squared

    # gamma stored in EulerEquation type
    gamma = eqn.gamma

#     println("pressure: ",pressure)
#     println("gamma: ",gamma)
#     println("rho: ",rho)
    # ideal gas law
    speed_sound = sqrt((gamma*pressure)/rho)

    # choice for edge stabilization constant: 
    #   refer to email from JH, 20150504:
    #   Anthony: there is little guidance in the literature regarding 
    #     gamma for the Euler equations.  I suggest starting with 
    #     gamma = 0.01.  If that fails (with a cfl of 1.0), then decrease 
    #     it by an order of magnitude at at time until RK is stable.  
    #     Once you find a value that works, try increasing it slowly.
    edge_stab_gamma = -0.01  # default
#     edge_stab_gamma = 0.0 
#     edge_stab_gamma = 0.00001

    # edge lengths component wise
    h_x = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
    h_y = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]

    # edge length
    h = sqrt(h_x^2 + h_y^2)

    # scaled velocity scalar
#     U = vel_x*(nrm[1]/h) + vel_y*(nrm[2]/h)
    U = vel_x*(h_x/h) + vel_y*(h_y/h)

#     return (U + speed_sound)*edge_stab_gamma*h^2
    return (abs(U) + speed_sound)*edge_stab_gamma*h^2

  end


function stabscale{Tmsh, Tsbp, Tsol}(u::AbstractArray{Tsol,1}, dxidx::AbstractArray{Tmsh,2}, nrm::AbstractArray{Tsbp,1}, eqn::EulerEquation{Tsol} )
# calculate stabscale for a single node

#     println("==== entering stabscale ====")

    # grabbing conserved variables
    rho = u[1]
    vel_x = u[2]/rho
    vel_y = u[3]/rho
    Energy = u[4]

    # from JC's code below, eqn should still be in scope
    pressure = calcPressure(u, eqn)

    # solved eqn for e: E = rho*e + (1/2)*rho*u^2
    vel_squared = vel_x^2 + vel_y^2
    energy = Energy/rho - (1/2)*vel_squared

    # gamma stored in EulerEquation type
    gamma = eqn.gamma

#     println("pressure: ",pressure)
#     println("gamma: ",gamma)
#     println("rho: ",rho)
    # ideal gas law
    speed_sound = sqrt((gamma*pressure)/rho)

    # choice for edge stabilization constant: 
    #   refer to email from JH, 20150504:
    #   Anthony: there is little guidance in the literature regarding 
    #     gamma for the Euler equations.  I suggest starting with 
    #     gamma = 0.01.  If that fails (with a cfl of 1.0), then decrease 
    #     it by an order of magnitude at at time until RK is stable.  
    #     Once you find a value that works, try increasing it slowly.
    edge_stab_gamma = -0.01  # default
#     edge_stab_gamma = 0.0 
#     edge_stab_gamma = 0.00001

    # edge lengths component wise
    h_x = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
    h_y = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]

    # edge length
    h = sqrt(h_x^2 + h_y^2)

    # scaled velocity scalar
#     U = vel_x*(nrm[1]/h) + vel_y*(nrm[2]/h)
    U = vel_x*(h_x/h) + vel_y*(h_y/h)

#     return (U + speed_sound)*edge_stab_gamma*h^2
    return (abs(U) + speed_sound)*edge_stab_gamma*h^2

  end




function stabscale{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator{Tsbp}, eqn::EulerEquation{Tsol})
# calculate stabscale for entire mesh

 nbrnodeindex = Array(sbp.numfacenodes:-1:1)

  for i=1:mesh.numInterfaces
    face_i = mesh.interfaces[i]
    for j=1:sbp.numfacenodes
      iL = sbp.facenodes[j, face_i.faceL]
      iR = sbp.facenodes[nbrnodeindex[j], face_i.faceR]
      q = unsafe_view(eqn.q, :, iL, face_i.elementL)
      dxidx = unsafe_view(mesh.dxidx, :, :, iL, face_i.elementL)
      nrm = unsafe_view(sbp.facenormal, :, face_i.faceL)
      
      eqn.stabscale[j,i] = stabscale(q, dxidx, nrm, eqn)
    end
  end

  return nothing
end
      

