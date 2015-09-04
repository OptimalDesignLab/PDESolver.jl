# Artificial Viscosity

function artificialViscosity{Tmsh,Tsol, Tdim}(mesh::AbstractMesh{Tmsh}, 
                             sbp::SBPOperator, eqn::EulerData{Tsol,Tdim})
  
  # Create the Artificial Viscosity flux matrix
  F = zeros(Tsol, mesh.numDofPerNode, sbp.numnodes, mesh.numEl)
  epsilonHat = 0.001
  hi = 1
  h = 1
  
  for i = 1:mesh.numEl
  	for j = 1:mesh.numNodesPerElement
  	  q_vals = view(eqn.q, :, j, i)
  	  calcArtViscosityFluxComp(eqn.params,q_vals,view(F,:,j,i))
  	end
  end

  # We now need to differentiate the above with xi and eta to get the actual
  # artificial viscosity fluxes
  
  # Intermediate variables which have been multiplied with mapping jacobian
  Fav = zeros(Tsol, mesh.numDofPerNode, sbp.numnodes, mesh.numEl, Tdim)
  phi = zeros(Fav)

  for i = 1:mesh.numEl
  	for j = 1:mesh.numNodesPerElement
      phi[:,j,i,1] = mesh.dxidx[1,1,j,i]*F[:,j,i]
      phi[:,j,i,2] = mesh.dxidx[1,2,j,i]*F[:,j,i]
    end
  end

  for k = 1:Tdim
    differentiate!(sbp,k,view(phi,:,:,:,k),view(Fav,:,:,:,1))
  end
  
  phi = zeros(Fav)
  for i = 1:mesh.numEl
    for j = 1:mesh.numNodesPerElement
      phi[:,j,i,1] = mesh.dxidx[2,1,j,i]*F[:,j,i]
      phi[:,j,i,2] = mesh.dxidx[2,2,j,i]*F[:,j,i]
    end
  end
  
  for k = 1:Tdim
    differentiate!(sbp,k,view(phi,:,:,:,k),view(Fav,:,:,:,2))
  end
  
  Fav = epsilonHat*(hi/h)*Fav

  for k = 1:Tdim
    weakdifferentiate!(sbp,k, view(Fav,:,:,:,k), eqn.res, trans=true)
  end


return nothing

end


function calcArtViscosityFluxComp{Tsol}(params::ParamType{2}, q::AbstractArray{Tsol,1},
                                        F::AbstractArray{Tsol,1})
  
  press = calcPressure(q, params)
  F[1] = q[1]
  F[2] = q[2]
  F[3] = q[3]
  F[4] = q[4] + press

return nothing
end

#=
function calcEpsilonHat{Tsol}(params::ParamType{2}, q::AbstractArray{Tsol,1}, 
                           epsilon::AbstractArray{Tsol,1})
  # This function works at the nodal level
  
  # Calculate the temperature
  T = (q[4] - 0.5*(q[2]*q[2] + q[3]*q[3])/q[1])*(1/(q[1]*params.cv))  
  c = sqrt(params.gamma*params.R*T)  # Speed of sound
  p = q[1]*params.R*T  # Pressure
  hBar = 1  # Arithmetic mean of element size metrics

  # Maximum Eigen value of the system
  lambdamax = sqrt((q[2]*q[2] + q[3]*q[3])/(q[1]^2)) + c
  
  thetah = lambdamax*hBar/p
  thetal = 0.01*thetah
  
  # Calculate the elsilonhat value
  if epsilon <= thetal
    epsilonHat = 0
  elseif epsilon > thetal && epsilon <= thetah
    epsilonHat = 0.5*thetah(1 + sin(pi*(epsilon - thetal)/(thetah - thetal) - 0.5))
  else
    epsilonHat = thetah
  end

return nothing
end

function EpsilonPDE{Tmsh,Tsol, Tdim}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator,
                                     eqn::EulerData{Tsol,Tdim})

  # epsilon = zeros(Tsol, 1, sbp.numnodes, mesh.numEl)
  epsilon = zeros(Tsol, sbp.numnodes, mesh.numEl)

return nothing
end

function shockIndicator{Tsol}(params::ParamType{2}, q::AbstractArray{Tsol,1})
# Operates at the Nodal Level

  p = calcPressure(q, params)  # Calculate Pressure

  deltaPsi = 0.5             # Empirical constants that determine when when the
  psi0 = -4 - 4.25*log10(p)  # shock indicator should take effect.

  thetaS = 1  # Maximum Value

  # Calculate F: Discontinuity sensor

  if F <= psi0 - deltaPsi
    Sk = 0
  elseif F >= psi0 + deltaPsi
    Sk = thetaS
  elseif abs(F - psi0) < deltaPsi
    Sk = 0.5*thetaS*(1 + sin(0.5*pi*(F-psi0)/deltaPsi))
  end

end

function boundingBox{Tmsh, Tdim}(coord::AbstractArray{Tmsh,2}, h::AbstractArray{Tmsh,1})
  # It works at the element level. Accepts the vertex coordinates of all the nodes
  # h = zeros(Tmsh,Tdim)  # Stores the dimensions of the bounding box
  for i = 1:Tdim
    h[i] = abs(maximum(coord[i,1:3]) - minimum(coord[i,1:3]))
  end

end

function AverageMeshDimen{Tmsh, Tdim}(mesh::AbstractMesh{Tmsh}, hBar::AbstractArray{Tmsh,1})
  # Operates at the element level
  # Reset the face iterator at the beginning of wherever this loop is called
  element = getFace() 
  (vertices, nvertices) = getDownward(mesh.m_ptr, element, 0)

  #hBar = zeros(Tmsh,nvertices)

  for i = 1:nvertices
    nadjElements = countAdjacent(mesh.m_ptr, vertices, 2) # Get number of adjacent elements
    adjElements = getAdjacent(nadjElements)  # Array of adjacent element pointers
    hArray = zeros(Tmsh, Tdim,nadjElements)  # Array of bounding box of all adjacent elements
    for j = 1:nadjElements  
      elemNum = getNumberJ(mesh.el_Nptr, adjElements, 0, 0) # Get the global element number
      coordinate = mesh.coords[:,:,elemNum]  # get element nodal coordinates
      boundingBox(coordinate, view(hArray,:,j)) # Get the bounding box
    end
    hBar[i] = sum(hArray)/(Tdim*nadjacencies) # Calculate the average mesh metric for the three vertices
  end

end

function AVSourceTerm{Tmsh, Tdim}(mesh::AbstractMesh{Tmsh}, 
                             sbp::SBPOperator, eqn::EulerData{Tsol,Tdim})

epsilon
end =#