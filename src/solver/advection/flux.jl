# this file contains the defitions of all the fluxes used for DG face integrals
@doc """
### AdvectionEquationMod.calcFaceFlux

  This function calculates the DG flux between a specified set of faces,
  using the solution data at the faces stored in eqn.q_face.
  Note that the flux is negated because the face integrals have a 
  negative sign in the weak form.

  Inputs:
    mesh
    sbp
    eqn
    functor: the functor that calculates the flux at a node
    interfaces: an array of type Interface that specifies which interfaces
                to calculate the flux for

  Inputs/Outputs:
    face_flux: array to store the flux in, numDofPerNode x nnodesPerFace
               x length(interfaces)

  The functor must have the signature:
  func( uL, qR, alpha_x, alpha_y, dxidx, nrm, params)

  where uL and uR are the solution values for a node on the left and right
  elements, alpha_x and alpha_y are the x and y advection velocities,
  dxidx is the scaled mapping jacobian for elementL, and nrm is the face
  normal in reference space.  params is eqn.params

"""->
function calcFaceFlux{Tmsh,  Tsol, Tres}( mesh::AbstractDGMesh{Tmsh}, 
                          sbp::AbstractSBP, eqn::AdvectionData{Tsol}, 
                          functor::FluxType, 
                          interfaces::AbstractArray{Interface,1}, 
                          face_flux::AbstractArray{Tres, 3})

  
  nfaces = length(interfaces)
  for i=1:nfaces  # loop over faces
    interface_i = interfaces[i]
    for j = 1:mesh.sbpface.numnodes
      eL = interface_i.elementL
      fL = interface_i.faceL

      # get components
      qL = eqn.q_face[1, 1, j, i]
      qR = eqn.q_face[1, 2, j, i]
      alpha_x = eqn.alpha_x
      alpha_y = eqn.alpha_y
      dxidx = view(mesh.dxidx_face, :, :, j, i)
      nrm = view(sbp.facenormal, :, fL)

      face_flux[1, j, i] = -functor(qL, qR, alpha_x, alpha_y, dxidx, nrm, eqn.params)
    end
  end

  return nothing
end



@doc """
### AdvectionEquationMod.avgFlux

  This flux function averages the two states to calculate the flux

  Inputs:
    uL: the left state
    uR: the right state
    alpha_x: the advection velocity in the x direction
    alpha_y: the advection velocity in the y direction
    dxidx: the scaled mapping jacobian for elementL
    nrm: the face normal vector for elementL in parametric space

  Outputs:
    the flux
"""->
type avgFlux <: FluxType
end

function call{Tmsh, Tsol}(obj::avgFlux, uL::Tsol, uR::Tsol,
              alpha_x, alpha_y, dxidx::AbstractArray{Tmsh,2}, 
              nrm::AbstractArray{Tmsh,1}, params::ParamType)

  alpha_xi = dxidx[1,1]*alpha_x + dxidx[1,2]*alpha_y
  alpha_eta = dxidx[2,1]*alpha_x + dxidx[2,2]*alpha_y
  alpha_n  = alpha_xi*nrm[1] + alpha_eta*nrm[2]

  u = alpha_n*(uL + uR)*0.5
  return u
end


@doc """
### AdvectionEquationMod.LFFlux

  The Lax-Friedrich flux, using a parameter alpha to control upwinding.
  alpha = 0 -> centered flux
  alpha = 1 -> completely upwinded
  The alpha parameter 

  Inputs:
    uL: the left state
    uR: the right state
    alpha_x: the advection velocity in the x direction
    alpha_y: the advection velocity in the y direction
    dxidx: the scaled mapping jacobian for elementL
    nrm: the face normal vector for elementL in parametric space

  Outputs:
    the flux
"""->
type LFFlux <: FluxType
end

function call{Tmsh, Tsol}(obj::LFFlux, uL::Tsol, uR::Tsol,
              alpha_x, alpha_y, dxidx::AbstractArray{Tmsh,2}, 
              nrm::AbstractArray{Tmsh,1}, params::ParamType)

  alpha_xi = dxidx[1,1]*alpha_x + dxidx[1,2]*alpha_y
  alpha_eta = dxidx[2,1]*alpha_x + dxidx[2,2]*alpha_y
  alpha_n  = alpha_xi*nrm[1] + alpha_eta*nrm[2]
  alpha_LF = params.LFalpha
  u = alpha_n*(uL + uR)*0.5 + absvalue(alpha_n)*(1 - alpha_LF)*0.5*(uL - uR)
  return u
end


global const FluxDict = Dict{ASCIIString, FluxType}(
"avgFlux" => avgFlux(),
"LFFlux" => LFFlux(),
)

function getFluxFunctors(mesh::AbstractDGMesh, sbp, eqn, opts)

  name = opts["Flux_name"]
  eqn.flux_func = FluxDict[name]
  return nothing
end
