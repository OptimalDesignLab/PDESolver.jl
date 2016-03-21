# functions for calculating the flux for interior face integrals

@doc """
### AdvectionEquationMod.calcFaceFlux

  This function calculates the DG flux between a specified set of faces,
  using the solution data at the faces stored in eqn.q_face.
  Note that the flux is negated because the face integrals have a 
  negative sign in the weak form.

  Conservative variables only!

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
function calcFaceFlux{Tmsh,  Tsol, Tres, Tdim}( mesh::AbstractDGMesh{Tmsh}, 
                          sbp::AbstractSBP, 
                          eqn::AdvectionData{Tsol, Tres, Tdim, :conservative}, 
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
      qL = view(eqn.q_face, :, 1, j, i)
      qR = view(eqn.q_face, :, 2, j, i)
      dxidx = view(mesh.dxidx_face, :, :, j, i)
      aux_vars = view(eqn.aux_vars_face, :, j, i)
      nrm = view(sbp.facenormal, :, fL)

      flux_j = view(face_flux, :, j, i)
      functor(qL, qR, aux_vars, dxidx, nrm, flux_j, eqn.params)
      # add the negative sign
      for k=1:mesh.numDofPerNode
        flux_j[k] = -flux_j[k]
      end
    end
  end

  return nothing
end


@doc """
### EulerEquationMod.interpolateFace

  This function interpolates the solution values from the internal nodes 
  to the face flux points of the elements

  Inputs:
    mesh: an AbstractDGMesh
    sbp
    eqn
    opts
    q: a 3D array of solution values at the nodes, numDofPerNode x
       numNodesPerElement x numEl

  Inputs/Outputs:
    q_face: a 4D array of solution values at each interface,
            numDofPerNode x 2 x numfacenodes x numInterface
            q_face[:, 1, j, i] stores the q values for elementL of interface
            i node j and q_face[:, 2, j, i] stores the values for elementR

    eqn.aux_vars_face is also populated
"""->
function interpolateFace(mesh::AbstractDGMesh, sbp, eqn, opts, q::Abstract3DArray, q_face::AbstractArray{Tsol, 4})

  # interpolate solution
  interiorfaceinterpolate!(mesh.sbpface, mesh.interface, q, q_face)

  # recalculte aux_vars
  for i=1:mesh.numInterface
    for j=1:mesh.sbpface.numnodes
      q_vals = view(q_face, :, 1, j, i) # always use elementL
      eqn.aux_vars_face[1, j, i] = calcPressure(q_vals, eqn.params)
    end
  end

  return nothing
end


type RoeFlux <: FluxType
end

#TODO: get rid of flux_parametric, aux_vars
function call{Tsol, Tres, Tmsh}(obj::RoeFlux, uL::AbstractArray{Tsol,1}, 
              uR::AbstractArray{Tsol,1}, 
              aux_vars, dxidx::AbstractArray{Tmsh, 2}, nrm::AbstractVector, 
              F::AbstractVector, params::ParamType)

  RoeSolver(uL, uR, aux_vars, dxidx, nrm, F, params)
end
