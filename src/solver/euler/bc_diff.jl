# functions for computing the jacobians of the boundary conditions calculated
# in bc.jl

"""
  Computes the jacobian of the boundary condition terms computed by
  [`getBCFluxes`](@ref) and assembles them into the Jacobian

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * assembler: AssembleElementData
"""
function getBCFluxes_diff(mesh::AbstractMesh, sbp::AbstractSBP, eqn::EulerData,
                          opts, assembler::AssembleElementData)
  #get all the fluxes for all the boundary conditions and save them in eqn.bndryflux

  #println("mesh.bndry_funcs = ", mesh.bndry_funcs)
  for i=1:mesh.numBC
  #  println("computing flux for boundary condition ", i)
    functor_i = mesh.bndry_funcs[i]
    start_index = mesh.bndry_offsets[i]
    end_index = mesh.bndry_offsets[i+1]
    idx_range = start_index:(end_index - 1)
    bndry_facenums_i = sview(mesh.bndryfaces, start_index:(end_index - 1))

    calcBoundaryFlux_nopre_diff(mesh, sbp, eqn, functor_i, idx_range, bndry_facenums_i, assembler)
  end

  return nothing
end


"""
  Computes the jacobian of the boundary condition terms for a given
  BC and assembles it into the Jacobian.  See [`calcBoundaryFlux_nopre`](@ref).
  Note that this function does not depend on any quatiites being precomputed.

  **Inputs**

   * mesh
   * sbp
   * eqn
   * functor: this is the regular BCType functor, not BCType_diff
   * idx_range: range of indices in mesh.bndryfaces that bndry_facenums
                came from
   * bndry_facenums: array of Boundary objects to compute the BC for
   * assembler: AssembleElementData needed to assemble the element level
                jacobian into the system jacobian


  
  **Implementation Notes**

  Currently this computes the jacobian of the flux with respect to the
  face via complex step (using eqn.params_complex).  This should be replaced
  with dual numbers eventually
"""
function calcBoundaryFlux_nopre_diff{Tsol1, Tres1, Tmsh}(mesh::AbstractDGMesh{Tmsh},
                          sbp::AbstractSBP, eqn::EulerData{Tsol1, Tres1},
                          functor::BCType, idx_range::UnitRange,
                          bndry_facenums::AbstractArray{Boundary,1},
                          assembler::AssembleElementData)
  # calculate the jacobian of the boundary flux for the boundary condition
  # evaluated by the functor
  # because we are complex stepping the boundary condition and Tsol is usually
  # a real number, the Tsol and Tres of the eqn object are not the right
  # Tsol and Tres for this function

  nfaces = length(bndry_facenums)
  params = eqn.params_complex
  q_face = params.q_faceL
  Tsol = eltype(q_face)
  Tres = promote_type(Tsol, Tmsh)

#  flux_face = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)
  aux_vars = Array(Tres, 1)
  flux_k = zeros(Tres, mesh.numDofPerNode)
  flux_jac = zeros(Tres1, mesh.numDofPerNode, mesh.numDofPerNode, mesh.numNodesPerFace)
  res_jac = zeros(Tres1, mesh.numDofPerNode, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numNodesPerElement)

  h = 1e-20
  pert = Tsol(0, h)
  for i=1:nfaces  # loop over faces with this BC
    bndry_i = bndry_facenums[i]
    global_facenum = idx_range[i]

    # interpolate to face
    q_volume = sview(eqn.q, :, :, bndry_i.element)
    boundaryFaceInterpolate!(mesh.sbpface, bndry_i.face, q_volume, q_face)

    for j = 1:mesh.numNodesPerFace

      # get components
#      q = ro_sview(eqn.q_bndry, :, j, global_facenum)
      q_j = sview(q_face, :, j)
      # convert to conservative variables if needed
#      aux_vars = ro_sview(eqn.aux_vars_bndry, :, j, global_facenum)
      coords = ro_sview(mesh.coords_bndry, :, j, global_facenum)
      nrm_xy = ro_sview(mesh.nrm_bndry, :, j, global_facenum)
#      bndryflux_i = sview(flux_face, :, j)

      # compute the jacobian of the flux wrt q_face
      for k=1:mesh.numDofPerNode
        q_j[k] += pert
        aux_vars[1] = calcPressure(params, q_j)

        functor(params, q_j, aux_vars, coords, nrm_xy, flux_k)

        for p=1:mesh.numDofPerNode
          flux_jac[p, k, j] = imag(flux_k[p])/h
        end

        q_j[k] -= pert
      end  # end loop k
    end  # end looop j

    boundaryFaceIntegrate_jac!(mesh.sbpface, bndry_i.face, flux_jac, res_jac,
                               SummationByParts.Subtract())

    assembleBoundary(assembler, mesh, bndry_i, res_jac)
    fill!(res_jac, 0.0)

  end  # end loop i

  return nothing
end


