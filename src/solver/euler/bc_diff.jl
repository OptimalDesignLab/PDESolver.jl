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
function calcBoundaryFlux_nopre_diff(mesh::AbstractDGMesh{Tmsh},
                          sbp::AbstractSBP, eqn::EulerData{Tsol1, Tres1},
                          functor::BCType, idx_range::UnitRange,
                          bndry_facenums::AbstractArray{Boundary,1},
                          assembler::AssembleElementData) where {Tsol1, Tres1, Tmsh}
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
  aux_vars = Array{Tres}(1)
  flux_k = zeros(Tres, mesh.numDofPerNode)
  flux_jac = params.flux_dotL
  res_jac = params.res_jacLL

  #=
  flux_jac = zeros(Tres1, mesh.numDofPerNode, mesh.numDofPerNode, mesh.numNodesPerFace)
  res_jac = zeros(Tres1, mesh.numDofPerNode, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numNodesPerElement)
  =#

  # Setting up temporary arrays for storage of the real parts of q_j, coords, and nrm_xy.
  #   We need to pass in only the real parts of these to the complex step 
  #   derivative evaluation below, because there might be a complex 
  #   perturbation already present in these quantities.

  # Note: can't type tmpreal_q_j on Tsol1; Tsol1 might be Float64 even if jac_method == 2 for complex step
  #   Solution: just type it based on q_face[1,1]
  #   Before the addition of the tmpreal_q_j feature, q_j was perturbed directly.
  #   q_j was an sview of q_face; see in the loop over nodes below.
  tmpreal_q_j = zeros(typeof(q_face[1,1]), mesh.numDofPerNode)
  tmpreal_coords = zeros(Tmsh, mesh.dim)
  tmpreal_nrm_xy = zeros(Tmsh, mesh.dim)

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
      # q = ro_sview(eqn.q_bndry, :, j, global_facenum)
      q_j = sview(q_face, :, j)
      coords = ro_sview(mesh.coords_bndry, :, j, global_facenum)
      nrm_xy = ro_sview(mesh.nrm_bndry, :, j, global_facenum)

      # Store real part of q_j in tmpreal_q_j
      for ix = 1:mesh.numDofPerNode
        tmpreal_q_j[ix] = real(q_j[ix])     # should still be of complex type; required for complex pert below
      end

      # Store real part of coords & nrm_xy in their respective tmpreal's
      for ix = 1:mesh.dim
        tmpreal_coords[ix] = real(coords[ix])
        tmpreal_nrm_xy[ix] = real(nrm_xy[ix])
      end

      bndry_node = BoundaryNode(bndry_i, i, j)
      # compute the jacobian of the flux wrt q_face
      for k=1:mesh.numDofPerNode
        tmpreal_q_j[k] += pert
        aux_vars[1] = calcPressure(params, tmpreal_q_j)

        functor(params, tmpreal_q_j, aux_vars, tmpreal_coords, tmpreal_nrm_xy, flux_k, bndry_node)

        for p=1:mesh.numDofPerNode
          flux_jac[p, k, j] = imag(flux_k[p])/h
        end

        tmpreal_q_j[k] -= pert
      end  # end loop k

      # No need to replace any imaginary components of anything, since the original
      #   q_j, coords, and nrm_xy arrays were unaltered

    end  # end looop j

    boundaryFaceIntegrate_jac!(mesh.sbpface, bndry_i.face, flux_jac, res_jac,
                               SummationByParts.Subtract())

    
    if eqn.params.use_Minv == 1
      # multiply by Minv if needed
      for q=1:mesh.numNodesPerElement
        for p=1:mesh.numNodesPerElement
          val = mesh.jac[p, bndry_i.element]/sbp.w[p]  # entry in Minv
          @simd for m=1:mesh.numDofPerNode
            @simd for n=1:mesh.numDofPerNode
              res_jac[n, m, p, q] *= val
            end
          end
        end
      end
    end
    


    assembleBoundary(assembler, mesh.sbpface, mesh, bndry_i, res_jac)
    fill!(res_jac, 0.0)

  end  # end loop i

  return nothing
end


