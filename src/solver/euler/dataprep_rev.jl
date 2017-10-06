# dataprep_rev.jl
function getEulerFlux_revm{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh},
                                        sbp::AbstractSBP,
                                        eqn::EulerData{Tsol, Tres, Tdim}, opts)
# calculate Euler flux in parametric coordinate directions, stores it in eqn.flux_parametric

  nrm = zeros(Tmsh, Tdim)
  nrm_bar = zeros(nrm)
  fill!(mesh.dxidx_bar, 0.0)
  for i=1:mesh.numEl  # loop over elements
    for j=1:mesh.numNodesPerElement  # loop over nodes on current element
      q_vals = sview(eqn.q, :, j, i)
      aux_vars = sview(eqn.aux_vars, :, j, i)

      for k=1:Tdim  # loop over dimensions
        for p=1:Tdim
          nrm[p] = mesh.dxidx[k, p, j, i]
          nrm_bar[p] = zero(Tmsh)
        end

        flux_bar = sview(eqn.flux_parametric_bar, :, j, i, k)

      	# this will dispatch to the proper calcEulerFlux
        calcEulerFlux_revm(eqn.params, q_vals, aux_vars, nrm, flux_bar, nrm_bar)

        for p=1:Tdim
          # nrm[p] = mesh.dxidx[k, p, j, i]
          mesh.dxidx_bar[k,p,j,i] += nrm_bar[p]
        end
      end
    end
  end


  # writeFlux(mesh, sbp, eqn, opts)

  return nothing
end

function getBCFluxes_revm(mesh::AbstractMesh, sbp::AbstractSBP, eqn::EulerData, opts)
  #get all the fluxes for all the boundary conditions and save them in eqn.bndryflux

  fill!(mesh.dxidx_bndry_bar, 0.0)
  for i=1:mesh.numBC
    functor_i = mesh.bndry_funcs_revm[i]
    start_index = mesh.bndry_offsets[i]
    end_index = mesh.bndry_offsets[i+1]
    idx_range = start_index:end_index  # TODO: should this be start_index:(end_index - 1) ?
    bndry_facenums_i = sview(mesh.bndryfaces, start_index:(end_index - 1))
    bndryflux_i = sview(eqn.bndryflux_bar, :, :, start_index:(end_index - 1))

    # call the function that calculates the flux for this boundary condition
    # passing the functor into another function avoid type instability
    calcBoundaryFlux_revm(mesh, sbp, eqn, functor_i, idx_range, bndry_facenums_i, bndryflux_i)
  end

  writeBoundary(mesh, sbp, eqn, opts)

  return nothing
end


function calcBoundaryFlux_revm{Tmsh,  Tsol, Tres}( mesh::AbstractDGMesh{Tmsh},
                          sbp::AbstractSBP, eqn::EulerData{Tsol},
                          functor_bar::BCType_revm, idx_range::UnitRange,
                          bndry_facenums::AbstractArray{Boundary,1},
                          bndryflux_bar::AbstractArray{Tres, 3})
  # calculate the boundary flux for the boundary condition evaluated by the
  # functor

  nfaces = length(bndry_facenums)
  q2 = zeros(Tsol, mesh.numDofPerNode)
  nrm = zeros(Tmsh, size(mesh.sbpface.normal,1))
  for i=1:nfaces  # loop over faces with this BC
    bndry_i = bndry_facenums[i]
    global_facenum = idx_range[i]
    for j = 1:mesh.numNodesPerFace

      # get components
      q = sview(eqn.q_bndry, :, j, global_facenum)
      convertToConservative(eqn.params, q, q2)
      aux_vars = sview(eqn.aux_vars_bndry, :, j, global_facenum)
      x = sview(mesh.coords_bndry, :, j, global_facenum)
      dxidx = sview(mesh.dxidx_bndry, :, :, j, global_facenum)
      dxidx_bar = sview(mesh.dxidx_bndry_bar, :, :, j, global_facenum)
      nrm[:] = mesh.sbpface.normal[:,bndry_i.face]
      bndryflux_i = sview(bndryflux_bar, :, j, i)

      functor_bar(q2, aux_vars, x, dxidx, dxidx_bar, nrm, bndryflux_i, eqn.params)
    end
  end

  return nothing
end

function calcFaceFlux_revm{Tmsh,  Tsol, Tres, Tdim}( mesh::AbstractDGMesh{Tmsh},
                          sbp::AbstractSBP,
                          eqn::EulerData{Tsol, Tres, Tdim, :conservative},
                          functor_bar::FluxType_revm,
                          interfaces::AbstractArray{Interface,1},
                          flux_face_bar::AbstractArray{Tres, 3})

  fill!(mesh.dxidx_face_bar, 0.0)
  nfaces = length(interfaces)
  nrm = zeros(Tmsh, size(mesh.sbpface.normal,1))
  for i=1:nfaces  # loop over faces
    interface_i = interfaces[i]
    for j = 1:mesh.numNodesPerFace
      eL = interface_i.elementL
      fL = interface_i.faceL

      # get components
      qL = sview(eqn.q_face, :, 1, j, i)
      qR = sview(eqn.q_face, :, 2, j, i)
      dxidx = sview(mesh.dxidx_face, :, :, j, i)
      dxidx_bar = sview(mesh.dxidx_face_bar, :, :, j, i)
      aux_vars = sview(eqn.aux_vars_face, :, j, i)
      nrm[:] = mesh.sbpface.normal[:,fL]

      #flux_j = sview(flux_face_bar, :, j, i)
      #functor(eqn.params, qL, qR, aux_vars, dxidx, nrm, flux_j)

      flux_j_bar = sview(flux_face_bar, :, j, i)
      functor_bar(eqn.params, qL, qR, aux_vars, dxidx, nrm, flux_j_bar, dxidx_bar)
    end
  end

  return nothing
end
