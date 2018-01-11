# differentiated version of functions in flux.jl

function calcFaceIntegral_nopre_diff{Tmsh, Tsol, Tres, Tdim}(
                                mesh::AbstractDGMesh{Tmsh},
                                sbp::AbstractSBP,
                                eqn::EulerData{Tsol, Tres, Tdim, :conservative},
                                opts,
                                functor::FluxType_diff,
                                interfaces::AbstractArray{Interface, 1},
                                assembler::AssembleElementData)


  nfaces = length(interfaces)
  params = eqn.params

  q_faceL = params.q_faceL
  q_faceR = params.q_faceR
  flux_dotL = params.flux_dotL
  flux_dotR = params.flux_dotR
  res_jacLL = params.res_jacLL
  res_jacLR = params.res_jacLR
  res_jacRL = params.res_jacRL
  res_jacRR = params.res_jacRR

  #=
  q_faceL = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
  q_faceR = zeros(q_faceL)

  flux_dotL = zeros(Tres, mesh.numDofPerNode, mesh.numDofPerNode, mesh.numNodesPerFace)
  flux_dotR = zeros(Tres, mesh.numDofPerNode, mesh.numDofPerNode, mesh.numNodesPerFace)
  
  res_jacLL = zeros(Tres, mesh.numDofPerNode, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numNodesPerElement)
  res_jacLR = zeros(res_jacLL)
  res_jacRL = zeros(res_jacLL)
  res_jacRR = zeros(res_jacLL)
  =#
#  flux_face = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)

  for i=1:nfaces
    iface_i = interfaces[i]

    qL = ro_sview(eqn.q, :, :, iface_i.elementL)
    qR = ro_sview(eqn.q, :, :, iface_i.elementR)
    interiorFaceInterpolate!(mesh.sbpface, iface_i, qL, qR, q_faceL,
                             q_faceR)

    # compute dF/dqface
    for j=1:mesh.numNodesPerFace
      qL_j = ro_sview(q_faceL, :, j)
      qR_j = ro_sview(q_faceR, :, j)

      eqn.aux_vars_face[1, j, i] = calcPressure(params, qL_j)
      aux_vars = ro_sview(eqn.aux_vars_face, :, j, i)

      nrm_xy = ro_sview(mesh.nrm_face, :, j, i)
#      flux_j = sview(flux_face, :, j)
      flux_dotL_j = sview(flux_dotL, :, :, j)
      flux_dotR_j = sview(flux_dotR, :, :, j)

      functor(params, qL_j, qR_j, aux_vars, nrm_xy, flux_dotL_j, flux_dotR_j)
    end  # end loop j

    # compute dR/dq
    interiorFaceIntegrate_jac!(mesh.sbpface, iface_i, flux_dotL, flux_dotR,
                             res_jacLL, res_jacLR, res_jacRL, res_jacRR,
                             SummationByParts.Subtract())

    
    # multiply by Minv if needed
    if params.use_Minv == 1
      for q=1:mesh.numNodesPerElement
        for p=1:mesh.numNodesPerElement
          valL = mesh.jac[p, iface_i.elementL]/sbp.w[p]  # entry in Minv
          valR = mesh.jac[p, iface_i.elementR]/sbp.w[p]
          @simd for m=1:mesh.numDofPerNode
            @simd for n=1:mesh.numDofPerNode
              res_jacLL[n, m, p, q] *= valL
              res_jacLR[n, m, p, q] *= valL
              res_jacRL[n, m, p, q] *= valR
              res_jacRR[n, m, p, q] *= valR
            end
          end
        end
      end
    end
    


    # assemble into the Jacobian
    assembleInterface(assembler, mesh.sbpface, mesh, iface_i, res_jacLL, res_jacLR,
                                                res_jacRL, res_jacRR)

    #TODO: for sparse faces, only zero out the needed entries, to avoid loading
    #      the entire array into cache
    # Tests show this doesn't make a difference
    fill!(res_jacLL, 0.0)
    fill!(res_jacLR, 0.0)
    fill!(res_jacRL, 0.0)
    fill!(res_jacRR, 0.0)
  end  # end loop i

  return nothing
end


#------------------------------------------------------------------------------
# Shared face integrals

"""
  Differentiated version of
  [`calcSharedFaceIntegrals_element_inner`](@ref).
  It presents the interface required by [`finishExchangeData`](@ref)
"""
function calcSharedFaceIntegrals_element_diff{Tmsh, Tsol}(
                            mesh::AbstractDGMesh{Tmsh},
                            sbp::AbstractSBP, eqn::EulerData{Tsol},
                            opts, data::SharedFaceData)

    calcSharedFaceIntegrals_nopre_element_inner_diff(mesh, sbp, eqn, opts, data, eqn.flux_func_diff, eqn.assembler)

  return nothing
end




"""
  Differentiated version of [`calcSharedFaceIntegrals_nopre_element_inner`](@ref),
"""
function calcSharedFaceIntegrals_nopre_element_inner_diff{Tmsh, Tsol, Tres}(
                            mesh::AbstractDGMesh{Tmsh},
                            sbp::AbstractSBP, eqn::EulerData{Tsol, Tres},
                            opts, data::SharedFaceData, functor::FluxType_diff,
                            assembler::AssembleElementData)

  params = eqn.params

  q_faceL = params.q_faceL
  q_faceR = params.q_faceR
  flux_dotL = params.flux_dotL
  flux_dotR = params.flux_dotR
  res_jacLL = params.res_jacLL
  res_jacLR = params.res_jacLR
  res_jacRL = params.res_jacRL
  res_jacRR = params.res_jacRR
  #=
  # TODO: make these fields of params
  q_faceL = Array(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
  q_faceR = Array(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
  flux_dotL = zeros(Tres, mesh.numDofPerNode, mesh.numDofPerNode, mesh.numNodesPerFace)
  flux_dotR = zeros(flux_dotL)

  res_jacLL = zeros(Tres, mesh.numDofPerNode, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numNodesPerElement)
  res_jacLR = zeros(res_jacLL)
  res_jacRL = zeros(res_jacLL)  # TODO: create a NoOp array for these
  res_jacRR = zeros(res_jacLL)
  =#

  # get data
  idx = data.peeridx
  interfaces = data.interfaces
  bndries_local = data.bndries_local
  bndries_remote = data.bndries_remote
#    qL_arr = eqn.q_face_send[i]
  qR_arr = data.q_recv
  nrm_arr = mesh.nrm_sharedface[idx]
  aux_vars_arr = eqn.aux_vars_sharedface[idx]
#  flux_face = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)

  start_elnum = mesh.shared_element_offsets[idx]

  for j=1:length(interfaces)
    iface_j = interfaces[j]
    bndryL_j = bndries_local[j]
    bndryR_j = bndries_remote[j]
    fL = bndryL_j.face

    # interpolate to face
    qL = ro_sview(eqn.q, :, :, iface_j.elementL)
    el_r = iface_j.elementR - start_elnum + 1
    qR = ro_sview(qR_arr, :, :, el_r)
    interiorFaceInterpolate!(mesh.sbpface, iface_j, qL, qR, q_faceL, q_faceR)

    # calculate flux jacobian
    for k=1:mesh.numNodesPerFace
      qL_k = ro_sview(q_faceL, :, k)
      qR_k = ro_sview(q_faceR, :, k)
      aux_vars = ro_sview(aux_vars_arr, :, k, j)
      nrm_xy = ro_sview(nrm_arr, :, k, j)
#      flux_k = sview(flux_face, :, k)
      flux_dotL_k = sview(flux_dotL, :, :, k)
      flux_dotR_k = sview(flux_dotR, :, :, k)

      parent(aux_vars)[1] = calcPressure(params, qL_k)

      functor(params, qL_k, qR_k, aux_vars, nrm_xy, flux_dotL_k, flux_dotR_k)
     end

     # this is excessive because we don't need jacRL, jacRR, but
     # boundaryFaceIntegrate_jac can't handle jacLR
     interiorFaceIntegrate_jac!(mesh.sbpface, iface_j, flux_dotL, flux_dotR,
                                res_jacLL, res_jacLR, res_jacRL, res_jacRR,
                                SummationByParts.Subtract())

    
    # multiply by Minv if needed
    if params.use_Minv == 1
      for q=1:mesh.numNodesPerElement
        for p=1:mesh.numNodesPerElement
          valL = mesh.jac[p, iface_j.elementL]/sbp.w[p]  # entry in Minv
          @simd for m=1:mesh.numDofPerNode
            @simd for n=1:mesh.numDofPerNode
              res_jacLL[n, m, p, q] *= valL
              res_jacLR[n, m, p, q] *= valL
            end
          end
        end
      end
    end
    


     assembleSharedFace(assembler, mesh.sbpface, mesh, iface_j, res_jacLL, res_jacLR)
     fill!(res_jacLL, 0.0)
     fill!(res_jacLR, 0.0)

   end  # end loop over interfaces

   fill!(res_jacRL, 0.0)
   fill!(res_jacRR, 0.0)

  return nothing
end




#------------------------------------------------------------------------------
# Functors

"""
  Calls the [`RoeSolver_diff`](@ref)
"""

type RoeFlux_diff <: FluxType_diff
end

function call{Tsol, Tres, Tmsh}(obj::RoeFlux_diff, params::ParamType,
              uL::AbstractArray{Tsol,1},
              uR::AbstractArray{Tsol,1},
              aux_vars::AbstractVector{Tres},
              nrm::AbstractVector{Tmsh},
              F_dotL::AbstractMatrix{Tres},
              F_dotR::AbstractMatrix{Tres})

  RoeSolver_diff(params, uL, uR, aux_vars, nrm, F_dotL, F_dotR)
  return nothing
end

type ErrorFlux_diff <: FluxType_diff
end

function call{Tsol, Tres, Tmsh}(obj::ErrorFlux_diff, params::ParamType,
              uL::AbstractArray{Tsol,1},
              uR::AbstractArray{Tsol,1},
              aux_vars::AbstractVector{Tres},
              nrm::AbstractVector{Tmsh},
              F_dotL::AbstractMatrix{Tres},
              F_dotR::AbstractMatrix{Tres})


  error("ErrorFlux() called, something unexpected has occured")

end


global const FluxDict_diff = Dict{ASCIIString, FluxType_diff}(
"RoeFlux" => RoeFlux_diff(),
"ErrorFlux" => ErrorFlux_diff(),
)

"""
  Gets the functor for the differentiated version of the flux function and
  stores is to `eqn.flux_func_diff`

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts

  **Options Keys**

  Uses the `Flux_name` key to get the functor (same as for the non-differentiated
  version). 

  This function only gets the functor if the jacobian will be computed explicitly
  as specified by `calc_jac_explicit`, otherwise it uses a functor that
  will throw an error if used.
"""
function getFluxFunctors_diff(mesh::AbstractDGMesh, sbp, eqn, opts)

  if opts["calc_jac_explicit"]
    name = opts["Flux_name"]
  else
    name = "ErrorFlux"
  end

  eqn.flux_func_diff = FluxDict_diff[name]

  return nothing
end



