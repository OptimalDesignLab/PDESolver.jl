# differentiated version of functions in flux.jl

function calcFaceIntegral_nopre_diff(
                                mesh::AbstractDGMesh{Tmsh},
                                sbp::AbstractSBP,
                                eqn::EulerData{Tsol, Tres, Tdim, :conservative},
                                opts,
                                functor::FluxType_diff,
                                interfaces::AbstractArray{Interface, 1},
                                assembler::AssembleElementData) where {Tmsh, Tsol, Tres, Tdim}


  nfaces = length(interfaces)
  params = eqn.params

  data = params.calc_face_integrals_data
  @unpack data q_faceL q_faceR flux_dotL flux_dotR res_jacLL res_jacLR res_jacRL res_jacRR

#  flux_face = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)

  for i=1:nfaces
    iface_i = interfaces[i]
    fill!(flux_dotL, 0)
    fill!(flux_dotR, 0)

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
    assembleInterface(assembler, mesh.sbpface, mesh, iface_i, res_jacLL,
                      res_jacLR, res_jacRL, res_jacRR)

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


"""
  Reverse mode wrt metrics of [`calcFaceIntegral_nopre`](@ref)

  **Inputs**

   * mesh: bar fields are updated
   * sbp
   * eqn: res_bar fields in input
   * opts
   * functor_revm: reverse mode wrt metrics functor (`FluxType_revm`)
   * interfaces
"""
function calcFaceIntegral_nopre_revm(
        mesh::AbstractDGMesh{Tmsh},
        sbp::AbstractSBP,
        eqn::EulerData{Tsol, Tres, Tdim, :conservative},
        opts,
        functor_revm::FluxType_revm,
        interfaces::AbstractArray{Interface, 1}) where {Tmsh, Tsol, Tres, Tdim}

  nfaces = length(interfaces)
  params = eqn.params

  q_faceL = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
  q_faceR = zeros(q_faceL)

  flux_face_bar = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)

  # reverse sweep
  for i=1:nfaces
    iface_i = interfaces[i]

    qL = ro_sview(eqn.q, :, :, iface_i.elementL)
    qR = ro_sview(eqn.q, :, :, iface_i.elementR)
    interiorFaceInterpolate!(mesh.sbpface, iface_i, qL, qR, q_faceL,
                             q_faceR)

    fill!(flux_face_bar, 0)
    resL_bar = sview(eqn.res_bar, :, :, iface_i.elementL)
    resR_bar = sview(eqn.res_bar, :, :, iface_i.elementR)
    interiorFaceIntegrate_rev!(mesh.sbpface, iface_i, flux_face_bar, resL_bar,
                               resR_bar, SummationByParts.Subtract())

    for j=1:mesh.numNodesPerFace
      qL_j = ro_sview(q_faceL, :, j)
      qR_j = ro_sview(q_faceR, :, j)

      eqn.aux_vars_face[1, j, i] = calcPressure(params, qL_j)
      aux_vars = ro_sview(eqn.aux_vars_face, :, j, i)

      nrm_xy = sview(mesh.nrm_face, :, j, i)
      nrm_bar = sview(mesh.nrm_face_bar, :, j, i)
      flux_bar_j = sview(flux_face_bar, :, j)

      functor_revm(params, qL_j, qR_j, aux_vars, nrm_xy, nrm_bar, flux_bar_j)
    end  # end loop j

  end  # end loop i

  return nothing
end




function getFaceElementIntegral_diff(
                           mesh::AbstractDGMesh{Tmsh},
                           sbp::AbstractSBP, eqn::EulerData{Tsol, Tres, Tdim},
                           face_integral_functor::FaceElementIntegralType,
                           flux_functor::FluxType_diff,
                           sbpface::AbstractFace,
                           interfaces::AbstractArray{Interface, 1},
                           assembler::AssembleElementData) where {Tmsh, Tsol, Tres, Tdim}

  params = eqn.params
#  sbpface = mesh.sbpface
  nfaces = length(interfaces)
#  resL2 = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement)
#  resR2 = zeros(resL2)

  data = params.calc_face_integrals_data
  @unpack data res_jacLL res_jacLR res_jacRL res_jacRR


  for i=1:nfaces
    iface = interfaces[i]
    elL = iface.elementL
    elR = iface.elementR
    qL = ro_sview(eqn.q, :, :, elL)
    qR = ro_sview(eqn.q, :, :, elR)
    aux_vars = ro_sview(eqn.aux_vars, :, :, elL)
    nrm_face = ro_sview(mesh.nrm_face, :, :, i)
    resL = sview(eqn.res, :, :, elL)
    resR = sview(eqn.res, :, :, elR)

    fill!(res_jacLL, 0); fill!(res_jacLR, 0)
    fill!(res_jacRL, 0); fill!(res_jacRR, 0)
    calcFaceElementIntegral_diff(face_integral_functor, params, sbpface, iface,
                       qL, qR, aux_vars,
                       nrm_face, flux_functor, res_jacLL, res_jacLR, res_jacRL, res_jacRR)

    assembleInterface(assembler, sbpface, mesh, iface, res_jacLL,
                      res_jacLR, res_jacRL, res_jacRR)

    # multiply by Minv if needed
    if params.use_Minv == 1
      for q=1:mesh.numNodesPerElement
        for p=1:mesh.numNodesPerElement
          valL = mesh.jac[p, iface.elementL]/sbp.w[p]  # entry in Minv
          valR = mesh.jac[p, iface.elementR]/sbp.w[p]
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
 
  end

  fill!(res_jacLL, 0); fill!(res_jacLR, 0)
  fill!(res_jacRL, 0); fill!(res_jacRR, 0)

  return nothing
end


"""
  Reverse mode wrt metrics of [`getFaceElementIntegral`](@ref)
"""
function getFaceElementIntegral_revm(
                           mesh::AbstractDGMesh{Tmsh},
                           sbp::AbstractSBP, eqn::EulerData{Tsol, Tres, Tdim},
                           face_integral_functor::FaceElementIntegralType,
                           flux_functor::FluxType,
                           sbpface::AbstractFace,
                           interfaces::AbstractArray{Interface, 1}) where {Tmsh, Tsol, Tres, Tdim}

  params = eqn.params
  nfaces = length(interfaces)
  for i=1:nfaces
    iface = interfaces[i]
    elL = iface.elementL
    elR = iface.elementR
    qL = ro_sview(eqn.q, :, :, elL)
    qR = ro_sview(eqn.q, :, :, elR)
    aux_vars = ro_sview(eqn.aux_vars, :, :, elL)
    nrm_face = ro_sview(mesh.nrm_face, :, :, i)
    nrm_face_bar = sview(mesh.nrm_face_bar, :, :, i)
    resL_bar = ro_sview(eqn.res_bar, :, :, elL)
    resR_bar = ro_sview(eqn.res_bar, :, :, elR)

    calcFaceElementIntegral_revm(face_integral_functor, params, sbpface, iface,
                        qL, qR, aux_vars, nrm_face, nrm_face_bar, flux_functor,
                        resL_bar, resR_bar)

  end

  return nothing
end





#------------------------------------------------------------------------------
# Shared face integrals

"""
  Differentiated version of
  [`calcSharedFaceIntegrals_element_inner`](@ref).
  It presents the interface required by [`finishExchangeData`](@ref)
"""
function calcSharedFaceIntegrals_element_diff(
                            mesh::AbstractDGMesh{Tmsh},
                            sbp::AbstractSBP, eqn::EulerData{Tsol},
                            opts, data::SharedFaceData) where {Tmsh, Tsol}

    calcSharedFaceIntegrals_nopre_element_inner_diff(mesh, sbp, eqn, opts, data, eqn.flux_func_diff, eqn.assembler)

  return nothing
end




"""
  Differentiated version of [`calcSharedFaceIntegrals_nopre_element_inner`](@ref),
"""
function calcSharedFaceIntegrals_nopre_element_inner_diff(
                            mesh::AbstractDGMesh{Tmsh},
                            sbp::AbstractSBP, eqn::EulerData{Tsol, Tres},
                            opts, data::SharedFaceData, functor::FluxType_diff,
                            assembler::AssembleElementData) where {Tmsh, Tsol, Tres}

  params = eqn.params
  fdata = params.calc_face_integrals_data
  @unpack fdata q_faceL q_faceR flux_dotL flux_dotR res_jacLL res_jacLR res_jacRL res_jacRR

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

    fill!(flux_dotL, 0); fill!(flux_dotR, 0)

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

"""
  Differentiated version of [`calcSharedFaceElementIntegrals_element`](@ref)
"""
function calcSharedFaceElementIntegrals_element_diff(
                            mesh::AbstractDGMesh{Tmsh},
                            sbp::AbstractSBP, eqn::EulerData{Tsol, Tres},
                            opts, data::SharedFaceData) where {Tmsh, Tsol, Tres}

  if opts["use_staggered_grid"]
    error("explicit computaton of jacobian not supported for staggered grids")

  else
    calcSharedFaceElementIntegrals_element_inner_diff(mesh, sbp, eqn, opts,
                            data, eqn.face_element_integral_func,  eqn.flux_func_diff, eqn.assembler)
  end

  return nothing
end

"""
  This function loops over given set of shared faces and computes a face
  integral that
  uses data from all volume nodes.  See [`FaceElementIntegralType`](@ref)
  for details on the integral performed.

  **Inputs**:

   * mesh
   * sbp
   * eqn
   * opts
   * data: a SharedFaceData specifying which shared faces to compute
   * face_integral_functor
   * flux_functor

"""
function calcSharedFaceElementIntegrals_element_inner_diff(
                            mesh::AbstractDGMesh{Tmsh},
                            sbp::AbstractSBP, eqn::EulerData{Tsol, Tres},
                            opts, data::SharedFaceData,
                            face_integral_functor::FaceElementIntegralType,
                            flux_functor::FluxType_diff,
                            assembler::AssembleElementData) where {Tmsh, Tsol, Tres}

  if opts["parallel_data"] != "element"
    throw(ErrorException("cannot use calcSharedFaceIntegrals_element_diff without parallel element data"))
  end

  #TODO: consider writing specialized 1-sided flux functions that only do
  #      the local half.
  fdata = eqn.params.calc_face_integrals_data
  @unpack fdata res_jacLL res_jacLR res_jacRL res_jacRR


  q = eqn.q
  params = eqn.params
  # we don't care about elementR here, so use this throwaway array
  resR = Array{Tres}(mesh.numDofPerNode, mesh.numNodesPerElement)

  # get the data for the parallel interface
  idx = data.peeridx
  interfaces = data.interfaces
  bndries_local = data.bndries_local
  bndries_remote = data.bndries_remote
#    qL_arr = eqn.q_face_send[i]
  qR_arr = data.q_recv
  nrm_face_arr = mesh.nrm_sharedface[idx]

  start_elnum = mesh.shared_element_offsets[idx]

  # compute the integrals
  for j=1:length(interfaces)
    iface_j = interfaces[j]
    elL = iface_j.elementL
    elR = iface_j.elementR - start_elnum + 1  # is this always equal to j?
    qL = ro_sview(q, :, :, elL)
    qR = ro_sview(qR_arr, :, :, elR)
    aux_vars = ro_sview(eqn.aux_vars, :, :, elL)
    nrm_face = ro_sview(nrm_face_arr, :, :, j)

    fill!(res_jacLL, 0); fill!(res_jacLR, 0)
    fill!(res_jacRL, 0); fill!(res_jacRR, 0)
    calcFaceElementIntegral_diff(face_integral_functor, eqn.params, mesh.sbpface,
                            iface_j, qL, qR, aux_vars,
                            nrm_face, flux_functor, res_jacLL, res_jacLR, res_jacRL, res_jacRR)

    assembleSharedFace(assembler, mesh.sbpface, mesh, iface_j, res_jacLL, res_jacLR)
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
 
  end  # end loop j

  fill!(res_jacLL, 0); fill!(res_jacLR, 0)
  fill!(res_jacRL, 0); fill!(res_jacRR, 0)

  return nothing
end



#------------------------------------------------------------------------------
# Functors

"""
  Calls the [`RoeSolver_diff`](@ref)
"""
mutable struct RoeFlux_diff <: FluxType_diff
end

function (obj::RoeFlux_diff)(params::ParamType,
              uL::AbstractArray{Tsol,1},
              uR::AbstractArray{Tsol,1},
              aux_vars::AbstractVector{Tres},
              nrm::AbstractVector{Tmsh},
              F_dotL::AbstractArray{Tres},
              F_dotR::AbstractArray{Tres}) where {Tsol, Tres, Tmsh}

  RoeSolver_diff(params, uL, uR, aux_vars, nrm, F_dotL, F_dotR)
  return nothing
end


mutable struct LFFlux_diff <: FluxType_diff
end

function (obj::LFFlux_diff)(params::ParamType,
              uL::AbstractArray{Tsol,1},
              uR::AbstractArray{Tsol,1},
              aux_vars::AbstractVector{Tres},
              nrm::AbstractVector{Tmsh},
              F_dotL::AbstractArray{Tres},
              F_dotR::AbstractArray{Tres}) where {Tsol, Tres, Tmsh}

  calcLFFlux_diff(params, uL, uR, aux_vars, nrm, F_dotL, F_dotR)

  return nothing
end

"""
  Calls the [`calcEulerFlux_IR_diff`](@ref)
"""
mutable struct IRFlux_diff <: FluxType_diff
end

function (obj::IRFlux_diff)(params::ParamType,
              uL::AbstractArray{Tsol,1},
              uR::AbstractArray{Tsol,1},
              aux_vars::AbstractVector{Tres},
              nrm::AbstractArray{Tmsh},
              F_dotL::AbstractArray{Tres},
              F_dotR::AbstractArray{Tres}) where {Tsol, Tres, Tmsh}

  calcEulerFlux_IR_diff(params, uL, uR, aux_vars, nrm, F_dotL, F_dotR)

  return nothing
end

mutable struct StandardFlux_diff <: FluxType_diff
end

function (obj::StandardFlux_diff)(params::ParamType,
              uL::AbstractArray{Tsol,1},
              uR::AbstractArray{Tsol,1},
              aux_vars::AbstractVector{Tres},
              nrm::AbstractVector{Tmsh},
              F_dotL::AbstractArray{Tres},
              F_dotR::AbstractArray{Tres}) where {Tsol, Tres, Tmsh}

  # this should be implemented eventually, but for now this method needs to
  # exist as a stub
  error("StandardFlux() called, something unexpected has occured")

end


mutable struct IRSLFFlux_diff <: FluxType_diff
end

function (obj::IRSLFFlux_diff)(params::ParamType,
              uL::AbstractArray{Tsol,1},
              uR::AbstractArray{Tsol,1},
              aux_vars::AbstractVector{Tres},
              nrm::AbstractVector{Tmsh},
              F_dotL::AbstractArray{Tres},
              F_dotR::AbstractArray{Tres}) where {Tsol, Tres, Tmsh}

  calcEulerFlux_IRSLF_diff(params, uL, uR, aux_vars, nrm, F_dotL, F_dotR)
end



mutable struct ErrorFlux_diff <: FluxType_diff
end

function (obj::ErrorFlux_diff)(params::ParamType,
              uL::AbstractArray{Tsol,1},
              uR::AbstractArray{Tsol,1},
              aux_vars::AbstractVector{Tres},
              nrm::AbstractVector{Tmsh},
              F_dotL::AbstractArray{Tres},
              F_dotR::AbstractArray{Tres}) where {Tsol, Tres, Tmsh}


  error("ErrorFlux() called, something unexpected has occured")

end



"""
  Container for all differentiated flux functors.  Maps name to object.
  The names are exactly the same as the non-differentiated functor.
"""
global const FluxDict_diff = Dict{String, FluxType_diff}(
"RoeFlux" => RoeFlux_diff(),
"LFFlux" => LFFlux_diff(),
"IRFlux" => IRFlux_diff(),
"IRSLFFlux" => IRSLFFlux_diff(),
"StandardFlux" => StandardFlux_diff(),
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
    name2 = opts["Volume_flux_name"]
  else
    name = "ErrorFlux"
    name2 = "ErrorFlux"
  end

  eqn.flux_func_diff = FluxDict_diff[name]
  eqn.volume_flux_func_diff = FluxDict_diff[name2]

  return nothing
end



