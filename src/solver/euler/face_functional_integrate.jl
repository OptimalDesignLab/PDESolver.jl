# functions for integrating over faces

"""
  Does an integral of the form:

  \sum \int f dGamma
  
  over the interior faces.  This is different than the regular face integrals
  used in the weak form, which are conservative (resL -= integral,
  resR += integral, so the sum is zero).  This is useful for computing
  functionals.  Up to `numDofPerNode` functionals can be computed at the
  same time.

  This function requires parallel communication to already have been completed
  before it is called.

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * flux_func: computes the integrand `f` at each face node, `FluxType`

  **Inputs/Outputs

   * val: vector of length `numDofPerNode`, one entry for each functional,
          to be overwritten
"""
function integrateFaceQuantity(mesh, sbp, eqn::EulerData{Tsol, Tres}, opts,
                               flux_func::FluxType, val::Vector) where {Tsol, Tres}

  q_faceL = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
  q_faceR = zeros(q_faceL)

  flux_face = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)
  interfaces = mesh.interfaces
  fill!(val, 0)

  for i=1:length(interfaces)
    iface_i = interfaces[i]

    qL = ro_sview(eqn.q, :, :, iface_i.elementL)
    qR = ro_sview(eqn.q, :, :, iface_i.elementR)
    interiorFaceInterpolate!(mesh.sbpface, iface_i, qL, qR, q_faceL,
                             q_faceR)

    for j=1:mesh.numNodesPerFace
      qL_j = ro_sview(q_faceL, :, j)
      qR_j = ro_sview(q_faceR, :, j)

      eqn.aux_vars_face[1, j, i] = calcPressure(eqn.params, qL_j)
      aux_vars = ro_sview(eqn.aux_vars_face, :, j, i)

      nrm_xy = ro_sview(mesh.nrm_face, :, j, i)
      flux_j = sview(flux_face, :, j)

      flux_func(eqn.params, qL_j, qR_j, aux_vars, nrm_xy, flux_j)

      for k=1:mesh.numDofPerNode
        val[k] += mesh.sbpface.wface[j]*flux_j[k]
      end
    end  # end loop j
  end

  if opts["parallel_data"] == PARALLEL_DATA_FACE
    integrateSharedFaceQuantity(mesh, sbp, eqn, opts, flux_func, val)
  elseif opts["parallel_data"] == PARALLEL_DATA_ELEMENT
    integrateSharedFaceQuantity_element(mesh, sbp, eqn, opts, flux_func, val)
  else
    error("unrecognized parallel data: $(opts["parallel_data"])")
  end

  return nothing
end


"""
  Integrate a quantity over the shared faces and sum.  Assumes parallel
  communication has already been done
"""
function integrateSharedFaceQuantity(mesh, sbp, eqn::EulerData{Tsol, Tres}, opts,
                          flux_func::FluxType, val::Vector) where {Tsol, Tres}


  params = eqn.params

  for idx=1:length(eqn.shared_data)
    data = eqn.shared_data[idx]

    if getParallelData(data) != PARALLEL_DATA_FACE
      throw(ErrorException("cannot use integrateSharedFaceQuantity without parallel face data"))
    end

    # only one process should do the integral for a given face.
    # Arbitrarily choose the lower ranked process
    if data.peernum < mesh.myrank
      continue
    end

    # calculate the flux
    idx = data.peeridx
    interfaces = data.interfaces
    qL_arr = data.q_send
    qR_arr = data.q_recv
    aux_vars_arr = eqn.aux_vars_sharedface[idx]
    nrm_arr = mesh.nrm_sharedface[idx]
    flux_j = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)

    for j=1:length(interfaces)
      interface_i = interfaces[j]
      eL = interface_i.elementL
      fL = interface_i.faceL

      # compute the flux
      for k=1:mesh.numNodesPerFace
        qL = ro_sview(qL_arr, :, k, j)
        qR = ro_sview(qR_arr, :, k, j)
        aux_vars = ro_sview(aux_vars_arr, :, k, j)
        parent(aux_vars)[1] = calcPressure(params, qL)
        nrm_xy = ro_sview(nrm_arr, :, k, j)
        flux_k = sview(flux_j, :, k)
        flux_func(params, qL, qR, aux_vars, nrm_xy, flux_k)

        for p=1:mesh.numDofPerNode
          val[p] += mesh.sbpface.wface[k]*flux_k[p]
        end
      end
    end
  end

  return nothing
end
   

function integrateSharedFaceQuantity_element(mesh, sbp, eqn::EulerData{Tsol, Tres}, opts,
                          flux_func::FluxType, val::Vector) where {Tsol, Tres}
  q_faceL = Array{Tsol}(mesh.numDofPerNode, mesh.numNodesPerFace)
  q_faceR = Array{Tsol}(mesh.numDofPerNode, mesh.numNodesPerFace)

  for idx=1:length(eqn.shared_data)

    data = eqn.shared_data[idx]

    if getParallelData(data) != PARALLEL_DATA_ELEMENT
      throw(ErrorException("cannot use integrateSharedFaceQuantity without parallel element data"))
    end


    # only one process should do the integral for a given face.
    # Arbitrarily choose the lower ranked process
    if data.peernum < mesh.myrank
      continue
    end

    interfaces = data.interfaces
    qR_arr = data.q_recv               
    nrm_arr = mesh.nrm_sharedface[idx]
    aux_vars_arr = eqn.aux_vars_sharedface[idx]
    flux_face = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)

    start_elnum = mesh.shared_element_offsets[idx]

    for j=1:length(interfaces)
      iface_j = interfaces[j]
      fL = iface_j.faceL

      # interpolate to face
      qL = ro_sview(eqn.q, :, :, iface_j.elementL)
      el_r = iface_j.elementR - start_elnum + 1
      qR = ro_sview(qR_arr, :, :, el_r)
      interiorFaceInterpolate!(mesh.sbpface, iface_j, qL, qR, q_faceL, q_faceR)

      # calculate flux
      for k=1:mesh.numNodesPerFace
        qL_k = ro_sview(q_faceL, :, k)
        qR_k = ro_sview(q_faceR, :, k)
        aux_vars = ro_sview(aux_vars_arr, :, k, j)
        nrm_xy = ro_sview(nrm_arr, :, k, j)
        flux_k = sview(flux_face, :, k)

        parent(aux_vars)[1] = calcPressure(eqn.params, qL_k)

        flux_func(eqn.params, qL_k, qR_k, aux_vars, nrm_xy, flux_k)

        for p=1:mesh.numDofPerNode
          val[p] += mesh.sbpface.wface[k]*flux_k[p]
        end
      end  # end k
    end  # end j
  end # end idx

  return nothing
end


#------------------------------------------------------------------------------
# revq


function integrateFaceQuantity_revq(mesh, sbp, eqn::EulerData{Tsol, Tres}, opts,
                               flux_func::FluxType_revq, val_bar::Vector) where {Tsol, Tres}

  q_faceL = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
  q_faceR = zeros(q_faceL)

  q_faceL_bar = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)
  q_faceR_bar = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)

  flux_face_bar = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)
  interfaces = mesh.interfaces

  for i=1:length(interfaces)
    iface_i = interfaces[i]
    fill!(q_faceL_bar, 0); fill!(q_faceR_bar, 0)

    qL = ro_sview(eqn.q, :, :, iface_i.elementL)
    qR = ro_sview(eqn.q, :, :, iface_i.elementR)
    interiorFaceInterpolate!(mesh.sbpface, iface_i, qL, qR, q_faceL,
                             q_faceR)

    for j=1:mesh.numNodesPerFace
      qL_j = ro_sview(q_faceL, :, j); qR_j = ro_sview(q_faceR, :, j)
      qL_bar_j = sview(q_faceL_bar, :, j); qR_bar_j = sview(q_faceR_bar, :, j)
  

      eqn.aux_vars_face[1, j, i] = calcPressure(eqn.params, qL_j)
      aux_vars = ro_sview(eqn.aux_vars_face, :, j, i)

      nrm_xy = ro_sview(mesh.nrm_face, :, j, i)
      flux_j_bar = sview(flux_face_bar, :, j)

      for k=1:mesh.numDofPerNode
        flux_j_bar[k] = mesh.sbpface.wface[j]*val_bar[k]
      end

      flux_func(eqn.params, qL_j, qL_bar_j, qR_j, qR_bar_j,  aux_vars, nrm_xy, flux_j_bar)

    end  # end loop j

    qL_bar = sview(eqn.q_bar, :, :, iface_i.elementL)
    qR_bar = sview(eqn.q_bar, :, :, iface_i.elementR)
    interiorFaceInterpolate_rev!(mesh.sbpface, iface_i, qL_bar, qR_bar,
                                 q_faceL_bar, q_faceR_bar)

  end

  if opts["parallel_data"] == PARALLEL_DATA_FACE
    error("parallel face data not supported")
  elseif opts["parallel_data"] == PARALLEL_DATA_ELEMENT
    func = (mesh, sbp, eqn, opts, data, data_bar) -> integrateSharedFaceQuantity_element_revq(mesh, sbp, eqn, opts, data, data_bar, flux_func, val_bar)

  else
    error("unrecognized parallel data: $(opts["parallel_data"])")
  end

  # do each set of shared faces and send q_bar back to the owning process
  exchangeData_rev(mesh, sbp, eqn, opts, eqn.shared_data, eqn.shared_data_bar,
                    func)

  # wait for communication to finish
  finishSolutionBarExchange(mesh, sbp, eqn, opts)

  return nothing
end


function integrateSharedFaceQuantity_element_revq(mesh, sbp,
                  eqn::EulerData{Tsol, Tres}, opts, data::SharedFaceData,
                  data_bar::SharedFaceData,
                  flux_func::FluxType_revq,
                  val_bar::Vector) where {Tsol, Tres}

  if getParallelData(data) != PARALLEL_DATA_ELEMENT
    throw(ErrorException("cannot use integrateSharedFaceQuantity without parallel element data"))
  end


  # only one process should do the integral for a given face.
  # Arbitrarily choose the lower ranked process
  if data.peernum < mesh.myrank
    fill!(data_bar.q_recv, 0)
    return nothing
  end


  q_faceL = Array{Tsol}(mesh.numDofPerNode, mesh.numNodesPerFace)
  q_faceR = Array{Tsol}(mesh.numDofPerNode, mesh.numNodesPerFace)
  q_faceL_bar = Array{Tsol}(mesh.numDofPerNode, mesh.numNodesPerFace)
  q_faceR_bar = Array{Tsol}(mesh.numDofPerNode, mesh.numNodesPerFace)
 
  idx = data_bar.peeridx

  interfaces = data.interfaces
  qR_arr = data.q_recv
  qR_arr_bar = data_bar.q_recv
  nrm_arr = mesh.nrm_sharedface[idx]
  aux_vars_arr = eqn.aux_vars_sharedface[idx]

  flux_face_bar = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)

  start_elnum = mesh.shared_element_offsets[idx]

  for j=1:length(interfaces)
    iface_j = interfaces[j]
    fill!(q_faceL_bar, 0); fill!(q_faceR_bar, 0)
    fL = iface_j.faceL

    # interpolate to face
    qL = ro_sview(eqn.q, :, :, iface_j.elementL)
    el_r = iface_j.elementR - start_elnum + 1
    qR = ro_sview(qR_arr, :, :, el_r)
    interiorFaceInterpolate!(mesh.sbpface, iface_j, qL, qR, q_faceL, q_faceR)

    # calculate flux
    for k=1:mesh.numNodesPerFace
      qL_k = ro_sview(q_faceL, :, k); qR_k = ro_sview(q_faceR, :, k)
      qL_k_bar = sview(q_faceL_bar, :, k); qR_k_bar = sview(q_faceR_bar, :, k)
      aux_vars = ro_sview(aux_vars_arr, :, k, j)
      nrm_xy = ro_sview(nrm_arr, :, k, j)
      flux_k_bar = sview(flux_face_bar, :, k)

      parent(aux_vars)[1] = calcPressure(eqn.params, qL_k)

      for p=1:mesh.numDofPerNode
        flux_k_bar[p] = mesh.sbpface.wface[k]*val_bar[p]
      end

      flux_func(eqn.params, qL_k, qL_k_bar, qR_k, qR_k_bar, aux_vars, nrm_xy, flux_k_bar)
    end  # end k

    qL_bar = sview(eqn.q_bar, :, :, iface_j.elementL)
    qR_bar = sview(qR_arr_bar, :, :, el_r)
    interiorFaceInterpolate_rev!(mesh.sbpface, iface_j, qL_bar, qR_bar,
                                 q_faceL_bar, q_faceR_bar)
  end  # end j

  return nothing
end


#------------------------------------------------------------------------------
# revm

function integrateFaceQuantity_revm(mesh, sbp, eqn::EulerData{Tsol, Tres}, opts,
                               flux_func::FluxType_revm, val_bar::Vector) where {Tsol, Tres}

  q_faceL = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
  q_faceR = zeros(q_faceL)

  flux_face_bar = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)
  interfaces = mesh.interfaces

  for i=1:length(interfaces)
    iface_i = interfaces[i]

    qL = ro_sview(eqn.q, :, :, iface_i.elementL)
    qR = ro_sview(eqn.q, :, :, iface_i.elementR)
    interiorFaceInterpolate!(mesh.sbpface, iface_i, qL, qR, q_faceL,
                             q_faceR)

    for j=1:mesh.numNodesPerFace
      qL_j = ro_sview(q_faceL, :, j); qR_j = ro_sview(q_faceR, :, j)
  

      eqn.aux_vars_face[1, j, i] = calcPressure(eqn.params, qL_j)
      aux_vars = ro_sview(eqn.aux_vars_face, :, j, i)

      nrm_xy = ro_sview(mesh.nrm_face, :, j, i)
      nrm_bar = sview(mesh.nrm_face_bar, :, j, i)
      flux_j_bar = sview(flux_face_bar, :, j)

      for k=1:mesh.numDofPerNode
        flux_j_bar[k] = mesh.sbpface.wface[j]*val_bar[k]
      end

      flux_func(eqn.params, qL_j, qR_j,  aux_vars, nrm_xy, nrm_bar, flux_j_bar)

    end  # end loop j
  end

  # no need to do parallel communication here because we update
  # mesh.nrm_sharedface_bar, and the mesh reverse mode is responsible for its
  # reverse mode
  if opts["parallel_data"] == PARALLEL_DATA_FACE
    error("parallel face data not supported")
  elseif opts["parallel_data"] == PARALLEL_DATA_ELEMENT
    integrateSharedFaceQuantity_element_revm(mesh, sbp, eqn, opts, flux_func,
                                             val_bar)
  else
    error("unrecognized parallel data: $(opts["parallel_data"])")
  end

  return nothing
end


function integrateSharedFaceQuantity_element_revm(mesh, sbp,
                  eqn::EulerData{Tsol, Tres}, opts,
                  flux_func::FluxType_revm,
                  val_bar::Vector) where {Tsol, Tres}

  q_faceL = Array{Tsol}(mesh.numDofPerNode, mesh.numNodesPerFace)
  q_faceR = Array{Tsol}(mesh.numDofPerNode, mesh.numNodesPerFace)
 

  for idx=1:length(eqn.shared_data)
    data = eqn.shared_data[idx]

    if getParallelData(data) != PARALLEL_DATA_ELEMENT
      throw(ErrorException("cannot use integrateSharedFaceQuantity without parallel element data"))
    end

    # only one process should do the integral for a given face.
    # Arbitrarily choose the lower ranked process
    if data.peernum < mesh.myrank
      continue
    end

    interfaces = data.interfaces
    qR_arr = data.q_recv
    nrm_arr = mesh.nrm_sharedface[idx]
    nrm_arr_bar = mesh.nrm_sharedface_bar[idx]
    aux_vars_arr = eqn.aux_vars_sharedface[idx]

    flux_face_bar = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)

    start_elnum = mesh.shared_element_offsets[idx]

    for j=1:length(interfaces)
      iface_j = interfaces[j]
      fL = iface_j.faceL

      # interpolate to face
      qL = ro_sview(eqn.q, :, :, iface_j.elementL)
      el_r = iface_j.elementR - start_elnum + 1
      qR = ro_sview(qR_arr, :, :, el_r)
      interiorFaceInterpolate!(mesh.sbpface, iface_j, qL, qR, q_faceL, q_faceR)

      # calculate flux
      for k=1:mesh.numNodesPerFace
        qL_k = ro_sview(q_faceL, :, k); qR_k = ro_sview(q_faceR, :, k)
        aux_vars = ro_sview(aux_vars_arr, :, k, j)
        nrm_xy = ro_sview(nrm_arr, :, k, j)
        nrm_bar = sview(nrm_arr_bar, :, k, j)
        flux_k_bar = sview(flux_face_bar, :, k)

        parent(aux_vars)[1] = calcPressure(eqn.params, qL_k)

        for p=1:mesh.numDofPerNode
          flux_k_bar[p] = mesh.sbpface.wface[k]*val_bar[p]
        end

        flux_func(eqn.params, qL_k, qR_k, aux_vars, nrm_xy, nrm_bar, flux_k_bar)
      end  # end k
    end  # end j
  end  # end idx

  return nothing
end


