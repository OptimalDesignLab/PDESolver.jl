# BoundaryEntropyDissipation functional

#------------------------------------------------------------------------------
# Boundary dissipation
function calcBndryFunctional(mesh::AbstractDGMesh{Tmsh},
            sbp::AbstractOperator, eqn::EulerData{Tsol, Tres}, opts,
            func::BoundaryEntropyDiss{Topt}) where {Tmsh, Tsol, Tres, Topt}

  val = zero(Topt)
  for i=1:mesh.numBC
    if i in func.bcnums
      functor_i = mesh.bndry_funcs[i]
      start_index = mesh.bndry_offsets[i]
      end_index = mesh.bndry_offsets[i+1] - 1
      idx_range = start_index:end_index
      bndry_facenums_i = sview(mesh.bndryfaces, idx_range)

      val += _calcBndryEntropyDiss(mesh, sbp, eqn, opts, func, functor_i, idx_range, bndry_facenums_i)
    end
  end

  return MPI.Allreduce(val, MPI.SUM, eqn.comm)
end

function calcBndryFunctional(mesh::AbstractDGMesh{Tmsh},
            sbp::AbstractOperator, eqn::EulerData{Tsol, Tres}, opts,
            func::NegBoundaryEntropyDiss{Topt}) where {Tmsh, Tsol, Tres, Topt}

  return -calcBndryFunctional(mesh, sbp, eqn, opts, func.func)
end

function calcBndryFunctional_revm(mesh::AbstractDGMesh{Tmsh},
     sbp::AbstractOperator, eqn::EulerData{Tsol, Tres, Tdim},
     opts, func::BoundaryEntropyDiss, _val_bar::Number=1,
     ) where {Tmsh, Tsol, Tres, Tdim}

  for i=1:mesh.numBC
    if i in func.bcnums
      functor_i = mesh.bndry_funcs_revm[i]
      start_index = mesh.bndry_offsets[i]
      end_index = mesh.bndry_offsets[i+1] - 1
      idx_range = start_index:end_index
      bndry_facenums_i = sview(mesh.bndryfaces, idx_range)

      _calcBndryEntropyDiss_revm(mesh, sbp, eqn, opts, func, functor_i, idx_range, bndry_facenums_i, _val_bar)
    end
  end

  return nothing
end

function calcBndryFunctional_revm(mesh::AbstractDGMesh{Tmsh},
     sbp::AbstractOperator, eqn::EulerData{Tsol, Tres, Tdim},
     opts, func::NegBoundaryEntropyDiss, _val_bar::Number=1,
     ) where {Tmsh, Tsol, Tres, Tdim}

  calcBndryFunctional_revm(mesh, sbp, eqn, opts, func.func, -_val_bar)
end


function calcFunctionalDeriv(mesh::AbstractDGMesh{Tmsh}, 
                           sbp::AbstractOperator,
                           eqn::EulerData{Tsol}, opts,
                           func::BoundaryEntropyDiss,
                           func_deriv_arr::Abstract3DArray,
                           _val_bar::Number=1) where {Tmsh, Tsol}

  for i=1:mesh.numBC
    if i in func.bcnums
      functor_i = mesh.bndry_funcs_revq[i]
      start_index = mesh.bndry_offsets[i]
      end_index = mesh.bndry_offsets[i+1] - 1
      idx_range = start_index:end_index
      bndry_facenums_i = sview(mesh.bndryfaces, idx_range)

      _calcBndryEntropyDiss_revq(mesh, sbp, eqn, opts, func, functor_i, idx_range, bndry_facenums_i, func_deriv_arr, _val_bar)
    end
  end

  return nothing
end


function calcFunctionalDeriv(mesh::AbstractDGMesh{Tmsh}, 
                           sbp::AbstractOperator,
                           eqn::EulerData{Tsol}, opts,
                           func::NegBoundaryEntropyDiss,
                           func_deriv_arr::Abstract3DArray,
                           _val_bar::Number=1) where {Tmsh, Tsol}

  calcFunctionalDeriv(mesh, sbp, eqn, opts, func.func, func_deriv_arr, -_val_bar)
end



function _calcBndryEntropyDiss(mesh::AbstractMesh{Tmsh},
            sbp::AbstractOperator, eqn::EulerData{Tsol, Tres}, opts,
            func::BoundaryEntropyDiss{Topt},
            functor::BCType,
            global_nums::UnitRange, faces::AbstractArray{Boundary}) where {Tmsh, Tsol, Tres, Topt}

  qface = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
  w_face_j = zeros(Tsol, mesh.numDofPerNode)
  q_bndry_j = zeros(Tsol, mesh.numDofPerNode)
  w_bndry_j = zeros(Tsol, mesh.numDofPerNode)
  delta_w = zeros(Tsol, mesh.numDofPerNode)
  q_avg = zeros(Tsol, mesh.numDofPerNode)
  aux_vars = Tres[0]
  kernel = eqn.params.entropy_lf_kernel
  flux = zeros(Tres, mesh.numDofPerNode)

  val = zero(Tres)
  for i=1:length(faces)  # loop over faces with this BC
    bndry_i = faces[i]
    global_facenum = global_nums[i]
    q_el = sview(eqn.q, :, :, bndry_i.element)
    boundaryFaceInterpolate!(mesh.sbpface, bndry_i.face, q_el, qface)

    val_i = zero(Tres)
    for j=1:mesh.numNodesPerFace

      # get components
      # convert to conservative variables if needed
      q_j = sview(qface, :, j)
      aux_vars[1] = calcPressure(eqn.params, q_j)
      coords = ro_sview(mesh.coords_bndry, :, j, global_facenum)
      nrm_xy = ro_sview(mesh.nrm_bndry, :, j, global_facenum)     
      bndry_node = BoundaryNode(bndry_i, i, j)

      getDirichletState(functor, eqn.params, q_j, aux_vars, coords, nrm_xy,
                        q_bndry_j)

      # compute Lambda(delta_w)
      convertToIR(eqn.params, q_j, w_face_j)
      convertToIR(eqn.params, q_bndry_j, w_bndry_j)
      for k=1:mesh.numDofPerNode
        q_avg[k]   = 0.5*(q_j[k] + q_bndry_j[k])
        delta_w[k] = w_face_j[k] - w_bndry_j[k]
      end

      applyEntropyKernel(kernel, eqn.params, q_avg, delta_w, nrm_xy, flux)

      # compute delta_w*B*Lambda(delta_w)
      for k=1:mesh.numDofPerNode
        val_i += delta_w[k]*mesh.sbpface.wface[j]*flux[k]
      end
    end  # end loop j
    println("contribution of boundary ", bndry_i, " = ", val_i)
    val += val_i
  end  # end loop i


  return -val
end


function _calcBndryEntropyDiss_revm(mesh::AbstractMesh{Tmsh},
            sbp::AbstractOperator, eqn::EulerData{Tsol, Tres}, opts,
            func::BoundaryEntropyDiss{Topt},
            functor::BCType_revm,
            global_nums::UnitRange, faces::AbstractArray{Boundary},
            val_bar::Number) where {Tmsh, Tsol, Tres, Topt}

  qface = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
  w_face_j = zeros(Tsol, mesh.numDofPerNode)
  q_bndry_j = zeros(Tsol, mesh.numDofPerNode)
  w_bndry_j = zeros(Tsol, mesh.numDofPerNode)
  delta_w = zeros(Tsol, mesh.numDofPerNode)
  q_avg = zeros(Tsol, mesh.numDofPerNode)
  aux_vars = Tres[0]
  kernel = eqn.params.entropy_lf_kernel
  flux = zeros(Tres, mesh.numDofPerNode)

  q_avg_bar = zeros(Tres, mesh.numDofPerNode)
  delta_w_bar = zeros(Tres, mesh.numDofPerNode)
  flux_bar = zeros(Tres, mesh.numDofPerNode)
  w_face_j_bar = zeros(Tres, mesh.numDofPerNode)
  w_bndry_j_bar = zeros(Tres, mesh.numDofPerNode)
  q_bndry_j_bar = zeros(Tres, mesh.numDofPerNode)
  qface_bar = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)
  A0inv = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode)

  val_bar = -val_bar

  for i=1:length(faces)  # loop over faces with this BC
    bndry_i = faces[i]
    global_facenum = global_nums[i]
    q_el = sview(eqn.q, :, :, bndry_i.element)
    boundaryFaceInterpolate!(mesh.sbpface, bndry_i.face, q_el, qface)

    for j=1:mesh.numNodesPerFace

      # get components
      # convert to conservative variables if needed
      q_j = sview(qface, :, j)
      q_j_bar = sview(qface_bar, :, j)
      aux_vars[1] = calcPressure(eqn.params, q_j)
      coords = ro_sview(mesh.coords_bndry, :, j, global_facenum)
      coords_bar = sview(mesh.coords_bndry_bar, :, j, global_facenum)
      nrm_xy = ro_sview(mesh.nrm_bndry, :, j, global_facenum)
      nrm_xy_bar = sview(mesh.nrm_bndry_bar, :, j, global_facenum)
      bndry_node = BoundaryNode(bndry_i, i, j)

      getDirichletState(functor, eqn.params, q_j, aux_vars, coords,
                        nrm_xy, q_bndry_j, bndry_node)

      # compute Lambda(delta_w)
      convertToIR(eqn.params, q_j, w_face_j)
      convertToIR(eqn.params, q_bndry_j, w_bndry_j)
      for k=1:mesh.numDofPerNode
        q_avg[k]   = 0.5*(q_j[k] + q_bndry_j[k])
        delta_w[k] = w_face_j[k] - w_bndry_j[k]
      end

      applyEntropyKernel(kernel, eqn.params, q_avg, delta_w, nrm_xy, flux)

      ## compute delta_w*B*Lambda(delta_w)
      #for k=1:mesh.numDofPerNode
      #  val += delta_w[k]*mesh.sbpface.wface[j]*flux[k]
      #end

      # reverse sweep
      fill!(delta_w_bar, 0); fill!(flux_bar, 0); fill!(q_avg_bar, 0)
      fill!(w_face_j_bar, 0); fill!(w_bndry_j_bar, 0); fill!(q_bndry_j_bar, 0)
      fill!(q_j_bar, 0)
      for k=1:mesh.numDofPerNode
        delta_w_bar[k] += mesh.sbpface.wface[j]*flux[k]*val_bar
        flux_bar[k]    += delta_w[k]*mesh.sbpface.wface[j]*val_bar
      end

      applyEntropyKernel_revm(kernel, eqn.params, q_avg, delta_w, nrm_xy,
                              nrm_xy_bar, flux, flux_bar)
      applyEntropyKernel_revq(kernel, eqn.params, q_avg, q_avg_bar, delta_w,
                              delta_w_bar, nrm_xy, flux, flux_bar)

      for k=1:mesh.numDofPerNode
        #w_face_j_bar[k] += delta_w_bar[k]
        w_bndry_j_bar[k] -= delta_w_bar[k]

        #q_j_bar[k]      += 0.5*q_avg_bar[k]
        q_bndry_j_bar[k] += 0.5*q_avg_bar[k]
      end

      # reverse of convertToIR
      #getIRA0inv(eqn.params, q_j, A0inv)
      #smallmatTvec_kernel!(A0inv, w_face_j_bar, q_j_bar, 1, 1)
      getIRA0inv(eqn.params, q_bndry_j, A0inv)
      smallmatTvec_kernel!(A0inv, w_bndry_j_bar, q_bndry_j_bar, 1, 1)

      getDirichletState_revm(functor, eqn.params, q_j, aux_vars,
                             coords, coords_bar, nrm_xy, nrm_xy_bar,
                             q_bndry_j_bar, bndry_node)

    end  # end loop j
  end  # end loop i


  return nothing
end




function _calcBndryEntropyDiss_revq(mesh::AbstractMesh{Tmsh},
            sbp::AbstractOperator, eqn::EulerData{Tsol, Tres}, opts,
            func::BoundaryEntropyDiss{Topt},
            functor::BCType_revq,
            global_nums::UnitRange, faces::AbstractArray{Boundary},
            func_deriv_arr::Abstract3DArray,
            val_bar::Number) where {Tmsh, Tsol, Tres, Topt}

  qface = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
  w_face_j = zeros(Tsol, mesh.numDofPerNode)
  q_bndry_j = zeros(Tsol, mesh.numDofPerNode)
  w_bndry_j = zeros(Tsol, mesh.numDofPerNode)
  delta_w = zeros(Tsol, mesh.numDofPerNode)
  q_avg = zeros(Tsol, mesh.numDofPerNode)
  aux_vars = Tres[0]
  kernel = eqn.params.entropy_lf_kernel
  flux = zeros(Tres, mesh.numDofPerNode)

  q_avg_bar = zeros(Tres, mesh.numDofPerNode)
  delta_w_bar = zeros(Tres, mesh.numDofPerNode)
  flux_bar = zeros(Tres, mesh.numDofPerNode)
  w_face_j_bar = zeros(Tres, mesh.numDofPerNode)
  w_bndry_j_bar = zeros(Tres, mesh.numDofPerNode)
  q_bndry_j_bar = zeros(Tres, mesh.numDofPerNode)
  qface_bar = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)
  A0inv = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode)

  val_bar = -val_bar
  for i=1:length(faces)  # loop over faces with this BC
    bndry_i = faces[i]
    global_facenum = global_nums[i]
    q_el = sview(eqn.q, :, :, bndry_i.element)
    boundaryFaceInterpolate!(mesh.sbpface, bndry_i.face, q_el, qface)

    for j=1:mesh.numNodesPerFace

      # get components
      # convert to conservative variables if needed
      q_j = sview(qface, :, j)
      q_j_bar = sview(qface_bar, :, j)
      aux_vars[1] = calcPressure(eqn.params, q_j)
      coords = ro_sview(mesh.coords_bndry, :, j, global_facenum)
      nrm_xy = ro_sview(mesh.nrm_bndry, :, j, global_facenum)
      bndry_node = BoundaryNode(bndry_i, i, j)

      getDirichletState(functor, eqn.params, q_j, aux_vars, coords, nrm_xy,
                        q_bndry_j, bndry_node)

      # compute Lambda(delta_w)
      convertToIR(eqn.params, q_j, w_face_j)
      convertToIR(eqn.params, q_bndry_j, w_bndry_j)
      for k=1:mesh.numDofPerNode
        q_avg[k]   = 0.5*(q_j[k] + q_bndry_j[k])
        delta_w[k] = w_face_j[k] - w_bndry_j[k]
      end

      applyEntropyKernel(kernel, eqn.params, q_avg, delta_w, nrm_xy, flux)

      ## compute delta_w*B*Lambda(delta_w)
      #for k=1:mesh.numDofPerNode
      #  val += delta_w[k]*mesh.sbpface.wface[j]*flux[k]
      #end

      # reverse sweep
      fill!(delta_w_bar, 0); fill!(flux_bar, 0); fill!(q_avg_bar, 0)
      fill!(w_face_j_bar, 0); fill!(w_bndry_j_bar, 0); fill!(q_bndry_j_bar, 0)
      fill!(q_j_bar, 0)
      for k=1:mesh.numDofPerNode
        delta_w_bar[k] += mesh.sbpface.wface[j]*flux[k]*val_bar
        flux_bar[k]    += delta_w[k]*mesh.sbpface.wface[j]*val_bar
      end

      applyEntropyKernel_revq(kernel, eqn.params, q_avg, q_avg_bar, delta_w,
                              delta_w_bar, nrm_xy, flux, flux_bar)

      for k=1:mesh.numDofPerNode
        w_face_j_bar[k]  += delta_w_bar[k]
        w_bndry_j_bar[k] -= delta_w_bar[k]

        q_j_bar[k]       += 0.5*q_avg_bar[k]
        q_bndry_j_bar[k] += 0.5*q_avg_bar[k]
      end

      # reverse of convertToIR
      getIRA0inv(eqn.params, q_j, A0inv)
      smallmatTvec_kernel!(A0inv, w_face_j_bar, q_j_bar, 1, 1)
      getIRA0inv(eqn.params, q_bndry_j, A0inv)
      smallmatTvec_kernel!(A0inv, w_bndry_j_bar, q_bndry_j_bar, 1, 1)

      getDirichletState_revq(functor, eqn.params, q_j, q_j_bar, aux_vars,
                             coords, nrm_xy, q_bndry_j_bar, bndry_node)

    end  # end loop j

    q_el_bar = sview(func_deriv_arr, :, :, bndry_i.element)
    boundaryFaceInterpolate_rev!(mesh.sbpface, bndry_i.face, q_el_bar, qface_bar)

  end  # end loop i


  return nothing
end


