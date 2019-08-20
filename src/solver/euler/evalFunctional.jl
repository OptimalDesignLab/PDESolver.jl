# evaluating EntropyPenaltyFunctionals

#------------------------------------------------------------------------------
# API functions

"""
  Evaluates [`EntropyPenaltyFunctional`](@ref)s

  Note that this function may overwrite eqn.res.
"""
function _evalFunctional(mesh::AbstractMesh{Tmsh},
            sbp::AbstractOperator, eqn::EulerData{Tsol, Tres}, opts,
            func::EntropyPenaltyFunctional) where {Tmsh, Tsol, Tres}

  val = calcFunctional(mesh, sbp, eqn, opts, func)
  val = MPI.Allreduce(val, MPI.SUM, eqn.comm)

  return val
end


#------------------------------------------------------------------------------
# derivative wrt q

#TODO: EntropyPenaltyFunctional should be EntropyDissipationData?
# derivative of functional wrt q
function _evalFunctionalDeriv_q(mesh::AbstractDGMesh{Tmsh}, 
                           sbp::AbstractOperator,
                           eqn::EulerData{Tsol, Tres}, opts,
                           func::EntropyPenaltyFunctional,
                           func_deriv_arr::Abstract3DArray) where {Tmsh, Tsol, Tres}

  # this is a trick to populate eqn.res (the forward sweep of reverse mode)
  calcFunctional(mesh, sbp, eqn, opts, func)

  # compute reverse mode of the contraction, take val_bar = 1
  w_j = zeros(Tsol, mesh.numDofPerNode)
  w_bar_j = zeros(Tres, mesh.numDofPerNode)
  fill!(eqn.res_bar, 0); fill!(eqn.q_bar, 0)
  A0inv = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode)

  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      q_j = sview(eqn.q, :, j, i)
      q_bar_j = sview(eqn.q_bar, :, j, i)
      fill!(w_bar_j, 0)
      convertToIR(eqn.params, q_j, w_j)
      for k=1:mesh.numDofPerNode
        #val += w_j[k]*eqn.res[k, j, i]
        #-----------------------------
        # reverse sweep
        w_bar_j[k] += eqn.res[k, j, i]
        eqn.res_bar[k, j, i] += w_j[k]
      end

      getIRA0inv(eqn.params, q_j, A0inv)

      smallmatTvec_kernel!(A0inv, w_bar_j, q_bar_j, 1, 1)
    end
  end
    
  if func.do_face_term
    setParallelData(eqn.shared_data_bar, PARALLEL_DATA_ELEMENT)
    startSolutionExchange_rev2(mesh, sbp, eqn, opts, send_q=false)


    # do reverse mode of the face integrals
    if typeof(mesh.sbpface) <: DenseFace
      @assert opts["parallel_data"] == PARALLEL_DATA_ELEMENT

      # local part
      face_integral_functor = func.func
      flux_functor = ErrorFlux_revq()  # not used, but required by the interface
      getFaceElementIntegral_revq(mesh, sbp, eqn, face_integral_functor, flux_functor, mesh.sbpface, mesh.interfaces)

      # define anonymous function for parallel faces
      calc_func = (mesh, sbp, eqn, opts, data, data_bar) ->
        calcSharedFaceElementIntegrals_element_inner_revq(mesh, sbp, eqn, opts,
                              data, data_bar, face_integral_functor, flux_functor)
    else # SparseFace
      flux_functor_revq = func.func_sparseface_revq
      calcFaceIntegral_nopre_revq(mesh, sbp, eqn, opts, flux_functor_revq,
                                  mesh.interfaces)
      # define anonymous functions for parallel faces
      calc_func = (mesh, sbp, eqn, opts, data, data_bar) ->
        calcSharedFaceIntegrals_nopre_element_inner_revq(mesh, sbp, eqn, opts,
          data, data_bar, flux_functor_revq)

    end
  end

  if func.do_lps && opts["addStabilization"]
    addStabilization_revq(mesh, sbp, eqn, opts)
  end

  if func.do_sc && opts["addShockCapturing"]
    evalShockCapturing_revq(mesh, sbp, eqn, opts)
  end

  # do the parallel face calculations
  if func.do_face_term
    finishExchangeData_rev2(mesh, sbp, eqn, opts, eqn.shared_data, eqn.shared_data_bar, calc_func)
  end

  copy!(func_deriv_arr, eqn.q_bar)

  return nothing
end


# wrapper method for NegEntropyDissipationData
function _evalFunctionalDeriv_q(mesh::AbstractDGMesh{Tmsh}, 
                           sbp::AbstractOperator,
                           eqn::EulerData{Tsol, Tres}, opts,
                           func::NegEntropyDissipations,
                           func_deriv_arr::Abstract3DArray) where {Tmsh, Tsol, Tres}

  _evalFunctionalDeriv_q(mesh, sbp, eqn, opts, func.func, func_deriv_arr)
  scale!(func_deriv_arr, -1)

  return nothing
end


# method for EntropyDissipaiton2Data
function _evalFunctionalDeriv_q(mesh::AbstractDGMesh{Tmsh}, 
                           sbp::AbstractOperator,
                           eqn::EulerData{Tsol, Tres}, opts,
                           func::EntropyDissipation2Data,
                           func_deriv_arr::Abstract3DArray) where {Tmsh, Tsol, Tres}


  _evalFunctionalDeriv_q(mesh, sbp, eqn, opts, func.flux_functional,
                        func_deriv_arr)

  fill!(eqn.q_bar, 0)
  flux_func = func.potential_flux_revq
  val_bar = zeros(Tres, mesh.numDofPerNode); val_bar[1] = 1
  integrateFaceQuantity_revq(mesh, sbp, eqn, opts, flux_func, val_bar)

  for i=1:length(func_deriv_arr)
    func_deriv_arr[i] += eqn.q_bar[i]
  end

  return nothing
end


function _evalFunctionalDeriv_q(mesh::AbstractDGMesh{Tmsh}, 
                           sbp::AbstractOperator,
                           eqn::EulerData{Tsol, Tres}, opts,
                           func::TotalEntropyDissipationData,
                           func_deriv_arr::Abstract3DArray) where {Tmsh, Tsol, Tres}

  fill!(eqn.res, 0)
  bcs_orig = setZeroBCs(mesh, sbp, eqn, opts, func)
  evalResidual(mesh, sbp, eqn, opts)

  # compute reverse mode of the contraction, take val_bar = 1
  w_j = zeros(Tsol, mesh.numDofPerNode)
  w_bar_j = zeros(Tres, mesh.numDofPerNode)
  fill!(eqn.res_bar, 0); fill!(eqn.q_bar, 0)
  A0inv = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode)

  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      q_j = sview(eqn.q, :, j, i)
      q_bar_j = sview(eqn.q_bar, :, j, i)
      fill!(w_bar_j, 0)
      convertToIR(eqn.params, q_j, w_j)
      for k=1:mesh.numDofPerNode
        #val += w_j[k]*eqn.res[k, j, i]
        #-----------------------------
        # reverse sweep
        w_bar_j[k] += eqn.res[k, j, i]
        eqn.res_bar[k, j, i] += w_j[k]
      end

      getIRA0inv(eqn.params, q_j, A0inv)

      smallmatTvec_kernel!(A0inv, w_bar_j, q_bar_j, 1, 1)
    end
  end

  # back propigate res_bar -> q_bar
  # This partially duplicates evalResidual_revq, but that function
  # produces a vector, whereas this function requires a 3D array as output.
  # That was probably a bad design decision, but we are stuck with it now
  setParallelData(eqn.shared_data_bar, PARALLEL_DATA_ELEMENT)
  startSolutionExchange_rev2(mesh, sbp, eqn, opts, send_q=false)

  evalResidual_revq(mesh, sbp, eqn, opts, 0.0)
  assertReceivesWaited(eqn.shared_data_bar)

  copy!(func_deriv_arr, eqn.q_bar)

  resetBCs(mesh, sbp, eqn, opts, bcs_orig)

  return nothing
end



#------------------------------------------------------------------------------
# derivative wrt metrics

function _evalFunctionalDeriv_m(mesh::AbstractDGMesh{Tmsh}, 
                           sbp::AbstractOperator,
                           eqn::AbstractSolutionData{Tsol}, opts,
                           func::EntropyPenaltyFunctional,
                           val_bar::Number=1,
                           ) where {Tmsh, Tsol}

  # compute reverse mode of the contraction, take val_bar = 1
  w_j = zeros(Tsol, mesh.numDofPerNode)
  fill!(eqn.res_bar, 0);

  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      q_j = sview(eqn.q, :, j, i)
      convertToIR(eqn.params, q_j, w_j)
      for k=1:mesh.numDofPerNode
        #val += w_j[k]*eqn.res[k, j, i]
        #-----------------------------
        # reverse sweep
        eqn.res_bar[k, j, i] += w_j[k]*val_bar
      end
    end
  end


  # do reverse mode of the face integrals
  if func.do_face_term
    if typeof(mesh.sbpface) <: DenseFace
      @assert opts["parallel_data"] == PARALLEL_DATA_ELEMENT

      face_integral_functor = func.func
      flux_functor = ErrorFlux()  # not used, but required by the interface

      # local part
      getFaceElementIntegral_revm(mesh, sbp, eqn, face_integral_functor, flux_functor, mesh.sbpface, mesh.interfaces)
   
      # define anonymous function for shared faces
      calc_func = (mesh, sbp, eqn, opts, data) ->
        calcSharedFaceElementIntegrals_element_inner_revm(mesh, sbp, eqn, opts,
                              data, face_integral_functor, flux_functor)
    else # SparseFace
      flux_functor_revm = func.func_sparseface_revm
      calcFaceIntegral_nopre_revm(mesh, sbp, eqn, opts, flux_functor_revm,
                                  mesh.interfaces)

      # anonymous function for shared faces
      calc_func = (mesh, sbp, eqn, opts, data) ->
        calcSharedFaceIntegrals_nopre_element_inner_revm(mesh, sbp, eqn, opts,
          data, flux_functor_revm)
    end
  end

  if func.do_lps && opts["addStabilization"]
    addStabilization_revm(mesh, sbp, eqn, opts)
  end

  if func.do_sc && opts["addShockCapturing"]
    evalShockCapturing_revm(mesh, sbp, eqn, opts)
  end


  # do the parallel part of the computation
  if func.do_face_term
    finishExchangeData(mesh, sbp, eqn, opts, eqn.shared_data, calc_func)
  end

  return nothing
end


# NegEntropyDissipationData method
function _evalFunctionalDeriv_m(mesh::AbstractDGMesh{Tmsh}, 
                           sbp::AbstractOperator,
                           eqn::AbstractSolutionData{Tsol}, opts,
                           func::NegEntropyDissipations,
                           val_bar::Number=1,
                           ) where {Tmsh, Tsol}

  _evalFunctionalDeriv_m(mesh, sbp, eqn, opts, func.func, -val_bar)

  return nothing
end


# method for EntropyDissipaiton2Data
function _evalFunctionalDeriv_m(mesh::AbstractDGMesh{Tmsh}, 
                           sbp::AbstractOperator,
                           eqn::EulerData{Tsol, Tres}, opts,
                           func::EntropyDissipation2Data,
                           val_bar::Number=1) where {Tmsh, Tsol, Tres}


  # back propigate val1 contribution
  _evalFunctionalDeriv_m(mesh, sbp, eqn, opts, func.flux_functional,
                        val_bar)

  # back propigate val2 contribution
  flux_func = func.potential_flux_revm
  vals_bar = zeros(Tres, mesh.numDofPerNode); vals_bar[1] = val_bar
  integrateFaceQuantity_revm(mesh, sbp, eqn, opts, flux_func, vals_bar)

  return nothing
end


function _evalFunctionalDeriv_m(mesh::AbstractDGMesh{Tmsh}, 
                           sbp::AbstractOperator,
                           eqn::AbstractSolutionData{Tsol}, opts,
                           func::TotalEntropyDissipationData,
                           val_bar::Number=1,
                           ) where {Tmsh, Tsol}

  # compute reverse mode of the contraction, take val_bar = 1
  w_j = zeros(Tsol, mesh.numDofPerNode)
  fill!(eqn.res_bar, 0);

  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      q_j = sview(eqn.q, :, j, i)
      convertToIR(eqn.params, q_j, w_j)
      for k=1:mesh.numDofPerNode
        #val += w_j[k]*eqn.res[k, j, i]
        #-----------------------------
        # reverse sweep
        eqn.res_bar[k, j, i] += w_j[k]*val_bar
      end
    end
  end

  bcs_orig = setZeroBCs(mesh, sbp, eqn, opts, func)
  evalResidual_revm(mesh, sbp, eqn, opts, 0.0)
  resetBCs(mesh, sbp, eqn, opts, bcs_orig)

  return nothing
end



#------------------------------------------------------------------------------
# Implementation for each functional

function calcFunctional(mesh::AbstractMesh{Tmsh},
            sbp::AbstractOperator, eqn::EulerData{Tsol, Tres}, opts,
            func::EntropyPenaltyFunctional) where {Tmsh, Tsol, Tres}


  fill!(eqn.res, 0.0)
  # use functions called from evalResidual to compute the dissipation, then
  # contract it with the entropy variables afterwards

  if func.do_face_term
    if typeof(mesh.sbpface) <: DenseFace
      @assert opts["parallel_data"] == PARALLEL_DATA_ELEMENT

      # local part
      face_integral_functor = func.func
      flux_functor = ErrorFlux()  # not used, but required by the interface
      getFaceElementIntegral(mesh, sbp, eqn, face_integral_functor, flux_functor, mesh.sbpface, mesh.interfaces)

      # parallel part

      # temporarily change the face integral functor
      face_integral_functor_orig = eqn.face_element_integral_func
      eqn.face_element_integral_func = face_integral_functor

      # finish data exchange/do the shared face integrals
      # this works even if the receives have already been waited on
      finishExchangeData(mesh, sbp, eqn, opts, eqn.shared_data, calcSharedFaceElementIntegrals_element)

      # restore eqn object to original state
      eqn.face_element_integral_func = face_integral_functor_orig
    else  # SparseFace
      flux_functor = func.func_sparseface
      calcFaceIntegral_nopre(mesh, sbp, eqn, opts, flux_functor, mesh.interfaces)

      # parallel part
      doSparseFaceSharedFaceIntegrals(mesh, sbp, eqn, opts, flux_functor)
    end
  end

  if func.do_lps && opts["addStabilization"]
    addStabilization(mesh, sbp, eqn, opts)
  end

  if func.do_sc && opts["addShockCapturing"]
    evalShockCapturing(mesh, sbp, eqn, opts)
  end

  # compute the contraction
  val = zero(Tres)
  w_j = zeros(Tsol, mesh.numDofPerNode)
  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      q_j = sview(eqn.q, :, j, i)
      convertToIR(eqn.params, q_j, w_j)
      for k=1:mesh.numDofPerNode
        val += w_j[k]*eqn.res[k, j, i]

      end
    end
  end

  return val
end

function calcFunctional(mesh::AbstractMesh{Tmsh},
            sbp::AbstractOperator, eqn::EulerData{Tsol, Tres}, opts,
            func::NegEntropyDissipations) where {Tmsh, Tsol, Tres}

  return -calcFunctional(mesh, sbp, eqn, opts, func.func)
end


function calcFunctional(mesh::AbstractMesh{Tmsh},
            sbp::AbstractOperator, eqn::EulerData{Tsol, Tres}, opts,
            func::EntropyDissipation2Data) where {Tmsh, Tsol, Tres}

  # get -\int (wL - wR)^T * F(uL, uR)  (negative sign because the face term
  # has a negative sign in the weak form.
  val1 = calcFunctional(mesh, sbp, eqn, opts, func.flux_functional)

  # now compute \int (psi_L - psi_R)
  # We can't use the regular face integral machinery for this, because it
  # is conservative: it computes resL += flux, resR -= flux, so the sum is
  # always zero.  Instead we need to integral and sum a quantity at the
  # face nodes
  # this only need to work for diagonal E

  flux_func = func.potential_flux
  val = zeros(Tres, mesh.numDofPerNode)
  integrateFaceQuantity(mesh, sbp, eqn, opts, flux_func, val)

  val2 = val[1]
  return val2 + val1
end


function calcFunctional(mesh::AbstractMesh{Tmsh},
            sbp::AbstractOperator, eqn::EulerData{Tsol, Tres}, opts,
            func::TotalEntropyDissipationData) where {Tmsh, Tsol, Tres}

  fill!(eqn.res, 0)
  # zero out all BCs except those specified
  bcs_orig = setZeroBCs(mesh, sbp, eqn, opts, func)

  evalResidual(mesh, sbp, eqn, opts)

  # reset original BCs
  resetBCs(mesh, sbp, eqn, opts, bcs_orig)

  # compute the contraction
  val = zero(Tres)
  w_j = zeros(Tsol, mesh.numDofPerNode)
  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      q_j = sview(eqn.q, :, j, i)
      convertToIR(eqn.params, q_j, w_j)
      for k=1:mesh.numDofPerNode
        val += w_j[k]*eqn.res[k, j, i]

      end
    end
  end

  return val
end

struct BCSet
  bndry_funcs::Vector{BCType}
  bndry_funcs_revm::Vector{BCType_revm}
  bndry_funcs_revq::Vector{BCType_revq}
end

function BCSet(n::Integer)

  bndry_funcs = Array{BCType}(n)
  bndry_funcs_revm = Array{BCType_revm}(n)
  bndry_funcs_revq = Array{BCType_revq}(n)

  return BCSet(bndry_funcs, bndry_funcs_revm, bndry_funcs_revq)
end

"""
  Sets the boundary flux functors to be ZeroFlux for all boundaries except
  those in func.bcnums.  Also does the revq and revm ones.
"""
function setZeroBCs(mesh, sbp, eqn, opts, func::TotalEntropyDissipationData)

  zfunc = ZeroFluxBC(mesh, eqn)
  zfunc_revq = ZeroFluxBC_revq(mesh, eqn)
  zfunc_revm = ZeroFluxBC_revm(mesh, eqn)

  nbcs = length(mesh.bndry_funcs)
  bcs_orig = BCSet(nbcs)


  for i=1:nbcs
    bcs_orig.bndry_funcs[i] = mesh.bndry_funcs[i]
    bcs_orig.bndry_funcs_revm[i] = mesh.bndry_funcs_revm[i]
    bcs_orig.bndry_funcs_revq[i] = mesh.bndry_funcs_revq[i]
    if !(i in func.bcnums)
      mesh.bndry_funcs[i] = zfunc
      mesh.bndry_funcs_revm[i] = zfunc_revm
      mesh.bndry_funcs_revq[i] = zfunc_revq
    end
  end

  return bcs_orig
end

function resetBCs(mesh, sbp, eqn, opts, bcs_orig::BCSet)

  nbcs = length(mesh.bndry_funcs)
  for i=1:nbcs
    mesh.bndry_funcs[i] = bcs_orig.bndry_funcs[i]
    mesh.bndry_funcs_revm[i] = bcs_orig.bndry_funcs_revm[i]
    mesh.bndry_funcs_revq[i] = bcs_orig.bndry_funcs_revq[i]
  end

  return nothing
end






"""
  Do the shared face integrals with the given flux function
"""
function doSparseFaceSharedFaceIntegrals(mesh, sbp, eqn, opts, flux_functor::FluxType)

  # figure out which parallel function to call
  if opts["parallel_data"] == PARALLEL_DATA_FACE
    pfunc = (mesh, sbp, eqn, opts, data) -> calcSharedFaceIntegrals_nopre_inner(mesh, sbp, eqn, opts, data, flux_functor)
  elseif opts["parallel_data"] == PARALLEL_DATA_ELEMENT
    pfunc = (mesh, sbp, eqn, opts, data) -> calcSharedFaceIntegrals_nopre_element_inner(mesh, sbp, eqn, opts, data, flux_functor)
  else
    error("unregonized parallel data: $(opts["parallel_data"])")
  end

  # do shared face integrals
  finishExchangeData(mesh, sbp, eqn, opts, eqn.shared_data, pfunc)

  return nothing
end
