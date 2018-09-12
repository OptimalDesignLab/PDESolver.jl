# evaluating AbstractFunctionals (those that are not AbstractIntegralFunctionals)

#------------------------------------------------------------------------------
# API functions

"""
  Evaluates [`AbstractFunctional`](@ref)s that are one of the more specific
  types of functional.

  Note that this function may overwrite eqn.res.


  **Keyword Arguments**

   * start_comm: start parallel communication, default true.
"""
function evalFunctional(mesh::AbstractMesh{Tmsh},
            sbp::AbstractSBP, eqn::EulerData{Tsol, Tres}, opts,
            functionalData::EntropyDissipationData; start_comm=true) where {Tmsh, Tsol, Tres}

  #TODO: figure out what the generalization is here

  array1DTo3D(mesh, sbp, eqn, opts, eqn.q_vec, eqn.q)
  if start_comm
    # TODO: wait=false
    startSolutionExchange(mesh, sbp, eqn, opts, wait=true)
  end

  val = calcFunctional(mesh, sbp, eqn, opts, functionalData)

  return val
end

"""
  Derivative for [`AbstractFunctional`](@ref)s that are not one of the more
  specific types.

  Currently this requires that the eqn object was creates with Tsol = Complex128
"""
function evalFunctionalDeriv(mesh::AbstractDGMesh{Tmsh}, 
                           sbp::AbstractSBP,
                           eqn::EulerData{Tsol, Tres}, opts,
                           functionalData::EntropyDissipationData,
                           func_deriv_arr::Abstract3DArray) where {Tmsh, Tsol, Tres}

  @assert size(func_deriv_arr, 1) == mesh.numDofPerNode
  @assert size(func_deriv_arr, 2) == mesh.numNodesPerElement
  @assert size(func_deriv_arr, 3) == mesh.numEl

  array1DTo3D(mesh, sbp, eqn, opts, eqn.q_vec, eqn.q)

  #TODO: use coloring at least
  # complex step
  h = 1e-20
  pert = Complex128(0, h)
  for i=1:length(eqn.q)
    eqn.q[i] += pert
    J_i = evalFunctional(mesh, sbp, eqn, opts, functionalData)
    func_deriv_arr[i] = imag(J_i)/h
    eqn.q[i] -= pert
  end

  return nothing
end



#------------------------------------------------------------------------------
# Implementation for each functional

function calcFunctional(mesh::AbstractMesh{Tmsh},
            sbp::AbstractSBP, eqn::EulerData{Tsol, Tres}, opts,
            functionalData::EntropyDissipationData) where {Tmsh, Tsol, Tres}


  fill!(eqn.res, 0.0)

  # local part
  face_integral_functor = ELFPenaltyFaceIntegral()
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
