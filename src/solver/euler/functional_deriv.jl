# functional derivative calculation

import PDESolver.evalFunctionalDeriv

"""
  Computes a 3D array of hte derivative of a functional wrt eqn.q.

  The derivative is evaluated at the state in eqn.q_vec.

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * functionalData: AbstractIntegralFunctional to evaluate

  **Inputs/Outputs**

   * func_deriv_arr: array to hold derivative of function wrt eqn.q, same
                     size as equation.q

  **Options Keys**

  This funciton is not compatible with `precompute_q_bndry` = false
"""
function evalFunctionalDeriv(mesh::AbstractDGMesh{Tmsh}, 
                           sbp::AbstractSBP,
                           eqn::EulerData{Tsol}, opts,
                           functionalData::AbstractIntegralFunctional,
                           func_deriv_arr::Abstract3DArray) where {Tmsh, Tsol}

  array1DTo3D(mesh, sbp, eqn, opts, eqn.q_vec, eqn.q)
  if mesh.isDG
    boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, eqn.q, eqn.q_bndry)
  end

  calcFunctionalDeriv(mesh, sbp, eqn, opts, functionalData, func_deriv_arr)

  return nothing
end


@doc """
### EulerEquationMod. calcFunctionalDeriv

Computes a 3D array of the derivative of a functional w.r.t eqn.q on all
mesh nodes.

**Inputs**

*  `mesh` : Abstract DG mesh type
*  `sbp`  : Summation-By-parts operator
*  `eqn`  : Euler equation object
*  `opts` : Options dictionary
*  `functionalData` : Functional object of super-type AbstractFunctional
                      that is needed for computing the adjoint vector.
                      Depending on the functional being computed, a different
                      method based on functional type may be needed to be
                      defined.
*  `func_deriv_arr` : 3D array that stores the derivative of the functional
                      w.r.t. eqn.q. The array is the same size as eqn.q

**Outputs**

*  None

"""->
function calcFunctionalDeriv(mesh::AbstractDGMesh{Tmsh}, 
                           sbp::AbstractSBP,
                           eqn::EulerData{Tsol}, opts,
                           functionalData::AbstractIntegralFunctional,
                           func_deriv_arr::Abstract3DArray) where {Tmsh, Tsol}

  integrand = zeros(eqn.q_bndry)

  # Populate integrand
  for itr = 1:length(functionalData.bcnums)
    bcnum = functionalData.bcnums[itr]

    start_index = mesh.bndry_offsets[bcnum]
    end_index = mesh.bndry_offsets[bcnum+1]
    idx_range = start_index:(end_index-1)
    bndry_facenums = sview(mesh.bndryfaces, idx_range) # faces on geometric edge i

    nfaces = length(bndry_facenums)
#    q2 = zeros(Tsol, mesh.numDofPerNode)
    for i = 1:nfaces
      bndry_i = bndry_facenums[i]
      global_facenum = idx_range[i]
      for j = 1:mesh.sbpface.numnodes
        vtx_arr = mesh.topo.face_verts[:,bndry_i.face]
        q = sview(eqn.q_bndry, :, j, global_facenum)
#        convertToConservative(eqn.params, q, q2)
        aux_vars = ro_sview(eqn.aux_vars_bndry, :, j, global_facenum)
        x = ro_sview(mesh.coords_bndry, :, j, global_facenum)
        nrm = ro_sview(mesh.nrm_bndry, :, j, global_facenum)
        node_info = Int[itr,j,i]
        integrand_i = sview(integrand, :, j, global_facenum)

        calcIntegrandDeriv(opts, eqn.params, q, aux_vars, nrm, integrand_i, node_info,
                           functionalData)
      end  # End for j = 1:mesh.sbpface.numnodes
    end    # End for i = 1:nfaces
  end      # End for itr = 1:length(functional_edges)

  boundaryintegrate!(mesh.sbpface, mesh.bndryfaces, integrand, func_deriv_arr)

  return nothing
end  # End function calcFunctionalDeriv


"""
  Method for LiftCoefficient, which delegates most of the computation to the
  Lift functional.
"""
function calcFunctionalDeriv(mesh::AbstractDGMesh{Tmsh}, 
                           sbp::AbstractSBP,
                           eqn::EulerData{Tsol}, opts,
                           functionalData::LiftCoefficient,
                           func_deriv_arr::Abstract3DArray) where {Tmsh, Tsol}

  calcFunctionalDeriv(mesh, sbp, eqn, opts, functionalData.lift, func_deriv_arr)

  Ma = eqn.params.Ma
  fac = 0.5*eqn.params.rho_free*Ma*Ma

  scale!(func_deriv_arr, 1./fac)

  return nothing
end



@doc """
### EulerEquationMod.calcIntegrandDeriv

Compute the derivative of the functional Integrand at a node w.r.t all the
degrees of freedom at the node.

**Inputs**

*  `opts`   : Options dictionary
*  `params` : parameter type
*  `q`      : Solution variable at a node
*  `aux_vars` : Auxiliary variables
*  `nrm`    : normal vector in the physical space
*  `integrand_deriv` : Derivative of the integrand at that particular node
*  `node_info` : Tuple containing information about the node
*  `functionalData` : Functional object that is a subtype of AbstractFunctional.

**Outputs**

*  None

"""->

function calcIntegrandDeriv(opts, params::ParamType{2},
          q::AbstractArray{Tsol,1},
          aux_vars::AbstractArray{Tres, 1}, nrm::AbstractArray{Tmsh},
          integrand_deriv::AbstractArray{Tsol, 1}, node_info,
          functionalData::BoundaryForceData{Tsol,:lift}) where {Tsol, Tres, Tmsh}

  pert = complex(0, 1e-20)
  aoa = params.aoa
  momentum = zeros(Tsol,2)

  for i = 1:length(q)
    q[i] += pert
    calcBoundaryFunctionalIntegrand(params, q, aux_vars, nrm, node_info, functionalData, momentum)
    val = -momentum[1]*sin(aoa) + momentum[2]*cos(aoa)
    integrand_deriv[i] = imag(val)/norm(pert)
    q[i] -= pert
  end # End for i = 1:length(q)

  return nothing
end

function calcIntegrandDeriv(opts, params::ParamType{2},
          q::AbstractArray{Tsol,1},
          aux_vars::AbstractArray{Tres, 1}, nrm::AbstractArray{Tmsh},
          integrand_deriv::AbstractArray{Tsol, 1}, node_info,
          functionalData::BoundaryForceData{Tsol,:drag}) where {Tsol, Tres, Tmsh}

  pert = complex(0, 1e-20)
  aoa = params.aoa
  momentum = zeros(Tsol,2)

  for i = 1:length(q)
    q[i] += pert
    calcBoundaryFunctionalIntegrand(params, q, aux_vars, nrm, node_info, functionalData, momentum)
    val = momentum[1]*cos(aoa) + momentum[2]*sin(aoa)
    integrand_deriv[i] = imag(val)/norm(pert)
    q[i] -= pert
  end # End for i = 1:length(q)

  return nothing
end

function calcIntegrandDeriv(opts, params::ParamType{2},
          q::AbstractArray{Tsol,1},
          aux_vars::AbstractArray{Tres, 1},
          nrm::AbstractArray{Tmsh},
          integrand_deriv::AbstractArray{Tsol, 1}, node_info,
          functionalData::MassFlowData) where {Tsol, Tres, Tmsh}

  node_info = [1, 2, 3]  #TODO: what is this doing??? node_info is an argument
  h = 1e-20
  pert = Complex128(0, h)
  val = zeros(Complex128, 1)
  for i=1:length(q)
    q[i] += pert
    calcBoundaryFunctionalIntegrand(params, q, aux_vars, nrm, node_info,
                                          functionalData, val)
    integrand_deriv[i] = imag(val[1])/h
    q[i] -= pert
  end
  # functional integrand rho*v
#  integrand_deriv[2] = 1*nrm[1]
#  integrand_deriv[3] = 1*nrm[2]

  return nothing
end


#------------------------------------------------------------------------------
# generic fallback: complex step it

function calcIntegrandDeriv(opts, params::ParamType{2},
          q::AbstractArray{Tsol,1},
          aux_vars::AbstractArray{Tres, 1},
          nrm::AbstractArray{Tmsh},
          integrand_deriv::AbstractArray{Tsol, 1}, node_info,
          functionalData::Tfunc) where {Tsol, Tres, Tmsh, Tfunc<:AbstractIntegralFunctional}

  h = 1e-20
  pert = Complex128(0, h)
  val = zeros(Complex128, 1)
  for i=1:length(q)
    q[i] += pert
    calcBoundaryFunctionalIntegrand(params, q, aux_vars, nrm, node_info,
                                          functionalData, val)
    integrand_deriv[i] = imag(val[1])/h
    q[i] -= pert
  end

  return nothing
end


