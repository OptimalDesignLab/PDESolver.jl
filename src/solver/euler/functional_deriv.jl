# functional derivative calculation

import PDESolver._evalFunctionalDeriv_q

"""
  Computes a 3D array of hte derivative of a functional wrt eqn.q.

  The derivative is evaluated at the state in eqn.q_vec.

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * func: AbstractBoundaryFunctional to evaluate

  **Inputs/Outputs**

   * func_deriv_arr: array to hold derivative of function wrt eqn.q, same
                     size as equation.q

  **Options Keys**

  This function is not compatible with `precompute_q_bndry` = false
"""
function _evalFunctionalDeriv_q(mesh::AbstractDGMesh{Tmsh}, 
                           sbp::AbstractOperator,
                           eqn::EulerData{Tsol}, opts,
                           func::AbstractBoundaryFunctional,
                           func_deriv_arr::Abstract3DArray) where {Tmsh, Tsol}

  if mesh.isDG
    boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, eqn.q, eqn.q_bndry)
  end

  calcFunctionalDeriv(mesh, sbp, eqn, opts, func, func_deriv_arr)

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
*  `func` : Functional object of super-type AbstractBoundaryFunctional
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
                           sbp::AbstractOperator,
                           eqn::EulerData{Tsol}, opts,
                           func::AbstractBoundaryFunctional,
                           func_deriv_arr::Abstract3DArray) where {Tmsh, Tsol}

  integrand = zeros(eqn.q_bndry)  #TODO: only allocate enough space for
                                  #      the boundaries the functional is using
  node_info = Array{Int}(3)

  # Populate integrand
  for itr = 1:length(func.bcnums)
    bcnum = func.bcnums[itr]

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
        node_info[1] = itr; node_info[2] = j; node_info[3] = i
        integrand_i = sview(integrand, :, j, global_facenum)

        calcIntegrandDeriv(opts, eqn.params, q, aux_vars, nrm, integrand_i, node_info,
                           func)
      end  # End for j = 1:mesh.sbpface.numnodes
    end    # End for i = 1:nfaces
  end      # End for itr = 1:length(functional_edges)

  fill!(func_deriv_arr, 0)
  boundaryintegrate!(mesh.sbpface, mesh.bndryfaces, integrand, func_deriv_arr)

  return nothing
end  # End function calcFunctionalDeriv


"""
  Method for LiftCoefficient, which delegates most of the computation to the
  Lift functional.
"""
function calcFunctionalDeriv(mesh::AbstractDGMesh{Tmsh}, 
                           sbp::AbstractOperator,
                           eqn::EulerData{Tsol}, opts,
                           func::AeroCoefficients,
                           func_deriv_arr::Abstract3DArray) where {Tmsh, Tsol}

  calcFunctionalDeriv(mesh, sbp, eqn, opts, func.func, func_deriv_arr)

  Ma = eqn.params.Ma
  fac = 0.5*eqn.params.rho_free*Ma*Ma

  scale!(func_deriv_arr, 1./fac)

  return nothing
end


"""
  Method for SolutionDeviation
"""
function calcFunctionalDeriv(mesh::AbstractDGMesh{Tmsh}, 
                           sbp::AbstractOperator,
                           eqn::EulerData{Tsol}, opts,
                           func::SolutionDeviation,
                           func_deriv_arr::Abstract3DArray) where {Tmsh, Tsol}


  val_bar = 1
  u_bar = zeros(Tsol, mesh.numDofPerNode)
  u_bar_bar = zeros(Tsol, mesh.numDofPerNode)
  for i=1:mesh.numEl

    # compute u_bar
    vol = zero(Tmsh)  # area/volume of the element
    fill!(u_bar, 0)
    for j=1:mesh.numNodesPerElement
      w_j = sbp.w[j]/mesh.jac[j, i]
      vol += w_j
      for k=1:mesh.numDofPerNode
        u_bar[k] += eqn.q[k, j, i]*w_j
      end
    end

    for k=1:mesh.numDofPerNode
      u_bar[k] /= vol
    end
#=
    # compute \int || u - u_bar ||
    for j=1:mesh.numNodesPerElement
      w_j = sbp.w[j]/mesh.jac[j, i]
      for k=1:mesh.numDofPerNode
        delta_u = eqn.q[k, j, i] - u_bar[k]
        val += delta_u*w_j*delta_u
      end
    end
=#
    #----------------------------------
    # reverse sweep
    fill!(u_bar_bar, 0)
    for j=1:mesh.numNodesPerElement
      w_j = sbp.w[j]/mesh.jac[j, i]
      for k=1:mesh.numDofPerNode
        delta_u = eqn.q[k, j, i] - u_bar[k]
        #val += delta_u*w_j*delta_u

        delta_u_bar = 2*delta_u*w_j*val_bar
        func_deriv_arr[k, j, i] += delta_u_bar
        u_bar_bar[k]       -= delta_u_bar
      end
    end

    vol_bar = zero(Tmsh)
    for k=1:mesh.numDofPerNode
      u_bar[k] *= vol  # restore primal value
      u_bar_bar[k] /= vol
    end

    for j=1:mesh.numNodesPerElement
      w_j = sbp.w[j]/mesh.jac[j, i]
      for k=1:mesh.numDofPerNode
        func_deriv_arr[k, j, i] += u_bar_bar[k]*w_j
      end
    end

  end  # end i

  return nothing
end




#------------------------------------------------------------------------------
# generic fallback: complex step it

function calcIntegrandDeriv(opts, params::ParamType{2},
          q::AbstractArray{Tsol,1},
          aux_vars::AbstractArray{Tres, 1},
          nrm::AbstractArray{Tmsh},
          integrand_deriv::AbstractArray{Tsol, 1}, node_info,
          func::Tfunc) where {Tsol, Tres, Tmsh, Tfunc<:AbstractBoundaryFunctional}

  h = 1e-20
  pert = Complex128(0, h)
  for i=1:length(q)
    q[i] += pert
    val = calcBoundaryFunctionalIntegrand(params, q, aux_vars, nrm, node_info,
                                          func)
    integrand_deriv[i] = imag(val)/h
    q[i] -= pert
  end

  return nothing
end


