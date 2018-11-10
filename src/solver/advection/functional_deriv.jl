# functional derivative computation

import PDESolver.evalFunctionalDeriv_q

"""
  evalFunctionalDeriv_q for Advection
"""
function evalFunctionalDeriv_q(mesh::AbstractDGMesh{Tmsh},
                           sbp::AbstractSBP,
                           eqn::AdvectionData{Tsol}, opts,
                           functionalData::AbstractIntegralFunctional,
                           func_deriv_arr::Abstract3DArray) where {Tmsh, Tsol}

  array1DTo3D(mesh, sbp, eqn, opts, eqn.q_vec, eqn.q)
  if mesh.isDG
    boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, eqn.q, eqn.q_bndry)
  end

  setupFunctional(mesh, sbp, eqn, opts, functionalData)

  calcFunctionalDeriv(mesh, sbp, eqn, opts, functionalData, func_deriv_arr)

  return nothing
end


@doc """
### AdvectionEquationMod.calcFunctionalDeriv

Computes a 3D array of the derivative of a functional w.r.t eqn.q on all
mesh nodes.

**Inputs**

*  `mesh`  : Abstract mesh object
*  `sbp`   : Summation-By-Parts operator
*  `eqn`   : Advection equation object
*  `opts`  : Options dictionary
*  `functionalData` : Object of subtype of AbstractFunctional. This is
                      the type associated with the adjoint of the functional
                      being computed and holds all the necessary data.
*  `func_deriv_arr` : 3D array that stors the derivative of functional w.r.t
                      eqn.q. It has a structure [1, numnodes_per_element, numEl]

**Outputs**

*  None

"""->
#=
function calcFunctionalDeriv{Tmsh, Tsol}(mesh::AbstractCGMesh{Tmsh}, sbp::AbstractSBP,
                             eqn ::AdvectionData{Tsol}, opts, functor, functional_edges,
                             functionalData, func_deriv_arr)

  alpha_x = eqn.params.alpha_x
  alpha_y = eqn.params.alpha_y

  # Obtain the derivative of the integrand at all meh.bndry
  integrand = zeros(eqn.bndryflux)
  # println("size of integrand = ", size(integrand))

  for itr = 1:length(functional_edges)
    g_edge_number = functional_edges[itr] # Extract geometric edge number
    # get the boundary array associated with the geometric edge
    itr2 = 0
    for itr2 = 1:mesh.numBC
      if findfirst(mesh.bndry_geo_nums[itr2],g_edge_number) > 0
        break
      end
    end

    start_index = mesh.bndry_offsets[g_edge_number]
    end_index = mesh.bndry_offsets[g_edge_number + 1]
    idx_range = start_index:(end_index-1)  # Index range
    bndry_facenums = sview(mesh.bndryfaces, idx_range) # faces on geometric edge i

    nfaces = length(bndry_facenums)

    for i = 1:nfaces
      bndry_i = bndry_facenums[i]
      global_facenum = idx_range[i]
      for j = 1:mesh.numNodesPerFace
        k = mesh.facenodes[j, bndry_i.face]
        q = eqn.q[1,k,bndry_i.element]
        x = ro_sview(mesh.coords, :, k, bndry_i.element)
        dxidx = ro_sview(mesh.dxidx, :, :, k, bndry_i.element)
        nrm = ro_sview(mesh.sbpface.normal, :, bndry_i.face)
        nx = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
        ny = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]
        integrand[1,j,global_facenum] = calcIntegrandDeriv(opts, functor, eqn.params,
                                        nx, ny, q, functionalData)
      end  # End for j = 1:mesh.numNodesPerFace
    end    # End for i = 1:nfaces
  end      # End for itr = 1:length(functional_edges)


  # println("mesh.bndryfaces = \n", mesh.bndryfaces)

  boundaryintegrate!(mesh.sbpface, mesh.bndryfaces, integrand, func_deriv_arr)

  return nothing
end
=#
# DG Version
function calcFunctionalDeriv(mesh::AbstractDGMesh{Tmsh}, sbp::AbstractSBP,
                 eqn::AdvectionData{Tsol}, opts,
                 functionalData::AbstractIntegralFunctional,
                 func_deriv_arr) where {Tmsh, Tsol}

  alpha_x = eqn.params.alpha_x
  alpha_y = eqn.params.alpha_y

  # Obtain the derivative of the integrand at all mesh.bndry
  integrand = zeros(eqn.q_bndry)
  # Populate integrand
  for itr = 1:length(functionalData.bcnums)
    bcnum = functionalData.bcnums[itr]

    start_index = mesh.bndry_offsets[bcnum]
    end_index = mesh.bndry_offsets[bcnum+1]
    idx_range = start_index:(end_index-1)
    bndry_facenums = sview(mesh.bndryfaces, idx_range) # faces on geometric edge i
    nfaces = length(bndry_facenums)

    for i = 1:nfaces
      bndry_i = bndry_facenums[i]
      global_facenum = idx_range[i]
      for j = 1:mesh.sbpface.numnodes
        q = eqn.q_bndry[ 1, j, global_facenum]
        coords = ro_sview(mesh.coords_bndry, :, j, global_facenum)
        nrm = ro_sview(mesh.nrm_bndry, :, j, global_facenum)
        nx = nrm[1]
        ny = nrm[2]
        integrand[1,j,global_facenum] = calcIntegrandDeriv(opts, eqn.params,
                                        nx, ny, q, functionalData)
      end  # End for j = 1:mesh.sbpface.numnodes
    end    # End for i = 1:nfaces
  end

  boundaryintegrate!(mesh.sbpface, mesh.bndryfaces, integrand, func_deriv_arr)

  return nothing
end    # End function calcFunctionalDeriv


@doc """
### AdvectionEquationMod.calcIntegrandDeriv

Compute the derivative of the integrand at a point. It presently uses complex
step to compute the derivative

**Inputs**

*  `opts`    : Input dictionary
*  `params`  : the ParamType for the equation
*  `nx` & `ny` : Normal vectors
*  `q`       : Solution variable
*  `functionalData` : Functional object

**Outputs**

*  `integrand_deriv` : derivative of the functor w.r.t q

"""->

function calcIntegrandDeriv(opts, params::ParamType2, nx, ny, q,
                            functionalData::AbstractIntegralFunctional)

  pert = complex(0, 1e-20)  # complex perturbation
  q += pert
  val = calcBoundaryFunctionalIntegrand(params, nx, ny, q, functionalData)
  integrand_deriv = imag(val)/norm(pert)

  return integrand_deriv
end
