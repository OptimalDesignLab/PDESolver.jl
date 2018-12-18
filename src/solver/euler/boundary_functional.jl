# boundary functional evaluation and derivatives

@doc """
### EulerEquationMod.evalFunctional

High level function that evaluates all the functionals specified over
various edges. This function is agnostic to the type of the functional being
computed and calls a mid level functional-type specific function for the actual
evaluation.

**Arguments**

*  `mesh` :  Abstract mesh object
*  `sbp`  : Summation-By-Parts operator
*  `eqn`  : Euler equation object
*  `opts` : Options dictionary
*  `func` : Object of type AbstractBoundaryFunctional. This is type is associated
                      with the functional being computed and holds all the
                      relevant data.

**Outputs**

 * val: functional value
"""->
function _evalFunctional(mesh::AbstractMesh{Tmsh},
            sbp::AbstractOperator, eqn::EulerData{Tsol}, opts,
            func::AbstractBoundaryFunctional) where {Tmsh, Tsol}

  if mesh.isDG
    boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, eqn.q, eqn.q_bndry)
  end

  # Calculate functional over edges
  val = calcBndryFunctional(mesh, sbp, eqn, opts, func)

  return val
end


@doc """
### EulerEquationMod.evalFunctional_revm

Reverse mode of EulerEquationMod.evalFunctional, It takes in functional value
and return `mesh.nrm_bndry_bar`. Different functionals will need to be added
to the if statement to further extend this function.

**Arguments**

*  `mesh` :  Abstract mesh object
*  `sbp`  : Summation-By-Parts operator
*  `eqn`  : Euler equation object
*  `opts` : Options dictionary
*  `func` : Object of type AbstractBoundaryFunctional. This is type is associated
                      with the functional being computed and holds all the
                      relevant data.
*  `functionalName` : Name of the functional being evaluated.

"""->
function _evalFunctionalDeriv_m(mesh::AbstractDGMesh{Tmsh}, 
                           sbp::AbstractOperator,
                           eqn::AbstractSolutionData{Tsol}, opts,
                           func::AbstractBoundaryFunctional,
                           val_bar::Number=1,
                           ) where {Tmsh, Tsol}

  if mesh.isDG
    boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, eqn.q, eqn.q_bndry)
  end

  calcBndryFunctional_revm(mesh, sbp, eqn, opts, func, val_bar)
  return nothing
end


"""
  This is the general case function,  It calls
  [`calcBoundaryFunctionalIntegrand`](@ref) to get the functional integrand at
  each boundary node and computes the integral.
"""
function calcBndryFunctional(mesh::AbstractDGMesh{Tmsh},
                      sbp::AbstractOperator, eqn::EulerData{Tsol, Tres, Tdim},
                      opts, func::AbstractBoundaryFunctional,
                      ) where {Tmsh, Tsol, Tres, Tdim}

  functional_val = zeros(Tres, 1)
  node_info = zeros(Int, 3)

  # loop over boundary conditions that have this functional
  for itr = 1:length(func.bcnums)
    bcnum = func.bcnums[itr]

    start_index = mesh.bndry_offsets[bcnum]
    end_index = mesh.bndry_offsets[bcnum+1]
    idx_range = start_index:(end_index-1)
    bndry_facenums = sview(mesh.bndryfaces, idx_range) # faces on geometric edge i


    nfaces = length(bndry_facenums)
    integrand = zeros(Tsol, 1, mesh.sbpface.numnodes, nfaces)
 
    for i = 1:nfaces
      bndry_i = bndry_facenums[i]
      global_facenum = idx_range[i]
      for j = 1:mesh.sbpface.numnodes
        q = ro_sview(eqn.q_bndry, :, j, global_facenum)
#        convertToConservative(eqn.params, q, q2)
        aux_vars = ro_sview(eqn.aux_vars_bndry, :, j, global_facenum)
        x = ro_sview(mesh.coords_bndry, :, j, global_facenum)
        phys_nrm = ro_sview(mesh.nrm_bndry, :, j, global_facenum)
        node_info[1] = itr; node_info[2] = j; node_info[3] = i

        integrand[1, j, i] = calcBoundaryFunctionalIntegrand(eqn.params, q,
                                           aux_vars, phys_nrm, node_info, func)
      end  # End for j = 1:mesh.sbpface.numnodes
    end    # End for i = 1:nfaces

    bndry_facenums_arr = sview(mesh.bndryfaces, idx_range)
    integratefunctional!(mesh.sbpface, bndry_facenums_arr, integrand, functional_val)
  end  # end loop itr

  # global sum
  val = MPI.allreduce(functional_val, MPI.SUM, eqn.comm)[1]

  return val
end  # function calcBndryFunctional


function calcBndryFunctional_revm(mesh::AbstractDGMesh{Tmsh},
     sbp::AbstractOperator, eqn::EulerData{Tsol, Tres, Tdim},
     opts, func::AbstractBoundaryFunctional, _val_bar::Number=1,
     ) where {Tmsh, Tsol, Tres, Tdim}

  functional_val = zeros(Tres, 1)
  node_info = zeros(Int, 3)

  val_bar = Array{Tres}(1)
  val_bar[1] = _val_bar
#  val_bar = Tres[_val_bar]  # value is implicitly known, no need to reverse allreduce

  # loop over boundary conditions that have this functional
  for itr = 1:length(func.bcnums)
    bcnum = func.bcnums[itr]

    start_index = mesh.bndry_offsets[bcnum]
    end_index = mesh.bndry_offsets[bcnum+1]
    idx_range = start_index:(end_index-1)
    bndry_facenums = sview(mesh.bndryfaces, idx_range) # faces on geometric edge i

    nfaces = length(bndry_facenums)
    integrand_bar = zeros(Tsol, 1, mesh.sbpface.numnodes, nfaces)
 
    bndry_facenums_arr = sview(mesh.bndryfaces, idx_range)
    integratefunctional_rev!(mesh.sbpface, bndry_facenums_arr, integrand_bar, val_bar)

    for i = 1:nfaces
      bndry_i = bndry_facenums[i]
      global_facenum = idx_range[i]
      for j = 1:mesh.sbpface.numnodes
        q = ro_sview(eqn.q_bndry, :, j, global_facenum)
#        convertToConservative(eqn.params, q, q2)
        aux_vars = ro_sview(eqn.aux_vars_bndry, :, j, global_facenum)
        x = ro_sview(mesh.coords_bndry, :, j, global_facenum)
        phys_nrm = ro_sview(mesh.nrm_bndry, :, j, global_facenum)
        nrm_bar = sview(mesh.nrm_bndry_bar, :, j, global_facenum)
        node_info[1] = itr; node_info[2] = j; node_info[3] = i

        calcBoundaryFunctionalIntegrand_revm(eqn.params, q,  aux_vars, phys_nrm,
                              node_info, func, nrm_bar, integrand_bar[1, j, i])
      end  # End for j = 1:mesh.sbpface.numnodes
    end    # End for i = 1:nfaces

  end  # end loop itr

  return
end  # function calcBndryFunctional


function calcBndryFunctional(mesh::AbstractDGMesh{Tmsh},
     sbp::AbstractOperator, eqn::EulerData{Tsol, Tres, Tdim},
     opts, func::LiftCoefficient) where {Tmsh, Tsol, Tres, Tdim}

  val = calcBndryFunctional(mesh, sbp, eqn, opts, func.lift)
  fac = 0.5*eqn.params.rho_free*eqn.params.Ma*eqn.params.Ma

  val = val/fac

  return val
end

#TODO _revm


@doc """
### EulerEquationMod.calcBoundaryFunctionalIntegrand

Computes the integrand for boundary functional at a surface SBP node. Every
functional of a different type may need a corresponding method to compute the
integrand. The type of the functional object, which is a subtype of
`AbstractBoundaryFunctional`.

**Arguments**

*  `params` : eqn.params object
*  `q` : Nodal solution
*  `aux_vars` : Auxiliary variables
*  `nrm` : Face normal vector in the physical space
*  `node_info` : Information about the SBP node
*  `objective` : Functional data type
*  `val` : Function output value

"""->
function calcBoundaryFunctionalIntegrand(params::ParamType{2},
                       q::AbstractArray{Tsol,1},
                       aux_vars::AbstractArray{Tres, 1},
                       nrm::AbstractArray{Tmsh},
                       node_info::AbstractArray{Int},
                       objective::BoundaryForceData,
                       ) where {Tsol, Tres, Tmsh}

  # Compute the numerical flux for the euler equation and extract the X & Y
  # momentum values. The normal vector supplied has already been converted
  # to the physical space from the parametric space.

  euler_flux = objective.euler_flux # Reuse existing memory

  fac = 1.0/(sqrt(nrm[1]*nrm[1] + nrm[2]*nrm[2]))
  # normalize normal vector
  nx = nrm[1]*fac
  ny = nrm[2]*fac

  normal_momentum = nx*q[2] + ny*q[3]

  qg = objective.qg
  for i=1:length(q)
    qg[i] = q[i]
  end
  qg[2] -= nx*normal_momentum
  qg[3] -= ny*normal_momentum

  calcEulerFlux(params, qg, aux_vars, nrm, euler_flux)

  return euler_flux[2]*objective.facx + euler_flux[3]*objective.facy
end # End calcBoundaryFunctionalIntegrand 2D

function calcBoundaryFunctionalIntegrand(params::ParamType{3},
                       q::AbstractArray{Tsol,1},
                       aux_vars::AbstractArray{Tres, 1},
                       nrm::AbstractArray{Tmsh},
                       node_info::AbstractArray{Int},
                       objective::BoundaryForceData,
                       ) where {Tsol, Tres, Tmsh}

  fac = 1.0/(sqrt(nrm[1]*nrm[1] + nrm[2]*nrm[2] + nrm[3]*nrm[3]))
  # normalize normal vector
  nx = nrm[1]*fac
  ny = nrm[2]*fac
  nz = nrm[3]*fac

  normal_momentum = nx*q[2] + ny*q[3] + nz*q[4]
  qg = objective.qg
  for i=1:length(q)
    qg[i] = q[i]
  end
  qg[2] -= nx*normal_momentum
  qg[3] -= ny*normal_momentum
  qg[4] -= nz*normal_momentum

  euler_flux = objective.euler_flux # Reuse existing memory
  calcEulerFlux(params, qg, aux_vars, nrm, euler_flux)

  val = euler_flux[2]*objective.facx + euler_flux[3]*objective.facy +
        euler_flux[4]*objective.facz
  return val
end # End calcBoundaryFunctionalIntegrand 3D

function calcBoundaryFunctionalIntegrand(params::ParamType{2},
                       q::AbstractArray{Tsol,1},
                       aux_vars::AbstractArray{Tres, 1},
                       nrm::AbstractArray{Tmsh},
                       node_info::AbstractArray{Int},
                       objective::MassFlowData,
                       ) where {Tsol, Tres, Tmsh}

  val = q[2]*nrm[1] + q[3]*nrm[2]

  return val
end

function calcBoundaryFunctionalIntegrand(params::ParamType{Tdim},
                       q::AbstractArray{Tsol,1},
                       aux_vars::AbstractArray{Tres, 1},
                       nrm::AbstractArray{Tmsh},
                       node_info::AbstractArray{Int},
                       objective::EntropyFluxData,
                       ) where {Tsol, Tres, Tmsh, Tdim}

  # compute u dot n
  unet = zero(Tres)
  for i=1:Tdim
    unet += q[i+1]*nrm[i]
  end
  unet/q[1]

  p = calcPressure(params, q)
  s = log(p) - params.gamma*log(q[1])
  U = -q[1]*s/params.gamma_1

  return unet*U
end



@doc """
### EulerEquationMod. calcBoundaryFunctionalIntegrand_revm

Reverse mode for boundary functional integrand w.r.t. nrm. Takes in input
val_bar and return nrm_bar for further reverse propagation.

**Arguments**

*  `params` : eqn.params object
*  `q` : Nodal solution
*  `aux_vars` : Auxiliary variables
*  `nrm` : Face normal vector in the physical space
*  `node_info` : Information about the SBP node
*  `objective` : Functional data type
*  `nrm_bar` : Resulting vector
*  `val_bar` : Nodal portion of the seeding vector

"""->

function calcBoundaryFunctionalIntegrand_revm(params::ParamType{2},
                                         q::AbstractArray{Tsol,1},
                                         aux_vars::AbstractArray{Tres, 1},
                                         nrm::AbstractArray{Tmsh},
                                         node_info::AbstractArray{Int},
                                         objective::BoundaryForceData,
                                         nrm_bar::AbstractArray{Tmsh,1},
                                         val_bar::Number) where {Tsol, Tres, Tmsh}

  #---- Forward sweep
  fac = 1.0/(sqrt(nrm[1]*nrm[1] + nrm[2]*nrm[2]))
  nx = nrm[1]*fac # Normalized unit vectors
  ny = nrm[2]*fac #
  normal_momentum = nx*q[2] + ny*q[3]
  qg = objective.qg
  for i=1:length(q)
    qg[i] = q[i]
  end
  qg[2] -= nx*normal_momentum
  qg[3] -= ny*normal_momentum

  #---- Reverse Sweep
  euler_flux_bar = zeros(Tsol, 4) # For 2D
  qg_bar = zeros(Tsol, 4)
  q_bar = zeros(Tsol,4)

  # Reverse diff val = euler_flux[2]*objective.facx + euler_flux[3]*objective.facy
  euler_flux_bar[2] = val_bar*objective.facx
  euler_flux_bar[3] = val_bar*objective.facy

  # Reverse diff calcEulerFlux
  calcEulerFlux_revm(params, qg, aux_vars, nrm, nrm_bar, euler_flux_bar)
  calcEulerFlux_revq(params, qg, qg_bar, aux_vars, nrm, euler_flux_bar)
  ny_bar = zero(Tres)               # Initialize
  nx_bar = zero(Tres)               #
  normal_momentum_bar = zero(Tsol)  #

  # Reverse diff qg[3] -= ny*normal_momentum
  ny_bar -= qg_bar[3]*normal_momentum
  normal_momentum_bar -= qg_bar[3]*ny
  qg_bar[3] += qg_bar[3]

  # Reverse diff qg[2] -= nx*normal_momentum
  nx_bar -= qg_bar[2]*normal_momentum
  normal_momentum_bar -= qg_bar[2]*nx
  qg_bar[2] += qg_bar[2]

  # Reverse diff qg[:] = q[:]
  for i=1:length(q_bar)
    q_bar[i] += qg_bar[i]
  end

  # Reverse diff normal_momentum = nx*q[2] + ny*q[3]
  nx_bar += normal_momentum_bar*q[2]
  ny_bar += normal_momentum_bar*q[3]
  q_bar[2] += normal_momentum_bar*nx
  q_bar[3] += normal_momentum_bar*ny

  # Reverse diff ny = nrm[2]*fac
  fac_bar = zero(Tsol)
  nrm_bar[2] += ny_bar*fac
  fac_bar += ny_bar*nrm[2]

  # Reverse diff nx = nrm[1]*fac
  nrm_bar[1] += nx_bar*fac
  fac_bar += nx_bar*nrm[1]

  # Reverse diff fac = 1.0/(sqrt(nrm[1]*nrm[1] + nrm[2]*nrm[2]))
  nrm_bar[1] -= fac_bar*((nrm[1]*nrm[1] + nrm[2]*nrm[2])^(-1.5))*nrm[1]
  nrm_bar[2] -= fac_bar*((nrm[1]*nrm[1] + nrm[2]*nrm[2])^(-1.5))*nrm[2]

  return nothing
end # End calcBoundaryFunctionalIntegrand_revm 2D

function calcBoundaryFunctionalIntegrand_revm(params::ParamType{3},
                                         q::AbstractArray{Tsol,1},
                                         aux_vars::AbstractArray{Tres, 1},
                                         nrm::AbstractArray{Tmsh},
                                         node_info::AbstractArray{Int},
                                         objective::BoundaryForceData,
                                         nrm_bar::AbstractArray{Tmsh,1},
                                         val_bar::Number) where {Tsol, Tres, Tmsh}

  # Forward Sweep
  fac = 1.0/(sqrt(nrm[1]*nrm[1] + nrm[2]*nrm[2] + nrm[3]*nrm[3]))
  nx = nrm[1]*fac
  ny = nrm[2]*fac
  nz = nrm[3]*fac

  normal_momentum = nx*q[2] + ny*q[3] + nz*q[4]
  qg = objective.qg
  for i=1:length(q)
    qg[i] = q[i]
  end
  qg[2] -= nx*normal_momentum
  qg[3] -= ny*normal_momentum
  qg[4] -= nz*normal_momentum

  # Reverse Sweep
  #TODO: preallocate these
  euler_flux_bar = zeros(Tsol, 5) # For 2D
  qg_bar = zeros(Tsol, 5)
  q_bar = zeros(Tsol,5)

  # Reverse diff val[:] = euler_flux[2:4]
  euler_flux_bar[2] = objective.facx*val_bar
  euler_flux_bar[3] = objective.facy*val_bar
  euler_flux_bar[4] = objective.facz*val_bar

  # Reverse diff calcEulerFlux
  calcEulerFlux_revm(params, qg, aux_vars, nrm, nrm_bar, euler_flux_bar)
  calcEulerFlux_revq(params, qg, qg_bar, aux_vars, nrm, euler_flux_bar)
  nz_bar = zero(Tsol)               #
  ny_bar = zero(Tsol)               # Initialize
  nx_bar = zero(Tsol)               #
  normal_momentum_bar = zero(Tsol)  #

  # qg[4] -= nz*normal_momentum
  nz_bar -= qg_bar[4]*normal_momentum
  normal_momentum_bar -= qg_bar[4]*nz
  qg_bar[4] += qg_bar[4]

  # Reverse diff qg[3] -= ny*normal_momentum
  ny_bar -= qg_bar[3]*normal_momentum
  normal_momentum_bar -= qg_bar[3]*ny
  qg_bar[3] += qg_bar[3]

  # Reverse diff qg[2] -= nx*normal_momentum
  nx_bar -= qg_bar[2]*normal_momentum
  normal_momentum_bar -= qg_bar[2]*nx
  qg_bar[2] += qg_bar[2]

  # Reverse diff qg[:] = q[:]a
  for i=1:length(q_bar)
    q_bar[i] += qg_bar[i]
  end

  # normal_momentum = nx*q[2] + ny*q[3] + nz*q[4]
  nx_bar += normal_momentum_bar*q[2]
  ny_bar += normal_momentum_bar*q[3]
  nz_bar += normal_momentum_bar*q[4]

  # nz = nrm[3]*fac
  nrm_bar[3] += nz_bar*fac
  fac_bar = nz_bar*nrm[3]

  # Reverse diff ny = nrm[2]*fac
  nrm_bar[2] += ny_bar*fac
  fac_bar += ny_bar*nrm[2]

  # Reverse diff nx = nrm[1]*fac
  nrm_bar[1] += nx_bar*fac
  fac_bar += nx_bar*nrm[1]

  # fac = 1.0/(sqrt(nrm[1]*nrm[1] + nrm[2]*nrm[2] + nrm[3]*nrm[3]))
  nrm_bar[1] -= fac_bar*((nrm[1]*nrm[1] + nrm[2]*nrm[2] + nrm[3]*nrm[3])^(-1.5))*nrm[1]
  nrm_bar[2] -= fac_bar*((nrm[1]*nrm[1] + nrm[2]*nrm[2] + nrm[3]*nrm[3])^(-1.5))*nrm[2]
  nrm_bar[3] -= fac_bar*((nrm[1]*nrm[1] + nrm[2]*nrm[2] + nrm[3]*nrm[3])^(-1.5))*nrm[3]

  return nothing
end # End calcBoundaryFunctionalIntegrand_revm 3D
