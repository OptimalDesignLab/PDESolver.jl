
import PDESolver.evalFunctional

@doc """
### EulerEquationMod.evalFunctional

Hight level function that evaluates all the functionals specified over
various edges. This function is agnostic to the type of the functional being
computed and calls a mid level functional-type specific function for the actual
evaluation.

**Arguments**

*  `mesh` :  Abstract mesh object
*  `sbp`  : Summation-By-Parts operator
*  `eqn`  : Euler equation object
*  `opts` : Options dictionary
*  `functionalData` : Object of type AbstractFunctional. This is type is associated
                      with the functional being computed and holds all the
                      relevant data.

**Outputs**

 * val: functional value
"""->
function evalFunctional{Tmsh, Tsol}(mesh::AbstractMesh{Tmsh},
                        sbp::AbstractSBP, eqn::EulerData{Tsol}, opts,
                        functionalData::AbstractFunctional)
#=
  if opts["parallel_type"] == 1
    startSolutionExchange(mesh, sbp, eqn, opts, wait=true)
  end
=#
  disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
  if mesh.isDG
    boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, eqn.q, eqn.q_bndry)
  end

  # Calculate functional over edges
  #TODO: have this return the value, get rid of functionalData.val field
  val = calcBndryFunctional(mesh, sbp, eqn, opts, functionalData)

  return val
  # return functionalData.val  # this doesn't work because val does not
                                   # always exist
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
*  `functionalData` : Object of type AbstractFunctional. This is type is associated
                      with the functional being computed and holds all the
                      relevant data.
*  `functionalName` : Name of the functional being evaluated.

"""->

function evalFunctional_revm{Tmsh, Tsol}(mesh::AbstractMesh{Tmsh},
                        sbp::AbstractSBP, eqn::EulerData{Tsol}, opts,
                        functionalData::AbstractFunctional,
                        functionalName::ASCIIString)


  if opts["parallel_type"] == 1
    startSolutionExchange(mesh, sbp, eqn, opts, wait=true)

  end

  disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
  if mesh.isDG
    boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, eqn.q, eqn.q_bndry)
  end

  # Calculate functional over edges
  if functionalName == "lift"

    bndry_force_bar = zeros(Tsol, mesh.dim)
    if mesh.dim == 2
      bndry_force_bar[1] -= sin(eqn.params.aoa)
      bndry_force_bar[2] += cos(eqn.params.aoa)
    else
      bndry_force_bar[1] -= sin(eqn.params.aoa)
      bndry_force_bar[3] += cos(eqn.params.aoa)
    end
    calcBndryFunctional_revm(mesh, sbp, eqn, opts, functionalData, bndry_force_bar)

  elseif functionalName == "drag"

    bndry_force_bar = zeros(Tsol, mesh.dim)
    if mesh.dim == 2
      bndry_force_bar[1] = cos(eqn.params.aoa)
      bndry_force_bar[2] = sin(eqn.params.aoa)
    else
      bndry_force_bar[1] = cos(eqn.params.aoa)
      bndry_force_bar[3] = sin(eqn.params.aoa)
    end
    calcBndryFunctional_revm(mesh, sbp, eqn, opts, functionalData, bndry_force_bar)

  else
    error("reverse mode of functional $functionalName not defined")
  end

  return nothing
end


@doc """
### EulerEquationMod.eval_dJdaoa

Compute the complete derivative of a functional w.r.t angle of attack

**Inputs**

* `mesh` : Abstract mesh object
* `sbp`  : Summation-By-Parts operator
* `eqn`  : Euler equation object
* `opts` : Options dictionary
* `functionalData` : Object of type AbstractFunctional. This is type is associated
                     with the functional being computed and holds all the
                     relevant data.
* `functionalName` : Name of the functional being evaluated
* `adjoint_vec` : Local portion of the adjoint vector owned by an MPI rank

**Output**

* `dJdaoa` : Complete derivative of the functional w.r.t angle of attack
             This is a scalar value that is the same across all MPI ranks

"""->

function eval_dJdaoa{Tmsh, Tsol}(mesh::AbstractMesh{Tmsh}, sbp::AbstractSBP,
                                 eqn::EulerData{Tsol}, opts,
                                 functionalData::BoundaryForceData,
                                 functionalName::ASCIIString,
                                 adjoint_vec::AbstractArray{Tsol,1})

  if functionalName == "lift"
    ∂J∂aoa = functionalData.dLiftdaoa
  elseif functionalName == "drag"
    ∂J∂aoa = functionalData.dDragdaoa
  end

  pert = 1e-20im
  eqn.params.aoa += pert # Imaginary perturbation
  fill!(eqn.res_vec, 0.0)
  fill!(eqn.res, 0.0)
  res_norm = physicsRhs(mesh, sbp, eqn, opts, eqn.res_vec, (evalResidual,))
  ∂R∂aoa = imag(eqn.res_vec)/imag(pert)
  eqn.params.aoa -= pert # Remove perturbation

  # Get the contribution from all MPI ranks
  local_ψT∂R∂aoa = dot(adjoint_vec, ∂R∂aoa)
  ψT∂R∂aoa = MPI.Allreduce(local_ψT∂R∂aoa, MPI.SUM, eqn.comm)
  dJdaoa = ∂J∂aoa + ψT∂R∂aoa

  return dJdaoa
end

function calcBndryFunctional{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractDGMesh{Tmsh},
                             sbp::AbstractSBP, eqn::EulerData{Tsol, Tres, Tdim},
                             opts, functionalData::MassFlowData)

  functional_val = zeros(Tres, 1)

  # loop over boundary conditions that have this functional
  for itr = 1:length(functionalData.bcnums)
    bcnum = functionalData.bcnums[itr]

    start_index = mesh.bndry_offsets[bcnum]
    end_index = mesh.bndry_offsets[bcnum+1]
    idx_range = start_index:(end_index-1)
    bndry_facenums = sview(mesh.bndryfaces, idx_range) # faces on geometric edge i


    nfaces = length(bndry_facenums)
    boundary_integrand = zeros(Tsol, functionalData.ndof, mesh.sbpface.numnodes, nfaces)
 
    for i = 1:nfaces
      bndry_i = bndry_facenums[i]
      global_facenum = idx_range[i]
      for j = 1:mesh.sbpface.numnodes
        q = ro_sview(eqn.q_bndry, :, j, global_facenum)
#        convertToConservative(eqn.params, q, q2)
        aux_vars = ro_sview(eqn.aux_vars_bndry, :, j, global_facenum)
        x = ro_sview(mesh.coords_bndry, :, j, global_facenum)
        phys_nrm = ro_sview(mesh.nrm_bndry, :, j, global_facenum)
        node_info = Int[itr,j,i]
        b_integrand_ji = sview(boundary_integrand,:,j,i)
        calcBoundaryFunctionalIntegrand(eqn.params, q, aux_vars, phys_nrm,
                                        node_info, functionalData, b_integrand_ji)
      end  # End for j = 1:mesh.sbpface.numnodes
    end    # End for i = 1:nfaces

    bndry_facenums_arr = mesh.bndryfaces[idx_range]
    integratefunctional!(mesh.sbpface, bndry_facenums_arr, boundary_integrand, functional_val)
  end  # end loop itr

  # global sum
  functionalData.val = MPI.allreduce(functional_val, MPI.SUM, eqn.comm)[1]

  return functionalData.val
end  # function calcBndryFunctional

@doc """
### EulerEquationMod.calcBndryFunctional

This function calculates a functional on a geometric boundary of a the
computational space. This is a mid level function that should not be called from
outside the module. Depending on the functional being computed, it may be
necessary to define another method for this function based on a different
boundary functional type or parameters.

**Inputs**

*  `mesh` :  Abstract mesh object
*  `sbp`  : Summation-By-Parts operator
*  `eqn`  : Euler equation object
*  `opts` : Options dictionary
*  `functionalData` : Object which is a subtype of Abstract OptimizationData.
                      This is type is associated with the functional being
                      computed and holds all the relevant data.

"""->
function calcBndryFunctional{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractDGMesh{Tmsh},
                             sbp::AbstractSBP, eqn::EulerData{Tsol, Tres, Tdim},
                             opts, functionalData::BoundaryForceData)

  local_functional_val = zeros(Tsol, functionalData.ndof) # Local processor share
  bndry_force = functionalData.bndry_force
  fill!(bndry_force, 0.0)
#  phys_nrm = zeros(Tmsh, Tdim)

  # loop over boundary conditions that have this functional
  for itr = 1:length(functionalData.bcnums)
    bcnum = functionalData.bcnums[itr]

    start_index = mesh.bndry_offsets[bcnum]
    end_index = mesh.bndry_offsets[bcnum+1]
    idx_range = start_index:(end_index-1)
    bndry_facenums = sview(mesh.bndryfaces, idx_range) # faces on geometric edge i

    nfaces = length(bndry_facenums)
    boundary_integrand = zeros(Tsol, functionalData.ndof, mesh.sbpface.numnodes, nfaces)
    q2 = zeros(Tsol, mesh.numDofPerNode)

    for i = 1:nfaces
      bndry_i = bndry_facenums[i]
      global_facenum = idx_range[i]
      for j = 1:mesh.sbpface.numnodes
        q = ro_sview(eqn.q_bndry, :, j, global_facenum)
        convertToConservative(eqn.params, q, q2)
        aux_vars = ro_sview(eqn.aux_vars_bndry, :, j, global_facenum)
        x = ro_sview(mesh.coords_bndry, :, j, global_facenum)
        phys_nrm = ro_sview(mesh.nrm_bndry, :, j, global_facenum)
        node_info = Int[itr,j,i]
        b_integrand_ji = sview(boundary_integrand,:,j,i)
        calcBoundaryFunctionalIntegrand(eqn.params, q2, aux_vars, phys_nrm,
                                        node_info, functionalData, b_integrand_ji)
      end  # End for j = 1:mesh.sbpface.numnodes
    end    # End for i = 1:nfaces

    val_per_geom_edge = zeros(Tsol, functionalData.ndof)

    integratefunctional!(mesh.sbpface, mesh.bndryfaces[idx_range],
                           boundary_integrand, val_per_geom_edge)

    local_functional_val[:] += val_per_geom_edge[:]

  end # End for itr = 1:length(functional_edges)

  for i = 1:functionalData.ndof
    bndry_force[i] = MPI.Allreduce(local_functional_val[i], MPI.SUM, eqn.comm)
  end

  # Compute lift, drag and their corresponding derivatives w.r.t alpha
  aoa = eqn.params.aoa # Angle of attack
  if mesh.dim == 2 # 2D Flow
    functionalData.lift_val = -bndry_force[1]*sin(aoa) + bndry_force[2]*cos(aoa)
    functionalData.drag_val = bndry_force[1]*cos(aoa) + bndry_force[2]*sin(aoa)
    functionalData.dLiftdaoa = -bndry_force[1]*cos(aoa) - bndry_force[2]*sin(aoa)
    functionalData.dDragdaoa = -bndry_force[1]*sin(aoa) + bndry_force[2]*cos(aoa)
  else # 3D Flow
    functionalData.lift_val = -bndry_force[1]*sin(aoa) + bndry_force[3]*cos(aoa)
    functionalData.drag_val = bndry_force[1]*cos(aoa) + bndry_force[3]*sin(aoa)
    functionalData.dLiftdaoa = -bndry_force[1]*cos(aoa) - bndry_force[3]*sin(aoa)
    functionalData.dDragdaoa = -bndry_force[1]*sin(aoa) + bndry_force[3]*cos(aoa)
  end

  if functionalData.isLift
    return functionalData.lift_val
  else
    return functionalData.drag_val
  end
end

@doc """
### EulerEquationMod.calcBndryFunctional_revm

Reverse mode of `calcBndryFunctional` that actually does the work. The commented
lines indicate the line being reverse diffed.

**Arguments**

*  `mesh` :  Abstract mesh object
*  `sbp`  : Summation-By-Parts operator
*  `eqn`  : Euler equation object
*  `opts` : Options dictionary
*  `functionalData` : Object of type AbstractFunctional. This is type is associated
                      with the functional being computed and holds all the
                      relevant data.
*  `bndry_force_bar`: Seed for the reverse mode. This is typically the adjoint
                      vector in our case but could be any vector that is being
                      multiplied to the partial derivative w.r.t mesh metrics
                      being evaluated here.

"""

function calcBndryFunctional_revm{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractDGMesh{Tmsh},
                                       sbp::AbstractSBP, eqn::EulerData{Tsol, Tres, Tdim},
                                       opts, functionalData::BoundaryForceData,
                                       bndry_force_bar::AbstractArray{Tsol, 1})

#  phys_nrm = zeros(Tmsh, Tdim)
  aoa = eqn.params.aoa # Angle of attack

  lift_bar = one(Tsol)
#  nxny_bar = zeros(Tmsh, functionalData.ndof)

  # TODO: Figure out the reverse of MPI.Allreduce. Is it even necessary
  local_functional_val_bar = zeros(Tsol, functionalData.ndof)
  # for i = 1:functionalData.ndof
  #   local_function_val_bar[i] = MPI.bcast(bndry_force_bar[i], 0, eqn.comm)
  # end
  local_functional_val_bar[:] += bndry_force_bar[:]

  # loop over boundary conditions that have this functional
  for itr = 1:length(functionalData.bcnums)
    bcnum = functionalData.bcnums[itr]

    start_index = mesh.bndry_offsets[bcnum]
    end_index = mesh.bndry_offsets[bcnum+1]
    idx_range = start_index:(end_index-1)
    bndry_facenums = sview(mesh.bndryfaces, idx_range) # faces on geometric edge i


    nfaces = length(bndry_facenums)
    boundary_integrand_bar = zeros(Tsol, functionalData.ndof, mesh.sbpface.numnodes, nfaces)
    q2 = zeros(Tsol, mesh.numDofPerNode)

    # local_functional_val[:] += val_per_geom_edge[:]
    val_per_geom_face_bar = zeros(Tsol, functionalData.ndof)
    val_per_geom_face_bar[:] += local_functional_val_bar[:]
    local_functional_val_bar[:] += local_functional_val_bar[:]
    integratefunctional_rev!(mesh.sbpface, mesh.bndryfaces[idx_range],
                             boundary_integrand_bar, val_per_geom_face_bar)


    for i = 1:nfaces
      bndry_i = bndry_facenums[i]
      global_facenum = idx_range[i]
      for j = 1:mesh.sbpface.numnodes
        q = ro_sview(eqn.q_bndry, :, j, global_facenum)
        convertToConservative(eqn.params, q, q2)
        aux_vars = ro_sview(eqn.aux_vars_bndry, :, j, global_facenum)
        x = ro_sview(mesh.coords_bndry, :, j, global_facenum)
        phys_nrm = ro_sview(mesh.nrm_bndry, :, j, global_facenum)
        phys_nrm_bar = sview(mesh.nrm_bndry_bar, :, j, global_facenum)
        node_info = Int[itr,j,i]
        b_integrand_ji_bar = sview(boundary_integrand_bar, :, j, i)
        # calcBoundaryFunctionalIntegrand(eqn.params, q2, aux_vars, phys_nrm,
        #                                node_info, functionalData, b_integrand_ji)
        # fill!(nxny_bar, 0.0)
        calcBoundaryFunctionalIntegrand_revm(eqn.params, q2, aux_vars, phys_nrm,
                                             node_info, functionalData,
                                             phys_nrm_bar, b_integrand_ji_bar)
      end  # End for j = 1:mesh.sbpface.numnodes
    end    # End for i = 1:nfaces

  end # End for itr = 1:length(functional_faces)

  return nothing
end

@doc """
### EulerEquationMod.calcBoundaryFunctionalIntegrand

Computes the integrand for boundary functional at a surface SBP node. Every
functional of a different type may need a corresponding method to compute the
integrand. The type of the functional object, which is a subtype of
`AbstractFunctional`.

**Arguments**

*  `params` : eqn.params object
*  `q` : Nodal solution
*  `aux_vars` : Auxiliary variables
*  `nrm` : Face normal vector in the physical space
*  `node_info` : Information about the SBP node
*  `objective` : Functional data type
*  `val` : Function output value

"""->
function calcBoundaryFunctionalIntegrand{Tsol, Tres, Tmsh}(params::ParamType{2},
                                         q::AbstractArray{Tsol,1},
                                         aux_vars::AbstractArray{Tres, 1},
                                         nrm::AbstractArray{Tmsh},
                                         node_info::AbstractArray{Int},
                                         objective::BoundaryForceData,
                                         val::AbstractArray{Tsol,1})

  # Compute the numerical flux for the euler equation and extract the X & Y
  # momentum values. The normal vector supplied has already been converted
  # to the physical space from the parametric space.

  euler_flux = params.flux_vals1 # Reuse existing memory

  fac = 1.0/(sqrt(nrm[1]*nrm[1] + nrm[2]*nrm[2]))
  # normalize normal vector
  nx = nrm[1]*fac
  ny = nrm[2]*fac

  normal_momentum = nx*q[2] + ny*q[3]

  qg = params.qg
  for i=1:length(q)
    qg[i] = q[i]
  end
  qg[2] -= nx*normal_momentum
  qg[3] -= ny*normal_momentum

  calcEulerFlux(params, qg, aux_vars, nrm, euler_flux)
  val[:] = euler_flux[2:3]      # grab just the momentum terms

  return nothing
end # End calcBoundaryFunctionalIntegrand 2D

@doc """
### EulerEquationMod.calcBoundaryFunctionalIntegrand_diff

Forward mode derivative of calcBoundaryFunctionalIntegrand.
Calculates the derivative of the integrand a boundary functional
wrt the dof's of q. This occurs at a surface SBP node.

**Arguments**

*  `params` : eqn.params object
*  `q` : Nodal solution
*  `aux_vars` : Auxiliary variables
*  `nrm` : Face normal vector in the physical space
*  `node_info` : Information about the SBP node
*  `objective` : Functional data type
*  `val_diff` : Function output value. Must be of size (2,4),
                for 2 components of integrand by four dofs of q that the
                derivative is taken with respect to.

"""->
function calcBoundaryFunctionalIntegrand_diff{Tsol, Tres, Tmsh}(params::ParamType{2},
                                         q::AbstractArray{Tsol,1},
                                         aux_vars::AbstractArray{Tres, 1},
                                         nrm::AbstractArray{Tmsh},
                                         node_info::AbstractArray{Int},
                                         objective::BoundaryForceData,
                                         val_diff::AbstractArray{Tsol,2})

  # forward mode differentiation for dJdu
  # Note: this happens at a node.

  @assert( length(q) == 4 )

  euler_flux = params.flux_vals1 # Reuse existing memory

  fac = 1.0/(sqrt(nrm[1]*nrm[1] + nrm[2]*nrm[2]))
  # normalize normal vector
  nx = nrm[1]*fac
  ny = nrm[2]*fac

  normal_momentum = nx*q[2] + ny*q[3]

  normal_momentum_dot1 = 0.0  # wrt q1
  normal_momentum_dot2 = nx   # wrt q2
  normal_momentum_dot3 = ny   # wrt q3
  normal_momentum_dot4 = 0.0  # wrt q4

  qg = params.qg
  qg_intermediate = params.tmp_qg_intermediate    # now alloc'd in params, was 'qg_temp = zeros(params.qg)'
  # qg_temp_diff = zeros(Tsol, length(qg), length(qg))
  for i=1:length(q)
    # qg_temp[i] = q[i]
    qg_intermediate[i] = q[i]

    # qg_temp_diff[i,i] = 1.0       # qg_temp_diff is actually not used; its values are hardcoded below.
                                    # See the comment above setting qg_diff.
  end

  #=
  qg[1] = qg_temp[1]
  qg[2] = qg_temp[2] - nx*normal_momentum
  qg[3] = qg_temp[3] - ny*normal_momentum
  qg[4] = qg_temp[4]
  =#
  qg[1] = qg_intermediate[1]
  qg[2] = qg_intermediate[2] - nx*normal_momentum
  qg[3] = qg_intermediate[3] - ny*normal_momentum
  qg[4] = qg_intermediate[4]

  qg_diff = params.tmp_qg_diff    # now allocated in params, was qg_diff = zeros(Tsol, length(qg), length(qg))
  fill!(qg_diff, 0.0)

  # This 4x4 array is d(qg)/d(q)
  # 1st dimension is qg element
  # 2nd dimension is the q element that deriv is with respect to
  # Note: this is skipping the use of qg_temp_diff. It's like d(gqt2)/d(q3) = qg_temp_diff[2,3], if you wanted to use it
  qg_diff[1,1] = 1.0      # = d(qg1)/d(qgt1)*d(qgt1)/d(q1) , and d(qgt1)/dq1 = 1.0    (note, qgt1 is qg_temp[1])
  # qg_diff[1,2] = 0.0      # = d(qg1)/d(qgt1)*d(qgt1)/d(q2) , and d(qgt1)/dq2 = 0.0    (note, qgt1 is qg_temp[1])
  # qg_diff[1,3] = 0.0      # = d(qg1)/d(qgt1)*d(qgt1)/d(q3) , and d(qgt1)/dq3 = 0.0    (note, qgt1 is qg_temp[1])
  # qg_diff[1,4] = 0.0      # = d(qg1)/d(qgt1)*d(qgt1)/d(q4) , and d(qgt1)/dq4 = 0.0    (note, qgt1 is qg_temp[1])

  # qg_diff[2,1] = 0.0      # = d(qg2)/d(qgt2)*d(qgt2)/d(q1) , and d(qgt2)/d(q1) = 0.0
  qg_diff[2,2] = 1.0 - normal_momentum_dot2*nx
                          # = d(qg2)/d(qgt2)*d(qgt2)/d(q2) - d(nx*normmom)/d(q2)
                          # = d(qg2)/d(qgt2)*d(qgt2)/d(q2) - [ d(nx)/d(q2)*(normmom) + d(normmom)/d(q2)*(nx) ]
                          #    and d(qg2)/d(qgt2) = 1.0
                          #    and d(qgt2)/d(q2) = 1.0
                          #    and d(nx)/d(q2) = 0.0   (from def)
                          #    and d(normmom)/d(q2) = normal_momentum_dot2
  qg_diff[2,3] = 0.0 - normal_momentum_dot3*nx
                          # = d(qg2)/d(qgt2)*d(qgt2)/d(q3) - d(nx*normmom)/d(q3)
                          # = d(qg2)/d(qgt2)*d(qgt2)/d(q3) - [ d(nx)/d(q3)*(normmom) + d(normmom)/d(q3)*(nx) ]
                          #    and d(qgt2)/d(q3) = 0.0
                          #    and d(nx)/d(q3) = 0.0   (from def)
                          #    and d(normmom)/d(q3) = normal_momentum_dot3
  # qg_diff[2,4] = 0.0      # = d(qg2)/d(qgt2)*d(qgt2)/d(q4) , and d(qgt2)/d(q4) = 0.0

  # qg_diff[3,1] = 0.0      # = d(qg3)/d(qgt3)*d(qgt3)/d(q1) , and d(qgt3)/d(q1) = 0.0
  qg_diff[3,2] = 0.0 - normal_momentum_dot2*ny
                          # = d(qg3)/d(qgt3)*d(qgt3)/d(q2) - d(ny*normmom)/d(q2)
                          # = d(qg3)/d(qgt3)*d(qgt3)/d(q2) - [ d(ny)/d(q2)*(normmom) + d(normmom)/d(q2)*(nx) ]
                          #    and d(qgt3)/d(q2) = 0.0
                          #    and d(ny)/d(q2) = 0.0   (from def)
                          #    and d(normmom)/d(q2) = normal_momentum_dot2
  qg_diff[3,3] = 1.0 - normal_momentum_dot3*ny
                          # = d(qg3)/d(qgt3)*d(qgt3)/d(q3) - d(nx*normmom)/d(q3)
                          # = d(qg3)/d(qgt3)*d(qgt3)/d(q3) - [ d(ny)/d(q3)*(normmom) + d(normmom)/d(q3)*(ny) ]
                          #    and d(qg3)/d(qgt3) = 1.0
                          #    and d(qgt3)/d(q3) = 1.0
                          #    and d(ny)/d(q3) = 0.0   (from def)
                          #    and d(normmom)/d(q3) = normal_momentum_dot3
  # qg_diff[3,4] = 0.0      # = d(qg3)/d(qgt3)*d(qgt3)/d(q4) , and d(qgt3)/d(q4) = 0.0

  # qg_diff[4,1] = 0.0      # = d(qg4)/d(qgt4)*d(qgt4)/d(q1) , and d(qgt4)/dq1 = 0.0
  # qg_diff[4,2] = 0.0      # = d(qg4)/d(qgt4)*d(qgt4)/d(q2) , and d(qgt4)/dq2 = 0.0
  # qg_diff[4,3] = 0.0      # = d(qg4)/d(qgt4)*d(qgt4)/d(q3) , and d(qgt4)/dq3 = 0.0
  qg_diff[4,4] = 1.0      # = d(qg4)/d(qgt4)*d(qgt4)/d(q4) , and d(qgt4)/dq4 = 1.0

  # euler_flux_Jac is the derivative of the euler flux wrt qg. so 4x4
  euler_flux_Jac = params.tmp_euler_flux_Jac    # now allocated in params, was
                                                # 'euler_flux_Jac = zeros(Tres, length(euler_flux), length(euler_flux))'
  fill!(euler_flux_Jac, 0.0)

  calcEulerFlux_diff(params, qg, aux_vars, nrm, euler_flux_Jac)

  # val[:] = euler_flux[2:3]

  # val_diff is passed in. Size: zeros(Tres, 2, 4)

  # EF: Euler Flux
  # now using smallmatvec from ODLCommonTools because it's faster
  # TODO: probably can be further optimized by doing smallmatmat here. need to do some math on it though
  # TODO: slices to sviews when setting val_diff
  # TODO: don't even need dEF_dq* values, they're just intermediate. later
  # dEF_dq1 = euler_flux_Jac*qg_diff[:,1]       # deriv wrt q1, so qg_diff is indexed as row 1
                                              # dEF_dq1 is a 4x1, 4 elements for each of the EF, 1 for d wrt q1
  # dEF_dq1 = smallmatvec(euler_flux_Jac, qg_diff[:,1])
  dEF_dq1 = smallmatvec(euler_flux_Jac, sview(qg_diff, :, 1))

  val_diff[:,1] = dEF_dq1[2:3]                # We want to pick only the momentum terms off the deriv of the EF wrt q1

  # dEF_dq2 = euler_flux_Jac*qg_diff[:,2]       # deriv wrt q2, so qg_diff is indexed as row 2
                                              # dEF_dq2 is a 4x1, 4 elements for each of the EF, 1 for d wrt q2
  # dEF_dq2 = smallmatvec(euler_flux_Jac, qg_diff[:,2])
  dEF_dq2 = smallmatvec(euler_flux_Jac, sview(qg_diff, :, 2))

  val_diff[:,2] = dEF_dq2[2:3]                # We want to pick only the momentum terms off the deriv of the EF wrt q2

  # dEF_dq3 = euler_flux_Jac*qg_diff[:,3]       # deriv wrt q3, so qg_diff is indexed as row 3
                                              # dEF_dq3 is a 4x1, 4 elements for each of the EF, 1 for d wrt q3
  # dEF_dq3 = smallmatvec(euler_flux_Jac, qg_diff[:,3])
  dEF_dq3 = smallmatvec(euler_flux_Jac, sview(qg_diff, :, 3))

  val_diff[:,3] = dEF_dq3[2:3]                # We want to pick only the momentum terms off the deriv of the EF wrt q3

  # dEF_dq4 = euler_flux_Jac*qg_diff[:,4]       # deriv wrt q4, so qg_diff is indexed as row 4
                                              # dEF_dq4 is a 4x1, 4 elements for each of the EF, 1 for d wrt q4
  # dEF_dq4 = smallmatvec(euler_flux_Jac, qg_diff[:,4])
  dEF_dq4 = smallmatvec(euler_flux_Jac, sview(qg_diff, :, 4))

  val_diff[:,4] = dEF_dq4[2:3]                # We want to pick only the momentum terms off the deriv of the EF wrt q4

  # Summary:
  #   val_diff[1,2] contains d/dq2 of val1, the x-component of the boundary functional
  #   val_diff[2,2] contains d/dq2 of val2, the y-component of the boundary functional
  #   val_diff[1,3] contains d/dq3 of val1, the x-component of the boundary functional
  #   val_diff[2,3] contains d/dq3 of val2, the y-component of the boundary functional
  #   same for val_diff[:,1] and val_diff[:,4] for d/dq1 and d/dq4 respectively, but those are zero

  return nothing

end


# 3D version
function calcBoundaryFunctionalIntegrand{Tsol, Tres, Tmsh}(params::ParamType{3},
                                         q::AbstractArray{Tsol,1},
                                         aux_vars::AbstractArray{Tres, 1},
                                         nrm::AbstractArray{Tmsh},
                                         node_info::AbstractArray{Int},
                                         objective::BoundaryForceData,
                                         val::AbstractArray{Tsol,1})

  fac = 1.0/(sqrt(nrm[1]*nrm[1] + nrm[2]*nrm[2] + nrm[3]*nrm[3]))
  # normalize normal vector
  nx = nrm[1]*fac
  ny = nrm[2]*fac
  nz = nrm[3]*fac

  normal_momentum = nx*q[2] + ny*q[3] + nz*q[4]
  qg = params.qg
  for i=1:length(q)
    qg[i] = q[i]
  end
  qg[2] -= nx*normal_momentum
  qg[3] -= ny*normal_momentum
  qg[4] -= nz*normal_momentum

  euler_flux = params.flux_vals1 # Reuse existing memory
  calcEulerFlux(params, qg, aux_vars, nrm, euler_flux)
  val[:] = euler_flux[2:4]

  return nothing
end # End calcBoundaryFunctionalIntegrand 3D

# MassFlow version
function calcBoundaryFunctionalIntegrand{Tsol, Tres, Tmsh}(params::ParamType{2},
                                         q::AbstractArray{Tsol,1},
                                         aux_vars::AbstractArray{Tres, 1},
                                         nrm::AbstractArray{Tmsh},
                                         node_info::AbstractArray{Int},
                                         objective::MassFlowData,
                                         val::AbstractArray{Tsol,1})

  val[1] = q[2]*nrm[1] + q[3]*nrm[2]

  return nothing
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

function calcBoundaryFunctionalIntegrand_revm{Tsol, Tres, Tmsh}(params::ParamType{2},
                                         q::AbstractArray{Tsol,1},
                                         aux_vars::AbstractArray{Tres, 1},
                                         nrm::AbstractArray{Tmsh},
                                         node_info::AbstractArray{Int},
                                         objective::BoundaryForceData,
                                         nrm_bar::AbstractArray{Tmsh,1},
                                         val_bar::AbstractArray{Tres, 1})

  #---- Forward sweep
  fac = 1.0/(sqrt(nrm[1]*nrm[1] + nrm[2]*nrm[2]))
  nx = nrm[1]*fac # Normalized unit vectors
  ny = nrm[2]*fac #
  normal_momentum = nx*q[2] + ny*q[3]
  qg = params.qg
  for i=1:length(q)
    qg[i] = q[i]
  end
  qg[2] -= nx*normal_momentum
  qg[3] -= ny*normal_momentum

  #---- Reverse Sweep
  euler_flux_bar = zeros(Tsol, 4) # For 2D
  qg_bar = zeros(Tsol, 4)
  q_bar = zeros(Tsol,4)

  # Reverse diff val[:] = euler_flux[2:3]
  euler_flux_bar[2:3] += val_bar[:]

  # Reverse diff calcEulerFlux
  calcEulerFlux_revm(params, qg, aux_vars, nrm, euler_flux_bar, nrm_bar)
  calcEulerFlux_revq(params, qg, aux_vars, nrm, euler_flux_bar, qg_bar)
  ny_bar = zero(Tsol)               # Initialize
  nx_bar = zero(Tsol)               #
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
  q_bar[:] += qg_bar[:]

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

function calcBoundaryFunctionalIntegrand_revm{Tsol, Tres, Tmsh}(params::ParamType{3},
                                         q::AbstractArray{Tsol,1},
                                         aux_vars::AbstractArray{Tres, 1},
                                         nrm::AbstractArray{Tmsh},
                                         node_info::AbstractArray{Int},
                                         objective::BoundaryForceData,
                                         nrm_bar::AbstractArray{Tmsh,1},
                                         val_bar::AbstractArray{Tres, 1})

  # Forward Sweep
  fac = 1.0/(sqrt(nrm[1]*nrm[1] + nrm[2]*nrm[2] + nrm[3]*nrm[3]))
  nx = nrm[1]*fac
  ny = nrm[2]*fac
  nz = nrm[3]*fac

  normal_momentum = nx*q[2] + ny*q[3] + nz*q[4]
  qg = params.qg
  for i=1:length(q)
    qg[i] = q[i]
  end
  qg[2] -= nx*normal_momentum
  qg[3] -= ny*normal_momentum
  qg[4] -= nz*normal_momentum

  # Reverse Sweep
  euler_flux_bar = zeros(Tsol, 5) # For 2D
  qg_bar = zeros(Tsol, 5)
  q_bar = zeros(Tsol,5)

  # Reverse diff val[:] = euler_flux[2:4]
  euler_flux_bar[2:4] += val_bar[:]

  # Reverse diff calcEulerFlux
  calcEulerFlux_revm(params, qg, aux_vars, nrm, euler_flux_bar, nrm_bar)
  calcEulerFlux_revq(params, qg, aux_vars, nrm, euler_flux_bar, qg_bar)
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

  # Reverse diff qg[:] = q[:]
  q_bar[:] += qg_bar[:]

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
