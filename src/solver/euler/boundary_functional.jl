export evalFunctional, calcBndryFunctional, eval_dJdaoa

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
*  `functionalData` : Object of type AbstractOptimizationData. This is type is associated
                      with the functional being computed and holds all the
                      relevant data.
*  `functional_number` : A number identifying which functional is being computed.
                         This is important when multiple functions, that aren't
                         objective functions are being evaluated. Default value
                         is 1.
"""->
function evalFunctional{Tmsh, Tsol}(mesh::AbstractMesh{Tmsh},
                        sbp::AbstractSBP, eqn::EulerData{Tsol}, opts,
                        functionalData::AbstractOptimizationData;
                        functional_number::Int=1)

  if opts["parallel_type"] == 1
    startSolutionExchange(mesh, sbp, eqn, opts, wait=true)
  end

  disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
  if mesh.isDG
    boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, eqn.q, eqn.q_bndry)
  end

  # Calculate functional over edges
  calcBndryFunctional(mesh, sbp, eqn, opts, functionalData)

  #TODO: return functional value as well
  return nothing
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
*  `functionalData` : Object of type AbstractOptimizationData. This is type is associated
                      with the functional being computed and holds all the
                      relevant data.
*  `functionalName` : Name of the functional being evaluated.

"""->

function evalFunctional_revm{Tmsh, Tsol}(mesh::AbstractMesh{Tmsh},
                        sbp::AbstractSBP, eqn::EulerData{Tsol}, opts,
                        functionalData::AbstractOptimizationData,
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
* `functionalData` : Object of type AbstractOptimizationData. This is type is associated
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
                                 functionalData::AbstractOptimizationData,
                                 functionalName::ASCIIString,
                                 adjoint_vec::AbstractArray{Tsol,1})

  if functionalName == "lift"
    ∂J∂aoa = functionalData.dLiftdaoa
  elseif functionalName == "drag"
    ∂J∂aoa = functionalData.dDragdaoa
  end

  println("norm(adjoint_vec) = ", norm(adjoint_vec))
  println("partial J partial alpha = ", ∂J∂aoa)


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

  println("norm(partial R partial alpha) = ", norm(∂R∂aoa))
  println("adjoint contribution = ", ψT∂R∂aoa)


  dJdaoa = ∂J∂aoa + ψT∂R∂aoa

  return dJdaoa
end

function calcBndryFunctional{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractDGMesh{Tmsh},
                             sbp::AbstractSBP, eqn::EulerData{Tsol, Tres, Tdim},
                             opts, functionalData::MassFlowData)

  functional_val = zeros(Tres, 1)
  # Get bndry_offsets for the functional edge concerned
  #!!!: no, this gives wrong results if there is more than one geometric edge in
  #     the group containing the specified edge

  functional_edges = functionalData.geom_faces_functional
  for itr = 1:length(functional_edges)
    g_edge_number = functional_edges[itr] # Extract geometric edge number
    # get the boundary array associated with the geometric edge
    itr2 = 0
    for itr2 = 1:mesh.numBC
      if findfirst(mesh.bndry_geo_nums[itr2],g_edge_number) > 0
        break
      end
    end

    start_index = mesh.bndry_offsets[itr2]
    end_index = mesh.bndry_offsets[itr2+1]
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

  return nothing
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
  functional_edges = functionalData.geom_faces_functional
#  phys_nrm = zeros(Tmsh, Tdim)

  # Get bndry_offsets for the functional edge concerned
  #!!!: no, this gives wrong results if there is more than one geometric edge in
  #     the group containing the specified edge
  for itr = 1:length(functional_edges)
    g_edge_number = functional_edges[itr] # Extract geometric edge number
    # get the boundary array associated with the geometric edge
    itr2 = 0
    for itr2 = 1:mesh.numBC
      if findfirst(mesh.bndry_geo_nums[itr2],g_edge_number) > 0
        break
      end
    end

    start_index = mesh.bndry_offsets[itr2]
    end_index = mesh.bndry_offsets[itr2+1]
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

  return nothing
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
*  `functionalData` : Object of type AbstractOptimizationData. This is type is associated
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

  functional_faces = functionalData.geom_faces_functional
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

  # Loop over geometrical functional faces
  for itr = 1:length(functional_faces)

    g_face_number = functional_faces[itr] # Extract geometric edge number
    # get the boundary array associated with the geometric edge
    itr2 = 0
    for itr2 = 1:mesh.numBC
      if findfirst(mesh.bndry_geo_nums[itr2],g_face_number) > 0
        break
      end
    end

    start_index = mesh.bndry_offsets[itr2]
    end_index = mesh.bndry_offsets[itr2+1]
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
`AbstractOptimizationData`.

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
  val[:] = euler_flux[2:3]

  return nothing
end # End calcBoundaryFunctionalIntegrand 2D

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
