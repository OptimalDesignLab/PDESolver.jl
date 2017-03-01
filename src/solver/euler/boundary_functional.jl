export evalFunctional, calcBndryFunctional, getFunctionalName

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

    startDataExchange(mesh, opts, eqn.q, eqn.q_face_send, eqn.q_face_recv,
                      params.f, wait=true)
    @debug1 println(params.f, "-----entered if statement around startDataExchange -----")

  end

  eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
  if mesh.isDG
    boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, eqn.q, eqn.q_bndry)
  end

  # Calculate functional over edges
  calcBndryFunctional(mesh, sbp, eqn, opts, functionalData)

  return nothing
end


@doc """
### EulerEquationMod.calcBndryFunctional

This function calculates a functional on a geometric boundary of a the
computational space. This is a mid level function that should not be called from
outside the module. DEpending on the functional being computd, it may be
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

function calcBndryFunctional{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractDGMesh{Tmsh},sbp::AbstractSBP,
                         eqn::EulerData{Tsol, Tres, Tdim}, opts, functionalData::BoundaryForceData)

  local_functional_val = zeros(Tsol, functionalData.ndof) # Local processor share
  bndry_force = functionalData.bndry_force
  fill!(bndry_force, 0.0)
  functional_edges = functionalData.geom_faces_functional
  phys_nrm = zeros(Tmsh, Tdim)

  # Get bndry_offsets for the functional edge concerned
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
        q = sview(eqn.q_bndry, :, j, global_facenum)
        convertToConservative(eqn.params, q, q2)
        aux_vars = sview(eqn.aux_vars_bndry, :, j, global_facenum)
        x = sview(mesh.coords_bndry, :, j, global_facenum)
        dxidx = sview(mesh.dxidx_bndry, :, :, j, global_facenum)
        nrm = sview(sbp.facenormal, :, bndry_i.face)
        for k = 1:Tdim
            # nx = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
            # ny = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]
            phys_nrm[k] = dxidx[1,k]*nrm[1] + dxidx[2,k]*nrm[2]
          end # End for k = 1:Tdim
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
  functionalData.lift_val = -bndry_force[1]*sin(aoa) + bndry_force[2]*cos(aoa)
  functionalData.drag_val = bndry_force[1]*cos(aoa) + bndry_force[2]*sin(aoa)
  functionalData.dLiftdAlpha = -bndry_force[1]*cos(aoa) - bndry_force[2]*sin(aoa)
  functionalData.dDragdAlpha = -bndry_force[1]*sin(aoa) + bndry_force[2]*cos(aoa)

  return nothing
end

@doc """
###EulerEquationMod.calcBoundaryFunctionalIntegrand

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
function calcBoundaryFunctionalIntegrand{Tsol, Tres, Tmsh}(params,
                                         q::AbstractArray{Tsol,1},
                                         aux_vars::AbstractArray{Tres, 1},
                                         nrm::AbstractArray{Tmsh},
                                         node_info::AbstractArray{Int},
                                         objective::BoundaryForceData,
                                         val::AbstractArray{Tsol,1})

  # Compute the numerical flux for the euler equation and extract the X & Y
  # momentum values

  aoa = params.aoa # Angle of attack
  euler_flux = params.flux_vals1 # Reuse existing memory
  nx = nrm[1]
  ny = nrm[2]

  fac = 1.0/(sqrt(nx*nx + ny*ny))
  # normalize normal vector
  nx *= fac
  ny *= fac

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
end

@doc """
calcBoundaryFunctionalIntegrand_revm

Reverse mode for boundary functional integrand w.r.t. nrm. Takes in input
val_bar and return nrm_bar for further reverse propagation.

"""->

function calcBoundaryFunctionalIntegrand_revm{Tsol, Tres, Tmsh}(params,
                                         q::AbstractArray{Tsol,1},
                                         aux_vars::AbstractArray{Tres, 1},
                                         nrm::AbstractArray{Tmsh},
                                         node_info::AbstractArray{Int},
                                         objective::BoundaryForceData,
                                         nrm_bar, val_bar)

  aoa = params.aoa # Angle of attack
  nx = nrm[1]
  ny = nrm[2]
  fac = 1.0/(sqrt(nx*nx + ny*ny))
  nx *= fac # Normalize nx & ny
  ny *= fac
  normal_momentum = nx*q[2] + ny*q[3]

  euler_flux_bar = zeros(Tsol, 4) # For 2D
  qg_bar = zeros(Tsol, 4)
  q_bar = zeros(Tsol,4)
  euler_flux_bar[2:3] += val_bar[:]

  calcEulerFlux_revm(params, q, aux_vars, nrm, euler_flux_bar, qg_bar, nrm_bar)
  normal_momentum_bar = zero(Tsol)

  # Reverse diff qg[3] -= ny*normal_momentum
  ny_bar = zero(Tsol)
  ny_bar -= qg_bar[3]*normal_momentum
  normal_momentum_bar -= qg_bar[3]*ny
  # Reverse diff qg[2] -= nx*normal_momentum
  nx_bar = zero(Tsol)
  nx_bar -= qg_bar[2]*normal_momentum
  normal_momentum_bar -= qg_bar[2]*nx

  q_bar[:] += qg_bar[:]

  # Reverse diff normal_momentum = nx*q[2] + ny*q[3]
  q_bar[2] += normal_momentum_bar*nx
  nx_bar += normal_momentum_bar*q[2]
  q_bar[3] += normal_momentum_bar*ny
  ny_bar += normal_momentum_bar*q[3]

  # Reverse diff ny *= fac
  fac_bar = zero(Tsol)
  fac_bar += ny_bar*ny
  ny_bar += ny_bar*fac

  # Reverse diff nx *= fac
  fac_bar += nx_bar*nx
  nx_bar += nx_bar*fac

  # Reverse diff fac = 1.0/(sqrt(nx*nx + ny*ny))
  nx_bar += -fac_bar*((nx*nx + ny*ny)^(-1.5))*nx
  ny_bar += -fac_bar*((nx*nx + ny*ny)^(-1.5))*ny

  nrm_bar[1] += nx_bar
  nrm_bar[2] += ny_bar

  return nothing
end

#=
@doc """
### EulerEquationMod.drag

Computes the integrand at a node for computing the drag force on a boundary
surface/edge. Note that the drag is always tangential to the free-stream
velocity.

**Inputs**

*  `params` : Parameter type
*  `q`      : Solution at a node
*  `aux_vars` : Vector of auxiliary variables
*  `nrm`    : Normal vector in the physical space
*  `node_info` : 1D, 3 element array containing information about the node.
                 node_info[1] = geometric edge number
                 node_info[2] = sbpface node number
                 node_info[3] = element face number on the geometric edge

**Outputs**

*  `val`    : Momentum derivative in the X-direction

"""->

type drag <: FunctionalType
end

function call{Tsol, Tres, Tmsh}(obj::drag, params, q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, nrm::AbstractArray{Tmsh},
              node_info::AbstractArray{Int},
              objective::AbstractOptimizationData, val::AbstractArray{Tsol,1})

  aoa = params.aoa # Angle of attack
  euler_flux = params.flux_vals1 # Reuse existing memory
  nx = nrm[1]
  ny = nrm[2]

  fac = 1.0/(sqrt(nx*nx + ny*ny))
  # normalize normal vector
  nx *= fac
  ny *= fac

  normal_momentum = nx*q[2] + ny*q[3]

  qg = params.qg
  for i=1:length(q)
    qg[i] = q[i]
  end
  qg[2] -= nx*normal_momentum
  qg[3] -= ny*normal_momentum

  calcEulerFlux(params, qg, aux_vars, nrm, euler_flux)

  val[1] = euler_flux[2]*cos(aoa) + euler_flux[3]*sin(aoa)

  return nothing
end

@doc """
###EulerEquationMod.dDragdALpha

Compute the integrand for computing \frac{\partial Drag}/{\partial Alpha}

"""->
type dDragdAlpha <: FunctionalType
end

function call{Tsol, Tres, Tmsh}(obj::dDragdAlpha, params, q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, nrm::AbstractArray{Tmsh},
              node_info::AbstractArray{Int},
              objective::AbstractOptimizationData, val::AbstractArray{Tsol,1})

  aoa = params.aoa # Angle of attack
  euler_flux = params.flux_vals1 # Reuse existing memory
  nx = nrm[1]
  ny = nrm[2]

  fac = 1.0/(sqrt(nx*nx + ny*ny))
  # normalize normal vector
  nx *= fac
  ny *= fac

  normal_momentum = nx*q[2] + ny*q[3]

  qg = params.qg
  for i=1:length(q)
    qg[i] = q[i]
  end
  qg[2] -= nx*normal_momentum
  qg[3] -= ny*normal_momentum
  calcEulerFlux(params, qg, aux_vars, nrm, euler_flux)

  # calcEulerFlux(params, q, aux_vars, nrm, euler_flux)

  val[1] = -euler_flux[2]*sin(aoa) + euler_flux[3]*cos(aoa)

  return nothing
end

@doc """
### EulerEquationMod.lift

Computes the lift force. Note lift is always perpendicular to the free-stream
velocity.

**Inputs**

*  `params` : Parameter type
*  `q`      : Solution at a node
*  `aux_vars` : Vector of auxiliary variables
*  `nrm`    : Normal vector in the physical space

**Outputs**

*  `val`    : Momentum derivative in the Y-direction

"""->

type lift <: FunctionalType
end

function call{Tsol, Tres, Tmsh}(obj::lift, params, q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, nrm::AbstractArray{Tmsh},
              node_info::AbstractArray{Int},
              objective::AbstractOptimizationData, val::AbstractArray{Tsol,1})

  aoa = params.aoa # Angle of attack
  euler_flux = params.flux_vals1 # Reuse existing memory
  nx = nrm[1]
  ny = nrm[2]

  fac = 1.0/(sqrt(nx*nx + ny*ny))
  # normalize normal vector
  nx *= fac
  ny *= fac

  normal_momentum = nx*q[2] + ny*q[3]

  qg = params.qg
  for i=1:length(q)
    qg[i] = q[i]
  end
  qg[2] -= nx*normal_momentum
  qg[3] -= ny*normal_momentum
  calcEulerFlux(params, qg, aux_vars, nrm, euler_flux)
  # calcEulerFlux(params, q, aux_vars, nrm, euler_flux)
  val[1] = -euler_flux[2]*cos(aoa) + euler_flux[3]*sin(aoa)

  return nothing
end

@doc """
### EulerEquationMod.dLiftdALpha

Compute the integrand for computing \frac{\partial lift}/{\partial alpha}

"""->
type dLiftdAlpha <: FunctionalType
end

function call{Tsol, Tres, Tmsh}(obj::dLiftdAlpha, params, q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, nrm::AbstractArray{Tmsh},
              node_info::AbstractArray{Int},
              objective::AbstractOptimizationData, val::AbstractArray{Tsol,1})

  aoa = params.aoa # Angle of attack
  euler_flux = params.flux_vals1 # Reuse existing memory
  nx = nrm[1]
  ny = nrm[2]

  fac = 1.0/(sqrt(nx*nx + ny*ny))
  # normalize normal vector
  nx *= fac
  ny *= fac

  normal_momentum = nx*q[2] + ny*q[3]

  qg = params.qg
  for i=1:length(q)
    qg[i] = q[i]
  end
  qg[2] -= nx*normal_momentum
  qg[3] -= ny*normal_momentum
  calcEulerFlux(params, qg, aux_vars, nrm, euler_flux)
  val[1] = euler_flux[2]*sin(aoa) + euler_flux[3]*cos(aoa)

  return nothing
end

@doc """
### EulerEquationMod.targetCp

"""

type targetCp <: FunctionalType
end

function call{Tsol, Tres, Tmsh}(obj::targetCp, params, q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, nrm::AbstractArray{Tmsh},
              node_info::AbstractArray{Int},
              objective::AbstractOptimizationData, val::AbstractArray{Tsol,1})

  cp_node = calcPressureCoeff(params, q)
  g_face = node_info[1]
  node = node_info[2]
  face = node_info[3]
  cp_target = objective.pressCoeff_obj.targetCp_arr[g_face][node, face]

  val[1] = 0.5*((cp_node - cp_target).^2)

  return nothing
end


@doc """
### EulerEquationMod.FunctionalDict

It stores the names of all possible functional options that can be computed.
Whenever a new functional is created, it should be added to FunctionalDict.

"""->
global const FunctionalDict = Dict{ASCIIString, FunctionalType} (
"drag" => drag(),
"lift" => lift(),
"targetCp" => targetCp(),
"dLiftdAlpha" => dLiftdAlpha(),
"dDragdAlpha" => dDragdAlpha(),
"boundaryForce" => boundaryForce()
)


@doc """
### EulerEquationMod.getFunctionalName

Gets the name of the functional that needs to be computed at a particular point

**Inputs**

*  `opts`     : Input dictionary
*  `f_number` : Number of the functional in the input dictionary

**Outputs**

*  `functional` : Returns the functional name in the dictionary. It is of type
                  `FucntionalType`,

"""->
function getFunctionalName(opts, f_number;is_objective_fn=false)

  key = string("functional_name", f_number)
  val = opts[key]

  return functional = FunctionalDict[val]
end

function getnFaces(mesh::AbstractDGMesh, g_face::Int)

  i = 0
  for i = 1:mesh.numBC
    if findfirst(mesh.bndry_geo_nums[i],g_face) > 0
      break
    end
  end

  start_index = mesh.bndry_offsets[i]
  end_index = mesh.bndry_offsets[i+1]
  idx_range = start_index:(end_index-1)
  bndry_facenums = sview(mesh.bndryfaces, idx_range) # faces on geometric edge i
  nfaces = length(bndry_facenums)

  return nfaces
end
=#

#=

function calcPhysicalEulerFlux{Tsol}(params::ParamType{2}, q::AbstractArray{Tsol,1},
                               F::AbstractArray{Tsol, 2})

  u = q[2]/q[1]
  v = q[3]/q[1]
  p = calcPressure(params, q)

  # Calculate Euler Flux in X-direction
  F[1,1] = q[2]
  F[2,1] = q[2]*u + p
  F[3,1] = q[2]*v
  F[4,1] = u*(q[4] + p)

  # Calculate Euler Flux in Y-direction

  F[1,2] = q[3]
  F[2,2] = q[3]*u
  F[3,2] = q[3]*v + p
  F[4,2] = v*(q[4] + p)

  return nothing
end

function calcBndryfunctional{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractCGMesh{Tmsh},
                            sbp::AbstractSBP, eqn::EulerData{Tsol, Tres, Tdim},
                            opts, g_edge_number)

  # Specify the boundary conditions for the edge on which the force needs to be computed
  # separately in the input dictionary. Use that boundary number to access the boundary
  # offset array. Then proceed the same as bndryflux to get the forces using
  # boundaryintegrate!


  # g_edge_number = 1 # Geometric boundary edge on which the force needs to be computed
  start_index = mesh.bndry_offsets[g_edge_number]
  end_index = mesh.bndry_offsets[g_edge_number+1]
  bndry_facenums = sview(mesh.bndryfaces, start_index:(end_index - 1)) # faces on geometric edge i
  # println("bndry_facenums = ", bndry_facenums)

  nfaces = length(bndry_facenums)
  boundary_press = zeros(Tsol, Tdim, sbp.numfacenodes, nfaces)
  boundary_force = zeros(Tsol, Tdim, sbp.numnodes, mesh.numEl)
  q2 = zeros(Tsol, mesh.numDofPerNode)
  # analytical_force = zeros(Tsol, sbp.numfacenodes, nfaces)


  for i = 1:nfaces
    bndry_i = bndry_facenums[i]
    for j = 1:sbp.numfacenodes
      k = sbp.facenodes[j, bndry_i.face]
      q = sview(eqn.q, :, k, bndry_i.element)
      convertToConservative(eqn.params, q, q2)
      aux_vars = sview(eqn.aux_vars, :, k, bndry_i.element)
      x = sview(mesh.coords, :, k, bndry_i.element)
      dxidx = sview(mesh.dxidx, :, :, k, bndry_i.element)
      nrm = sview(sbp.facenormal, :, bndry_i.face)

      # analytical_force[k,bndry_i.element] = calc_analytical_forces(mesh, eqn.params, x)
      nx = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
      ny = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]

      # Calculate euler flux for the current iteration
      euler_flux = zeros(Tsol, mesh.numDofPerNode)
      calcEulerFlux(eqn.params, q2, aux_vars, [nx, ny], euler_flux)

      # Boundary pressure in "ndimensions" direcion
      boundary_press[:,j,i] =  euler_flux[2:3]
    end # end for j = 1:sbp.numfacenodes
  end   # end for i = 1:nfaces
  boundaryintegrate!(mesh.sbpface, mesh.bndryfaces[start_index:(end_index - 1)],
                     boundary_press, boundary_force)

  functional_val = zeros(Tsol,2)

  for (bindex, bndry) in enumerate(mesh.bndryfaces[start_index:(end_index - 1)])
    for i = 1:sbp.numfacenodes
      k = sbp.facenodes[i, bndry.face]
      functional_val[:] += boundary_force[:,k,bndry.element]
    end
  end  # end enumerate


  return functional_val
end
=#
