export evalFunctional, calcBndryFunctional, getFunctionalName

@doc """
### EulerEquationMod.evalFunctional

Hight level function that evaluates all the functionals specified over
various edges. At the moment, it can only handle one functional.

**Arguments**

*  `mesh` :  Abstract mesh object
*  `sbp`  : Summation-By-Parts operator
*  `eqn`  : Euler equation object
*  `opts` : Options dictionary

"""->
function evalFunctional{Tmsh, Tsol}(mesh::AbstractMesh{Tmsh},
                        sbp::AbstractSBP, eqn::EulerData{Tsol}, opts,
                        functional::AbstractOptimizationData;
                        functional_number::Int=1)


  eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
  if mesh.isDG
    boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, eqn.q, eqn.q_bndry)
  end

  # Calculate functional over edges
  if functional.is_objective_fn == true
    # The function to be evaluated is an objective function
    functional_edges = opts["geom_faces_objective"]
    functional_name = FunctionalDict[opts["objective_function"]]
    functional.val = calcBndryFunctional(mesh, sbp, eqn, opts, functional,
                     functional_name, functional_edges)
  else
    # Geometric edge at which the functional needs to be integrated
    key = string("geom_edges_functional", functional_number)
    functional_edges = opts[key]
    functional_name = getFunctionalName(opts, functional_number)

    functional.val = calcBndryFunctional(mesh, sbp, eqn, opts, functional,
                     functional_name, functional_edges)

    # Print statements
    if MPI.Comm_rank(eqn.comm) == 0 # If rank is master
      if opts["functional_error"]
        println("\nNumerical functional value on geometric edges ",
                    functional_edges, " = ", functional.val)
        analytical_functional_val = opts["analytical_functional_val"]
        println("analytical_functional_val = ", analytical_functional_val)

        absolute_functional_error = norm((functional.val -
                                         analytical_functional_val), 2)
        relative_functional_error = absolute_functional_error/
                                    norm(analytical_functional_val, 2)

        mesh_metric = 1/sqrt(mesh.numEl/2)  # TODO: Find a suitable mesh metric
        # write functional error to file
        outname = string(opts["functional_error_outfname"], functional_number, ".dat")
        println("printed relative functional error = ",
                relative_functional_error, " to file ", outname, '\n')
        f = open(outname, "w")
        println(f, relative_functional_error, " ", mesh_metric)
        close(f)
      end  # End if opts["functional_error"]
    end    # End @mpi_master
  end # End if is_objective_fn == true

  return nothing
end


@doc """
### EulerEquationMod.calcBndryFunctional

This function calculates a functional on a geometric boundary of a the
computational space. There is no need to call this function withing the
nonlinear solve while computing eqn.q

**Inputs**

*  `mesh` :  Abstract mesh object
*  `sbp`  : Summation-By-Parts operator
*  `eqn`  : Euler equation object
*  `opts` : Options dictionary
*  `g_edge_number` : Geometric edge number

**Outputs**

*  `functional_val` : computed numerical force at the boundary.

"""->

function calcBndryFunctional{Tmsh, Tsol}(mesh::AbstractDGMesh{Tmsh},sbp::AbstractSBP,
                         eqn::EulerData{Tsol}, opts, objective::AbstractOptimizationData,
                         functor::FunctionalType, functional_edges::AbstractArray{Int,1})


  local_functional_val = zero(Tsol)

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
    boundary_integrand = zeros(Tsol, 1, mesh.sbpface.numnodes, nfaces)
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
        nx = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
        ny = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]
        node_info = Int[itr,j,i]
        boundary_integrand[1,j,i] = functor(eqn.params, q2, aux_vars, [nx, ny],
                                            node_info, objective)

      end  # End for j = 1:mesh.sbpface.numnodes
    end    # End for i = 1:nfaces

    val_per_geom_edge = zeros(Tsol, 1)

    integratefunctional!(mesh.sbpface, mesh.bndryfaces[idx_range],
                           boundary_integrand, val_per_geom_edge)

    local_functional_val += val_per_geom_edge[1]

  end # End for itr = 1:length(functional_edges)

  functional_val = zero(Tsol)
  functional_val = MPI.Allreduce(local_functional_val, MPI.SUM, eqn.comm)

  return functional_val
end

@doc """
### EulerEquationMod.drag

Computes the force in the X-direction.

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
              node_info::AbstractArray{Int}, objective::AbstractOptimizationData)

  euler_flux = zeros(Tsol, length(q))
  calcEulerFlux(params, q, aux_vars, nrm, euler_flux)
  val = euler_flux[2]

  return val
end

@doc """
### EulerEquationMod.lift

Computes the force in the Y-direction.

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
              node_info::AbstractArray{Int}, objective::AbstractOptimizationData)

  euler_flux = zeros(Tsol, length(q))
  calcEulerFlux(params, q, aux_vars, nrm, euler_flux)
  val = euler_flux[3]

  return val
end

@doc """
### EulerEquationMod.targetCp

"""

type targetCp <: FunctionalType
end

function call{Tsol, Tres, Tmsh}(obj::targetCp, params, q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, nrm::AbstractArray{Tmsh},
              node_info::AbstractArray{Int}, objective::AbstractOptimizationData)

  cp_node = calcPressureCoeff(params, q)
  g_face = node_info[1]
  node = node_info[2]
  face = node_info[3]
  cp_target = objective.pressCoeff_obj.targetCp_arr[g_face][node, face]

  val = 0.5*((cp_node - cp_target).^2)

  return val
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


#=
function calcForceError()

  #--- Getting nodal errors
  numfacenodes, nfaces, ndimensions = size(boundary_force)
  force_error_mag = zeros(analytical_force)

  for i = 1:nfaces
    for j = 1:numfacenodes
      boundary_force_mag = sqrt(boundary_force[j,i,1]*boundary_force[j,i,1] +
                                boundary_force[j,i,2]*boundary_force[j,i,2])
      force_error_mag[j,i] = analytical_force[j,i] - boundary_force_mag
    end
  end


  #--- Calculate the integral norm of the force error

  # Get a 1D force error array


  M = Array(Tmsh, 1)
  # Get the corresponding mass matrix
  for i = 1:nfaces
    for j = 1:numfacenodes
      for k = 1:mesh.numDofPerNode
        dofnum_k = mesh.dofs[k,j,i]
        M[dofnum_k] += (sbp.w[j]/mesh.jac[j,i])
      end
    end
  end

  return nothing
end

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

=#

#=
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
