# Calculate boundary "forces" in advection
export evalFunctional, calcBndryfunctional, getFunctionalName

@doc """
### AdvectionEquationMod.evalFunctional

Hight level function that evaluates functionals specified in the options
dictionary. This function is agnostic which type of a functional is being
computed and calls a mid level type specific function for the actual functional
evaluation.

**Arguments**

*  `mesh` :  Abstract mesh object
*  `sbp`  : Summation-By-Parts operator
*  `eqn`  : Euler equation object
*  `opts` : Options dictionary
*  `functionalData` : Object of the functional being computed.
*  `functional_number` : Optional argument. This needs to be specified for all
                         non-objective functionals being computed, if there are
                         more than 1 of them.

"""->
function evalFunctional{Tmsh, Tsol}(mesh::AbstractMesh{Tmsh},
                        sbp::AbstractSBP, eqn::AdvectionData{Tsol}, opts,
                        functionalData::AbstractOptimizationData;
                        functional_number::Int=1)

  if opts["parallel_type"] == 1

    startSolutionExchange(mesh, sbp, eqn, opts, wait=true)
    @debug1 println(params.f, "-----entered if statement around startDataExchange -----")

  end

  eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
  if mesh.isDG
    boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, eqn.q, eqn.q_bndry)
  end

  calcBndryFunctional(mesh, sbp, eqn, opts, functionalData)

  return nothing
end


@doc """
AdvectionEquationMod.calcBndryfunctional

This function calculates the functional on a geometric boundary of a the
computational space. This is a mid level function that should not be called
from outside the module. Depending on the functional being computed, it
may be necessary to define another method for this function based on a
different boundary functional type.

**Arguments**

*  `mesh` :  Abstract mesh object
*  `sbp`  : Summation-By-Parts operator
*  `eqn`  : Advection equation object
*  `opts` : Options dictionary
*  `functionalData` : Object of the functional being computed

"""->
#=  TODO: uncomment and run test when mesh.bndry_geo_nums gets added to CG meshes

function calcBndryFunctional{Tmsh, Tsol}(mesh::AbstractCGMesh{Tmsh},sbp::AbstractSBP,
                         eqn::AdvectionData{Tsol}, opts, functor, functional_edges)

  # Specify the boundary conditions for the edge on which the force needs to be
  # computed separately. Use that boundary number to access the boundary
  # offset array. Then proceed the same as bndryflux to get the forces using
  # boundaryintegrate!

  functional_val = zero(Tsol)
  alpha_x = eqn.params.alpha_x
  alpha_y = eqn.params.alpha_y

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
    end_index = mesh.bndry_offsets[itr2 + 1]
    idx_range = start_index:(end_index-1)  # Index range
    bndry_facenums = sview(mesh.bndryfaces, idx_range) # faces on geometric edge i

    nfaces = length(bndry_facenums)
    boundary_integrand = zeros(Tsol, 1, sbp.numfacenodes, nfaces)
    boundary_functional = zeros(Tsol, 1, sbp.numnodes, mesh.numEl)

    for i = 1:nfaces
    	bndry_i = bndry_facenums[i]
    	for j = 1:sbp.numfacenodes
        k = sbp.facenodes[j, bndry_i.face]
        q = eqn.q[1,k,bndry_i.element]
        x = sview(mesh.coords, :, k, bndry_i.element)
        dxidx = sview(mesh.dxidx, :, :, k, bndry_i.element)
        nrm = sview(sbp.facenormal, :, bndry_i.face)
        nx = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
        ny = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]
        boundary_integrand[1,j,i] = functor(eqn.params, nx, ny, q) # Boundary Flux
    	end
    end

    boundaryintegrate!(mesh.sbpface, mesh.bndryfaces[idx_range], boundary_integrand, boundary_functional)
    # Add all boundary_force nodal values along the edge to get the nodal force value
    edge_functional_val = zero(Tsol) # functional value over a geometric edge
    for (bindex, bndry) in enumerate(mesh.bndryfaces[idx_range])
      for i = 1:sbp.numfacenodes
        k = sbp.facenodes[i, bndry.face]
        edge_functional_val += boundary_functional[1,k,bndry.element]
      end  # end for i = 1:sbp.numfacenodes
    end    # end enumerate

    functional_val += edge_functional_val
  end      # for itr = 1:length(functional_edges)

  return functional_val
end
=#

function calcBndryFunctional{Tmsh, Tsol, Topt}(mesh::AbstractDGMesh{Tmsh},sbp::AbstractSBP,
                            eqn::AdvectionData{Tsol}, opts,
                            functionalData::QfluxData{Topt})

  # Specify the boundary conditions for the edge on which the force needs to be
  # computed separately. Use that boundary number to access the boundary
  # offset array. Then proceed the same as bndryflux to get the integrand. Finally
  # use integratefunctional! to get the solution.

  functionalData.val = zero(Topt)
  local_functional_val = zero(Tsol)
  alpha_x = eqn.params.alpha_x
  alpha_y = eqn.params.alpha_y

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
    boundary_integrand = zeros(Tsol, 1, mesh.sbpface.numnodes, nfaces)

    for i = 1:nfaces
      bndry_i = bndry_facenums[i]
      global_facenum = idx_range[i]
      for j = 1:mesh.sbpface.numnodes
        q = eqn.q_bndry[ 1, j, global_facenum]
        coords = sview(mesh.coords_bndry, :, j, global_facenum)
        dxidx = sview(mesh.dxidx_bndry, :, :, j, global_facenum)
        nrm = sview(sbp.facenormal, :, bndry_i.face)
        nx = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
        ny = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]
        boundary_integrand[1,j,i] = calcBoundaryFunctionalIntegrand(eqn.params, nx, ny, q,
                                    functionalData) # Boundary Flux
      end
    end

    val_per_geom_edge = zeros(Tsol, 1)
    integratefunctional!(mesh.sbpface, mesh.bndryfaces[idx_range],
                         boundary_integrand, val_per_geom_edge)

    # Add contributions of multiple geometric edges to the final functional value
    local_functional_val += val_per_geom_edge[1]

  end   # End for itr = 1:length(functional_edges)

  functionalData.val = MPI.Allreduce(local_functional_val, MPI.SUM, eqn.comm)

  return nothing
end


#-----    Computing Functional Integrand    -----#
@doc """
### AdvectionEquationMod.calcBoundaryFunctionalIntegrand

Computes the integrand for boundary functional at a surface SBP node. Every
functional needs to have its own method and the functional type determines
which method is called.

**Inputs**

*  `params` : eqn.params object
*  `nx` : X component of face normal vector
*  `ny` : Y component of face normal vector
*  `q`  : Nodal solution variable
*  `functionalData` : Object of the functional being computed

**Outputs**

* `functional_integrand` : Computed integrand at the surface node

"""->

function calcBoundaryFunctionalIntegrand(params::ParamType2, nx, ny, q,
                                         functionalData::QfluxData)

  alpha_x = params.alpha_x
  alpha_y = params.alpha_y

  return functional_integrand = (alpha_x*nx + alpha_y*ny)*q
end


#=
type qflux <: FunctionalType
end

function call(obj::qflux, params::ParamType2, nx, ny, q,
              functionalData::AbstractOptimizationData)

  alpha_x = params.alpha_x
  alpha_y = params.alpha_y
  return functional_integrand = (alpha_x*nx + alpha_y*ny)*q
end

@doc """
### AdvectionEquationMod.FunctionalDict

It stores the names of all possible functional options that can be computed.
Whenever a new functional is created, it should be added to FunctionalDict.

"""->
global const FunctionalDict = Dict{ASCIIString, FunctionalType} (
"qflux" => qflux(),
)

@doc """
### AdvectionEquationMod.getFunctionalName

Gets the name of the functional that needs to be computed at a particular point

**Inputs**

*  `opts`     : Input dictionary
*  `f_number` : Number of the functional in the input dictionary

**Outputs**

*  `functional` : Returns the functional name in the dictionary

"""->
function getFunctionalName(opts, f_number)

  key = string("functional_name", f_number)
  val = opts[key]

  return functional = FunctionalDict[val]
end
=#
