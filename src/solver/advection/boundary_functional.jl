# Calculate boundary "forces" in advection
export evalFunctional, calcBndryfunctional, getFunctionalName

@doc """
### EulerEquationMod.evalFunctional

Hight level function that evaluates all the functionals specified over 
various edges 

**Arguments**

*  `mesh` :  Abstract mesh object
*  `sbp`  : Summation-By-Parts operator
*  `eqn`  : Euler equation object
*  `opts` : Options dictionary

"""->
function evalFunctional{Tmsh, Tsol}(mesh::AbstractMesh{Tmsh},
                        sbp::AbstractSBP, eqn::AdvectionData{Tsol}, opts)

  
  eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
  if mesh.isDG
    boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, eqn.q, eqn.q_bndry)
  end

  # Calculate functional over edges
  num_functionals = opts["num_functionals"]
  for j = 1:num_functionals
    # Geometric edge at which the functional needs to be integrated
    key_j = string("geom_edges_functional", j)
    functional_edges = opts[key_j]
    functional_name = getFunctionalName(opts, j)

    functional_val = zero(Tsol)
    functional_val = calcBndryFunctional(mesh, sbp, eqn, opts, 
                     functional_name, functional_edges)

    # Print statements
    MPI.Barrier(eqn.comm)
    if MPI.Comm_rank(eqn.comm) == 0 # If rank is master
      if opts["functional_error"]
        println("\nNumerical functional value on geometric edges $functional_edges = $functional_val")
        
        analytical_functional_val = opts["analytical_functional_val"]
        println("analytical_functional_val = $analytical_functional_val")

        absolute_functional_error = norm((functional_val - 
                                         analytical_functional_val), 2)
        relative_functional_error = absolute_functional_error/
                                    norm(analytical_functional_val, 2)

        mesh_metric = 1/sqrt(mesh.numEl/2)  # TODO: Find a suitable mesh metric
        
        # write functional error to file
        outname = string(opts["functional_error_outfname"], j, ".dat")
        println("printed relative functional error = $relative_functional_error to file $outname\n") 
        f = open(outname, "w")
        println(f, relative_functional_error, " ", mesh_metric)
        close(f)
      end  # End if opts["functional_error"]
    end    # End @mpi_master

  end  # End for i = 1:num_functionals

  return nothing
end


@doc """
AdvectionEquationMod.calcBndryfunctional

This function calculates the forces on a geometric boundary of a the 
computational space. There is no need to call this function withing the 
nonlinear solve while computing eqn.q

**Inputs**

*  `mesh` :  Abstract mesh object
*  `sbp`  : Summation-By-Parts operator
*  `eqn`  : Advection equation object
*  `opts` : Options dictionary
*  `g_edge_number` : Geometric edge number

**Outputs**

*  `functional_val` : computed numerical functional at the boundary.

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

function calcBndryFunctional{Tmsh, Tsol}(mesh::AbstractDGMesh{Tmsh},sbp::AbstractSBP,
                         eqn::AdvectionData{Tsol}, opts, functor, functional_edges)

  # Specify the boundary conditions for the edge on which the force needs to be
  # computed separately. Use that boundary number to access the boundary 
  # offset array. Then proceed the same as bndryflux to get the integrand. Finally
  # use integratefunctional! to get the solution.

  functional_val = zero(Tsol)
  local_functional_val = zero(Tsol)
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
        boundary_integrand[1,j,i] = functor(eqn.params, nx, ny, q) # Boundary Flux
      end
    end

    val_per_geom_edge = zeros(Tsol, 1)
    integratefunctional!(mesh.sbpface, mesh.bndryfaces[idx_range], 
                         boundary_integrand, val_per_geom_edge)

    # Add contributions of multiple geometric edges to the final functional value
    local_functional_val += val_per_geom_edge[1] 

  end   # End for itr = 1:length(functional_edges)

  functional_val = MPI.Allreduce(local_functional_val, MPI.SUM, eqn.comm)
  
  return functional_val
end


#-----    Computing Functional Integrand    -----#

@doc """
### AdvectionEquationMod.qflux

Computes the flux direction and multiplies it with eqn.q_bndry. This is nodal
level operation

"""->

type qflux <: FunctionalType
end

function call(obj::qflux, params::ParamType2, nx, ny, q)
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

#=
@doc """
###AdvectionEquationMod.calcFunctionalIntegrand

Calculates the functional integrand for a particular q.

**Inputs**

*  `alpha_x` : eqn.params.alpha_x
*  `alpha_y` : eqn.params.alpha_y
*  `nx`      : x-component of normal vector
*  `ny`      : y-component of normal vector
*  `q`       : eqn.q or eqn.q_bndry at the particular node

**Outputs**

*  `functional_integrand` : value of the functional integrand at that node 

"""->

function calcFunctionalIntegrand(alpha_x, alpha_y, nx, ny, q)
  
  return functional_integrand = (alpha_x*nx + alpha_y*ny)*q

end
=# 
