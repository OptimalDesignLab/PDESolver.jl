# Calculate boundary "forces" in advection
export calcBndryfunctional, getFunctionalName

@doc """
AdvectionEquationMod.calcBndryforces

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

function calcBndryfunctional{Tmsh, Tsol}(mesh::AbstractCGMesh{Tmsh},sbp::AbstractSBP,
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
    start_index = mesh.bndry_offsets[g_edge_number]
    end_index = mesh.bndry_offsets[g_edge_number + 1]
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


function calcBndryfunctional{Tmsh, Tsol}(mesh::AbstractDGMesh{Tmsh},sbp::AbstractSBP,
                         eqn::AdvectionData{Tsol}, opts, functor, functional_edges)

  # Specify the boundary conditions for the edge on which the force needs to be
  # computed separately. Use that boundary number to access the boundary 
  # offset array. Then proceed the same as bndryflux to get the integrand. Finally
  # use integratefunctional! to get the solution.

  functional_val = zero(Tsol)
  alpha_x = eqn.params.alpha_x
  alpha_y = eqn.params.alpha_y

  for itr = 1:length(functional_edges)
    g_edge_number = functional_edges[itr] # Extract geometric edge number
    start_index = mesh.bndry_offsets[g_edge_number]
    end_index = mesh.bndry_offsets[g_edge_number+1]
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
    functional_val += val_per_geom_edge[1] 

  end   # End for itr = 1:length(functional_edges)
  
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
