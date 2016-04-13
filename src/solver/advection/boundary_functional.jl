# Calculate boundary "forces" in advection
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
                         eqn::AdvectionData{Tsol}, opts, functional_edges)

  # Specify the boundary conditions for the edge on which the force needs to be
  # computed separately. Use that boundary number to access the boundary 
  # offset array. Then proceed the same as bndryflux to get the forces using 
  # boundaryintegrate!

  functional_val = zero(Tsol)
  for itr = 1:length(functional_edges)
    g_edge_number = functional_edges[itr] # Extract geometric edge number
    start_index = mesh.bndry_offsets[g_edge_number]
    end_index = mesh.bndry_offsets[g_edge_number+1]
    idx_range = start_index:(end_index-1)  # Index range
    bndry_facenums = view(mesh.bndryfaces, idx_range) # faces on geometric edge i

    nfaces = length(bndry_facenums)
    boundary_integrand = zeros(Tsol, 1, sbp.numfacenodes, nfaces)
    boundary_functional = zeros(Tsol, 1, sbp.numnodes, mesh.numEl)
    alpha_x = eqn.alpha_x
    alpha_y = eqn.alpha_y

    for i = 1:nfaces
    	bndry_i = bndry_facenums[i]
    	for j = 1:sbp.numfacenodes
        k = sbp.facenodes[j, bndry_i.face]
        q = eqn.q[1,k,bndry_i.element]
        x = view(mesh.coords, :, k, bndry_i.element)
        dxidx = view(mesh.dxidx, :, :, k, bndry_i.element)
        nrm = view(sbp.facenormal, :, bndry_i.face)
        nx = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
        ny = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]
        boundary_integrand[1,j,i] = (alpha_x*nx + alpha_y*ny)*q # Boundary Flux
    	end
    end

    boundaryintegrate!(sbp, mesh.bndryfaces[idx_range], boundary_integrand, boundary_functional)
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
                         eqn::AdvectionData{Tsol}, opts, functional_edges)

  # Specify the boundary conditions for the edge on which the force needs to be
  # computed separately. Use that boundary number to access the boundary 
  # offset array. Then proceed the same as bndryflux to get the integrand. Finally
  # use integratefunctional! to get the solution.

  functional_val = zero(Tsol)

  for itr = 1:length(functional_edges)
    g_edge_number = functional_edges[itr] # Extract geometric edge number
    start_index = mesh.bndry_offsets[g_edge_number]
    end_index = mesh.bndry_offsets[g_edge_number+1]
    idx_range = start_index:(end_index-1)
    bndry_facenums = view(mesh.bndryfaces, idx_range) # faces on geometric edge i

    nfaces = length(bndry_facenums)
    boundary_integrand = zeros(Tsol, 1, mesh.sbpface.numnodes, nfaces)
    # boundary_force = zeros(Tsol, 1, sbp.numnodes, mesh.numEl)
    alpha_x = eqn.alpha_x
    alpha_y = eqn.alpha_y

    for i = 1:nfaces
      bndry_i = bndry_facenums[i]
      global_facenum = idx_range[i]
      for j = 1:mesh.sbpface.numnodes
        q = eqn.q_bndry[ 1, j, global_facenum]
        coords = view(mesh.coords_bndry, :, j, global_facenum)
        dxidx = view(mesh.dxidx_bndry, :, :, j, global_facenum)
        nrm = view(sbp.facenormal, :, bndry_i.face)
        nx = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
        ny = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]
        # println("[nx, ny] = ", [nx, ny])
        boundary_integrand[1,j,i] = calcFunctionalIntegrand(alpha_x, alpha_y, nx, ny, q) # Boundary Flux
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

@doc """
###AdvectionEquationMod.calcFunctionalIntegrand

Calculates the functional integrand for a particular q.

**Inputs**

*  `alpha_x` : eqn.alpha_x
*  `alpha_y` : eqn.alpha_y
*  `nx`      : x-component of normal vector
*  `ny`      : y-component of normal vector
*  `q`       : eqn.q or eqn.q_bndry at the particular node

**Outputs**

*  `functional_integrand` : value of the functional integrand at that node 

"""->

function calcFunctionalIntegrand(alpha_x, alpha_y, nx, ny, q)
  
  return functional_integrand = (alpha_x*nx + alpha_y*ny)*q

end

function getFunctionalBoundaryq()

  n_functional_faces = 0  # Total length of the interpolated q values across all geometric functional edges
  for i = 1:length(functional_edges)
    nfaces = 0
    g_edge_number = functional_edges[i]
    start_index = mesh.bndry_offsets[g_edge_number]
    end_index = mesh.bndry_offsets[g_edge_number+1]
    idx_range = start_index:(end_index-1)  # Index range
    nfaces = length(mesh.bndryfaces[idx_range])
    n_functional_faces += nfaces
  end  # End for i = 1:length(functional_edges)

  fq_bndry_vec = zeros(Tsol,n_functional_faces)

  # Populate fq_bndry
  for i = 1:length(functional_edges)
    g_edge_number = functional_edges[itr] # Extract geometric edge number
    start_index = mesh.bndry_offsets[g_edge_number]
    end_index = mesh.bndry_offsets[g_edge_number+1]
    idx_range = start_index:(end_index-1)
    bndry_facenums = view(mesh.bndryfaces, idx_range) # faces on geometric edge i
    nfaces = length(bndry_facenums)

    for j = 1:nfaces
      bndry_j = bndry_facenums[j]
      global_facenum = idx_range[j]
      for k  = 1:mesh.sbpface.numnodes
        # indexing to jump nfaces per geometric edge and populate q_values
        fq_bndry_vec[(i-1)*nfaces + k] = eqn.q_bndry[1, k, global_facenum]
      end  # End for j = 1:nfaces
    end    # End for j = 1:nfaces
  end      # End for i = 1:length(functional_edges)

  return nothing
end


function functionalIntegration()

  # disassemble fq_bndry_vec
  for i = 1:length(functional_edges)


  return nothing
end