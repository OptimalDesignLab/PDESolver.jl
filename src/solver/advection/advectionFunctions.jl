# advectionFunctions.jl

@doc """
### evalAdvection

This function evaluates the Advection equation and preps it for the RK4 solver.
Pass this function as an input argument to the RK4 solver just like evalEuler.

"""->

function evalAdvection(mesh::AbstractMesh, sbp::SBPOperator, eqn::AdvectionData, opts, t=0.0)

  # u_i_1 = zeros(mesh.numDof)
  # eqn.res_vec = fill!(eqn.res_vec, 0.0)
  evalSCResidual(mesh, sbp, eqn.res_vec, eqn.q_vec, alpha_x, alpha_y)
  evalBndry(mesh, sbp, eqn.res_vec, eqn.q_vec, alpha_x, alpha_y)

  eqn.res_vec = mass_matrix\eqn.res_vec
  # dissamble to eqn.u
  eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.res_vec)

  eqn.res = copy(eqn.u) # transfer eqn.u to eqn.res for assembly in RK
  # return u_i_1
  return nothing
end

@doc """
### evalSCResidual

Evaluate the residual using summation by parts (not including boundary 
integrals) this only works for triangular meshes, where are elements are same

* mesh : mesh type
* operator: what operator to use (SBP, lagrange polynomial ...)
* u : solution vector (must be ndof by 1) to be populated
* u0 : solution vector at previous timestep
* alpha_x and alpha_y : advection velocities in x and y directions

"""->

function evalSCResidual(mesh::AbstractMesh, operator::SBPOperator, 
	                      u::AbstractVector, u0, alpha_x::FloatingPoint,
	                      alpha_y::FloatingPoint)

# not clear how to do source term using SBP operators

	# initialize

  ndof = mesh.numDof  # Total number of dofs
  numEl = mesh.numEl  # Total number of elements
  nnodes = mesh.numNodesPerElement # count the number of nodes per element 
  ub = zeros(nnodes) # holds original solution at for an element
  fluxes = zeros(nnodes, 2)  # jacobian term times advection velocity divided
                             # by jac

  # Q_xi = operator.Q[:,:,1]
  # Q_eta = operator.Q[:,:,2]

  for i=1:numEl  # loop over element

    dofnums_i = getGlobalNodeNumbers(mesh, i)
    ub = u0[dofnums_i]  # get original solution values
    vert_coords = getElementVertCoords(mesh, [i])
    vert_coords_extract = vert_coords[1:2, :, :]
    dxi_dx = zeros(2, 2, nnodes, 1)
    jac = zeros(nnodes, 1)
    mappingjacobian!(operator, vert_coords_extract, dxi_dx, jac)
    dxi_dxu = zeros(nnodes,2)

    for j=1:nnodes
      dxi_dxu[j,1] += (dxi_dx[1,1,j,1]*alpha_x + dxi_dx[1,2,j,1]*alpha_y)*ub[j]
      dxi_dxu[j,2] += (dxi_dx[2,1,j,1]*alpha_x + dxi_dx[2,2,j,1]*alpha_y)*ub[j]
    end

    res = zeros(nnodes,1)
 
    # copy dxi_dxu into column matricies to satisfy SBP
    dxi_dxu_x = zeros(nnodes,1)
    dxi_dxu_x[:] = dxi_dxu[:,1]
    dxi_dxu_y = zeros(nnodes,1)
    dxi_dxu_y[:] = dxi_dxu[:,2]
    weakdifferentiate!(sbp, 1, dxi_dxu_x, res, trans=true)
    weakdifferentiate!(sbp, 2, dxi_dxu_y, res, trans=true)
    u[dofnums_i.'] += res
  end  # end loop over elements

  #println("after volume integral, u = \n", u)
end  # end function

@doc """
### evalBndry

Evaluate boundary integrals for advection equation

*  mesh : mesh type, must be PumiMesh (need adjacency information
*  operator : SBP operator
*  u : solution vector to be populated, has already has volume integrals done
*  u0 : solution vector at previous timestep
*  alpha_x and alpha_y : advection veocities in x and y directions

"""->

function evalBndry(mesh::PumiMesh2, operator::SBPOperator, u::AbstractVector,
                   u0, alpha_x::FloatingPoint, alpha_y::FloatingPoint)

  println("entered evalBndry")

  # get arguments needed for sbp boundaryintegrate!

  bndry_edges = mesh.bndryfaces
  num_bndry_edges = mesh.numBoundaryEdges

  if length(mesh.bndryfaces) != mesh.numBoundaryEdges
    println("Error with Boundary!!!!")
  end
  #=
  bndry_faces = Array(Boundary, num_bndry_edges)
  bndry_sign = zeros(Int, num_bndry_edges)
  bndry_orientation = zeros(Int, num_bndry_edges)

  for i=1:num_bndry_edges  # extract element and edge numbers
    edgenum_i = bndry_edges[i]
    edge_i = mesh.edges[edgenum_i]

    # get mesh face associated with edge
    countAdjacent(mesh.m_ptr, edge_i, 2)
    face = getAdjacent(1)[1]  # get the single face (not an array)
    facenum_i = getFaceNumber2(face) + 1  # convert to 1 based indexing
    (down_edges, numedges) = getDownward(mesh.m_ptr, face, 1)
    edgenum_local = 0
    for j = 1:numedges  # find which local edge is edge_i
      if down_edges[j] == edge_i
        edgenum_local = j
      end # end if
    end # end for j = 1:numedges

    println("after correction, edgenum_local = ", edgenum_local)
    bndry_faces[i] = Boundary(facenum_i, edgenum_local)
  end  # end for i=1:num_bndry_edges 
  
  for i=1:num_bndry_edges
    bndry_i = bndry_faces[i]
  end # end for i=1:num_bndry_edges
  
  # create new u formatted for SBP
  # also get coordinates needed for  mapping jacobian
  u_sbp = zeros(operator.numnodes, mesh.numEl)
  x = zeros(2, operator.numnodes, mesh.numEl)

  for i=1:mesh.numEl
    dofnums_i = getGlobalNodeNumbers(mesh, i)
    u_sbp[:,i] = u0[dofnums_i]
    x_i = getElementVertCoords(mesh, [i])
    x[:,:,i] = x_i[1:2,:]
  end

  dxi_dx = zeros(2,2,operator.numnodes, mesh.numEl)
  jac = zeros(operator.numnodes, mesh.numEl)
  mappingjacobian!(operator, x, dxi_dx, jac) # populate dxi_dx and jac
  res = zeros(operator.numnodes, mesh.numEl)
  =#

  # define the flux function
  u_bc = sin(-1)  # boundary condition (for all sides)
  cntr = 1  # count number of times flux function is called
  
  # boundaryintegrate!(operator, bndry_faces, u_sbp, dxi_dx, flux1, res)
  boundaryintegrate!(sbp, mesh.bndryfaces, eqn.bndryflux, eqn.res)
  bndry_contrib = zeros(mesh.numDof)

  for i=1:mesh.numEl
    dofnums_i = getGlobalNodeNumbers(mesh, i)
    vals_i = res[:,i]
    u[dofnums_i.'] -= vals_i
    bndry_contrib[dofnums_i.'] -= vals_i
  end # end for i=1:mesh.numEl

end # end function evalBndry
#=


"""-> =#
