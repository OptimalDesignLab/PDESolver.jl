# solver code for the Euler equation

function evalSCResidual(mesh::AbstractMesh, operator::AbstractSBP, u::AbstractVector, u0, alpha_x::FloatingPoint, alpha_y::FloatingPoint)
# evaluate the residual using summation by parts (not including boundary integrals)
# this only works for triangular meshes, where are elements are same
# mesh : mesh type
# operator: what operator to use (SBP, lagrange polynomial ...)
# u : solution vector (must be ndof by 1) to be populated
# u0 : solution vector at previous timestep
# alpha_x and alpha_y : advection velocities in x and y directions

# not clear how to do source term using SBP operators

# initialize
#println("entered evalResidual")
#println("u = \n", u)
#println("u0 =\n", u0)

ndof = getNumNodes(mesh)
numEl = getNumEl(mesh)

# count the number of nodes per element
nnodes = operator.numnodes  # BAD: should not directly access fields
                            # SBP needs to provide access this
ub = zeros(nnodes) # holds original solution at for an element
fluxes = zeros(nnodes, 2)  # jacobian term times advection velocity divided by jac

Q_xi = operator.Q[:,:,1]
Q_eta = operator.Q[:,:,2]
#println("Q_xi = \n", Q_xi)
#println("Q_eta = \n", Q_eta)

p#rintln("performing volume integrals")
for i=1:numEl  # loop over element
#  println("\nelement i = ", i)
  dofnums_i = getGlobalNodeNumbers(mesh, i)

  ub = u0[dofnums_i]  # get original solution values
#  println("ub = \n", ub)

  vert_coords = getElementVertCoords(mesh, [i])
#  println("vert_coords = ", vert_coords)
#  println("size(vert_coords) = ", size(vert_coords))
  vert_coords_extract = vert_coords[1:2, :, :]
#  println("vert_coords_extract = ", vert_coords_extract)
#  println("type of vert_coords_extract = ", typeof(vert_coords_extract))
#  println("size of vert_coords_extract = ", size(vert_coords_extract))

  dxi_dx = zeros(2, 2, nnodes, 1)
#  println("size(dxi_dx) = ", size(dxi_dx))
  jac = zeros(nnodes, 1)
  mappingjacobian!(operator, vert_coords_extract, dxi_dx, jac)

#  println("jac = ", jac)
  dxi_dxu = zeros(nnodes,2)


  for j=1:nnodes
#    println("\nj = ", j)
#    println("jac_j = \n", dxi_dx[:,:,j,1])


    dxi_dxu[j,1] += (dxi_dx[1,1,j,1]*alpha_x + dxi_dx[1,2,j,1]*alpha_y)*ub[j]
    dxi_dxu[j,2] += (dxi_dx[2,1,j,1]*alpha_x + dxi_dx[2,2,j,1]*alpha_y)*ub[j]
  end

#  println("dxi_dxu = ", dxi_dxu)
  # multiply by stiffness matrix here because it is the same for all points in the element


  res = zeros(nnodes,1)
#  res_x = zeros(nnodes, 1)
#  res_y = zeros(nnodes, 1)
 
  # copy dxi_dxu into column matricies to satisfy SBP
  dxi_dxu_x = zeros(nnodes,1)
#  println("type of dxi_dxu_x = ", dxi_dxu_x)
  dxi_dxu_x[:] = dxi_dxu[:,1]
  dxi_dxu_y = zeros(nnodes,1)
  dxi_dxu_y[:] = dxi_dxu[:,2]
#  println("dxi_dxu_x = ", dxi_dxu_x)
#  println("dxi_dxu_y = ", dxi_dxu_y)
  weakdifferentiate!(sbp, 1, dxi_dxu_x, res, trans=true)
  weakdifferentiate!(sbp, 2, dxi_dxu_y, res, trans=true)

#  println("dxi_dxu_x = \n", dxi_dxu_x)
#  println("dxi_dxu_y = \n", dxi_dxu_y)
#  println("res_x = \n", res_x)
#  println("res_y = \n", res_y)

#  res = res_x + res_y

#  println("res = \n", res)
#  println("dofnums_i = ", dofnums_i)

  u[dofnums_i.'] += res
#  println("u = \n", u)
end  # end loop over elements

#println("after volume integral, u = \n", u)

end  # end function


function evalBndry(mesh::PumiMesh2, operator::AbstractSBP, u::AbstractVector, u0, alpha_x::FloatingPoint, alpha_y::FloatingPoint)
# evaluate boundary integrals for advection equation
# mesh : mesh type, must be PumiMesh (need adjacency information
# operator : SBP operator
# u : solution vector to be populated, has already has volume integrals done
# u0 : solution vector at previous timestep
# alpha_x and alpha_y : advection veocities in x and y directions


println("entered evalBndry")
#println("u = \n", u)
#println("u0 = \n", u0)

# get arguments needed for sbp boundaryintegrate!

bndry_edges = getBoundaryEdgeNums(mesh)
num_bndry_edges = length(bndry_edges)

bndry_faces = Array{Boundary}( num_bndry_edges)
bndry_sign = zeros(Int, num_bndry_edges)
bndry_orientation = zeros(Int, num_bndry_edges)
#println("getting boundary edge numbers")
for i=1:num_bndry_edges  # extract element and edge numbers
#  println("i = ", i)
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
    end
  end

#=
  # modify for orientation of parent element
  coords = getElementVertCoords(mesh, [facenum_i])
  println("coords = ", coords)
  println("before correction, edgenum = ", edgenum_local)
  if (coords[1,1] < coords[1,2]) && (coords[2,1] = coords[2,2] < 1e-12)  # first point is to left of second pointa
    println("correcting horizontally flipped element")
    if edgenum_local == 2  # reverse 2nd, 3rd edges
      edgenum_local = 3
      bndry_sign[i] = -1  # reverse sign of normal
    elseif edgenum_local == 3
      edgenum_local = 2
      bndry_sign[i] = -1  # reverse sign of normal
    end
  end


  if (coords[2,1] < coords[2,2]) && (coords[1,1] - coords [3,1] < 1e-12)
    println("correcting vertically flipped element")
    if edgenum_local == 1
      edgenum_local == 2  # reverse normal components here
      bndry_orientation[i] = -1  # reverse normal components
    elseif edgenum_local == 2
      edgenum_local = 1
      bndry_sign[i] = -1
    end
  end
=#
#  println("after correction, edgenum_local = ", edgenum_local)
  bndry_faces[i] = Boundary(facenum_i, edgenum_local)
end

#println("bndry_faces = \n", bndry_faces)

for i=1:num_bndry_edges
  bndry_i = bndry_faces[i]
#  println("i = ", i, " element ", bndry_i.element, "face ", bndry_i.face)
end

#println("size(bndry_faces) = ", size(bndry_faces))

# create new u formatted for SBP
# also get coordinates needed for  mapping jacobian
u_sbp = zeros(operator.numnodes, mesh.numEl)
x = zeros(2, operator.numnodes, mesh.numEl)

#println("u = ", u)
for i=1:mesh.numEl
  dofnums_i = getGlobalNodeNumbers(mesh, i)
  u_sbp[:,i] = u0[dofnums_i]

  x_i = getElementVertCoords(mesh, [i])
  x[:,:,i] = x_i[1:2,:]
end

dxi_dx = zeros(2,2,operator.numnodes, mesh.numEl)
jac = zeros(operator.numnodes, mesh.numEl)

#println("typeof dxi_dx = ", typeof(dxi_dx))

# populate dxi_dx and jac
mappingjacobian!(operator, x, dxi_dx, jac)

res = zeros(operator.numnodes, mesh.numEl)

# define the flux function
u_bc = sin(-1)  # boundary condition (for all sides)
#u_bc = 0.0

  cntr = 1  # count number of times flux function is called
  function flux1(u_sbp_, dxi_dx_, nrm)
    # u_sbp_ is the entry from u_sbp for this node
    # dxi_dx is the jacobian for this node
    # nrm is the normal vector


#    println("\nentered flux1")
#    println("cntr = ", cntr)
#    println("u_sbp_ = ", u_sbp_)
#    println("dxi_dx_ = \n", dxi_dx_)
#    println("u_sbp_ = ", u_sbp_)
#    println("u_bc = ", u_bc)
    u_xi = dxi_dx_[1,1]*alpha_x + dxi_dx_[1,2]*alpha_y
    u_eta = dxi_dx_[2,1]*alpha_x + dxi_dx_[2,2]*alpha_y

#    println("nrm = ", nrm)
#=
    # modify normal vector
    cnt_i = fld((cntr-1),2)
    cnt_i += 1

    println("cnt_i = ", cnt_i)
    if bndry_sign[cnt_i] == -1  # modify sign
      println("negating sign of normal")
      nrm = -1*nrm
    end

    if bndry_orientation[cnt_i] == -1  # reverse componentsa
      println("reversing normal components")
      tmp = nrm[1]
      nrm[1] = nrm[2]
      nrm[2] = tmp
    end
=#
#    println("nrm_mod = ", nrm)
#    println("u_xi = ", u_xi)
#    println("u_eta = ", u_eta)

    mag_x = u_xi*nrm[1]  # x flow magnitude
    mag_y = u_eta*nrm[2]  # y flow magnitude

    mag = mag_x + mag_y

#    println("mag_x = ", mag_x, " mag_y = ", mag_y)

    # because mag is scalar,not vector, can combine these if statements
    if mag < 0  # inflow condition, added factor of nrm[2]
#      println("xi inflow condition")
      flx_xi = u_xi*u_bc
    else
#      println("xi outflow condition")
      flx_xi = u_xi*u_sbp_
    end

    if mag < 0  # inflow condition, added factor of nrm[1]
#      println("eta inflow condition")
      flx_eta = u_eta*u_bc
    else
#      println("eta outflow condition")
      flx_eta = u_eta*u_sbp_
    end


#    println("flx_xi = ", flx_xi, " flx_eta = ", flx_eta)
    net_flux = flx_xi*nrm[1] + flx_eta*nrm[2]
#    println("net_flux = ", net_flux)

    cntr += 1
    return net_flux

  end

#  println("\nu_sbp = \n", u_sbp)
  boundaryintegrate!(operator, bndry_faces, u_sbp, dxi_dx, flux1, res)
#  println("finished performing boundaryintegrate!()\n")

#  println("before applying boundary integral, u = \n", u)
#  println("\nu_sbp = \n", u_sbp)

  bndry_contrib = zeros(mesh.numDof)
#println("assembling boundary terms")
for i=1:mesh.numEl
  dofnums_i = getGlobalNodeNumbers(mesh, i)
#  println("dofnums_i = ", dofnums_i)
  vals_i = res[:,i]

  u[dofnums_i.'] -= vals_i
  bndry_contrib[dofnums_i.'] -= vals_i

end

#  println("res = \n", res)
#  println("bndry_contrib = \n", bndry_contrib)
#  println("after boundary integral, u = \n", u)

end # end of outer function

function getMass(operator, mesh)
  # assemble mesh
  numnodes = mesh.numNodes  # number of dofs
  numdof = numnodes*mesh.numDofPerNode
  mass_matrix = zeros(numdof, numdof)

  for i=1:mesh.numEl
    dofnums_i = getGlobalNodeNumbers(mesh, i)
    nnodes = size(dofnums_i)[2]  # number of nodes
    for j=1:nnodes
      for k=1:mesh.numDofPerNode
        dofnum_k = dofnums_i[k,j]
        mass_matrix[dofnum_k, dofnum_k] += sbp.w[j]
      end
    end
  end

  return mass_matrix
end



  

