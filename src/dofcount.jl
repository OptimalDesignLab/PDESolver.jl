"""
  This function returns an array with the number of nodes classified on each
  mesh vertex, edge, face, and region (in 3d) for the different operator
  types and degrees.

  op_type    Description
    1        SBP Gamma, CG, 2d
    2        SBPGamma, DG, 2d
    3        SBPOmega, DG 2d
    4        SBPGamma, DG, 3d
    5        SBPOmega, DG, 3d

  Degree 1 - 4 is supported for all operators.
  An exception will be thrown if an unsupported op_type and degree
  combinatino is specified.

  Inputs:
    op_type: integer specifying which operator type
    p: degree of the operator (linear, quadratic, etc.)

  Outputs:
    nodecnts: array of length 3 in 2D or 4 in 3D specifying number of
              nodes classified on each dimension mesh entity
"""
function getNodeCounts(op_type, order)
# get the number of nodes classifies on each vert, edge, and element
# op type specifies the type of operator, order the order of accuracy
  if op_type == 1  # SBPGamma, CG, 2d
    if order == 1
      return nodecnts = [1, 0, 0]
    elseif order == 2
      return nodecnts = [1, 1, 1]
    elseif order == 3
      return nodecnts = [1, 2, 3]
    elseif order == 4
      return nodecnts = [1, 3, 6]
    else
      error("unsuported op_type/order")
    end
  elseif op_type == 2  # SBPGamma DG, 2D
    if order == 1
      return [0, 0, 3]
    elseif order == 2
      return [0, 0, 7]
    elseif order == 3
      return [0, 0, 12]
    elseif order == 4
      return [0, 0, 18]
    else
      error("unsupported op_type/order")
    end
  elseif op_type == 3  # SBPOmega DG, 2D
    if order == 1
      return [0, 0, 3]
    elseif order == 2
      return [0, 0, 6]
    elseif order == 3
      return [0, 0, 10]
    elseif order == 4
      return [0, 0, 15]
    else
      error("unsupported op_type/order")
    end
  elseif op_type == 4  # SBPGamma DG, 3D
    if order == 1
      return [0, 0, 0, 4]
    elseif order == 2
      return [0, 0, 0, 11]
    elseif order == 3
      return [0, 0, 0, 24]
    elseif order == 4
      return [0, 0, 0, 45]
    else
      error("unsupported op_type/order")
    end
  elseif op_type == 5  # SBPOmega DG, 3D
    if order == 1
      return [0, 0, 0, 4]
    elseif order == 2
      return [0, 0, 0, 10]
    elseif order == 3
      return [0, 0, 0, 20]
    elseif order == 4
      return [0, 0, 0, 38]
    else
      error("unsupported op_type/order")
    end
 
  end

end


function calcDofs{T <: Integer}(m, nodecnt::AbstractArray{T})
# calculate number of dofs in a rectangular domain with the same number of
# elements along both axes
# m = number of elements per side
# p = order of SBP operators

  if length(nodecnt) == 3
    n = m  # number of elements in horiztonal direction
    numV = (m+1)*(n+1)
    numEdge = m*(n+1) + n*(m+1) + m*n
    numEl = 2*m*n

  #  counts = [numV, numEdge, numEl]
    return numV*nodecnt[1] + numEdge*nodecnt[2] + numEl*nodecnt[3]
  else
    nx = m
    ny = nx
    nz = nx
    nvert = (nx+1)*(ny+1)*(nz+1)

    nedges = nx*(ny+1)*(nz+1) +  # cube dges
             ny*(nx+1)*(nz+1) +
             nz*(nx+1)*(ny+1) +  
             nx*nz*(ny+1) +  # cube face edges
             ny*nz*(nx+1) +
             nx*ny*(nz+1) +
             nx*ny*nz  # cube interior

    nfaces = nx*ny*(nz+1)*2 +
             ny*nz*(nx+1)*2 +
             nx*ny*(ny+1)*2 +
             nx*ny*nz*6

    nel = nx*ny*nz*6

    return nvert*nodecnt[1] + nedges*nodecnt[2] + nfaces*nodecnt[3] + nel*nodecnt[4]
  end


  return 0
#  numDof = dot(nodecnts, counts)

#  println("for a $m x $m mesh, there are $numDof dofs")
#  return numDof
end

function calcDofs3D()

  ny = nx
  nz = nx
  nvert = (nx+1)*(ny+1)*(nz+1)

  nedges = nx*(ny+1)*(nz+1) +  # cube dges
           ny*(nx+1)*(nz+1) +
           nz*(nx+1)*(ny+1) +  
           nx*nz*(ny+1) +  # cube face edges
           ny*nz*(nx+1) +
           nx*ny*(nz+1) +
           nx*ny*nz  # cube interior

  nfaces = nx*ny*(nz+1)*2 +
           ny*nz*(nx+1)*2 +
           nx*ny*(ny+1)*2 +
           nx*ny*nz*6

  nel = nx*ny*nz*6

  return nvert*nodecnt[1] + nedges*nodecnt[2] + nfaces*nodecnt[3] + nel*nodecnt[4]

end

"""
  Find the mesh with specified operator operator type and degree elements 
  with the closest number of degrees of freedom to the specified number

  Inputs:
    numdof: the number of dofs to try to match
    op_type: operator type (see getNumNodes)
    p: degree of operator

  Outputs:
    mesh_size_i: number of degree p element in each direction
    dof_mesh_size_i: total number of degrees of freedom on the mesh

"""
function findMeshSize(numdof::Integer, op_type, p::Integer)
# find the size of the square mesh with elements of order p 
# with the closest number of dofs to 
# the p1 mesh with p1_el elements per side,
# returns the number of elements per side, total number of dofs
  
#  p1dofs = calcDofs(p1_el, 1)
#  println("the p1 mesh has ", p1dofs, " dofs")

  # need to calculate an upper bound on the number of meshes to search
  # use the number of elements per side of a p1 mesh
  println("searching for p$p mesh with $numdof dofs")

  nodecnt = getNodeCounts(op_type, p)

  p1_el = round(Int, ceil(sqrt(numdof))) + 1 + 1

  mesh_size_i = 0  # declare in this scope
  dof_mesh_size_i = 0
    mesh_size = 1
    mesh_size_upper = 0
    for i=1:p1_el  # loop over mesh sizes
      dofs_i = calcDofs(i, nodecnt)
      if dofs_i > numdof  # break on first mesh larger than needed
        mesh_size_upper = i
        break
      end
    end
    # check whether the slightly larger of slightly smaller mesh is closest
    dofs_upper = calcDofs(mesh_size_upper, nodecnt)
    dofs_lower = calcDofs(mesh_size_upper - 1, nodecnt)

    if abs(dofs_upper - numdof) < abs(dofs_lower - numdof)
      mesh_size_i = mesh_size_upper
      dof_mesh_size_i = dofs_upper
    else
      mesh_size_i = mesh_size_upper - 1
      dof_mesh_size_i = dofs_lower
    end
    println("p $p mesh should have $mesh_size_i elements per side, and therefore $dof_mesh_size_i dofs")

  return mesh_size_i, dof_mesh_size_i
end

"""
  This function finds the number of element in each direction of a square
  (or cubic) mesh that has the closest number of degrees of freedom to
  a p=1 mesh with the specified number of elements

  Inputs:
    p1_el: number of element in each direction of the p=1 mesh
    op_type: operator type (see getNodeCounts)
    p: degree of operator

  Outputs:
    mesh_size_i: number of degree p element in each direction
    dof_mesh_size_i: total number of degrees of freedom on the mesh
"""
function findMeshSizeEl(p1_el, op_type, p)
# find the size of of mesh with the closest numbere of dofs to a p1 mesh
# of the specified size
# the mesh has the same number of element along the x and y axes
# returns the number of element per side and the total number of dofs

  nodecnt = getNodeCounts(op_type, 1)
  p1dofs = calcDofs(p1_el, nodecnt)
  println("the p1 mesh has ", p1dofs, " dofs")
  return findMeshSize(p1dofs, op_type, p)
end

"""
  This function calls findMeshSizeEl for a range of element sizes
  and writes the results to a file.
"""
function run_el()
  p1_mesh_sizes = collect(11:5:50)
  p = 4
  fname = "mesh_counts.dat"
  if isfile(fname)
    rm(fname)
  end
  f = open(fname, "w")
  for i=1:length(p1_mesh_sizes)
    els, dofs = findMeshSizeEl(p1_mesh_sizes[i], p)
    println(f, els)
  end

  close(f)

  return nothing
end

"""
  This function calls findMeshSize for a range of dof values and writes the
  resutls to a file.
"""
function run_dof()
  dof_spacing = collect(linspace(2751, 12865, 10))
  dof_counts = zeros(Int, length(dof_spacing))
  for i=1:length(dof_spacing)
    dof_counts[i] = round(dof_spacing[i])
  end
  p = 2
  fname = "mesh_counts.dat"
  if isfile(fname)
    rm(fname)
  end
  f = open(fname, "w")
  for i=1:length(dof_counts)
    els, dofs = findMeshSize(dof_counts[i], p)
    println(f, els)
  end

  close(f)

  return nothing
end


#run_dof()
op_type = 3
p = 4
numel = 18
numdof = 64^3
findMeshSizeEl(numel, op_type, p)
#findMeshSize(numdof, op_type, p)


