function getNodeCounts(op_type, order)
# get the number of nodes classifies on each vert, edge, and element
# op type specifies the type of operator, order the order of accuracy
  if op_type == 1  # SBPGamma, CG
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
  elseif op_type == 2  # SBPGamma DG
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
  elseif op_type == 3  # SBPOmega DG
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
  end

end


function calcDofs{T <: Integer}(m, nodecnt::AbstractArray{T})
# calculate number of dofs in a rectangular domain with the same number of
# elements along both axes
# m = number of elements per side
# p = order of SBP operators

  n = m  # number of elements in horiztonal direction
  numV = (m+1)*(n+1)
  numEdge = m*(n+1) + n*(m+1) + m*n
  numEl = 2*m*n

#  counts = [numV, numEdge, numEl]

  return numV*nodecnt[1] + numEdge*nodecnt[2] + numEl*nodecnt[3]
#  numDof = dot(nodecnts, counts)

#  println("for a $m x $m mesh, there are $numDof dofs")
#  return numDof
end

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
findMeshSizeEl(30, op_type, 2)
#findMeshSize(771, 4)


