function calcDofs(m, p)
# calculate number of dofs in a rectangular domain with the same number of
# elements along both axes
# m = number of elements per side
# p = order of SBP operators

 n = m  # number of elements in horiztonal direction
 numV = (m+1)*(n+1)
 numEdge = m*(n+1) + n*(m+1)
 numEl = m*n

 counts = [numV, numEdge, numEl]

  if p == 1
    nodecnts = [1., 0, 0]
  elseif p == 2
    nodecnts = [1., 1, 0]
  elseif p == 3
    nodecnts = [1., 2, 3]
  elseif p == 4
    nodecnts = [1., 3, 6]
  else
    println("Error: unsupported p value requested")
  end

  numDof = dot(nodecnts, counts)

#  println("for a $m x $m mesh, there are $numDof dofs")
  return numDof
end


function findMeshSize(p1_el)
# find the size of the mesh with the closest number of dofs to 
# the p1 mesh

  p1dofs = calcDofs(p1_el, 1)
  println("the p1 mesh has ", p1dofs, " dofs")
  for p=2:4
    mesh_size = 1
    mesh_size_upper = 0
    for i=1:p1_el  # loop over mesh sizes
      dofs_i = calcDofs(i, p)
      if dofs_i > p1dofs  # break on first mesh larger than needed
        mesh_size_upper = i
        break
      end
    end
      # check whether the slightly larger of slightly smaller mesh is closest
      dofs_upper = calcDofs(mesh_size_upper, p)
      dofs_lower = calcDofs(mesh_size_upper - 1, p)

      if abs(dofs_upper - p1dofs) < abs(dofs_lower - p1dofs)
        mesh_size_i = mesh_size_upper
        dof_mesh_size_i = dofs_upper
      else
        mesh_size_i = mesh_size_upper - 1
        dof_mesh_size_i = dofs_lower
      end

        println("p $p mesh should have $mesh_size_i elements per side, and therefore $dof_mesh_size_i dofs")
    end  # end loop over p values

  return nothing
end

findMeshSize(100)

