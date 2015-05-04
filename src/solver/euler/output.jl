# output functions

function printSolution(mesh::PumiMesh2, u::AbstractVector)
# print solution to file
# format = for each element, node, print u rho*u rho*v E

println("entered printSolution")
f = open("solution.dat", "w")

for i=1:mesh.numEl
  dofnums_i = getGlobalNodeNumbers(mesh, i)

  for j=1:3
    u_vals = u[dofnums_i[:, j]]
    str = @sprintf("%d %d %16.15e %16.15e %16.15e %16.15e \n", i, j, u_vals[1], u_vals[2], u_vals[3], u_vals[4] )  # print element number

    print(f, str)
  end
#  print(f, "\n")
end

close(f)
return nothing

end


function printCoordinates(mesh::PumiMesh2)
# print solution to file
# format = for each element, node, print u rho*u rho*v E

println("entered printCoordinates")
f = open("coords_output.dat", "w")

for i=1:mesh.numEl
  dofnums_i = getGlobalNodeNumbers(mesh, i)
  coords = getElementVertCoords(mesh, [i])

  for j=1:3
    coords_j = coords[:,j]
    str = @sprintf("%d %d %16.15e %16.15e %16.15e \n", i, j, coords_j[1], coords_j[2], coords_j[3] )  # print element number

    print(f, str)
  end
#  print(f, "\n")
end

close(f)
return nothing

end




