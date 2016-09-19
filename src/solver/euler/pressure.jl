# Store the coefficient of pressure in and compare it against the

function calcPressureCoeff(params, q)

  # Assumption of calcFeeStream is speed of sound, c = 1 and rho_inf = 1

  p_inf = 1/params.gamma # Free stream pressure. Refer Pulliam's notes
  m_inf = 0.5*params.Ma*params.Ma

  press = calcPressure(params, q)

  pressCoeff = (press - p_inf) #/m_inf

  return pressCoeff
end

@doc """
### writeSurfacePressureCoeff

Writes the pressure coefficient for multiple boundary edges to file. Every rank
writes for its portion of the geometric edge. the files are stored in a
separate directory called pressCoeff. The function deletes any existing
`pressCoeff` directory and creates a new one every time its called.

**Arguments**

*  `mesh` : Abstract mesh type
*  `sbp`  : Summation-By-Parts Operator
*  `eqn`  : Euler equation object
*  `opts` : Options dictionary
*  `g_edges` : geometric edges over which the coefficient of pressure is to be
               computed

"""

function writeSurfacePressureCoeff{Tmsh, Tsol}(mesh::AbstractDGMesh{Tmsh},
         sbp::AbstractSBP, eqn::EulerData{Tsol}, opts, g_edges, nodal_pressCoeff)

  my_rank = MPI.Comm_rank(eqn.comm)  # Get the rank of the process

  if my_rank == 0 # Delete existing directory & Create new directory
    if isdir("./pressCoeff")
      rm("./pressCoeff"; recursive=true)
    end
    mkdir("./pressCoeff")
  end

  MPI.Barrier(eqn.comm)

  for itr = 1:length(g_edges)
    g_edge_number = g_edges[itr] # Extract geometric edge number
    filename = string("pressCoeff_edge", g_edge_number, "_rank", my_rank,".dat")
    f = open(joinpath("./pressCoeff",filename), "w")
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
    q2 = zeros(Tsol, mesh.numDofPerNode)

    for i = 1:nfaces
      bndry_i = bndry_facenums[i]
      global_facenum = idx_range[i]
      for j = 1:mesh.sbpface.numnodes
        q = sview(eqn.q_bndry, :, j, global_facenum)
        convertToConservative(eqn.params, q, q2)
        pressCoeff = calcPressureCoeff(eqn.params, q2)
        println(f, pressCoeff)
        nodal_pressCoeff[itr][j,i] = pressCoeff
      end  # End for j = 1:mesh.sbpface.numnodes
    end    # End for i = 1:nfaces

    close(f)  # Close filestream

  end   # End for itr = 1:length(edges)


  return nothing
end # End evalSurfacePressureCoeff

function readSurfacePressureCoeff{Tmsh, Tsol}(mesh::AbstractDGMesh{Tmsh},
         sbp::AbstractSBP, eqn::EulerData{Tsol}, opts, g_edges, nodal_pressCoeff)

  my_rank = MPI.Comm_rank(eqn.comm)  # Get the rank of the process

  for itr = 1:length(g_edges)
    g_edge_number = g_edges[itr] # Extract geometric edge number

    # Read in the pressure coefficient from file
    filename = string("pressCoeff_edge", g_edge_number, "_rank", my_rank,".dat")
    f = open(joinpath("./pressCoeff",filename), "r")
    Cp = readdlm(f,'n',Float64)
    close(f)

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
    q2 = zeros(Tsol, mesh.numDofPerNode)

    itr3 = 1
    for i = 1:nfaces
      bndry_i = bndry_facenums[i]
      global_facenum = idx_range[i]
      for j = 1:mesh.sbpface.numnodes
        nodal_pressCoeff[itr][j,i] = Cp[itr3]
        itr3 += 1
      end
    end # fir i = 1:nfaces
  end # End for itr = 1:length(g_edges)

  return nothing
end
