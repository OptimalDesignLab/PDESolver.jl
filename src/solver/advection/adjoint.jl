# Adjoint Computation
@doc """
### AdvectionEquationMod.calcAdjoint

Calcualtes the adjoint vector for a single functional

**Inputs**

*  `mesh` : Abstract mesh type
*  `sbp`  : Summation-By-Parts operator
*  `eqn`  : Advection equation object
*  `opts` : Options dictionary
*  `functional_number` : The functional for which the adjoint vector is being
                         computed

**Outputs**

*  None

"""->
function calcAdjoint{Tmsh, Tsol}(mesh::AbstractDGMesh{Tmsh},sbp::AbstractSBP,
                                 eqn::AdvectionData{Tsol}, opts, functional_number)

  # Get information corresponding to functional_number
  key = string("geom_edges_functional", functional_number)
  functional_edges = opts[key]

  # Calculate the Jacobian of the residual
  newton_data = NonlinearSolvers.NewtonData(mesh, sbp, eqn, opts)
  res_jac = zeros(Tres, mesh.numDof, mesh.numDof)
  pert = complex(0, opts["epsilon"])
  NonlinearSolvers.calcJacobianComplex(newton_data, mesh, sbp, eqn, opts, evalAdvection, pert, res_jac)


  # calculate the derivative of the function w.r.t q_vec
  func_deriv = zeros(Tsol, mesh.numDof)
  adjoint_vec = zeros(Tsol, mesh.numDof)
  
  # Get a vector of interpolated q values along all element faces of the 
  # the functional edges
  n_functional_faces = 0  # Total length of the interpolated q values across all geometric functional edges
  for i = 1:length(functional_edges)
    g_edge_number = functional_edges[i]
    start_index = mesh.bndry_offsets[g_edge_number]
    end_index = mesh.bndry_offsets[g_edge_number+1]
    idx_range = start_index:(end_index-1)  # Index range
    nfaces = length(mesh.bndryfaces[idx_range])
    n_functional_faces += nfaces
  end  # End for i = 1:length(functional_edges)
  
  # Create a book-keeping array of tuples for storing information on 
  # 1. sbpface.numnodes index
  # 2. global facenum
  # 3. functional edge number

  fq_bndry_info = Array{Tuple{Int, Int, Int}}(n_functional_faces) # Creat
  starting_index = 0
  for itr = 1:length(functional_edges)
    g_edge_number = functional_edges[itr]
    start_index = mesh.bndry_offsets[g_edge_number]
    end_index = mesh.bndry_offsets[g_edge_number+1]
    idx_range = start_index:(end_index-1)
    bndry_facenums = view(mesh.bndryfaces, idx_range) # faces on geometric edge i
    nfaces = length(bndry_facenums)

    for i = 1:nfaces
      bndry_i = bndry_facenums[i]
      global_facenum = idx_range[i]
      for j = 1:mesh.sbpface.numnodes
        fq_bndry_info[starting_index + (i-1)*mesh.sbpface.numnodes + j] = 
                                               (j,global_facenum,g_edge_number)
      end
    endc
    nfaces_prev = nfaces
    starting_index += nfaces_prev*mesh.sbpface.numnodes
  end
end  # End function calcAdjoint

function calcFunctionalDeriv(mesh, sbp, eqn, opts, fq_bndry, functional_edges)
  
  # create a 1D array of length n_functional_faces to store df/dq_bndry
  df_dqbndry = zeros(n_functional_faces)

  pert = complex(0, opts["epsilon"])  # complex perturbation

  # populate df_dqbndry
  for i = 1:n_functional_faces
     
    # get q_bndry
    sbpface_index, global_facenum, g_edge_number = fq_bndry_info[i]

    # perturb the q_bndry corresponding to the index
    eqn.q_bndry[1,sbpface_index,global_facenum] += pert

    # Use the above perturbed value to get the functional value
    functional_val = calcBndryfunctional(mesh, ebp, eqn, opts, functional_edges)
    df_dqbndry[i] = imag(functional_val)/norm(pert)
    
  end  # end for i = 1:n_functional_faces

  # unperturb this eqn.q_bndry
  eqn.q_bndry[1,sbpface_index,global_facenum] -= pert

  return nothing
end