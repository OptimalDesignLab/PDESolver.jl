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
    nfaces = 0
    g_edge_number = functional_edges[i]
    start_index = mesh.bndry_offsets[g_edge_number]
    end_index = mesh.bndry_offsets[g_edge_number+1]
    idx_range = start_index:(end_index-1)  # Index range
    nfaces = length(mesh.bndryfaces[idx_range])
    n_functional_faces += nfaces
  end  # End for i = 1:length(functional_edges)
  
  # create an array to store all q_bndry values for that functional 
  fq_bndry = zeros(Tsol,n_functional_faces)

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
        fq_bndry[(i-1)*nfaces + k] = eqn.q_bndry[1, k, global_facenum]
      end  # End for j = 1:nfaces
    end    # End for j = 1:nfaces
  end      # End for i = 1:length(functional_edges)

  return nothing
end  # End function calcAdjoint

function calcFunctionalDeriv(mesh, sbp, eqn, opts, fq_bndry)

  entry_orig = zero(eltype(fq_bndry))

  for i = 1:length(fq_bndry)
    if i == 1
      entry_orig = fq_bndry[i]
      fq_bndry[i] += pert
    else
      fq_bndry[i-1] = entry_orig # Undo previous iteration perturbation
      entry_orig = fq_bndry[i]
      fq_bndry[i] += pert
    end

    # Compute the functional with perturbed q_bndry
    # - This would require disassembly of fq_bndry to q_bndry
    # - supply q_bndry to a function that computes the functional value
    # 
    # use that to compute the scomplex step derivative using calcJacCol


  end  # end for i = 1:length(fq_bndry)

  # Undo final perturbation
  fq_bndry[i] = entry_orig

  return nothing
end