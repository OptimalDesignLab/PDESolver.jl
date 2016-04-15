# Adjoint Computation
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/NonlinearSolvers"))
using NonlinearSolvers


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
function calcAdjoint{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractDGMesh{Tmsh},sbp::AbstractSBP,
                    eqn::AdvectionData{Tsol, Tres, Tdim}, opts, functional_number, adjoint_vec::Array{Tsol,1})

  # println("\n\nbefore Computation, adjoint_vec = \n", adjoint_vec)

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
  
  func_deriv_arr = zeros(eqn.q) # 3D array into which func_deriv_arr_bndry gets interpolated
  
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
  fq_bndry_info = Array{Tuple{Int, Int, Int}}(n_functional_faces*mesh.sbpface.numnodes)
  starting_index = 0
  for itr = 1:length(functional_edges)
    g_edge_number = functional_edges[itr]
    start_index = mesh.bndry_offsets[g_edge_number]
    end_index = mesh.bndry_offsets[g_edge_number+1]
    idx_range = start_index:(end_index-1)
    bndry_facenums = view(mesh.bndryfaces, idx_range) # faces on geometric edge i
    nfaces = length(mesh.bndryfaces[idx_range])

    for i = 1:nfaces
      bndry_i = bndry_facenums[i]
      global_facenum = idx_range[i]
      for j = 1:mesh.sbpface.numnodes
        #println("starting_index + (i-1)*mesh.sbpface.numnodes + j = ", 
        #  starting_index + (i-1)*mesh.sbpface.numnodes + j, " itr = ", itr, " i = ", i, " j = ", j)
        fq_bndry_info[starting_index + (i-1)*mesh.sbpface.numnodes + j] = 
                                               (j,global_facenum,g_edge_number)
      end   # End for j = 1:mesh.sbpface.numnodes
    end     # End for i = 1:nfaces
    nfaces_prev = nfaces
    starting_index += nfaces_prev*mesh.sbpface.numnodes
  end       # End for itr = 1:length(functional_edges)

  # Calculate df/dq_bndry on edges where the functional is calculated and put 
  # it back in func_deriv_arr
  calcFunctionalDeriv(mesh, sbp, eqn, opts, fq_bndry_info, functional_edges, 
                      func_deriv_arr)  # populate df_dq_bndry

  # Assemble func_deriv
  assembleArray(mesh, sbp, eqn, opts, func_deriv_arr, func_deriv)
  # println("\nAfter Assemble array in calcAdjoint")
  # println("func_deriv_arr = \n", func_deriv_arr)
  # println("\nfunc_deriv = \n", func_deriv)

  # Solve for adjoint vector
  # println("\n\n res_jac = \n", res_jac)

  # TODO: The following operation creates a temporary copy of adjoint_vec, does
  #       the '\' computation and then puts it back into adjoint_vec. This
  #       needs to change.
  adjoint_vec[:] = -(res_jac.')\func_deriv
  
  # println("\n\nadjoint_vec = \n", real(adjoint_vec))
  
  return nothing
end  # End function calcAdjoint

function calcFunctionalDeriv(mesh, sbp, eqn, opts, fq_bndry_info, functional_edges,
                             func_deriv_arr)
  
  # Create a 1D array of length n_functional_faces*sbpface.numnodes to store df/dq_bndry ONLY
  # along the the faces on which the particular functional is defined.
  n_functional_q = length(fq_bndry_info)
  df_dq_bndry = zeros(n_functional_q)

  pert = complex(0, opts["epsilon"])  # complex perturbation

  # populate df_dq_bndry using complex step method
  # TODO: This needs to change to support algorithmic differentiation (reverse
  #       mode) in the future.
  for i = 1:n_functional_q
     
    # get q_bndry
    sbpface_index, global_facenum, g_edge_number = fq_bndry_info[i]

    # perturb the q_bndry corresponding to the index
    eqn.q_bndry[1,sbpface_index,global_facenum] += pert

    # Use the above perturbed value to get the functional value
    functional_val = calcBndryfunctional(mesh, sbp, eqn, opts, functional_edges)
    df_dq_bndry[i] = imag(functional_val)/norm(pert)

    
    # unperturb this eqn.q_bndry
    eqn.q_bndry[1,sbpface_index,global_facenum] -= pert
  
  end  # end for i = 1:n_functional_q
  
  # eqn.q_bndry[1,:,:] += pert
  # functional_val = calcBndryfunctional(mesh, sbp, eqn, opts, functional_edges)

  # println("df_dq_bndry = ", df_dq_bndry)

  func_deriv_arr_bndry = zeros(eqn.q_bndry)  # 3D array that holds df/dq_bndry on ALL edges
  # disassemble df_dq_bndry into appropriate places in func_deriv_arr_bndry
  for i = 1:n_functional_q
    sbpface_index, global_facenum, g_edge_number = fq_bndry_info[i]
    func_deriv_arr_bndry[1,sbpface_index,global_facenum] = df_dq_bndry[i]
  end
  # println("\nfunc_deriv_arr_bndry = \n", func_deriv_arr_bndry)

  # Interpolate func_deriv_arr_bndry to func_deriv_arr
  boundaryToVolumeInterpolation(mesh.sbpface, mesh.bndryfaces, func_deriv_arr,
    func_deriv_arr_bndry)
  # println("\nfunc_deriv_arr = \n", func_deriv_arr)

  return nothing
end

function boundaryToVolumeInterpolation{Tsbp,Tsol}(sbpface::TriFace{Tsbp},
                                         bndryfaces::Array{Boundary},
                                         uvol::AbstractArray{Tsol,3},
                                         uface::AbstractArray{Tsol,3})

  for (bindex, bndry) in enumerate(bndryfaces)
    for i = 1:sbpface.numnodes
      uvol[1,:,bndry.element]
      for j = 1:sbpface.stencilsize
        uvol[1, sbpface.perm[j,bndry.face], bndry.element] += uface[1,i,bindex]*
                                                             sbpface.interp[j,i]
      end  # end for j = 1:sbpface.stencilsize
    end # for i = 1:sbpface.numnodes
  end # end enumerate

  return nothing
end