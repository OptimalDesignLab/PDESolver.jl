# Adjoint for Euler Equations

@doc """
### EulerEquationMod.calcAdjoint

Calculates the adjoint vector for a single functional

**Inputs**

*  `mesh` : Abstract DG mesh type
*  `sbp`  : Summation-By-parts operator
*  `eqn`  : Euler equation object
*  `functor` : functional to be evaluated
*  `functional_number` : Numerical identifier to obtain geometric edges on
                         which a functional acts
*  `adjoint_vec` : Resulting adjoint vector

**Outputs**

*  None

"""->

function calcAdjoint{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractDGMesh{Tmsh},
                                             sbp::AbstractSBP, 
                                             eqn::EulerData{Tsol, Tres, Tdim}, 
                                             opts,
                                             functor, 
                                             functional_number, 
                                             adjoint_vec::Array{Tsol, 1})

  # Get information corresponding to functional_number
  key = string("geom_edges_functional", functional_number)
  functional_edges = opts[key]

  # Calculate the Jacobian of the residual
  newton_data = NonlinearSolvers.NewtonData(mesh, sbp, eqn, opts)
  res_jac = zeros(Tres, mesh.numDof, mesh.numDof)
  pert = complex(0, opts["epsilon"])
  NonlinearSolvers.calcJacobianComplex(newton_data, mesh, sbp, eqn, opts,
                                       evalResidual, pert, res_jac)
  
  # Re-interpolate interior q to q_bndry. This is done because the above step
  # pollutes the existing eqn.q_bndry with complex values.
  boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, eqn.q, eqn.q_bndry)

  # calculate the derivative of the function w.r.t q_vec
  func_deriv = zeros(Tsol, mesh.numDof)
  
  # 3D array into which func_deriv_arr gets interpolated
  func_deriv_arr = zeros(eqn.q) 
  

  # Calculate df/dq_bndry on edges where the functional is calculated and put 
  # it back in func_deriv_arr
  calcFunctionalDeriv_adjoint_consistent(mesh, sbp, eqn, opts, functor, 
                                         functional_edges, func_deriv_arr)  # populate df_dq_bndry

  # Assemble func_deriv
  assembleArray(mesh, sbp, eqn, opts, func_deriv_arr, func_deriv)

  # Solve for adjoint vector

  # TODO: The following operation creates a temporary copy of adjoint_vec, does
  #       the '\' computation and then puts it back into adjoint_vec. This
  #       needs to change.
  adjoint_vec[:] = (res_jac.')\func_deriv # There is no negative sign because 
                                          # the weak residual is computed on
                                          # the right hand side

  return nothing
end

@doc """
### EulerEquationMod. calcFunctionalDeriv

Computes a 3D array of the derivative of a functional w.r.t eqn.q on all
mesh nodes.

**Inputs**

*  `mesh` : Abstract DG mesh type
*  `sbp`  : Summation-By-parts operator
*  `eqn`  : Euler equation object
*  `opts` : Options dictionary
*  `functor` : Functional name which is to be evaluated
*  `functional_edges` : Numerical identifier to obtain geometric edges on
                         which a functional acts

**Outputs**

*  None

"""->

function calcFunctionalDeriv{Tmsh, Tsol}(mesh::AbstractDGMesh{Tmsh}, 
                                         sbp::AbstractSBP,
                                         eqn::EulerData{Tsol}, 
                                         opts, 
                                         functor, 
                                         functional_edges,
                                         func_deriv_arr)

  integrand = zeros(eqn.q_bndry)

  # Populate integrand
  for itr = 1:length(functional_edges)
    g_edge_number = functional_edges[itr] # Extract geometric edge number
    start_index = mesh.bndry_offsets[g_edge_number]
    end_index = mesh.bndry_offsets[g_edge_number+1]
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
        aux_vars = sview(eqn.aux_vars_bndry, :, j, global_facenum)
        x = sview(mesh.coords_bndry, :, j, global_facenum)
        dxidx = sview(mesh.dxidx_bndry, :, :, j, global_facenum)
        nrm = sview(sbp.facenormal, :, bndry_i.face)
        nx = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
        ny = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]
        integrand_i = sview(integrand, :, j, global_facenum)

        calcIntegrandDeriv(opts, eqn.params, q, aux_vars, [nx, ny], integrand_i, functor)
      end  # End for j = 1:mesh.sbpface.numnodes
    end    # End for i = 1:nfaces
  end      # End for itr = 1:length(functional_edges)

  boundaryintegrate!(mesh.sbpface, mesh.bndryfaces, integrand, func_deriv_arr)

  return nothing
end  # End function calcFunctionalDeriv

function calcFunctionalDeriv_adjoint_consistency{Tmsh, Tsol}(mesh::AbstractDGMesh{Tmsh}, 
                                                             sbp::AbstractSBP,
                                                             eqn::EulerData{Tsol}, 
                                                             opts, 
                                                             functor, 
                                                             functional_edges,
                                                             func_deriv_arr)

  if haskey(opts, "epsilon")
    pert = complex(0, opts["epsilon"])
  else
    pert = complex(0, 1.0e-12)
  end

  # in order to use integratefunctional!, constract array with only one element
  faces_integrand = Array(Tsol, 1, mesh.numNodesPerFace, 1)
  faces = Array(Boundary, 1)

  # Populate integrand
  for itr = 1:length(functional_edges)
    g_edge_number = functional_edges[itr] # Extract geometric edge number
    start_index = mesh.bndry_offsets[g_edge_number]
    end_index = mesh.bndry_offsets[g_edge_number+1]
    idx_range = start_index:(end_index-1)
    bndry_facenums = sview(mesh.bndryfaces, idx_range) # faces on geometric edge i

    nfaces = length(bndry_facenums)
    q2 = zeros(Tsol, mesh.numDofPerNode)

    for i = 1:nfaces
      bndry_i = bndry_facenums[i]
      elem = bndry_i.element
      global_facenum = idx_range[i]
      q_elem = sview(eqn.q, :,:,elem)

      faces[1] = bndry_i

      for j = 1:mesh.numNodesPerElement
        for dof = 1 : mesh.numDofPerNode
          # perturb
          q_elem[dof, j] += pert
          q_face = slice(eqn.q_bndry[:,:,global_facenum])
          boundaryinterpolate(mesh.sbp, bndry_i, q_elem, q_face)
          q_face_bak = q_face[:]

          # Start to compute contribution to functional from this boundary face
          # first call functor to compute integrand
          integrand = slice(faces_integrand, 1, :, 1)
          functor(mesh, sbp,eqn, opts, global_facenum, integrand)
          face_integral = zeros(Tsol, 1)
          
          # integrate
          integratefunctional!(mesh.sbpface, faces, faces_integrand, face_integral)

          # complex finite difference
          func_deriv_arr[dof, i, elem] = imag(face_integral[1]/norm(pert))

          # un-perturb
          q_elem[dof, j] -= pert
          q_face = q_face_bak[:]
        end
      end
    end    # End for i = 1:nfaces
  end      # End for itr = 1:length(functional_edges)

  boundaryintegrate!(mesh.sbpface, mesh.bndryfaces, integrand, func_deriv_arr)

  return nothing
end  # End function calcFunctionalDeriv

@doc """
### EulerEquationMod.calcIntegrandDeriv

Compute the derivative of the functional Integrand at a node w.r.t all the 
degrees of freedom at the node.

**Inputs**

*  `opts`   : Options dictionary
*  `params` : parameter type
*  `q`      : Solution variable at a node
*  `aux_vars` : Auxiliary variables
*  `nrm`    : normal vector in the physical space
*  `integrand_deriv` : Derivative of the integrand at that particular node
*  `functor`: Functional that is to be evaluated

**Outputs**

*  None

"""->

function calcIntegrandDeriv{Tsol, Tres, Tmsh}(mesh::AbstractDGMesh{Tmsh},
                                              sbp::AbstractSBP,
                                              eqn::EulerData{Tsol},
                                              opts, 
                                              bndry_indx,
                                              functor,
	                                            integrand_deriv::AbstractArray{Tsol, 2})


  pert = complex(0, opts["epsilon"])

  num_facenodes = mesh.numNodesPerFace
  for n = 1 : num_facenodes
    q = sview(eqn.q_bndry, :, n, bndry_indx)
    for i = 1:length(q)
      q[i] += pert
      functor(mesh, sbp, eqn, opts, bndry_indx, face_integrand)
      val = functor(params, q, aux_vars, nrm)
      integrand_deriv[i,n] = imag(face_integrand[n])/norm(pert)
      q[i] -= pert
    end
  end

  return nothing
end  # End function calcIntegrandDeriv
