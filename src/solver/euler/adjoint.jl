# Adjoint for Euler Equations

# @doc """
# ### EulerEquationMod.calcAdjoint

# Calculates the adjoint vector for a single functional

# **Inputs**

# *  `mesh` : Abstract DG mesh type
# *  `sbp`  : Summation-By-parts operator
# *  `eqn`  : Euler equation object
# *  `functor` : functional to be evaluated
# *  `functional_number` : Numerical identifier to obtain geometric edges on
                         # which a functional acts
# *  `adjoint_vec` : Resulting adjoint vector

# **Outputs**

# *  None

# """->
function calcAdjoint{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractDGMesh{Tmsh},
                                             sbp::AbstractSBP, 
                                             eqn::EulerData{Tsol, Tres, Tdim}, 
                                             opts,
                                             # functional_number, 
                                             adjoint_vec::Array{Tsol, 2})

  # Get information corresponding to functional_number
  key = string("geom_edges_functional")
  functional_edges = opts[key]

  for itr = 1:length(functional_edges)
    g_edge_number = functional_edges[itr] # Extract geometric edge number
    start_index = mesh.bndry_offsets[g_edge_number]
    end_index = mesh.bndry_offsets[g_edge_number+1]
    idx_range = start_index:(end_index-1)
  end
  
  # Re-interpolate interior q to q_bndry. This is done because the above step
  # pollutes the existing eqn.q_bndry with complex values.
  boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, eqn.q, eqn.q_bndry)

  newton_data = NonlinearSolvers.NewtonData{Tsol, Tres}(mesh, sbp, eqn, opts)
  jac_type = opts["jac_type"]
  Tjac = typeof(real(eqn.res_vec[1]))  # type of jacobian, residual
  m = mesh.numDof

  if jac_type == 2  # sparse
    if typeof(mesh) <: AbstractCGMesh
      println("creating CG SparseMatrix")
      jac = SparseMatrixCSC(mesh.sparsity_bnds, Tjac)
    else
      println("Creating DG sparse matrix")
      jac = SparseMatrixCSC(mesh, Tjac)
    end
    ctx_newton = ()
  # elseif jac_type == 3 || jac_type == 4 # petsc
    # jac, jacp, x, b, ksp = createPetscData(mesh, pmesh, sbp, eqn, opts, newton_data, f)
    # ctx_newton = (jacp, x, b, ksp)
  else
    error("only jac_type = 2 supported in adjoint solver")
  end

  func = evalResidual
  res_dummy = Array(Float64, 0, 0, 0)  # not used, so don't allocation memory
  if haskey(opts, "epsilon")
    pert = complex(0, opts["epsilon"])
  else
    pert = complex(0, 1.0e-16)
  end
  NonlinearSolvers.calcJacobianSparse(newton_data, mesh, sbp, eqn, opts, func, 
                                      res_dummy, pert, jac, 0.0)
  
  # calculate the derivative of the function w.r.t q_vec
  func_deriv = zeros(Float64, mesh.numDof)
  
  # 3D array into which dJ_dq gets interpolated
  dCl_dq = zeros(Float64, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl) 
  dCd_dq = zeros(Float64, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl) 
  # Calculate functional_deric_arr = dJ/dq 
  calcFunctionalDeriv(mesh, sbp, eqn, opts, functional_edges, dCl_dq, dCd_dq)
  
  # solve the linear system
  lambda = zeros(Tjac, m)  # newton update
  if jac_type == 2
    my_transposeJac(jac)
    jac_f = factorize(jac)

    assembleSolution(mesh, sbp, eqn, opts, dCl_dq, func_deriv)
    adjoint_vec[:,1] = jac_f\(func_deriv) 

    assembleSolution(mesh, sbp, eqn, opts, dCd_dq, func_deriv)
    adjoint_vec[:,2] = jac_f\(func_deriv) 
  else
    error("only jac_type = 2 supported in adjoint solver")
  end

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
                                         functional_edges,
                                         dCl_dq,
                                         dCd_dq)

  if haskey(opts, "epsilon")
    pert = complex(0, opts["epsilon"])
  else
    pert = complex(0, 1.0e-16)
  end

  faces_integrand = Array(Tsol, 2, mesh.numNodesPerFace, 1)
  # in order to use integratefunctional!, constract array with only one element
  faces = Array(Boundary, 1)

  # q_face = Array(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
  # qf_bak = Array(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
  # Populate integrand
  for itr = 1:length(functional_edges)
    g_edge_number = functional_edges[itr] # Extract geometric edge number
    start_index = mesh.bndry_offsets[g_edge_number]
    end_index = mesh.bndry_offsets[g_edge_number+1]
    idx_range = start_index:(end_index-1)
    bndry_facenums = sview(mesh.bndryfaces, idx_range) # faces on geometric edge i

    nfaces = length(bndry_facenums)
    
    if eqn.params.isViscous
      functor = calcForceCoef_viscous
    else
      functor = calcForceCoef_inviscid
    end
    for i = 1:nfaces
      bndry_i = bndry_facenums[i]
      elem = bndry_i.element
      face_id = idx_range[i]
      q_elem = sview(eqn.q, :,:,elem)

      faces[1] = bndry_i

      for j = 1:mesh.numNodesPerElement
        for dof = 1 : mesh.numDofPerNode
          # perturbation
          q_elem[dof, j] += pert
          # q_face = slice(eqn.q_bndry[:, :, face_id])
          # q_face = eqn.q_bndry[:, :, face_id]
          # qf_bak[:,:] = q_face[:,:]

          # boundaryinterpolate(mesh.sbpface, bndry_i, q_elem, q_face)
          # eqn.q_bndry[:,:,face_id] = q_face[:,:]

          q_face = sview(eqn.q_bndry, :, :, face_id)
          qf_bak = q_face[:,:]
          boundaryinterpolate(mesh.sbpface, bndry_i, q_elem, q_face)
          # println(eqn.q_bndry[:,:,face_id])

          integrand = slice(faces_integrand, :, :, 1)

          functor(mesh, sbp, eqn, opts, face_id, integrand)
          face_integral = zeros(Tsol, 2)
          integratefunctional!(mesh.sbpface, faces, faces_integrand, face_integral)
          dCl_dq[dof, j, elem] = imag(face_integral[1])/imag(pert)
          dCd_dq[dof, j, elem] = imag(face_integral[2])/imag(pert)

          # undo perturbation
          q_elem[dof, j] -= pert
          eqn.q_bndry[:,:,face_id] = qf_bak[:,:]
        end
      end
    end    # End for i = 1:nfaces
  end      # End for itr = 1:length(functional_edges)

  return nothing
end  # End function calcFunctionalDeriv


@doc """
# Transpose jac. Given that jac is symmetric, 
# we don't need additional intermediate variables
"""->

function my_transposeJac(jac)
  for col = 1 : jac.m
    idx0 = jac.colptr[col]
    idx1 = jac.colptr[col + 1] - 1
    for idx = idx0 : idx1
      row = jac.rowval[idx] 
      if row <= col
        continue
      end
      j0 = jac.colptr[row]
      j1 = jac.colptr[row + 1] - 1
      for j = j0 : j1
        row1 = jac.rowval[j]
        if row1 != col
          continue
        end
        tmp = jac.nzval[idx]
        jac.nzval[idx] = jac.nzval[j]
        jac.nzval[j] = tmp;
      end
    end
  end
  return nothing
end
