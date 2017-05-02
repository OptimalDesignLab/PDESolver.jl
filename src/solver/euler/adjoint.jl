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

  
  # Re-interpolate interior q to q_bndry. This is done because the above step
  # pollutes the existing eqn.q_bndry with complex values.
  boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, eqn.q, eqn.q_bndry)

  # calculate the derivative of the function w.r.t q_vec
  func_deriv = zeros(Tsol, mesh.numDof)
  
  # 3D array into which dJ_dq gets interpolated
  dJ_dq = zeros(eqn.q) 
  

  # Calculate functional_deric_arr = dJ/dq 
  calcFunctionalDeriv(mesh, sbp, eqn, opts, functor, 
                                         functional_edges, dJ_dq)  

  # Assemble func_deriv
  assembleArray(mesh, sbp, eqn, opts, dJ_dq, func_deriv)

  # Solve for adjoint vector

  newton_data = NonlinearSolvers.NewtonData{Tsol, Tres}(mesh, sbp, eqn, opt)
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
  elseif jac_type == 3 || jac_type == 4 # petsc
    jac, jacp, x, b, ksp = createPetscData(mesh, pmesh, sbp, eqn, opts, newton_data, f)
    ctx_newton = (jacp, x, b, ksp)
  end

  func = evalResidual
  res_dummy = Array(Float64, 0, 0, 0)  # not used, so don't allocation memory
  calcJacobianSparse(newton_data, mesh, sbp, eqn, opts, func, res_dummy, 
                     pert, jac, 0.0)
  
  # solve the linear system
  lambda = zeros(Tjac, m)  # newton update
  if jac_type == 2
    transpose(jac)
    jac_f = factorize(jac)
    adjoint_vec[:] = jac_f\(func_deriv)  #  calculate Newton update
  else
    println("only jac_type = 2 supported in adjoint solver")
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
                                         functor, 
                                         functional_edges,
                                         dJ_dq)

  if haskey(opts, "epsilon")
    pert = complex(0, opts["epsilon"])
  else
    pert = complex(0, 1.0e-12)
  end

  faces_integrand = Array(Tsol, 1, mesh.numNodesPerFace, 1)
  # in order to use integratefunctional!, constract array with only one element
  faces = Array(Boundary, 1)

  # Populate integrand
  for itr = 1:length(functional_edges)
    g_edge_number = functional_edges[itr] # Extract geometric edge number
    start_index = mesh.bndry_offsets[g_edge_number]
    end_index = mesh.bndry_offsets[g_edge_number+1]
    idx_range = start_index:(end_index-1)
    bndry_facenums = sview(mesh.bndryfaces, idx_range) # faces on geometric edge i

    nfaces = length(bndry_facenums)

    for i = 1:nfaces
      bndry_i = bndry_facenums[i]
      elem = bndry_i.element
      face_gid = idx_range[i]
      q_elem = sview(eqn.q, :,:,elem)

      faces[1] = bndry_i

      for j = 1:mesh.numNodesPerElement
        for dof = 1 : mesh.numDofPerNode
          # perturbation
          q_elem[dof, j] += pert
          # q_face = slice(eqn.q_bndry[:, :, face_gid])
          q_face = sview(eqn.q_bndry[:, :, face_gid])
          boundaryinterpolate(mesh.sbpface, bndry_i, q_elem, q_face)
          q_face_bak = q_face[:]

          integrand = slice(faces_integrand, 1, :, 1)

          functor(mesh, sbp,eqn, opts, face_gid, integrand)
          face_integral = zeros(Tsol, 1)
          integratefunctional!(mesh.sbpface, faces, faces_integrand, face_integral)
          dJ_dq[dof, j, elem] = imag(face_integral)/norm(pert)
          # undo perturbation
          q_elem[dof, j] -= pert
          q_face = q_face_bak[:]
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

function transposeJac(jac)
  for col = 1 : jac.m
    idx0 = jac.colptr[col]
    idx1 = jac.colptr[col + 1] - 1
    for idx = idx0 : idx1
      row = jac.rowind[idx] 
      if row <= col_0
        continue
      end
      j0 = jac.colptr[row_0]
      j1 = jac.colptr[row_0 + 1] - 1
      for j = j0 : j1
        row1 = jac.rowind[j]
        if row1 != col
          continue
        end
        tmp = jac.val[idx]
        jac.val[idx] = jac.val[j]
        jac.val[j] = tmp;
      end
    end
  end
  return nothing
end
