# Adjoint for Euler Equations

push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/NonlinearSolvers"))
using NonlinearSolvers

export calcAdjoint

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
	                sbp::AbstractSBP, eqn::EulerData{Tsol, Tres, Tdim}, opts,
                  functionalData::AbstractOptimizationData,
                  adjoint_vec::Array{Tsol,1}; functional_number::Int=1)
                  #functor, functional_number, adjoint_vec::Array{Tsol, 1})

  # Get information corresponding to functional
  functional_edges = []
  if functionalData.is_objective_fn == true
    functional_edges = opts["geom_faces_objective"]
    functional_name = FunctionalDict[opts["objective_function"]]
  else
    key = string("geom_edges_functional", functional_number)
    functional_name = getFunctionalName(opts, functional_number)
    functional_edges = opts[key]
  end

  # Check if PETSc is initialized
  if PetscInitialized() == 0 # PETSc Not initialized before
    PetscInitialize(["-malloc", "-malloc_debug", "-ksp_monitor",  "-pc_type", "bjacobi", "-sub_pc_type", "ilu", "-sub_pc_factor_levels", "4", "ksp_gmres_modifiedgramschmidt", "-ksp_pc_side", "right", "-ksp_gmres_restart", "30" ])
  end
  # Calculate the Jacobian of the residual

  # res_jac = zeros(Tres, mesh.numDof, mesh.numDof)
  # pert = complex(0, opts["epsilon"])
  # NonlinearSolvers.calcJacobianComplex(newton_data, mesh, sbp, eqn, opts,
  #                                      evalEuler, pert, res_jac)
  res_jac, jacData = calcResidualJacobian(mesh, sbp, eqn, opts)
  println("typeof res_jac = ", typeof(res_jac))
  # Re-interpolate interior q to q_bndry. This is done because the above step
  # pollutes the existing eqn.q_bndry with complex values.
  boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, eqn.q, eqn.q_bndry)

  # calculate the derivative of the function w.r.t q_vec
  func_deriv = zeros(Tsol, mesh.numDof)

  # 3D array into which func_deriv_arr gets interpolated
  func_deriv_arr = zeros(eqn.q)


  # Calculate df/dq_bndry on edges where the functional is calculated and put
  # it back in func_deriv_arr
  calcFunctionalDeriv(mesh, sbp, eqn, opts, functional_name, functional_edges,
                      functionalData, func_deriv_arr)  # populate df_dq_bndry

  # Assemble func_deriv
  assembleArray(mesh, sbp, eqn, opts, func_deriv_arr, func_deriv)

  # Solve for adjoint vector. This depends on whether PETSc is used or not.

  # TODO: The following operation creates a temporary copy of adjoint_vec, does
  #       the '\' computation and then puts it back into adjoint_vec. This
  #       needs to change.

  if opts["jac_type"] == 1 || opts["jac_type"] == 2
    adjoint_vec[:] = (res_jac.')\func_deriv # There is no negative sign because
                                            # the weak residual is computed on
                                            # the right hand side
  elseif opts["jac_type"] == 3
    b = PetscVec(eqn.comm)
    PetscVecSetType(b, VECMPI)
    PetscVecSetSizes(b, PetscInt(mesh.numDof), PETSC_DECIDE)
    x = PetscVec(eqn.comm)
    PetscVecSetType(x, VECMPI)
    PetscVecSetSizes(x, PetscInt(mesh.numDof), PETSC_DECIDE)
    ksp = KSP(eqn.comm)
    KSPSetFromOptions(ksp)
    KSPSetOperators(ksp, res_jac, res_jac)  # this was A, Ap
    println("Before NonlinearSolvers.petscSolve")
    NonlinearSolvers.petscSolve(jacData, res_jac, res_jac, x, b, ksp, opts, func_deriv, adjoint_vec)
  end # End how to solve for adjoint_vec

  outname = string("adjoint_vec.dat")
  f = open(outname, "w")
  for i = 1:length(adjoint_vec)
    println(f, real(adjoint_vec[i]))
  end
  close(f)

  return nothing
end

@doc """
###EulerEquationMod.calcResidualJacobian

The function calculates the residual for computing the adjoint vector. The
function allows for jacobian to be computed depending on the jacobian type
specified in the options dictionary `jac_type`.

**Input**

* `mesh` : Abstract mesh object
* `sbp`  : Summation-By-parts operator
* `eqn`  : Euler equation object
* `opts` : options dictionary

**Output**

* `jac` : Jacobian matrix

"""->

function calcResidualJacobian{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh},
         sbp::AbstractSBP, eqn::EulerData{Tsol, Tres, Tdim}, opts)

  jac_type = opts["jac_type"]
  if jac_type == 4 # For now. Date: 28/11/2016
    error("jac_type = 4 not yet supported")
  end
  jacData = NonlinearSolvers.NewtonData{Tsol, Tres}(mesh, sbp, eqn, opts)

  # Initialize shape of Jacobian Matrix
  if jac_type == 1
    jac = zeros(Tres, mesh.numDof, mesh.numDof)
  elseif jac_type == 2
    if typeof(mesh) <: AbstractCGMesh
      jac = SparseMatrixCSC(mesh.sparsity_bnds, Tres)
    else
      jac = SparseMatrixCSC(mesh, Tres)
    end
  elseif jac_type == 3
    obj_size = PetscInt(mesh.numDof)
    jac = PetscMat(eqn.comm)
    PetscMatSetFromOptions(jac)
    PetscMatSetType(jac, PETSc.MATMPIAIJ)
    PetscMatSetSizes(jac, obj_size, obj_size, PETSC_DECIDE, PETSC_DECIDE)
    if mesh.isDG
      MatSetOption(jac, PETSc.MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_TRUE)
    end
    # Preallocate matrix jac
    dnnz = zeros(PetscInt, mesh.numDof)  # diagonal non zeros per row
    onnz = zeros(PetscInt, mesh.numDof)
    for i = 1:mesh.numDof
      dnnz[i] = mesh.sparsity_counts[1, i]
      onnz[i] = mesh.sparsity_counts[2, i]
    end
    PetscMatMPIAIJSetPreallocation(jac, PetscInt(0),  dnnz, PetscInt(0), onnz)
    MatSetOption(jac, PETSc.MAT_ROW_ORIENTED, PETSC_FALSE)
    PetscMatZeroEntries(jac)
    matinfo = PetscMatGetInfo(jac, Int32(1))
  end

  # Now function call for calculating Jacobian
  ctx_residual = (evalEuler,)
  NonlinearSolvers.physicsJac(jacData, mesh, sbp, eqn, opts, jac, ctx_residual)

  return jac, jacData
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

function calcFunctionalDeriv{Tmsh, Tsol}(mesh::AbstractDGMesh{Tmsh}, sbp::AbstractSBP,
	                         eqn::EulerData{Tsol}, opts, functor, functional_edges,
	                         functionalData, func_deriv_arr)

  integrand = zeros(eqn.q_bndry)

  # Populate integrand
  for itr = 1:length(functional_edges)
    g_edge_number = functional_edges[itr] # Extract geometric edge number
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
        aux_vars = sview(eqn.aux_vars_bndry, :, j, global_facenum)
        x = sview(mesh.coords_bndry, :, j, global_facenum)
        dxidx = sview(mesh.dxidx_bndry, :, :, j, global_facenum)
        nrm = sview(sbp.facenormal, :, bndry_i.face)
        nx = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
        ny = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]
        node_info = Int[itr,j,i]
        integrand_i = sview(integrand, :, j, global_facenum)

        calcIntegrandDeriv(opts, eqn.params, q2, aux_vars, [nx, ny], integrand_i,
                           node_info, functor, functionalData)
      end  # End for j = 1:mesh.sbpface.numnodes
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

function calcIntegrandDeriv{Tsol, Tres, Tmsh}(opts, params, q::AbstractArray{Tsol,1},
	                        aux_vars::AbstractArray{Tres, 1}, nrm::AbstractArray{Tmsh},
	                        integrand_deriv::AbstractArray{Tsol, 1}, node_info,
                          functor, functionalData)


  pert = complex(0, opts["epsilon"])

  for i = 1:length(q)
    q[i] += pert
    val = functor(params, q, aux_vars, nrm, node_info, functionalData)
    integrand_deriv[i] = imag(val)/norm(pert)
    q[i] -= pert
  end

  return nothing
end  # End function calcIntegrandDeriv
