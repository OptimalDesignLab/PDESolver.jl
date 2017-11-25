# Adjoint for Euler Equations

@doc """
### EulerEquationMod.calcAdjoint

Calculates the adjoint vector, ψ, for a single functional. The user must always
call this function in order to compute the adjoint vector. Currently only DG meshes
are supported. The function performs a direct solve using Julia's  `\\` operator.
For parallel meshes, a PETSc solve is done using ILU factorization.

**Inputs**

*  `mesh` : Abstract DG mesh type
*  `sbp`  : Summation-By-parts operator
*  `eqn`  : Euler equation object
*  `opts` : Options dictionary
*  `functionalData` : Object corresponding the boundary functional being
                      computed. It must be a subtype of `AbstractOptimizationData`
*  `adjoint_vec` : Resulting adjoint vector. In the parallel case, the adjoint
                   vector is distributed over the processors similar to
                   eqn.q_vec i.e. every rank has its
                   share of the adjoint vector corresponding to the dofs on the
                   rank.

*  `functional_number` : Numerical identifier to obtain geometric edges on
                         which a functional acts

**Outputs**

*  None

"""->

function calcAdjoint{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractDGMesh{Tmsh},
                  sbp::AbstractSBP, eqn::EulerData{Tsol, Tres, Tdim}, opts,
                  functionalData::AbstractOptimizationData,
                  adjoint_vec::Array{Tsol,1}; functional_number::Int=1)
                  #functor, functional_number, adjoint_vec::Array{Tsol, 1})
#=
  # Check if PETSc is initialized
  if PetscInitialized() == 0 # PETSc Not initialized before
    PetscInitialize(["-malloc", "-malloc_debug", "-ksp_monitor",  "-pc_type",
      "bjacobi", "-sub_pc_type", "ilu", "-sub_pc_factor_levels", "4",
      "ksp_gmres_modifiedgramschmidt", "-ksp_pc_side", "right",
      "-ksp_gmres_restart", "30" ])
  end
=#
  println("computing adjoint")
  println("vecnorm(eqn.q) = ", vecnorm(eqn.q))
  # This doesn't work in parallel?  paralle_type == 2 if calculating the
  # jacobian
  if opts["parallel_type"] == 1
    startSolutionExchange(mesh, sbp, eqn, opts, wait=true)
  end

  # Allocate space for adjoint solve
  pc, lo = NonlinearSolvers.getNewtonPCandLO(mesh, sbp, eqn, opts)
  ls = StandardLinearSolver(pc, lo, eqn.comm)
  ctx_residual = (evalResidual,)
  calcPCandLO(ls, mesh, sbp, eqn, opts, ctx_residual, 0.0)

#=
  # Allocate space for adjoint solve
  jacData, res_jac, rhs_vec = NonlinearSolvers.setupNewton(mesh, mesh, sbp, eqn,
                                             opts, alloc_rhs=true)

  # Get the residual jacobian
  ctx_residual = (evalResidual,)
  NonlinearSolvers.physicsJac(jacData, mesh, sbp, eqn, opts, res_jac, ctx_residual)
=#
  # Re-interpolate interior q to q_bndry. This is done because the above step
  # pollutes the existing eqn.q_bndry with complex values.
  boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, eqn.q, eqn.q_bndry)
  # calculate the derivative of the function w.r.t q_vec
  func_deriv = zeros(Tsol, mesh.numDof)

  # 3D array into which func_deriv_arr gets interpolated
  func_deriv_arr = zeros(eqn.q)

  # Calculate df/dq_bndry on edges where the functional is calculated and put
  # it back in func_deriv_arr
  calcFunctionalDeriv(mesh, sbp, eqn, opts, functionalData, func_deriv_arr)

  # Assemble func_deriv
  assembleSolution(mesh, sbp, eqn, opts, func_deriv_arr, func_deriv)
  func_deriv[:] = -func_deriv[:]

  # do transpose solve

  _adjoint_vec = zeros(real(Tsol), length(adjoint_vec))
 
  println("vecnorm(func_deriv) = ", vecnorm(real(func_deriv)))

  # debugging
  lo2 = getBaseLO(lo)
  jacT = lo2.A.'
  _adjoint_vec = jacT\real(func_deriv)
#  linearSolveTranspose(ls, real(func_deriv), _adjoint_vec)
  copy!(adjoint_vec, _adjoint_vec)

  println("norm(adjoint_vec) = ", norm(adjoint_vec))
#=
  # Solve for adjoint vector. residual jacobian needs to be transposed first.
  jac_type = typeof(res_jac)
  if jac_type <: Array || jac_type <: SparseMatrixCSC
    res_jac = res_jac.'
  elseif  jac_type <: PetscMat
    PetscMatAssemblyBegin(res_jac) # Assemble residual jacobian
    PetscMatAssemblyEnd(res_jac)
    res_jac = MatTranspose(res_jac, inplace=true)
  else
    error("Unsupported jacobian type")
  end
  step_norm = NonlinearSolvers.matrixSolve(jacData, eqn, mesh, opts, res_jac,
                                         adjoint_vec, real(func_deriv), BSTDOUT)

=#
  # Output/Visualization options for Adjoint
  if opts["write_adjoint"]
    outname = string("adjoint_vec_", mesh.myrank,".dat")
    f = open(outname, "w")
    for i = 1:length(adjoint_vec)
      println(f, real(adjoint_vec[i]))
    end
    close(f)
  end
  
  if opts["write_adjoint_vis"]
    saveSolutionToMesh(mesh, adjoint_vec)
    fname = "adjoint_field"
    writeVisFiles(mesh, fname)
  end

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
* `opts` : Options dictionary

**Output**

* `jac` : Jacobian matrix

"""->

function calcResidualJacobian{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh},
         sbp::AbstractSBP, eqn::EulerData{Tsol, Tres, Tdim}, opts)

  #TODO: get rid of this function, use NonlinearSolvers.physicsJac instead
  jac_type = opts["jac_type"]
  Tjac = typeof(real(eqn.res_vec[1]))  # type of jacobian, residual
  if jac_type == 4 # For now. Date: 28/11/2016
    error("jac_type = 4 not yet supported")
  end
  jacData = NonlinearSolvers.NewtonData{Tsol, Tres}(mesh, sbp, eqn, opts)

  # Initialize shape of Jacobian Matrix
  if jac_type == 1
    jac = zeros(Tjac, mesh.numDof, mesh.numDof)
  elseif jac_type == 2
    if typeof(mesh) <: AbstractCGMesh
      jac = SparseMatrixCSC(mesh.sparsity_bnds, Tjac)
    else
      jac = SparseMatrixCSC(mesh, Tjac)
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
  ctx_residual = (evalResidual,)
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
*  `functionalData` : Functional object of super-type AbstractOptimizationData
                      that is needed for computing the adjoint vector.
                      Depending on the functional being computed, a different
                      method based on functional type may be needed to be
                      defined.
*  `func_deriv_arr` : 3D array that stores the derivative of the functional
                      w.r.t. eqn.q. The array is the same size as eqn.q

**Outputs**

*  None

"""->

function calcFunctionalDeriv{Tmsh, Tsol}(mesh::AbstractDGMesh{Tmsh}, sbp::AbstractSBP,
                           eqn::EulerData{Tsol}, opts,
                           functionalData::AbstractIntegralOptimizationData, func_deriv_arr)

  integrand = zeros(eqn.q_bndry)
  functional_edges = functionalData.geom_faces_functional

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
#    q2 = zeros(Tsol, mesh.numDofPerNode)
    for i = 1:nfaces
      bndry_i = bndry_facenums[i]
      global_facenum = idx_range[i]
      for j = 1:mesh.sbpface.numnodes
        vtx_arr = mesh.topo.face_verts[:,bndry_i.face]
        q = sview(eqn.q_bndry, :, j, global_facenum)
#        convertToConservative(eqn.params, q, q2)
        aux_vars = ro_sview(eqn.aux_vars_bndry, :, j, global_facenum)
        x = ro_sview(mesh.coords_bndry, :, j, global_facenum)
        nrm = ro_sview(mesh.nrm_bndry, :, j, global_facenum)
        node_info = Int[itr,j,i]
        integrand_i = sview(integrand, :, j, global_facenum)

        calcIntegrandDeriv(opts, eqn.params, q, aux_vars, nrm, integrand_i, node_info,
                           functionalData)
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
*  `node_info` : Tuple containing information about the node
*  `functionalData` : Functional object that is a subtype of AbstractOptimizationData.

**Outputs**

*  None

"""->

function calcIntegrandDeriv{Tsol, Tres, Tmsh}(opts, params::ParamType{2},
                            q::AbstractArray{Tsol,1},
                            aux_vars::AbstractArray{Tres, 1}, nrm::AbstractArray{Tmsh},
                            integrand_deriv::AbstractArray{Tsol, 1}, node_info,
                            functionalData::BoundaryForceData{Tsol,:lift})

  pert = complex(0, 1e-20)
  aoa = params.aoa
  momentum = zeros(Tsol,2)

  for i = 1:length(q)
    q[i] += pert
    calcBoundaryFunctionalIntegrand(params, q, aux_vars, nrm, node_info, functionalData, momentum)
    val = -momentum[1]*sin(aoa) + momentum[2]*cos(aoa)
    integrand_deriv[i] = imag(val)/norm(pert)
    q[i] -= pert
  end # End for i = 1:length(q)

  return nothing
end

function calcIntegrandDeriv{Tsol, Tres, Tmsh}(opts, params::ParamType{2},
                            q::AbstractArray{Tsol,1},
                            aux_vars::AbstractArray{Tres, 1}, nrm::AbstractArray{Tmsh},
                            integrand_deriv::AbstractArray{Tsol, 1}, node_info,
                            functionalData::BoundaryForceData{Tsol,:drag})

  pert = complex(0, 1e-20)
  aoa = params.aoa
  momentum = zeros(Tsol,2)

  for i = 1:length(q)
    q[i] += pert
    calcBoundaryFunctionalIntegrand(params, q, aux_vars, nrm, node_info, functionalData, momentum)
    val = momentum[1]*cos(aoa) + momentum[2]*sin(aoa)
    integrand_deriv[i] = imag(val)/norm(pert)
    q[i] -= pert
  end # End for i = 1:length(q)

  return nothing
end

function calcIntegrandDeriv{Tsol, Tres, Tmsh}(opts, params::ParamType{2},
                            q::AbstractArray{Tsol,1},
                            aux_vars::AbstractArray{Tres, 1},
                            nrm::AbstractArray{Tmsh},
                            integrand_deriv::AbstractArray{Tsol, 1}, node_info,
                            functionalData::MassFlowData)

  node_info = [1, 2, 3]
  h = 1e-20
  pert = Complex128(0, h)
  val = zeros(Complex128, 1)
  for i=1:length(q)
    q[i] += pert
    calcBoundaryFunctionalIntegrand(params, q, aux_vars, nrm, node_info,
                                          functionalData, val)
    integrand_deriv[i] = imag(val[1])/h
    q[i] -= pert
  end
  # functional integrand rho*v
#  integrand_deriv[2] = 1*nrm[1]
#  integrand_deriv[3] = 1*nrm[2]

  return nothing
end
  


#=
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
=#
