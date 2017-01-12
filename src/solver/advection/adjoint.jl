# Adjoint Computation
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/NonlinearSolvers"))
using NonlinearSolvers

export calcAdjoint

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
*  `adjoint_vec` : Adjoint vector corresponding to the particular functional

**Outputs**

*  None

"""->
#=
function calcAdjoint{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractCGMesh{Tmsh},
                    sbp::AbstractSBP, eqn::AdvectionData{Tsol, Tres, Tdim},
                    opts, functor, functional_number, adjoint_vec::Array{Tsol,1})

  # Get information corresponding to functional_number
  key = string("geom_edges_functional", functional_number)
  functional_edges = opts[key]

  # Calculate the Jacobian of the residual
  newton_data = NonlinearSolvers.NewtonData(mesh, sbp, eqn, opts)
  res_jac = zeros(Tres, mesh.numDof, mesh.numDof)
  pert = complex(0, opts["epsilon"])
  NonlinearSolvers.calcJacobianComplex(newton_data, mesh, sbp, eqn, opts,
                                       evalAdvection, pert, res_jac)

  func_deriv = zeros(Tsol, mesh.numDof)
  func_deriv_arr = zeros(eqn.q)

  # Calculate df/dq_bndry on edges where the functional is calculated and put
  # it back in func_deriv_arr
  calcFunctionalDeriv(mesh, sbp, eqn, opts, functor, functional_edges,
                      func_deriv_arr)  # populate df_dq_bndry

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
=#

function calcAdjoint{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractDGMesh{Tmsh}, sbp::AbstractSBP,
                     eqn::AdvectionData{Tsol, Tres, Tdim}, opts,
                     functionalData::AbstractOptimizationData, adjoint_vec::Array{Tsol,1};
                     functional_number::Int=1)

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
  res_jac, jacData = calcResidualJacobian(mesh, sbp, eqn, opts)

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
  #=
  for i = 1:size(func_deriv_arr,3)
    for j = 1:size(func_deriv_arr,2)
      println("func_deriv_arr[1,$j,$i] = $(func_deriv_arr[1,j,i])")
    end
  end
  =#
  # Assemble func_deriv
  assembleArray(mesh, sbp, eqn, opts, func_deriv_arr, func_deriv)

  # Solve for adjoint vector. This depends on whether PETSc is used or not.
  if opts["jac_type"] == 1 || opts["jac_type"] == 2
    adjoint_vec[:] = -(res_jac.')\func_deriv 
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
    adjoint_vec = -adjoint_vec
  end # End how to solve for adjoint_vec

  outname = string("adjoint_vec.dat")
  f = open(outname, "w")
  for i = 1:length(adjoint_vec)
    println(f, real(adjoint_vec[i]))
  end
  close(f)


  return nothing
end # End function calcAdjoint

#= # Legacy
function calcAdjoint{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractDGMesh{Tmsh},
                    sbp::AbstractSBP, eqn::AdvectionData{Tsol, Tres, Tdim},
                    opts, functor, functional_number, adjoint_vec::Array{Tsol,1})

  # Get information corresponding to functional_number
  key = string("geom_edges_functional", functional_number)
  functional_edges = opts[key]

  # Calculate the Jacobian of the residual
  newton_data = NonlinearSolvers.NewtonData{Tsol, Tres}(mesh, sbp, eqn, opts)
  res_jac = zeros(Tres, mesh.numDof, mesh.numDof)
  pert = complex(0, opts["epsilon"])
  NonlinearSolvers.calcJacobianComplex(newton_data, mesh, sbp, eqn, opts,
                                       evalAdvection, pert, res_jac)

  # Re-interpolate interior q to q_bndry. This is done because the above step
  # pollutes the existing eqn.q_bndry with complex values.
  boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, eqn.q, eqn.q_bndry)

  # calculate the derivative of the function w.r.t q_vec
  func_deriv = zeros(Tsol, mesh.numDof)

  # 3D array into which func_deriv_arr_bndry gets interpolated
  func_deriv_arr = zeros(eqn.q)


  # Calculate df/dq_bndry on edges where the functional is calculated and put
  # it back in func_deriv_arr
  calcFunctionalDeriv(mesh, sbp, eqn, opts, functor, functional_edges,
                      func_deriv_arr)  # populate df_dq_bndry

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
end  # End function calcAdjoint
=#

@doc """
###AdvectionEquationMod.calcResidualJacobian

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
         sbp::AbstractSBP, eqn::AdvectionData{Tsol, Tres, Tdim}, opts)

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
  ctx_residual = (evalAdvection,)
  NonlinearSolvers.physicsJac(jacData, mesh, sbp, eqn, opts, jac, ctx_residual)

  return jac, jacData
end


@doc """
### AdvectionEquationMod.calcFunctionalDeriv

Computes a 3D array of the derivative of a functional w.r.t eqn.q on all
mesh nodes.

**Inputs**

*  `mesh`  : Abstract mesh object
*  `sbp`   : Summation-By-Parts operator
*  `eqn`   : Advection equation object
*  `opts`  : Options dictionary
*  `functional_edges` : Array of geometric edges over which a functional acts
*  `func_deriv_arr`   : 3D array that stors the derivative of functional w.r.t
                        eqn.q. It has a structure [1, numnodes_per_element, numEl]

**Outputs**

*  None

"""->

function calcFunctionalDeriv{Tmsh, Tsol}(mesh::AbstractCGMesh{Tmsh}, sbp::AbstractSBP,
                             eqn ::AdvectionData{Tsol}, opts, functor, functional_edges,
                             functionalData, func_deriv_arr)

  alpha_x = eqn.params.alpha_x
  alpha_y = eqn.params.alpha_y

  # Obtain the derivative of the integrand at all meh.bndry
  integrand = zeros(eqn.bndryflux)
  # println("size of integrand = ", size(integrand))

  for itr = 1:length(functional_edges)
    g_edge_number = functional_edges[itr] # Extract geometric edge number
    # get the boundary array associated with the geometric edge
    itr2 = 0
    for itr2 = 1:mesh.numBC
      if findfirst(mesh.bndry_geo_nums[itr2],g_edge_number) > 0
        break
      end
    end

    start_index = mesh.bndry_offsets[g_edge_number]
    end_index = mesh.bndry_offsets[g_edge_number + 1]
    idx_range = start_index:(end_index-1)  # Index range
    bndry_facenums = sview(mesh.bndryfaces, idx_range) # faces on geometric edge i

    nfaces = length(bndry_facenums)

    for i = 1:nfaces
      bndry_i = bndry_facenums[i]
      global_facenum = idx_range[i]
      for j = 1:sbp.numfacenodes
        k = sbp.facenodes[j, bndry_i.face]
        q = eqn.q[1,k,bndry_i.element]
        x = sview(mesh.coords, :, k, bndry_i.element)
        dxidx = sview(mesh.dxidx, :, :, k, bndry_i.element)
        nrm = sview(sbp.facenormal, :, bndry_i.face)
        nx = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
        ny = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]
        integrand[1,j,global_facenum] = calcIntegrandDeriv(opts, functor, eqn.params,
                                        nx, ny, q, functionalData)
      end  # End for j = 1:sbp.numfacenodes
    end    # End for i = 1:nfaces
  end      # End for itr = 1:length(functional_edges)


  # println("mesh.bndryfaces = \n", mesh.bndryfaces)

  boundaryintegrate!(mesh.sbpface, mesh.bndryfaces, integrand, func_deriv_arr)

  return nothing
end

# DG Version
function calcFunctionalDeriv{Tmsh, Tsol}(mesh::AbstractDGMesh{Tmsh}, sbp::AbstractSBP,
                             eqn::AdvectionData{Tsol}, opts, functor, functional_edges,
                             functionalData, func_deriv_arr)

  alpha_x = eqn.params.alpha_x
  alpha_y = eqn.params.alpha_y

  # Obtain the derivative of the integrand at all mesh.bndry
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

    for i = 1:nfaces
      bndry_i = bndry_facenums[i]
      global_facenum = idx_range[i]
      for j = 1:mesh.sbpface.numnodes
        q = eqn.q_bndry[ 1, j, global_facenum]
        coords = sview(mesh.coords_bndry, :, j, global_facenum)
        dxidx = sview(mesh.dxidx_bndry, :, :, j, global_facenum)
        nrm = sview(sbp.facenormal, :, bndry_i.face)
        nx = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
        ny = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]
        integrand[1,j,global_facenum] = calcIntegrandDeriv(opts, functor, eqn.params,
                                        nx, ny, q, functionalData)
      end  # End for j = 1:mesh.sbpface.numnodes
    end    # End for i = 1:nfaces
  end

  boundaryintegrate!(mesh.sbpface, mesh.bndryfaces, integrand, func_deriv_arr)

  return nothing
end    # End function calcFunctionalDeriv


@doc """
### AdvectionEquationMod.calcIntegrandDeriv

Compute the derivative of the integrand at a point. It presently uses complex
step to compute the derivative

**Inputs**

*  `opts`    : Input dictionary
*  `functor` : Functional name
*  `params`  : the ParamType for the equation
*  `nx` & `ny` : Normal vectors
*  `q`       : Solution variable

**Outputs**

*  `integrand_deriv` : derivative of the functor w.r.t q

"""->

function calcIntegrandDeriv(opts, functor, params::ParamType2, nx, ny, q, functionalData)

  pert = complex(0, opts["epsilon"])  # complex perturbation
  q += pert
  val = functor(params, nx, ny, q, functionalData)
  integrand_deriv = imag(val)/norm(pert)


  return integrand_deriv
end


@doc """
### AdvectionEquationMod.functionalBoundaryInfo

It creates a book-keeping array of tuples of information on all the mesh faces
over which a functional acts. The tuple contains
  1. sbpface.numnodes : number of nodes on an SBP face
  2. global_facenum : Face number on the global mesh
  3. functional_edge number : The geometric edge on which the functional acts
  4. sbpface_faceno : Face number in an SBP element

**Inputs**

*  `mesh` : Abstract mesh object
*  `sbp`  : Summation-By-Parts operator
*  `eqn`  : Advection equation object
*  `fq_bndry_info` : The book keeping array of tuples
*  `functional_edges` : The geometric edge on which the functional acts

**Outputs**

*  None

"""->

function functionalBoundaryInfo{Tmsh, Tsol}(mesh::AbstractDGMesh{Tmsh}, sbp::AbstractSBP,
         eqn::AdvectionData{Tsol}, fq_bndry_info::Array{Tuple{Int64,Int64,Int64,Int64},1},
         functional_edges)

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
  starting_index = 0
  for itr = 1:length(functional_edges)
    g_edge_number = functional_edges[itr]
    start_index = mesh.bndry_offsets[g_edge_number]
    end_index = mesh.bndry_offsets[g_edge_number+1]
    idx_range = start_index:(end_index-1)
    bndry_facenums = sview(mesh.bndryfaces, idx_range) # faces on geometric edge i
    nfaces = length(mesh.bndryfaces[idx_range])

    for i = 1:nfaces
      bndry_i = bndry_facenums[i]
      global_facenum = idx_range[i]

      for j = 1:mesh.sbpface.numnodes

        sbpface_faceno = bndry_i.face  # Face number of the edge (for 2D Tri
                                       # element: 1, 2, or 3)
        fq_bndry_info[starting_index + (i-1)*mesh.sbpface.numnodes + j] =
                               (j,global_facenum,g_edge_number, sbpface_faceno)

      end   # End for j = 1:mesh.sbpface.numnodes
    end     # End for i = 1:nfaces

    nfaces_prev = nfaces
    starting_index += nfaces_prev*mesh.sbpface.numnodes
  end       # End for itr = 1:length(functional_edges)

  return nothing
end


@doc """
###AdvectionEquationMod.boundaryToVolumeInterpolation

Interpolates values from a array based on Summation-By-Parts sbpface boundary
to the interior nodes of the SBP Omega operator.

**Inputs**

*  `sbpface` : Summation-By-Parts face operator
*  `bndryfaces` : Array of type mesh.bndryfaces
*  `uvol`    : Array into which the values have to be interpolated
*  `uface`   : Array which has to be interpolated.

**Outputs**

*  None

"""->
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
