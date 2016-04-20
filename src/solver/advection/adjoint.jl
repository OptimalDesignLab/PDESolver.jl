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
function calcAdjoint{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractDGMesh{Tmsh},
                    sbp::AbstractSBP, eqn::AdvectionData{Tsol, Tres, Tdim},
                    opts, functional_number, adjoint_vec::Array{Tsol,1})

  # Get information corresponding to functional_number
  key = string("geom_edges_functional", functional_number)
  functional_edges = opts[key]

  # Calculate the Jacobian of the residual
  newton_data = NonlinearSolvers.NewtonData(mesh, sbp, eqn, opts)
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
  calcFunctionalDeriv(mesh, sbp, eqn, opts, functional_edges, 
                      func_deriv_arr)  # populate df_dq_bndry

  # Assemble func_deriv
  assembleArray(mesh, sbp, eqn, opts, func_deriv_arr, func_deriv)

  # Solve for adjoint vector

  # TODO: The following operation creates a temporary copy of adjoint_vec, does
  #       the '\' computation and then puts it back into adjoint_vec. This
  #       needs to change.
  adjoint_vec[:] = -(res_jac.')\func_deriv
  
  return nothing
end  # End function calcAdjoint


function calcFunctionalDeriv{Tmsh, Tsol}(mesh::AbstractDGMesh{Tmsh}, sbp::AbstractSBP,
                             eqn ::AdvectionData{Tsol}, opts, functional_edges, 
                             func_deriv_arr)

  alpha_x = eqn.alpha_x
  alpha_y = eqn.alpha_y

  # Initialize array of derivative of functional integrand over all edges
  # func_deriv_arr = zeros(eqn.q)  # 3D array that holds df/dq_bndry on ALL edges

  # Obtain the derivative of the integrand at all meh.bndry
  integrand = zeros(eqn.q_bndry)

  # Populate integrand
  for itr = 1:length(functional_edges)
    g_edge_number = functional_edges[itr] # Extract geometric edge number
    start_index = mesh.bndry_offsets[g_edge_number]
    end_index = mesh.bndry_offsets[g_edge_number+1]
    idx_range = start_index:(end_index-1)
    bndry_facenums = view(mesh.bndryfaces, idx_range) # faces on geometric edge i

    nfaces = length(bndry_facenums)

    for i = 1:nfaces
      bndry_i = bndry_facenums[i]
      global_facenum = idx_range[i]
      for j = 1:mesh.sbpface.numnodes
        q = eqn.q_bndry[ 1, j, global_facenum]
        # println("eqn.q_bndry[1,$j,$global_facenum] = ", q)
        coords = view(mesh.coords_bndry, :, j, global_facenum)
        dxidx = view(mesh.dxidx_bndry, :, :, j, global_facenum)
        nrm = view(sbp.facenormal, :, bndry_i.face)
        nx = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
        ny = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]
        integrand[1,j,global_facenum] = calcIntegrandDeriv(opts, alpha_x, alpha_y, nx, ny, q)
      end  # End for j = 1:mesh.sbpface.numnodes
    end    # End for i = 1:nfaces
  end

  #=
  for (bindex, bndry) in enumerate (mesh.bndryfaces)
    for i = 1:mesh.sbpface.numnodes
      wflux = mesh.sbpface.wface[i]*integrand[1,i,bindex]
      for j = 1:mesh.sbpface.stencilsize
        func_deriv_arr[1, mesh.sbpface.perm[j,bndry.face], bndry.element] += 
                                            mesh.sbpface.interp[j,i]*wflux
      end
    end  # End for j = 1:mesh.sbpface.numnodes
  end    # End enumerate
  =#

  boundaryintegrate!(mesh.sbpface, mesh.bndryfaces, integrand, func_deriv_arr)

  return nothing
end    # End function calcFunctionalDeriv


@doc """
### AdvectionEquationMod.calcIntegrandDeriv

Compute the analytical derivative of the integrand at a point

"""->

function calcIntegrandDeriv(opts, alpha_x, alpha_y, nx, ny, q)

  pert = complex(0, opts["epsilon"])  # complex perturbation
  q += pert
  val = calcFunctionalIntegrand(alpha_x, alpha_y, nx, ny, q)
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
    bndry_facenums = view(mesh.bndryfaces, idx_range) # faces on geometric edge i
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