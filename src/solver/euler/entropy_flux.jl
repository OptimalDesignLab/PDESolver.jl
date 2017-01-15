@doc """
  This function calculates the potential flux psi from entropy stability 
  theory.  Methods are available for both 2 and 3D
"""

function getPsi(params::ParamType{2, :conservative}, q_vals::AbstractVector, nrm::AbstractVector)
  return nrm[1]*q_vals[2] + nrm[2]*q_vals[3]
end

function getPsi(params::ParamType{3, :conservative}, q_vals::AbstractVector, nrm::AbstractVector)
  return nrm[1]*q_vals[2] + nrm[2]*q_vals[3] +nrm[3]*q_vals[4]
end

function getPsi{Tsol}(params, qL::AbstractMatrix{Tsol}, qR::AbstractMatrix{Tsol}, nrm::AbstractVector)

  numDofPerNode, numNodesPerElement = size(qL)
  psiL = zeros(Tsol, numNodesPerElement)
  psiR = zeros(psiL)
  dim = length(nrm)
  # calculate psi at volume nodes
  for i=1:numNodesPerElement
    qL_i = qL[:, i]
    qR_i = qR[:, i]

    psi_vecL = psi_vec(params, qL_i)
    psi_vecR = psi_vec(params, qR_i)
    for d=1:dim
      psiL[i] += nrm[d]*psi_vecL[d]
      psiR[i] += nrm[d]*psi_vecR[d]
    end
  end

  return psiL, psiR
end

@doc """
  This function computes the integral of the potential flux over an interface
"""
function computeInterfacePotentialFlux{Tdim, Tsol, Tres}(
                params::ParamType{Tdim, :conservative, Tsol, Tres}, 
                iface::Interface, sbpface, dxidx_face, 
                qL::AbstractMatrix, qR::AbstractMatrix)
# compute the potential flux then compute the reduction with Eface

  rhs = zero(Tres)
  for dim=1:Tdim
    rhs += reduceEface(params, iface, sbpface, dxidx_face, dim, qL, qR)
  end

  return rhs
end

@doc """
  Computes 1^T(E_faceL * psiL + E_faceR * psiR), where E_face performs the 
  face integral over the specified face and psi is the potential flux 
  (from entropy stability theory).  The dir argument indicates whether 
  the normal vector used for the integral is x, y, or z, ie. in 

  integral psi dot n dGamma

  it determines the normal vector n
"""
function reduceEface{Tdim, Tsol, Tres, Tmsh}(params::ParamType{Tdim, :conservative, Tsol, Tres, Tmsh}, iface, sbpface, dxidx_face::Abstract3DArray{Tmsh}, dir::Integer, qL::AbstractMatrix, qR::AbstractMatrix)
  # compute Ex_gamma kappa * psiL + Ex_gamma_nu * psiR, where x is one 
  # of either x or y, as specified by dir

  RHS1 = zero(Tres)
  RHS2 = zero(Tres)

  flux_nrm = params.nrm
  fill!(flux_nrm, 0.0)
  flux_nrm[dir] = 1

  for i=1:sbpface.stencilsize
    for j=1:sbpface.stencilsize
      p_jL = sbpface.perm[j, iface.faceL]
      p_jR = sbpface.perm[j, iface.faceR]
      psiL = getPsi(params, sview(qL, :, p_jL), flux_nrm)
      psiR = getPsi(params, sview(qR, :, p_jR), flux_nrm)
      for k=1:sbpface.numnodes
        nrm_k = zero(Tmsh)
        for d=1:Tdim
          nrm_k += sbpface.normal[d, iface.faceL]*dxidx_face[d, dir, k]
        end
        val = sbpface.interp[i,k]*sbpface.interp[j,k]*sbpface.wface[k]*nrm_k
        RHS1 += val*psiL

        kR = sbpface.nbrperm[k, iface.orient]
        val = sbpface.interp[i, kR]*sbpface.interp[j, kR]*sbpface.wface[k]*nrm_k
        RHS1 -= val*psiR
      end
    end
  end

  return RHS1 + RHS2
end

#TODO: this can be made more efficient once SBP stores E
function computeVolumePotentialFlux{Tdim, Tsol, Tres}(params::ParamType{Tdim, :conservative, Tsol, Tres}, sbp, q_i::AbstractMatrix, dxidx)
  
  numDofPerNode, numNodesPerElement = size(q_i)

  nrm = params.nrm
  # calculate psi vector
  psi = zeros(numNodesPerElement, 2)
  for j=1:numNodesPerElement
    q_j = sview(q_i, :, j)
    for d=1:Tdim
      fill!(nrm, 0.0)
      nrm[d] = 1
      psi[j, d] = getPsi(params, q_j, nrm)
    end
  end

  rhs_reduced = zero(Tres)
  for d=1:Tdim
#    println("dimension = ", d)
    E_d = (sbp.Q[:, :, d] + sbp.Q[:, :, d].')
    psi_nrm = dxidx[d, 1, 1]*psi[:, 1] + dxidx[d, 2, 1]*psi[:, 2]
#    println("psi_nrm = ", psi_nrm)
    val = sum(E_d*psi_nrm)
#    println("rhs_reduced = ", val)

    rhs_reduced += val
  end

  return -rhs_reduced
end

@doc """
  Calculate the integral of entropy over the entire domain for the given 
  solution vector q_vec.  
  This performs an MPI blocking collective operation, so all processes must 
  call this function at the same time.

  The (scalar) value is returned.
"""
function calcEntropyIntegral{Tsol, Tres, Tmsh, Tdim}(mesh::AbstractMesh{Tmsh}, 
             sbp, eqn::EulerData{Tsol, Tres, Tdim}, opts, q_vec::AbstractVector)

  val = zero(Tsol)
  for i=1:mesh.numDofPerNode:mesh.numDof
    q_vals = sview(eqn.q_vec, i:(i+(mesh.numDofPerNode - 1)))
#      s = calcEntropy(eqn.params, q_vals)
    s = calcEntropyIR(eqn.params, q_vals)
#      val += real(q_vals[1]*s)*eqn.M[i]
    val += s*eqn.M[i]
  end

  val2 = MPI.Allreduce(val, MPI.SUM, eqn.comm)
  return val2
end

@doc """
  Compute w^T * res_vec, where w is the vector of entropy variables.
  This performs an MPI blocking collective operation, so all processes must 
  call this funciton at the same time.

  The (scalar) value is returned.
"""
function contractResEntropyVars{Tsol, Tres, Tmsh, Tdim}(
             mesh::AbstractDGMesh{Tmsh}, 
             sbp, eqn::EulerData{Tsol, Tres, Tdim}, opts, q_vec::AbstractVector,
             res_vec::AbstractVector)

  val = zero(Tres)
  w_vals = eqn.params.v_vals
  for i=1:mesh.numDofPerNode:mesh.numDof
    q_vals_i = sview(eqn.q_vec, i:(i+mesh.numDofPerNode - 1))
    convertToEntropy(eqn.params, q_vals_i, w_vals)
    scale!(w_vals, 1./eqn.params.gamma_1)  # the IR entropy variables are
                                           # scaled by 1/gamma compared to
                                           # Hugh's
    res_vals = sview(res_vec, i:(i+mesh.numDofPerNode - 1))
    for p=1:mesh.numDofPerNode
      val += w_vals[p]*res_vals[p]
    end
  end

  val = MPI.Allreduce(val, MPI.SUM, eqn.comm)
  return val
end

@doc """
  Compute the SBP approximation to integral q dOmega, ie. the mass matrix
  times the vector of conservative variables at each node in the mesh.  

  This function returns an array of length numDofPerNode containing the 
  result for each variable.
"""
function integrateQ{Tsol, Tres, Tmsh, Tdim}( mesh::AbstractDGMesh{Tmsh}, 
             sbp, eqn::EulerData{Tsol, Tres, Tdim}, opts, q_vec::AbstractVector)

  vals = zeros(Tsol, mesh.numDofPerNode)
  for i=1:mesh.numDofPerNode:mesh.numDof
    q_vals_i = sview(eqn.q_vec, i:(i+mesh.numDofPerNode - 1))
    w_val = eqn.M[i]
    for j=1:length(q_vals_i)
      vals[j] += w_val*q_vals_i[j]
    end
  end

  vals2 = zeros(vals)
  MPI.Allreduce(vals, vals2, MPI.SUM, eqn.comm)

  return vals2
end




@doc """
  Computes the net potential flux integral over all interfaces, where the 
  potential flux is calculated from q_arr.
  Does not work in parallel
"""
function calcInterfacePotentialFlux{Tsol, Tres, Tdim, Tmsh}(
                                   mesh::AbstractDGMesh{Tmsh}, sbp, 
                                   eqn::EulerData{Tsol, Tres, Tdim}, opts, 
                                   q_arr::Abstract3DArray)

  if mesh.commsize != 1
    throw(ErrorException("cannot perform reduction over interfaces in parallel"))
  end

  val = zero(Tres)
  for i=1:mesh.numInterfaces
    iface = mesh.interfaces[i]
    elL = iface.elementL
    elR = iface.elementR
    qL = sview(q_arr, :, :, elL)
    qR = sview(q_arr, :, :, elR)
    aux_vars = sview(eqn.aux_vars, :, :, elL)
    dxidx_face = sview(mesh.dxidx_face, :, :, :, i)

    bndry_potentialflux = -computeInterfacePotentialFlux(eqn.params, iface, mesh.sbpface, dxidx_face, qL, qR)
    val += bndry_potentialflux
  end

  return val
end

@doc """
  Computes the net potential flux from the volume terms, calulated from 
  q_arr.
  This performs a blocking MPI collective operation, so all processes must
  call this function at the same time.
"""
function calcVolumePotentialFlux{Tsol, Tres, Tmsh}(mesh::AbstractMesh{Tmsh}, 
                                 sbp, eqn::EulerData{Tsol, Tres}, opts, 
                                 q_arr::Abstract3DArray)

  val = zero(Tres)
  for i=1:mesh.numEl
    q_i = sview(q_arr, :, :, i)
    dxidx_i = sview(mesh.dxidx, :, :, :, i)
    volume_potentialflux = -computeVolumePotentialFlux(eqn.params, sbp, q_i, dxidx_i)
    val += volume_potentialflux
  end

  val = MPI.Allreduce(val, MPI.SUM, mesh.comm)

  return val
end


"""
  Calculates ( 1/(2V) )*integral(rho * omega dot omega dV), where V is the
  volume of the mesh omega is the vorticity vector, and rho is the density.
  3D, conservative variables only.  This should work for CG and DG, but
  has only been tested for the latter

  Inputs:
    mesh: an AbstractMesh
    sbp: an SBP operator
    eqn: an EulerData object
    opts: options dictionary
    q_arr: a 3D array of conservative variables 

  Outputs:
    val: the value of the integral (over the entire domain, not just the
         part owned by this process)

  Aliasing restrictions: see calcVorticity
"""
function calcEnstrophy{Tsol, Tres, Tmsh}(mesh::AbstractMesh{Tmsh}, sbp,
                                         eqn::EulerData{Tsol, Tres}, opts,
                                         q_arr::Abstract3DArray{Tsol})

  @assert mesh.dim == 3
  Tdim = 3

  val = zero(Tres)
  vorticity = zeros(Tres, Tdim, mesh.numNodesPerElement)
  for i=1:mesh.numEl
    q_i = sview(q_arr, :, :, i)
    dxidx_i = sview(mesh.dxidx, :, :, :, i)
    jac_i = sview(mesh.jac, :, i)

    calcVorticity(eqn.params, sbp, q_i, dxidx_i, jac_i, vorticity)

    for j=1:mesh.numNodesPerElement

      rho = q_i[1, j]

      # accumulate vorticity dot vorticity
      vorticity_mag = zero(Tres)
      for k=1:Tdim
        vorticity_mag += vorticity[k, j]*vorticity[k, j]
      end

      val += sbp.w[j]*rho*vorticity_mag/jac_i[j]
    end  # end loop j
  end  # end loop i

  val = MPI.Allreduce(val, MPI.SUM, mesh.comm)

  return 0.5*val/mesh.volume
end

"""
  This function calculates ( 1/(2*V) )*integral(rho * v dot v dV), where
  V is the volume of the mesh, v is the velocity vector, and rho is the density.
  This is the total kinetic energy normalized by the volume of the domain
  Conservative variables only.

  This function contains an MPI blocking collective operation.  It must be
  called by all processes at the same time.

  This function relies on the sequential numbering of dofs on the same node

  Inputs:
    mesh: a mesh
    sbp: an SBP Operator
    eqn: an EulerData object
    opts: options dictionary
    q_vec: the vector of conservative variables for the entire mesh

  Outputs:
    val: the value of the integral (over the entire domain, not just teh part
         owned by this procss)

  Aliasing restrictions: none
"""
function calcKineticEnergy{Tsol, Tres, Tdim, Tmsh}(mesh::AbstractMesh{Tmsh}, sbp, 
                           eqn::EulerData{Tsol, Tres, Tdim}, opts, 
                           q_vec::AbstractVector{Tsol})


  val = zero(Tsol)
  for i=1:mesh.numDofPerNode:mesh.numDof
    rho_i = q_vec[i]

    # compute v dot v
    v_magnitude = zero(Tsol)
    for j=1:Tdim
      v_j = q_vec[i+j]/rho_i
      v_magnitude += v_j*v_j
    end

    val += eqn.M[i]*rho_i*v_magnitude
  end

  val = MPI.Allreduce(val, MPI.SUM, mesh.comm)

  return 0.5*val/mesh.volume
end


"""
  This function calclates the time derivative of calcKineticEnergy.

  The idea is to expand the left hand side of d rho*u/dt = res using the 
  product rule and solve for du/dt, then use it to compute the integral.
  Inputs:
    mesh: a mesh
    sbp: a SBP operator
    eqn: an EulerData
    opts: options dictionary
    q_vec: vector of conservative variables for the entire mesh
    res_vec: residual vector (dq/dt) of entire mesh

  Aliasing restrictions: none
"""
function calcKineticEnergydt{Tsol, Tres, Tdim, Tmsh}(mesh::AbstractMesh{Tmsh},
                              sbp, eqn::EulerData{Tsol, Tres, Tdim}, opts, 
                              q_vec::AbstractVector{Tsol}, 
                              res_vec::AbstractVector{Tres})

  val = zero(Tres)
  for i=1:mesh.numDofPerNode:mesh.numDof
    rho_i = q_vec[i]
    drhodt = res_vec[i]  # time derivative of rho
    # accumulate v dot rho*dv/dt
    term_i = zero(Tres)
    for j=1:Tdim
      v_j = q_vec[i + j]/rho_i
      rho_dvdt = res_vec[i + j] - drhodt*v_j

      term_i += v_j*rho_dvdt
    end

    val += eqn.M[i]*term_i
  end

  val = MPI.Allreduce(val, MPI.SUM, mesh.comm)

  return val/mesh.volume
end


