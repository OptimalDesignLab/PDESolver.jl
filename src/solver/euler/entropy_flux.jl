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

function getPsi(params, qL::AbstractMatrix{Tsol}, qR::AbstractMatrix{Tsol}, nrm::AbstractVector) where Tsol

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

  Inputs:
    params: a ParamType
    iface: the Interface
    sbpface: an AbstractFace
    nrm_scaled: the scaled normal vector at the face nodes in x-y space
    qL: the solution at the volume nodes of elementL
    qR: the solution at the volume nodes of elementR

  Outputs:
    rhs: the value of the integral
"""
function computeInterfacePotentialFlux(
                params::ParamType{Tdim, :conservative, Tsol, Tres}, 
                iface::Interface, sbpface, nrm_scaled::AbstractMatrix, 
                qL::AbstractMatrix, qR::AbstractMatrix) where {Tdim, Tsol, Tres}
# compute the potential flux then compute the reduction with Eface

  
  rhs = reduceEface(params, iface, sbpface, nrm_scaled, dim, qL, qR)
  return rhs
end

@doc """
  Computes 1^T(E_faceL * psiL + E_faceR * psiR), where E_face performs the 
  face integral over the specified face and psi is the potential flux 
  (from entropy stability theory).  The dir argument indicates whether 
  the normal vector used for the integral is x, y, or z, ie. in 

  integral psi dot n dGamma

  it determines the normal vector n

  TODO: this is currently broken for SparseFace SBP operators
"""
function reduceEface(
params::ParamType{Tdim, :conservative, Tsol, Tres, Tmsh},
iface::Interface, sbpface, nrm_scaled::AbstractMatrix,
qL::AbstractMatrix, qR::AbstractMatrix) where {Tdim, Tsol, Tres, Tmsh}
  # compute Ex_gamma kappa * psiL + Ex_gamma_nu * psiR, where x is one 
  # of either x or y, as specified by dir

  # interpolate q to the faces
  qfaceL = zeros(Tsol, size(qL, 1), sbpface.numnodes)
  qfaceR = zeros(qfaceL)
  interiorFaceInterpolate!(sbpface, iface, qL, qR, qfaceL, qfaceR)

  psiL = zeros(Tres, sbpface.numnodes)
  psiR = zeros(Tres, sbpface.numnodes)
  # compute psi dot n
  for i=1:sbpface.numnodes
    nrm_i = sview(nrm_scaled, :, i)
    psiL[i] = getPsi(params, sview(qfaceL, :, i), nrm_i)
    psiR[i] = -getPsi(params, sview(qfaceR, :, i), nrm_i)  # negative sign for
                                                           # direction
  end

  val = integrateBoundaryFunctional!(sbpface, iface.faceL, psiL)
  val += integrateBoundaryFunctional!(sbpface, iface.faceR, psiR)

  return val
end

#TODO: this can be made more efficient once SBP stores E
function computeVolumePotentialFlux(params::ParamType{Tdim, :conservative, Tsol, Tres}, sbp, q_i::AbstractMatrix, dxidx) where {Tdim, Tsol, Tres}
  
  numDofPerNode, numNodesPerElement = size(q_i)

  nrm = params.nrm
  # calculate psi vector
  psi = zeros(numNodesPerElement, 2)
  for j=1:numNodesPerElement
    q_j = ro_sview(q_i, :, j)
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
function calcEntropyIntegral(mesh::AbstractMesh{Tmsh}, 
             sbp, eqn::EulerData{Tsol, Tres, Tdim}, opts, q_vec::AbstractVector) where {Tsol, Tres, Tmsh, Tdim}

  val = zero(Tsol)
  for i=1:mesh.numDofPerNode:mesh.numDof
    q_vals = ro_sview(eqn.q_vec, i:(i+(mesh.numDofPerNode - 1)))
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
function contractResEntropyVars(
             mesh::AbstractDGMesh{Tmsh}, 
             sbp, eqn::EulerData{Tsol, Tres, Tdim}, opts, q_vec::AbstractVector,
             res_vec::AbstractVector) where {Tsol, Tres, Tmsh, Tdim}

  val = zero(Tres)
  w_vals = eqn.params.contract_res_entropy_vars_data.w_vals
  for i=1:mesh.numDofPerNode:mesh.numDof
    q_vals_i = ro_sview(eqn.q_vec, i:(i+mesh.numDofPerNode - 1))
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

"""
  Like contractResEntropyVars, but uses a special summation technique
"""
function contractResEntropyVars2(
             mesh::AbstractDGMesh{Tmsh}, 
             sbp, eqn::EulerData{Tsol, Tres, Tdim}, opts, q_vec::AbstractVector,
             res_vec::AbstractVector) where {Tsol, Tres, Tmsh, Tdim}

  vals = zeros(Tres, mesh.numNodes)
  w_vals = eqn.params.contract_res_entropy_vars_data.w_vals
  idx = 1
  for i=1:mesh.numDofPerNode:mesh.numDof
    q_vals_i = ro_sview(eqn.q_vec, i:(i+mesh.numDofPerNode - 1))
    convertToEntropy(eqn.params, q_vals_i, w_vals)
    scale!(w_vals, 1./eqn.params.gamma_1)  # the IR entropy variables are
                                           # scaled by 1/gamma compared to
                                           # Hugh's
    res_vals = sview(res_vec, i:(i+mesh.numDofPerNode - 1))
    for p=1:mesh.numDofPerNode
      vals[idx] += w_vals[p]*res_vals[p]
    end
    idx += 1
  end

  # now do the sum in the order of decreasing absolute value
  pvec = sortperm(abs(vals))
  val = zero(Tres)
  for i=mesh.numNodes:-1:1  # decreasing absolute value
    val += vals[pvec[i]]
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
function integrateQ( mesh::AbstractDGMesh{Tmsh}, 
             sbp, eqn::EulerData{Tsol, Tres, Tdim}, opts, q_vec::AbstractVector) where {Tsol, Tres, Tmsh, Tdim}

  vals = zeros(Tsol, mesh.numDofPerNode)
  for i=1:mesh.numDofPerNode:mesh.numDof
    q_vals_i = ro_sview(eqn.q_vec, i:(i+mesh.numDofPerNode - 1))
    w_val = eqn.M[i]
    for j=1:length(q_vals_i)
      vals[j] += w_val*q_vals_i[j]
    end
  end

  vals2 = zeros(vals)
  MPI.Allreduce!(vals, vals2, MPI.SUM, eqn.comm)

  return vals2
end

@doc """
  Computes the net potential flux integral over all interfaces, where the 
  potential flux is calculated from q_arr.
  Does not work in parallel
"""
function calcInterfacePotentialFlux(
           mesh::AbstractDGMesh{Tmsh}, sbp, 
           eqn::EulerData{Tsol, Tres, Tdim}, opts, 
           q_arr::Abstract3DArray) where {Tsol, Tres, Tdim, Tmsh}

  if mesh.commsize != 1
    throw(ErrorException("cannot perform reduction over interfaces in parallel"))
  end

  val = zero(Tres)
  for i=1:mesh.numInterfaces
    iface = mesh.interfaces[i]
    elL = iface.elementL
    elR = iface.elementR
    qL = ro_sview(q_arr, :, :, elL)
    qR = ro_sview(q_arr, :, :, elR)
    nrm_scaled = ro_sview(mesh.nrm_face, :, :, i)

    bndry_potentialflux = -computeInterfacePotentialFlux(eqn.params, iface, mesh.sbpface, nrm_scaled, qL, qR)
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
function calcVolumePotentialFlux(mesh::AbstractMesh{Tmsh}, 
               sbp, eqn::EulerData{Tsol, Tres}, opts, 
               q_arr::Abstract3DArray) where {Tsol, Tres, Tmsh}

  val = zero(Tres)
  for i=1:mesh.numEl
    q_i = ro_sview(q_arr, :, :, i)
    dxidx_i = ro_sview(mesh.dxidx, :, :, :, i)
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
function calcEnstrophy(mesh::AbstractMesh{Tmsh}, sbp,
                       eqn::EulerData{Tsol, Tres}, opts,
                       q_arr::Abstract3DArray{Tsol}) where {Tsol, Tres, Tmsh}

  @assert mesh.dim == 3
  Tdim = 3

  val = zero(Tres)
  vorticity = zeros(Tres, Tdim, mesh.numNodesPerElement)
  for i=1:mesh.numEl
    q_i = ro_sview(q_arr, :, :, i)
    dxidx_i = ro_sview(mesh.dxidx, :, :, :, i)
    jac_i = ro_sview(mesh.jac, :, i)

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


function getVorticity(mesh::AbstractMesh{Tmsh}, sbp, eqn::EulerData{Tsol, Tres}, opts) where {Tmsh, Tsol, Tres}
# get an array of the vorticity at every node
# also computes enstrophy

  new_field = zeros(mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)
  vorticity = zeros(Tres, mesh.dim, mesh.numNodesPerElement)
  for i=1:mesh.numEl
    q_i = ro_sview(eqn.q, :, :, i)
    dxidx_i = ro_sview(mesh.dxidx, :, :, :, i)
    jac_i = ro_sview(mesh.jac, :, i)

    calcVorticity(eqn.params, sbp, q_i, dxidx_i, jac_i, vorticity)

    # copy into new_field
    for j=1:mesh.numNodesPerElement
      for k=1:3
        new_field[k, j, i] = real(vorticity[k, j])
      end
      vorticity_mag = real(sqrt(vorticity[1, j]^2 + vorticity[2, j]^2 + vorticity[3, j]^2))
      new_field[4, j, i] = real(vorticity_mag)
      new_field[5, j, i] = real(eqn.q[1, j, i]*vorticity_mag)

    end

  end

  myrank = mesh.myrank
  writedlm("q2.dat", new_field)

  return new_field
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
function calcKineticEnergy(mesh::AbstractMesh{Tmsh}, sbp, 
   eqn::EulerData{Tsol, Tres, Tdim}, opts, 
   q_vec::AbstractVector{Tsol}) where {Tsol, Tres, Tdim, Tmsh}


  val = zero(Tsol)
  for i=1:mesh.numDofPerNode:mesh.numDof
#    println(f, "node ", i)
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

  val = 0.5*val/mesh.volume

  return val
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
function calcKineticEnergydt(mesh::AbstractMesh{Tmsh},
      sbp, eqn::EulerData{Tsol, Tres, Tdim}, opts, 
      q_vec::AbstractVector{Tsol}, 
      res_vec::AbstractVector{Tres}) where {Tsol, Tres, Tdim, Tmsh}

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


