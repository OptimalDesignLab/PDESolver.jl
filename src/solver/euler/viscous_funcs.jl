abstract AbstractBoundaryValueType

@doc """

Compute derivative operators

Input:
mesh
sbp
dxidx : derivatives of mapping, i.e., Jacobian matrix
jac   : determinant of Jacobian

Output:
Dx    : (in/out)derivative operators in physical domain, incling Dx, Dy
"""->
function calcDx{Tmsh}(sbp::AbstractSBP,
                      dxidx::AbstractArray{Tmsh, 3},
                      jac::AbstractArray{Tmsh, 1},
                      Dx::AbstractArray{Tmsh, 3})
  @assert(size(Dx, 1) == sbp.numnodes)
  @assert(size(Dx, 1) == size(dxidx, 3))
  @assert(size(Dx, 2) == size(Dx, 1))
  @assert(size(Dx, 3) == size(dxidx, 1))

  Dx[:,:,:] = 0.0
  dim = size(Dx, 3)
  numNodes = sbp.numnodes

  for d=1:dim            # loop over direction in which derivative is computing
    for dd=1:dim
      for n1 = 1:numNodes
        for n2 = 1:numNodes
          # Since dxidx is scaled by 1/|J|, we need to get it back,
          # that's why jac is here
          Dx[n1, n2, d] += dxidx[dd, d, n1]*jac[n1]*sbp.Q[n1, n2, dd]
        end
      end
    end

    # Until here Dx stores Qx, we need to left multiply H^(-1)
    for n2=1:numNodes
      for n1 = 1:numNodes
        Dx[n1, n2, d] /= sbp.w[n1]
      end
    end
  end
  return nothing
end


@doc """

Compute derivative operators

Input:
mesh
sbp
elem    : index of element of which we are computing the derivatives

Output:
Dx        : (in/out) derivative operators in physical domain, incling Dx, Dy
"""->
function calcDx{Tmsh}(mesh::AbstractMesh{Tmsh},
                      sbp::AbstractSBP,
                      elem::Integer,
                      Dx::AbstractArray{Tmsh, 3})
  @assert(size(Dx, 1) == mesh.numNodesPerElement)
  @assert(size(Dx, 2) == mesh.numNodesPerElement)
  @assert(size(Dx, 3) == size(mesh.dxidx, 1))
  dxidx = sview(mesh.dxidx, :,:,:,elem) # (dim, dim, numNodesPerElement)
  jac = sview(mesh.jac, :, elem)
  dim = size(Dx, 3)

  for i = 1 : length(Dx)
    Dx[i] = 0.0
  end

  for d = 1 : dim            # loop over direction in which derivative is computing
    for dd = 1 : dim
      for n1 = 1 : mesh.numNodesPerElement
        for n2 = 1 : mesh.numNodesPerElement
          # Since dxidx is scaled by 1/|J|, we need to get it back,
          # that's why jac is here
          Dx[n1, n2, d] += dxidx[dd, d, n1] * jac[n1] * sbp.Q[n1, n2, dd]
        end
      end
    end

    # Until here Dx stores Qx, we need to left multiply H^(-1)
    for n2=1:mesh.numNodesPerElement
      for n1 = 1:mesh.numNodesPerElement
        Dx[n1, n2, d] /= sbp.w[n1]
      end
    end
  end

  return nothing
end


function calcQx{Tmsh}(mesh::AbstractMesh{Tmsh},
                      sbp::AbstractSBP,
                      elem::Integer,
                      Qx::AbstractArray{Tmsh, 3})
  @assert(size(Qx, 1) == mesh.numNodesPerElement)
  @assert(size(Qx, 2) == mesh.numNodesPerElement)
  @assert(size(Qx, 3) == size(mesh.dxidx, 1))

  dxidx = sview(mesh.dxidx, :,:,:,elem) # (dim, dim, numNodesPerElement)
  jac = mesh.jac[:, elem]
  Qx[:,:,:] = 0.0
  dim = size(Qx, 3)

  for d = 1 : dim
    for dd = 1 : dim
      for l=1:mesh.numNodesPerElement
        # Since dxidx is scaled by 1/|J|, we need to get it back,
        # that's why jac is here
        for m = 1:mesh.numNodesPerElement
          Qx[l,m,d] += dxidx[dd,d,l] * jac[l] * sbp.Q[l,m,dd]
        end
      end
    end
  end
end


@doc """

Given variables q at element nodes, compute corresponding gradients

Input:
sbp        : sbp operator
dxidx    : derivatives of mapping, i.e., jacobian matrix
jac        : determinant of jacobian
q        : element node value
q_grad    : gradient of q
Output:
nothing
"""->
function calcGradient{Tmsh, Tsol, Tsbp}(sbp::AbstractSBP{Tsbp},
                                        dxidx::AbstractArray{Tmsh, 3},
                                        jac::AbstractArray{Tmsh, 1},
                                        q::AbstractArray{Tsol, 2},
                                        q_grad::AbstractArray{Tsol, 3})
  @assert(size(q, 2) == sbp.numnodes)

  @assert(size(dxidx, 1) == size(dxidx, 2))
  @assert(size(q_grad, 1) == size(dxidx, 1))
  @assert(size(q_grad, 2) == size(q, 1))
  @assert(size(q_grad, 3) == sbp.numnodes)

  numNodes = sbp.numnodes
  numDofs = size(q, 1)
  dim = size(q_grad, 1)

  Dx = Array(Tsbp, numNodes, numNodes, 2)

  calcDx(sbp, dxidx, jac, Dx)

  for i = 1 : length(q_grad)
    q_grad[i] = 0.0
  end

  for n=1:numNodes
    for iDof=1:numDofs
      for d=1:dim
        for col=1:numNodes
          q_grad[d, iDof, n] += Dx[n,col,d] * q[iDof, col]
        end
      end
    end
  end
  return nothing
end


@doc """

Given variables q at element nodes, compute corresponding gradients

Input:
mesh:
sbp:
q      : element node value
elem   : index of element
Output :
q_grad : (in/out) gradient of q
"""->
function calcGradient{Tmsh, Tsol, Tsbp}(mesh::AbstractDGMesh{Tmsh},
                                        sbp::AbstractSBP{Tsbp},
                                        elem::Integer,
                                        q::AbstractArray{Tsol, 2},
                                        q_grad::AbstractArray{Tsol, 3})
  @assert(size(q, 2) == mesh.numNodesPerElement)
  @assert(size(q, 1) == mesh.numDofPerNode)

  @assert(size(q_grad, 1) == size(mesh.coords, 1))
  @assert(size(q_grad, 3) == mesh.numNodesPerElement)
  @assert(size(q_grad, 2) == mesh.numDofPerNode)

  numNodes = mesh.numNodesPerElement
  numDofs = mesh.numDofPerNode
  dim = size(q_grad, 1)

  Dx = Array(Tsbp, numNodes, numNodes, 2)
  # for e=1:numElems
  # First compute Dx for this element
  calcDx(mesh, sbp, elem, Dx)

  q_grad[:,:,:] = 0.0

  for n=1:numNodes
    for iDof=1:numDofs
      for d=1:dim
        for col=1:numNodes
          q_grad[d, iDof, n] += Dx[n,col,d]*q[iDof, col]
        end
      end
    end
  end
  return nothing
end


@doc """
Another(single face) version of interiorfaceintegrate.
Given Q on element L and element R, interpolate Q to interface shared by L and R

Input:
sbpface        : face SBP operator
face        : the interface which we are interpolating Q onto
qvolL        : Q at nodes in element L, qvolL(idof, iNode)    
qvolR        : Q at nodes in element R    
Output:
qface        : Q at nodes on interface, qface(idof, L/R, ifacenode)

function call examples
"""->
function interiorfaceinterpolate{Tsbp, Tsol}(sbpface::AbstractFace{Tsbp},
                                             face::Interface,
                                             qvolL::AbstractArray{Tsol, 2},
                                             qvolR::AbstractArray{Tsol, 2},
                                             qface::AbstractArray{Tsol, 3})
  @assert(size(qvolR, 1) == size(qvolL, 1))
  @assert(size(qvolR, 2) == size(qvolL, 2))
  @assert(size(qvolR, 1) == size(qface, 1))
  @assert(size(sbpface.interp, 1) <= size(qvolL, 2))
  @assert(size(sbpface.interp, 2) == size(qface, 3))

  numDofs = size(qvolL, 1)    

  for n = 1 : sbpface.numnodes
    iR = sbpface.nbrperm[n, face.orient]

    for dof = 1 : numDofs
      qface[dof, 1, n] = zero(Tsol)
      qface[dof, 2, n] = zero(Tsol)
    end

    for j = 1 : sbpface.stencilsize
      permL = sbpface.perm[j, face.faceL]
      permR = sbpface.perm[j, face.faceR]

      for dof = 1 : numDofs
        qface[dof, 1, n] += sbpface.interp[j, n] * qvolL[dof, permL]
        qface[dof, 2, n] += sbpface.interp[j, iR]* qvolR[dof, permR]
      end
    end
  end

  return nothing
end


@doc """
Another version of boundaryinterpolate.
Given Q on the parent element, interpolate Q to boundary face owned by parent element

Input:
sbpface        : face SBP operator
bndface        : the interface which we are interpolating Q onto
qvol        : Q at nodes in element L, qvolL(idof, iNode)    
Output:
qface        : Q at nodes on interface, qface(idof, L/R, ifacenode)

"""->
function boundaryinterpolate{Tsbp, Tsol}(sbpface::AbstractFace{Tsbp},
                                         bndface::Boundary,
                                         qvol::AbstractArray{Tsol, 2},
                                         qface::AbstractArray{Tsol, 2})

  @assert(size(qvol, 1) == size(qface, 1))            # dof
  @assert(size(qface,2) == sbpface.numnodes)
  @assert(size(sbpface.interp, 1) <= size(qvol, 2))    # stencilsize <= numNodesPerElem
  @assert(size(sbpface.interp, 2) == size(qface, 2))    # check number of face nodes

  numDofs = size(qvol, 1)    

  for n = 1 : sbpface.numnodes
    for dof = 1 : numDofs
      qface[dof, n] = zero(Tsol)
    end

    for j = 1 : sbpface.stencilsize
      perm = sbpface.perm[j, bndface.face]
      for dof = 1 : numDofs
        qface[dof, n] += sbpface.interp[j, n]*qvol[dof, perm]
      end
    end
  end
end

@doc """

Compute dynamic viscousity for an element using Sutherland's Law

Input:
temp        : temperature
Output:
rmu            : dynamic viscousity

"""->
function getMuK{Tsol}(temp::Tsol, rMuK::AbstractArray{Tsol, 1})
  @assert(length(rMuK) == 2)

  Tref = 460.0
  cstar = 198.6/Tref

  rMuK[1] = (1.0 + cstar)/(temp + cstar)*temp^1.5
  # rMuK[1] = 1.0
  rMuK[2] = 1.0

  return nothing
end
@doc """

Compute dynamic viscousity for an element using Sutherland's Law

Input:
temp        : temperature
Output:
rmu            : dynamic viscousity

"""->
function getMuK{Tsol}(temp::AbstractArray{Tsol, 1}, rMuK::AbstractArray{Tsol, 2})
  @assert(size(rMuK, 1) == 2)
  @assert(size(rMuK, 2) == length(temp))

  Tref = 460.0
  cstar = 198.6/Tref

  for n = 1:length(temp)
    rMuK[1, n] = (1.0 + cstar)/(temp[n] + cstar)*temp[n]^1.5
    # rMuK[1, n] = 1.0
    rMuK[2, n] = 1.0
  end

  return nothing
end


# @doc """

# Compute the viscous flux if we don't have derivatives yet
# Only work for 2D.
# Since derivatives are involved, it's difficult (or should say, expensive)
# to define node-based version.

# Input: 
# q        : conservative variable for a element
# dxidx    : derivatives of mapping
# jac        : determinant of mapping
# Output
# Fv        : viscous flux
# """->
#
function calcFvis_elem_1{Tsol, Tmsh}(sbp::AbstractSBP,
                                     q::AbstractArray{Tsol, 2},
                                     dxidx::AbstractArray{Tmsh, 3},
                                     jac::AbstractArray{Tmsh, 1},
                                     Fv::AbstractArray{Tsol, 3})
  @assert(size(q, 2) == sbp.numnodes)        # number of element nodes
  @assert(size(dxidx, 1) == size(dxidx, 2))  # dimension
  @assert(size(dxidx, 3) == size(q, 2))      # node
  @assert(size(Fv, 1) == size(dxidx, 2))     # dimension
  @assert(size(Fv, 2) == size(q, 1))         # dof
  @assert(size(Fv, 3) == size(q, 2))         # node

  Pr            = 0.72
  gamma        = 1.4
  gamma_1        = gamma - 1.0
  coef_nondim = 1.0/(Pr*gamma_1)
  two3rd = 2.0/3.0
  numNodes = size(dxidx, 3)
  dim = size(dxidx, 1)
  #
  # we need μ, κ, V, dVdx, dTdx
  #
  rMuK    = Array(Tsol, 2, numNodes)          # dynamic viscousity
  dVdx    = zeros(Tsol, dim, dim, numNodes)   # gradient of velocity, (velocity, dimension, node)
  dTdx    = zeros(Tsol, dim, numNodes)        # gradient of velocity, (dimension, node)
  tau = Array(Tsol, dim, dim, numNodes)       # stress    
  Dx  = Array(Tmsh, numNodes, numNodes, dim)  # derivative operators in physical domain

  #
  # primitive variables
  #
  T = Array(Tsol, numNodes)
  v = Array(Tsol, dim, numNodes)
  for n = 1:numNodes
    v[1, n] = q[2, n]/q[1, n]
    v[2, n] = q[3, n]/q[1, n]
    T[n] = gamma*gamma_1*(q[4, n]/q[1, n] - 0.5*(v[1, n]*v[1, n] + v[2, n]*v[2, n]))
  end

  # Compute viscousity and conductivity coefficients
  getMuK(T, rMuK)

  # First we compute derivatives, dudx, dudy, dTdx, dTdy
  calcDx(sbp, dxidx, jac, Dx)

  for d = 1:dim        # loop over dimensions, in which we are computing the derivatives
    for n = 1:numNodes    # loop over rows
      for col = 1:numNodes        # loop over columns
        dVdx[1, d, n] += Dx[n, col, d] * v[1, col]
        dVdx[2, d, n] += Dx[n, col, d] * v[2, col]
        dTdx[d, n]    += Dx[n, col, d] * T[col]
      end
    end
  end

  # compute viscous viscous stress, tau
  for n = 1:numNodes
    for d1 = 1:dim
      for d2 = 1:dim
        tau[d1, d2, n] = dVdx[d1, d2, n] + dVdx[d2, d1, n]
      end

      tau[d1, d1, n] -= two3rd*(dVdx[1, 1, n] + dVdx[2, 2, n])
    end
    tau[:, :, n] *= rMuK[1, n]
  end

  # start to compute viscous flux
  for n = 1:numNodes
    Fv[1, 1, n] = 0.0
    Fv[1, 2, n] = tau[1, 1, n]
    Fv[1, 3, n] = tau[2, 1, n]
    Fv[1, 4, n] = v[1, n]*tau[1, 1, n] + v[2, n]*tau[1, 2, n]
    Fv[1, 4, n] += rMuK[2, n]*dTdx[1, n]*coef_nondim  

    Fv[2, 1, n] = 0.0
    Fv[2, 2, n] = tau[1, 2, n]
    Fv[2, 3, n] = tau[2, 2, n]
    Fv[2, 4, n] = v[1, n]*tau[2, 1, n] + v[2, n]*tau[2, 2, n]
    Fv[2, 4, n] += rMuK[2, n]*dTdx[2, n]*coef_nondim 
  end

  return nothing
end

# @doc """

# Another version of calculating viscous flux. Fv = G(q)∇q
# Since    ∇q is available, we compute G(q) then do the multiplication

# Input:
# q        : conservative variable
# dqdx    : gradient of conservative variable
# Output:
# Fv    : viscous flux
# """->
function calcFvis{Tsol}(q::AbstractArray{Tsol, 2},
                        dqdx::AbstractArray{Tsol, 3},
                        Fv::AbstractArray{Tsol, 3})

  @assert(size(q, 2) == size(dqdx, 3)) # number of nodes per element
  @assert(size(q, 1) == size(dqdx, 2)) # number of dof per node
  @assert(size(dqdx, 1) == 2)

  dim      = size(dqdx, 1)
  numDofs  = size(q, 1)
  numNodes = size(q, 2)
  Gv = zeros(Tsol,  numDofs, numDofs, dim, dim, numNodes)

  calcDiffusionTensor(q, Gv)
  calcFvis(Gv, dqdx, Fv)

  return nothing
end


# @doc """

# Another version of calculating viscous flux. Fv = G(q)∇q
# In this one G(q) and ∇q are both available, so what we need to do 
# is just the multiplication

# Input:
# sbp        : sbp operator
# q        : viscousity tensor
# dqdx    : gradient of conservative variable
# Output:
# Fv        : viscous flux
# """->
function calcFvis{Tsol}(Gv::AbstractArray{Tsol, 5},
                        dqdx::AbstractArray{Tsol, 3},
                        Fv::AbstractArray{Tsol, 3})

  @assert(size(Gv, 5) == size(dqdx, 3)) # number of nodes per element
  @assert(size(Gv, 4) == size(dqdx, 1)) # spatial dimention 
  @assert(size(Gv, 3) == size(dqdx, 1)) # spatial dimention 
  @assert(size(Gv, 2) == size(dqdx, 2)) # number of dofs per element
  @assert(size(Gv, 1) == size(dqdx, 2)) # number of dofs per element

  dim      = size(dqdx, 1)
  numDofs  = size(Gv, 1)
  numNodes = size(Gv, 5)

  for n = 1 : numNodes
    # brutal loop
    # for d = 1 : dim
    # for d2 = 1 : dim
    # for row = 1 : numDofs
    # for col = 1 : numDofs
    # Fv[d, row, n] += Gv[row, col, d, d2, n]*dqdx[d2, col, n]
    # end
    # end
    # end
    # end

    Fv[1, 1, n] = 0.0
    Fv[1, 2, n] = (Gv[2, 1, 1, 1,n]*dqdx[1, 1, n] + Gv[2, 2, 1, 1,n]*dqdx[1, 2, n]
                   + Gv[2, 1, 1, 2,n]*dqdx[2, 1, n] + Gv[2, 3, 1, 2,n]*dqdx[2, 3, n] )
    Fv[1, 3, n] = (Gv[3, 1, 1, 1,n]*dqdx[1, 1, n] + Gv[3, 3, 1, 1,n]*dqdx[1, 3, n]
                   + Gv[3, 1, 1, 2,n]*dqdx[2, 1, n] + Gv[3, 2, 1, 2,n]*dqdx[2, 2, n] )
    Fv[1, 4, n] = (Gv[4, 1, 1, 1,n]*dqdx[1, 1, n] + Gv[4, 2, 1, 1,n]*dqdx[1, 2, n]
                   + Gv[4, 3, 1, 1,n]*dqdx[1, 3, n] + Gv[4, 4, 1, 1,n]*dqdx[1, 4, n] 
                   + Gv[4, 1, 1, 2,n]*dqdx[2, 1, n] + Gv[4, 2, 1, 2,n]*dqdx[2, 2, n] 
                   + Gv[4, 3, 1, 2,n]*dqdx[2, 3, n] )

    Fv[2, 1, n] = 0.0
    Fv[2, 2, n] = (Gv[2, 1, 2, 1,n]*dqdx[1, 1, n] + Gv[2, 3, 2, 1,n]*dqdx[1, 3, n]
                   + Gv[2, 1, 2, 2,n]*dqdx[2, 1, n] + Gv[2, 2, 2, 2,n]*dqdx[2, 2, n] )
    Fv[2, 3, n] = (Gv[3, 1, 2, 1,n]*dqdx[1, 1, n] + Gv[3, 2, 2, 1,n]*dqdx[1, 2, n]
                   + Gv[3, 1, 2, 2,n]*dqdx[2, 1, n] + Gv[3, 3, 2, 2,n]*dqdx[2, 3, n] )
    Fv[2, 4, n] = (Gv[4, 1, 2, 1,n]*dqdx[1, 1, n] + Gv[4, 2, 2, 1,n]*dqdx[1, 2, n]
                   + Gv[4, 3, 2, 1,n]*dqdx[1, 3, n] 
                   + Gv[4, 1, 2, 2,n]*dqdx[2, 1, n] + Gv[4, 2, 2, 2,n]*dqdx[2, 2, n] 
                   + Gv[4, 3, 2, 2,n]*dqdx[2, 3, n] + Gv[4, 4, 2, 2,n]*dqdx[2, 4, n] )
  end

  return nothing
end



@doc """

Compute the viscous flux. This is an element level function.

Input: 
sbp :
q : conservative variable at all element node
dxidx:
jac:
Fv : viscous flux
Output
jac = viscous flux jacobian at each node, dimension = (dim+2, dim+2, dim, dim, numNodes)
"""->
function calcFvis_elem{Tsol, Tmsh}(sbp::AbstractSBP,
                                   q::AbstractArray{Tsol, 2},
                                   dxidx::AbstractArray{Tmsh, 3},
                                   jac::AbstractArray{Tmsh, 1},
                                   Fv::AbstractArray{Tsol, 3})

  @assert(size(dxidx, 1) == size(dxidx, 2))    # dimension
  @assert(size(dxidx, 3) == size(q, 2))        # node
  @assert(length(jac) == sbp.numnodes)
  @assert(size(q, 2)  == sbp.numnodes)
  @assert(size(Fv, 1) == size(dxidx, 2))       # dimension
  @assert(size(Fv, 2) == size(q, 1))           # dof
  @assert(size(Fv, 3) == size(q, 2))           # node

  dim      = size(dxidx, 1)
  numNodes = size(dxidx, 3)
  numDofs  = size(q, 1)
  dqdx     = zeros(Tsol, dim, numDofs, numNodes)

  calcGradient(sbp, dxidx, jac, q, dqdx)

  calcFvis(q, dqdx, Fv)

  return nothing
end

@doc """

Compute the viscous diffusion matrix, only work for 2D

Input: 
q = conservative variable at a node
Output
jac = viscous flux jacobian at each node, dimension = (dim+2, dim+2, dim, dim, numNodes)
"""->
function calcDiffusionTensor{Tsol}(q::AbstractArray{Tsol, 2},
                                   Gv::AbstractArray{Tsol, 5})
  @assert(size(q, 2) == size(Gv, 5))
  @assert(size(Gv, 4) == 2)
  @assert(size(Gv, 3) == 2)
  @assert(size(Gv, 2) == 2+2)
  @assert(size(Gv, 1) == 2+2)
  numNodes = size(q, 2)
  gamma = 1.4
  gamma_1 = gamma - 1.0
  Pr = 0.72
  gamma_pr = gamma/Pr
  one3rd = 1.0/3.0
  two3rd = 2.0/3.0
  four3rd = 4.0/3.0
  E    = Array(Tsol, numNodes)
  T    = Array(Tsol, numNodes)
  v2   = Array(Tsol, numNodes)
  v1v2 = Array(Tsol, numNodes)
  v    = Array(Tsol, 2, numNodes)
  rmuk = Array(Tsol, 2, numNodes)

  for n = 1 : numNodes
    v[1, n] = q[2,n]/q[1,n]
    v[2, n] = q[3,n]/q[1,n]
    E[n]    = q[4,n]/q[1,n]
    v2[n]   = v[1,n]*v[1,n] + v[2,n]*v[2,n]
    v1v2[n] = v[1,n]*v[2,n]
    T[n]    = gamma*gamma_1 * (E[n] - 0.5*v2[n] )
  end

  getMuK(T, rmuk)

  for n = 1 : numNodes
    # Gv[1,1,1,1,n] = 0.0
    # Gv[1,2,1,1,n] = 0.0
    # Gv[1,3,1,1,n] = 0.0
    # Gv[1,4,1,1,n] = 0.0

    Gv[2,1,1,1,n] = -four3rd*v[1,n]
    Gv[2,2,1,1,n] = four3rd 
    # Gv[2,3,1,1,n] = 0.0
    # Gv[2,4,1,1,n] = 0.0

    Gv[3,1,1,1,n] = -v[2,n]
    # Gv[3,2,1,1,n] = 0.0
    Gv[3,3,1,1,n] = 1.0
    # Gv[3,4,1,1,n] = 0.0

    Gv[4,1,1,1,n] = -(four3rd*v[1,n]*v[1,n] + v[2,n]*v[2,n] + gamma_pr*(E[n] - v2[n]))
    Gv[4,2,1,1,n] = (four3rd - gamma_pr)*v[1,n]
    Gv[4,3,1,1,n] = (1.0 - gamma_pr)*v[2,n]
    Gv[4,4,1,1,n] = gamma_pr

    #
    # G12
    #
    # Gv[1,1,1,2,n] = 0.0
    # Gv[1,2,1,2,n] = 0.0
    # Gv[1,3,1,2,n] = 0.0
    # Gv[1,4,1,2,n] = 0.0

    Gv[2,1,1,2,n] = two3rd*v[2,n]
    # Gv[2,2,1,2,n] = 0.0 
    Gv[2,3,1,2,n] = -two3rd
    # Gv[2,4,1,2,n] = 0.0

    Gv[3,1,1,2,n] = -v[1,n]
    Gv[3,2,1,2,n] = 1.0
    # Gv[3,3,1,2,n] = 0.0
    # Gv[3,4,1,2,n] = 0.0

    Gv[4,1,1,2,n] = -one3rd*v1v2[n]
    Gv[4,2,1,2,n] = v[2,n]
    Gv[4,3,1,2,n] = -two3rd*v[1,n] 
    # Gv[4,4,1,2,n] = 0.0
    #
    # G21
    #
    # Gv[1,1,2,1,n] = 0.0
    # Gv[1,2,2,1,n] = 0.0
    # Gv[1,3,2,1,n] = 0.0
    # Gv[1,4,2,1,n] = 0.0

    Gv[2,1,2,1,n] = -v[2,n]
    # Gv[2,2,2,1,n] = 0.0 
    Gv[2,3,2,1,n] = 1.0
    # Gv[2,4,2,1,n] = 0.0

    Gv[3,1,2,1,n] = two3rd*v[1,n]
    Gv[3,2,2,1,n] = -two3rd
    # Gv[3,3,2,1,n] = 0.0
    # Gv[3,4,2,1,n] = 0.0

    Gv[4,1,2,1,n] = -one3rd*v1v2[n]
    Gv[4,2,2,1,n] = -two3rd*v[2,n]
    Gv[4,3,2,1,n] = v[1,n] 
    # Gv[4,4,2,1,n] = 0.0
    #
    # G22
    #
    # Gv[1,1,2,2,n] = 0.0
    # Gv[1,2,2,2,n] = 0.0
    # Gv[1,3,2,2,n] = 0.0
    # Gv[1,4,2,2,n] = 0.0

    Gv[2,1,2,2,n] = -v[1,n]
    Gv[2,2,2,2,n] = 1.0 
    # Gv[2,3,2,2,n] = 0.0
    # Gv[2,4,2,2,n] = 0.0

    Gv[3,1,2,2,n] = -four3rd*v[2,n]
    # Gv[3,2,2,2,n] = 0.0
    Gv[3,3,2,2,n] = four3rd
    # Gv[3,4,2,2,n] = 0.0

    Gv[4,1,2,2,n] = -(v[1,n]*v[1,n] + four3rd*v[2,n]*v[2,n] + gamma_pr*(E[n] - v2[n]))
    Gv[4,2,2,2,n] = (1.0 - gamma_pr)*v[1,n]
    Gv[4,3,2,2,n] = (four3rd - gamma_pr)*v[2,n]
    Gv[4,4,2,2,n] = gamma_pr

    coef = rmuk[1,n]/q[1,n]

    for iDof = 1: 4
      for jDof = 1 : 4
        Gv[iDof, jDof, 1, 1 ,n] *= coef 
        Gv[iDof, jDof, 1, 2 ,n] *= coef 
        Gv[iDof, jDof, 2, 1 ,n] *= coef 
        Gv[iDof, jDof, 2, 2 ,n] *= coef 
      end
    end
  end
  return nothing
end
@doc """

Compute the viscous diffusion matrix, only work for 2D

Input: 
q = conservative variable at a node
Output
jac = viscous flux jacobian at each node, dimension = (dim+2, dim+2, dim, dim, numNodes)
"""->
function calcDiffusionTensor_adiabaticWall{Tsol}(q::AbstractArray{Tsol, 2},
                                                 Gv::AbstractArray{Tsol, 5})
  @assert(size(q, 2) == size(Gv, 5))
  @assert(size(Gv, 4) == 2)
  @assert(size(Gv, 3) == 2)
  @assert(size(Gv, 2) == 2+2)
  @assert(size(Gv, 1) == 2+2)
  numNodes = size(q, 2)
  gamma = 1.4
  gamma_1 = gamma - 1.0
  Pr = 0.72
  gamma_pr = gamma/Pr
  #
  # set this coefficient to zero to enforce ∇T=0
  # TODO: remove everything in Gv directly.
  # 
  gamma_pr = 0.0

  one3rd = 1.0/3.0
  two3rd = 2.0/3.0
  four3rd = 4.0/3.0
  E    = Array(Tsol, numNodes)
  T    = Array(Tsol, numNodes)
  v2   = Array(Tsol, numNodes)
  v1v2 = Array(Tsol, numNodes)
  v    = Array(Tsol, 2, numNodes)
  rmuk = Array(Tsol, 2, numNodes)

  for n = 1 : numNodes
    v[1, n] = q[2,n]/q[1,n]
    v[2, n] = q[3,n]/q[1,n]
    E[n]    = q[4,n]/q[1,n]
    v2[n]   = v[1,n]*v[1,n] + v[2,n]*v[2,n]
    v1v2[n] = v[1,n]*v[2,n]
    T[n]    = gamma*gamma_1 * (E[n] - 0.5*v2[n] )
  end

  getMuK(T, rmuk)

  for n = 1 : numNodes
    # Gv[1,1,1,1,n] = 0.0
    # Gv[1,2,1,1,n] = 0.0
    # Gv[1,3,1,1,n] = 0.0
    # Gv[1,4,1,1,n] = 0.0

    Gv[2,1,1,1,n] = -four3rd*v[1,n]
    Gv[2,2,1,1,n] = four3rd 
    # Gv[2,3,1,1,n] = 0.0
    # Gv[2,4,1,1,n] = 0.0

    Gv[3,1,1,1,n] = -v[2,n]
    # Gv[3,2,1,1,n] = 0.0
    Gv[3,3,1,1,n] = 1.0
    # Gv[3,4,1,1,n] = 0.0

    Gv[4,1,1,1,n] = -(four3rd*v[1,n]*v[1,n] + v[2,n]*v[2,n] + gamma_pr*(E[n] - v2[n]))
    Gv[4,2,1,1,n] = (four3rd - gamma_pr)*v[1,n]
    Gv[4,3,1,1,n] = (1.0 - gamma_pr)*v[2,n]
    Gv[4,4,1,1,n] = gamma_pr

    #
    # G12
    #
    # Gv[1,1,1,2,n] = 0.0
    # Gv[1,2,1,2,n] = 0.0
    # Gv[1,3,1,2,n] = 0.0
    # Gv[1,4,1,2,n] = 0.0

    Gv[2,1,1,2,n] = two3rd*v[2,n]
    # Gv[2,2,1,2,n] = 0.0 
    Gv[2,3,1,2,n] = -two3rd
    # Gv[2,4,1,2,n] = 0.0

    Gv[3,1,1,2,n] = -v[1,n]
    Gv[3,2,1,2,n] = 1.0
    # Gv[3,3,1,2,n] = 0.0
    # Gv[3,4,1,2,n] = 0.0

    Gv[4,1,1,2,n] = -one3rd*v1v2[n]
    Gv[4,2,1,2,n] = v[2,n]
    Gv[4,3,1,2,n] = -two3rd*v[1,n] 
    # Gv[4,4,1,2,n] = 0.0
    #
    # G21
    #
    # Gv[1,1,2,1,n] = 0.0
    # Gv[1,2,2,1,n] = 0.0
    # Gv[1,3,2,1,n] = 0.0
    # Gv[1,4,2,1,n] = 0.0

    Gv[2,1,2,1,n] = -v[2,n]
    # Gv[2,2,2,1,n] = 0.0 
    Gv[2,3,2,1,n] = 1.0
    # Gv[2,4,2,1,n] = 0.0

    Gv[3,1,2,1,n] = two3rd*v[1,n]
    Gv[3,2,2,1,n] = -two3rd
    # Gv[3,3,2,1,n] = 0.0
    # Gv[3,4,2,1,n] = 0.0

    Gv[4,1,2,1,n] = -one3rd*v1v2[n]
    Gv[4,2,2,1,n] = -two3rd*v[2,n]
    Gv[4,3,2,1,n] = v[1,n] 
    # Gv[4,4,2,1,n] = 0.0
    #
    # G22
    #
    # Gv[1,1,2,2,n] = 0.0
    # Gv[1,2,2,2,n] = 0.0
    # Gv[1,3,2,2,n] = 0.0
    # Gv[1,4,2,2,n] = 0.0

    Gv[2,1,2,2,n] = -v[1,n]
    Gv[2,2,2,2,n] = 1.0 
    # Gv[2,3,2,2,n] = 0.0
    # Gv[2,4,2,2,n] = 0.0

    Gv[3,1,2,2,n] = -four3rd*v[2,n]
    # Gv[3,2,2,2,n] = 0.0
    Gv[3,3,2,2,n] = four3rd
    # Gv[3,4,2,2,n] = 0.0

    Gv[4,1,2,2,n] = -(v[1,n]*v[1,n] + four3rd*v[2,n]*v[2,n] + gamma_pr*(E[n] - v2[n]))
    Gv[4,2,2,2,n] = (1.0 - gamma_pr)*v[1,n]
    Gv[4,3,2,2,n] = (four3rd - gamma_pr)*v[2,n]
    Gv[4,4,2,2,n] = gamma_pr

    coef = rmuk[1,n]/q[1,n]
    for iDof = 1: 4
      for jDof = 1 : 4
        Gv[iDof, jDof, 1, 1 ,n] *= coef 
        Gv[iDof, jDof, 1, 2 ,n] *= coef 
        Gv[iDof, jDof, 2, 1 ,n] *= coef 
        Gv[iDof, jDof, 2, 2 ,n] *= coef 
      end
    end
  end
  return nothing
end


@doc """

compute the adiabatic wall boundary value

Input:
q_in    :: conservative variables on boundary face
norm    :: outward unit normal of the boundary face
Output:
q_bnd    :: the boundary value

"""->
type AdiabaticWall <: AbstractBoundaryValueType
end
function call{Tsol, Tmsh}(obj::AdiabaticWall, 
                          q_in::AbstractArray{Tsol, 2},
                          norm::AbstractArray{Tmsh, 2},
                          params::ParamType{2},
                          q_bnd::AbstractArray{Tsol, 2})
  @assert( size(q_in, 1) == size(q_bnd,  1))
  @assert( size(q_in, 2) == size(q_bnd,  2))
  @assert( size(q_in, 2) == size(norm,  2))

  dim = size(norm, 1)
  numNodes = size(q_in, 2)

  for n = 1 : numNodes
    q_bnd[1, n] = q_in[1, n]
    # noslip condition
    q_bnd[2:dim+1, n] = 0.0
    q_bnd[dim+2, n] = q_in[4, n]
    # rhoV2 = (q_in[2,n]*q_in[2,n] + q_in[3,n]*q_in[3,n])/q_in[1,n]
    # q_bnd[2:dim+1, n] = 0.0
    # q_bnd[dim+2, n] = q_in[4, n] - 0.5*rhoV2
  end

  return nothing
end

@doc """

compute the farfield boundary value

Input:
q_in    :: conservative variables on boundary face
norm    :: outward unit normal of the boundary face
Output:
q_bnd    :: the boundary value

"""->
type Farfield <: AbstractBoundaryValueType
end

function call{Tsol, Tmsh}(obj::Farfield,
                          q_in::AbstractArray{Tsol, 2},
                          norm::AbstractArray{Tmsh, 2},
                          params::ParamType{2},
                          q_bnd::AbstractArray{Tsol, 2})
  @assert( size(q_in, 1) == size(q_bnd,  1))
  @assert( size(q_in, 2) == size(q_bnd,  2))
  @assert( size(q_in, 2) == size(norm,  2))

  dim      = size(norm, 1)
  numNodes = size(q_in, 2)
  lambda   = zeros(Float64, 3) 
  gamma    = params.gamma
  gamma_1  = params.gamma_1
  aoa      = params.aoa
  gg_1     = (gamma*gamma_1)
  MaInf    = params.Ma    

  # freestream values
  qInf = zeros(Float64, dim + 2)
  qInf[1] = 1.0
  qInf[2] = qInf[1]*MaInf*cos(aoa)
  qInf[3] = qInf[1]*MaInf*sin(aoa)
  qInf[4] = qInf[1]*(1.0/gg_1 + 0.5*MaInf*MaInf)

  for n = 1 : numNodes
    #
    # contravariant velocity
    #
    vn = q_in[2, n]*norm[1, n] + q_in[3, n]*norm[2, n]
    vn = vn/q_in[1,n]
    v2 = q_in[2,n]*q_in[2,n] + q_in[3,n]*q_in[3,n]
    v2 /= q_in[1,n]*q_in[1,n]
    T = gg_1*(q_in[4,n]/q_in[1,n] - 0.5*v2)
    a = sqrt(T)
    #
    # eigenvalues
    # 
    lambda[1] = real(vn)
    lambda[2] = real(vn + a)
    lambda[3] = real(vn - a)

    if lambda[2] <= 0.0            # supersonic inflow
      q_bnd[:, n] = qInf[:]
    elseif lambda[3] >= 0.0        # supersonic outflow
      q_bnd[:, n] = q_in[:, n]
    elseif lambda[1] <= 0.0        # subsonic inflow
      p = q_in[1,n]*T/gamma 
      q_bnd[1:dim+1, n] = qInf[1:dim+1]
      q_bnd[dim+2, n] = p/gamma_1 + 0.5*MaInf*MaInf*qInf[1]
    else                        # subsonic outflow
      pInf = 1.0/gamma
      q_bnd[1:dim+1, n] = q_in[1:dim+1, n]
      rhoV2 = q_in[1, n] * v2
      q_bnd[dim+2, n] = pInf/gamma_1 + 0.5*rhoV2
    end    

    # q_bnd[:,n] = qInf[:]
  end

  return nothing
end

