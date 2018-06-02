@doc """

Compute dynamic viscosity for an element using Sutherland's Law

Input:
  temp        : temperature
Output:
  rmu         : dynamic viscosity

"""->
function getMuK(temp::Tsol, rMuK::AbstractArray{Tsol, 1}) where Tsol
  @assert(length(rMuK) == 2)

  Tref = 460.0
  cstar = 198.6/Tref

  # rMuK[1] = (1.0 + cstar)/(temp + cstar)*temp^1.5
  rMuK[1] = 1.0
  rMuK[2] = 1.0

  return nothing
end

@doc """
Compute dynamic viscosity for an element using Sutherland's Law

Input:
  temp        : temperature
Output:
  rmu         : dynamic viscosity

"""->
function getMuK(temp::AbstractArray{Tsol, 1}, rMuK::AbstractArray{Tsol, 2}) where Tsol
  @assert(size(rMuK, 1) == 2)
  @assert(size(rMuK, 2) == length(temp))

  Tref = 460.0
  cstar = 198.6/Tref

  for n = 1:length(temp)
    # rMuK[1, n] = (1.0 + cstar)/(temp[n] + cstar)*temp[n]^1.5
    rMuK[1, n] = 1.0
    rMuK[2, n] = 1.0
  end

  return nothing
end


# @doc """

# Compute the viscous flux if we don't have derivatives yet.
# Only work for 2D.
# Since derivatives are involved, it's difficult (or should say, expensive)
# to define node-based version.

# Input: 
#   q        : conservative variable for a element
#   dxidx    : derivatives of element mapping
#   jac      : determinant of element mapping
# Output
#   Fv       : viscous flux
# """->
#
function calcFvis_elem_direct(params::ParamType{2, :conservative},
                              sbp::AbstractSBP,
                              q::AbstractArray{Tsol, 2},
                              dxidx::AbstractArray{Tmsh, 3},
                              jac::AbstractArray{Tmsh, 1},
                              Fv::AbstractArray{Tsol, 3}) where {Tsol, Tmsh}
  @assert(size(q, 2) == sbp.numnodes)        # number of element nodes
  @assert(size(dxidx, 1) == size(dxidx, 2))  # dimension
  @assert(size(dxidx, 3) == size(q, 2))      # node
  @assert(size(dxidx, 3) == sbp.numnodes)      # node
  @assert(size(Fv, 1) == size(dxidx, 2))     # dimension
  @assert(size(Fv, 2) == size(q, 1))         # dof
  @assert(size(Fv, 3) == size(q, 2))         # node

  dim = 2
  Pr      = 0.72
  gamma   = params.gamma
  gamma_1 = gamma - 1.0
  coef_nondim = 1.0/(Pr*gamma_1)
  two3rd   = 2.0/3.0
  four3rd  = 4.0/3.0
  numNodes = size(dxidx, 3)
  #
  # we need μ, κ, V, dVdx, dTdx
  #
  rMuK = Array(Tsol, 2, numNodes)          # dynamic viscosity
  dVdx = zeros(Tsol, dim, dim, numNodes)   # gradient of velocity, (velocity, dimension, node)
  dTdx = zeros(Tsol, dim, numNodes)        # gradient of velocity, (dimension, node)
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

  # Compute viscosity and conductivity coefficients
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
    # for d1 = 1:dim
      # for d2 = 1 : dim
        # tau[d1, d2, n] = dVdx[d1, d2, n] + dVdx[d2, d1, n]
      # end

      # tau[d1, d1, n] -= two3rd*(dVdx[1, 1, n] + dVdx[2, 2, n])
      # for d2 = 1 : dim
        # tau[d1, d2, n] *= rMuK[1,n]
      # end
    # end

    tau[1,1,n] = rMuK[1,n] * two3rd *(2 * dVdx[1,1,n] - dVdx[2,2,n])
    tau[1,2,n] = rMuK[1,n] * (dVdx[1,2,n] + dVdx[2,1,n])
    tau[2,2,n] = rMuK[1,n] * two3rd *(2 * dVdx[2,2,n] - dVdx[1,1,n])
  end

  # start to compute viscous flux
  for n = 1:numNodes
    # Fv[1, 1, n] = 0.0
    Fv[1, 2, n] = tau[1, 1, n]
    Fv[1, 3, n] = tau[1, 2, n]
    Fv[1, 4, n] = v[1, n]*tau[1, 1, n] + v[2, n]*tau[1, 2, n]
    Fv[1, 4, n] += rMuK[2, n] * dTdx[1, n] * coef_nondim  

    # Fv[2, 1, n] = 0.0
    Fv[2, 2, n] = tau[1, 2, n]
    Fv[2, 3, n] = tau[2, 2, n]
    Fv[2, 4, n] = v[1, n]*tau[1, 2, n] + v[2, n]*tau[2, 2, n]
    Fv[2, 4, n] += rMuK[2, n] * dTdx[2, n] * coef_nondim 
  end

  return nothing
end

function calcFvis_elem_direct(params::ParamType{3, :conservative},
                              sbp::AbstractSBP,
                              q::AbstractArray{Tsol, 2},
                              dxidx::AbstractArray{Tmsh, 3},
                              jac::AbstractArray{Tmsh, 1},
                              Fv::AbstractArray{Tsol, 3}) where {Tsol, Tmsh}
  @assert(size(q, 2) == sbp.numnodes)        # number of element nodes
  @assert(size(dxidx, 1) == size(dxidx, 2))  # dimension
  @assert(size(dxidx, 3) == size(q, 2))      # node
  @assert(size(Fv, 1) == size(dxidx, 2))     # dimension
  @assert(size(Fv, 2) == size(q, 1))         # dof
  @assert(size(Fv, 3) == size(q, 2))         # node

  Pr          = 0.72
  gamma       = 1.4
  gamma_1     = gamma - 1.0
  coef_nondim = 1.0/(Pr*gamma_1)
  two3rd  = 2.0/3.0
  four3rd = 4.0/3.0
  numNodes = size(dxidx, 3)
  dim = 3
  #
  # we need μ, κ, V, dVdx, dTdx
  #
  rMuK = Array(Tsol, 2, numNodes)          # dynamic viscosity
  dVdx = zeros(Tsol, dim, dim, numNodes)   # gradient of velocity, (velocity, dimension, node)
  dTdx = zeros(Tsol, dim, numNodes)        # gradient of velocity, (dimension, node)
  tau  = Array(Tsol, dim, dim, numNodes)       # stress    
  Dx   = Array(Tmsh, numNodes, numNodes, dim)  # derivative operators in physical domain

  #
  # primitive variables
  #
  T = Array(Tsol, numNodes)
  v = Array(Tsol, dim, numNodes)
  for n = 1:numNodes
    v[1, n] = q[2, n]/q[1, n]
    v[2, n] = q[3, n]/q[1, n]
    v[3, n] = q[4, n]/q[1, n]
    vel2 = (v[1,n]*v[1,n] + v[2,n]*v[2,n] + v[3,n]*v[3,n])
    T[n] = gamma*gamma_1*(q[5, n]/q[1, n] - 0.5*vel2)
  end

  # Compute viscosity and conductivity coefficients
  getMuK(T, rMuK)

  # First we compute derivatives, dudx, dudy, dTdx, dTdy
  calcDx(sbp, dxidx, jac, Dx)

  for d = 1:dim        # loop over dimensions, in which we are computing the derivatives
    for n = 1:numNodes    # loop over rows
      for col = 1:numNodes        # loop over columns
        dVdx[1, d, n] += Dx[n, col, d] * v[1, col]
        dVdx[2, d, n] += Dx[n, col, d] * v[2, col]
        dVdx[3, d, n] += Dx[n, col, d] * v[3, col]
        dTdx[d, n]    += Dx[n, col, d] * T[col]
      end
    end
  end

  # compute viscous viscous stress, tau
  for n = 1:numNodes
    tau[1,1,n] = rMuK[1,n] * two3rd * (2 * dVdx[1,1,n] - dVdx[2,2,n] - dVdx[3,3,n])
    tau[2,2,n] = rMuK[1,n] * two3rd * (2 * dVdx[2,2,n] - dVdx[1,1,n] - dVdx[3,3,n])
    tau[3,3,n] = rMuK[1,n] * two3rd * (2 * dVdx[3,3,n] - dVdx[1,1,n] - dVdx[2,2,n])
    tau[1,2,n] = rMuK[1,n] * (dVdx[1,2,n] + dVdx[2,1,n])
    tau[1,3,n] = rMuK[1,n] * (dVdx[1,3,n] + dVdx[3,1,n])
    tau[2,3,n] = rMuK[1,n] * (dVdx[2,3,n] + dVdx[3,2,n])
  end

  # start to compute viscous flux
  for n = 1:numNodes
    # Fv[1, 1, n] = 0.0
    Fv[1, 2, n] = tau[1, 1, n]
    Fv[1, 3, n] = tau[1, 2, n]
    Fv[1, 4, n] = tau[1, 3, n]
    Fv[1, 5, n] = v[1, n]*tau[1, 1, n] + v[2, n]*tau[1, 2, n] + v[3,n]*tau[1,3,n]
    Fv[1, 5, n] += rMuK[2, n] * dTdx[1, n] * coef_nondim  

    # Fv[2, 1, n] = 0.0
    Fv[2, 2, n] = tau[1, 2, n]
    Fv[2, 3, n] = tau[2, 2, n]
    Fv[2, 4, n] = tau[2, 3, n]
    Fv[2, 5, n] = v[1, n]*tau[1, 2, n] + v[2, n]*tau[2, 2, n] + v[3,n]*tau[2,3,n]
    Fv[2, 5, n] += rMuK[2, n] * dTdx[2, n] * coef_nondim 

    # Fv[2, 1, n] = 0.0
    Fv[3, 2, n] = tau[1, 3, n]
    Fv[3, 3, n] = tau[2, 3, n]
    Fv[3, 4, n] = tau[3, 3, n]
    Fv[3, 5, n] = v[1, n]*tau[1, 3, n] + v[2, n]*tau[2, 3, n] + v[3,n]*tau[3,3,n]
    Fv[3, 5, n] += rMuK[2, n] * dTdx[3, n] * coef_nondim 
  end

  return nothing
end


"""
  Similar to [`calcBndryFvis`](@ref);
  Compute the viscous flux on an interface. There are possibly three ways:
    1). interpolate the solution and gradient operator from element to interface, and then compute the flux.
    2). compute the elememt flux directly, and then interpolate the flux instead of solution to interface.
    3). compute the elememt flux throught F = G∇q, and then interpolate the flux to interface.
  The experiments show subtle difference between three approaches. But the theoretical
  proof in the journal paper (Jianfeng's journal paper on SBP parabolic) uses the third one.

  **Input**
   * params
   * qL : 
   * qR : 
   * dxidxL    : derivatives of left element mapping
   * jacL      : determinant of left element mapping
   * dxidxR    : derivatives of right element mapping
   * jacR      : determinant of right element mapping
   * interface: interface

  **Input/Output**
   * FvL : viscous flux on left side of interface
   * FvR : viscous flux on left side of interface

"""
function calcFaceFvis(params::ParamType{Tdim, :conservative},
                      sbp::AbstractSBP,
                      sbpface::AbstractFace,
                      qL::AbstractArray{Tsol, 2},
                      qR::AbstractArray{Tsol, 2},
                      dxidxL::AbstractArray{Tmsh, 3},
                      jacL::AbstractArray{Tmsh, 1},
                      dxidxR::AbstractArray{Tmsh, 3},
                      jacR::AbstractArray{Tmsh, 1},
                      face::Interface,
                      Fv_face::AbstractArray{Tsol, 4}) where {Tsol, Tmsh, Tdim}

  #
  # method-1 NOT FINISHED!!!
  #
  # elemL = face.elementL
  # elemR = face.elementL
  # dqdx_eL = Array(Tsol, Tdim, numDofs, numNodes)
  # dqdx_eR = Array(Tsol, Tdim, numDofs, numNodes)
  # calcGradient(mesh, sbp, elemL, qL, dqdx_eL)
  # calcGradient(mesh, sbp, elemR, qR, dqdx_eR)
  # for d = 1 : Tdim
    # dqdxL = Base.view(dqdx_elemL, d, :, :)
    # dqdxR = Base.view(dqdx_elemR, d, :, :)
    # dqdx_f = Base.view(dqdx_face, d, :, :, :)
    # interiorfaceinterpolate(sbpface, face, dqdxL, dqdxR, dqdx_f)
  # end
  # # Now both G and dqdx are avaiable at face nodes  
  # dqdx_faceL = Base.view(dqdx_face, :, :, 1, :)
  # dqdx_faceR = Base.view(dqdx_face, :, :, 2, :)
  # Fv_faceL = Base.view(Fv_face, :, :, 1, :)
  # Fv_faceR = Base.view(Fv_face, :, :, 2, :)
  # calcFvis(params, GtL, dqdx_faceL, Fv_faceL)
  # calcFvis(params, GtR, dqdx_faceR, Fv_faceR)
  # return nothing

  numDofs = size(qL, 1)
  numNodes = size(qL, 2)
  Fv_eL = zeros(Tsol, Tdim, numDofs, numNodes)
  Fv_eR = zeros(Tsol, Tdim, numDofs, numNodes)
  #
  # method-2
  #
  # calcFvis_elem_direct(params, sbp, qL, dxidxL, jacL, Fv_eL)
  # calcFvis_elem_direct(params, sbp, qR, dxidxR, jacR, Fv_eR)
  # for d = 1 : Tdim
    # Fv_eL_d = Base.view(Fv_eL, d, :, :)
    # Fv_eR_d = Base.view(Fv_eR, d, :, :)
    # Fv_face_d = Base.view(Fv_face, d, :, :, :)
    # interiorfaceinterpolate(sbpface, face, Fv_eL_d, Fv_eR_d, Fv_face_d)
  # end

  #
  # method-3
  #
  calcFvis_elem(params, sbp, qL, dxidxL, jacL, Fv_eL)
  calcFvis_elem(params, sbp, qR, dxidxR, jacR, Fv_eR)

  #
  # TODO: we can combine the first 2 dimension as a single 
  # dimension, then we will not need slice here any more.
  #
  for d = 1 : Tdim
    Fv_eL_d = Base.view(Fv_eL, d, :, :)
    Fv_eR_d = Base.view(Fv_eR, d, :, :)
    Fv_face_d = Base.view(Fv_face, d, :, :, :)
    interiorfaceinterpolate(sbpface, face, Fv_eL_d, Fv_eR_d, Fv_face_d)
  end

  return nothing
end

# @doc """

# Another version of calculating viscous flux. Fv = G(q)∇q
# Since ∇q is available, we compute G(q) then do the multiplication

# Input:
#   q       : conservative variable
#   dqdx    : gradient of conservative variable
# Output:
#   Fv    : viscous flux
# """->
function calcFvis(params::ParamType{Tdim, :conservative},
                  q::AbstractArray{Tsol, 2},
                  dqdx::AbstractArray{Tsol, 3},
                  Fv::AbstractArray{Tsol, 3}) where {Tsol, Tdim}

  @assert(size(q, 2) == size(dqdx, 3)) # number of nodes per element
  @assert(size(q, 1) == size(dqdx, 2)) # number of dof per node
  @assert(size(dqdx, 1) == Tdim)

  numDofs  = size(q, 1)
  numNodes = size(q, 2)
  Gv = zeros(Tsol,  numDofs, numDofs, Tdim, Tdim, numNodes)

  calcDiffusionTensor(params, q, Gv)
  calcFvis(params, Gv, dqdx, Fv)

  return nothing
end


# @doc """

# Another version of calculating viscous flux. Fv = G(q)∇q
# In this one G(q) and ∇q are both available, so what we need to do 
# is just the multiplication

# Input:
#   sbp     : sbp operator
#   q       : viscosity tensor
#   dqdx    : gradient of conservative variable
# Output:
#   Fv      : viscous flux
# """->
function calcFvis(params::ParamType{2, :conservative},
                  Gv::AbstractArray{Tsol, 5},
                  dqdx::AbstractArray{Tsol, 3},
                  Fv::AbstractArray{Tsol, 3}) where Tsol

  @assert(size(Gv, 5) == size(dqdx, 3)) # number of nodes per element
  @assert(size(Gv, 4) == size(dqdx, 1)) # spatial dimension
  @assert(size(Gv, 3) == size(dqdx, 1)) # spatial dimension
  @assert(size(Gv, 2) == size(dqdx, 2)) # number of dofs per element
  @assert(size(Gv, 1) == size(dqdx, 2)) # number of dofs per element

  dim      = size(dqdx, 1)
  numDofs  = size(Gv, 1)
  numNodes = size(Gv, 5)

  for n = 1 : numNodes

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
                 + Gv[4, 3, 2, 1,n]*dqdx[1, 3, n] + Gv[4, 4, 2, 1,n]*dqdx[1, 4, n] 
                 + Gv[4, 1, 2, 2,n]*dqdx[2, 1, n] + Gv[4, 2, 2, 2,n]*dqdx[2, 2, n] 
                 + Gv[4, 3, 2, 2,n]*dqdx[2, 3, n] + Gv[4, 4, 2, 2,n]*dqdx[2, 4, n] )

    # DEBUG ONLY
    # brutal loop
    # Fv0 = zeros(Tsol, 2, 4)
    # for d = 1 : dim
      # for d2 = 1 : dim
        # for row = 1 : numDofs
          # for col = 1 : numDofs
            # Fv0[d, row] += Gv[row, col, d, d2, n]*dqdx[d2, col, n]
          # end
        # end
      # end
    # end
    # diff = maximum(abs(real(Fv0 - ro_sview(Fv, :,:,n))))
    # if diff > 1.e-14
      # println(Fv0)
      # println(Fv)
    # end
    # DEBUG END
  end

  return nothing
end
function calcFvis(params::ParamType{3, :conservative},
                  Gv::AbstractArray{Tsol, 5},
                  dqdx::AbstractArray{Tsol, 3},
                  Fv::AbstractArray{Tsol, 3}) where Tsol

  @assert(size(Gv, 5) == size(dqdx, 3)) # number of nodes per element
  @assert(size(Gv, 4) == size(dqdx, 1)) # spatial dimention 
  @assert(size(Gv, 3) == size(dqdx, 1)) # spatial dimention 
  @assert(size(Gv, 2) == size(dqdx, 2)) # number of dofs per element
  @assert(size(Gv, 1) == size(dqdx, 2)) # number of dofs per element

  dim      = size(dqdx, 1)
  numDofs  = size(Gv, 1)
  numNodes = size(Gv, 5)

  for n = 1 : numNodes

    Fv[1, 1, n] = 0.0

    Fv[1, 2, n] = (Gv[2, 1, 1, 1,n]*dqdx[1, 1, n] + Gv[2, 2, 1, 1,n]*dqdx[1, 2, n]
                 + Gv[2, 1, 1, 2,n]*dqdx[2, 1, n] + Gv[2, 3, 1, 2,n]*dqdx[2, 3, n]
                 + Gv[2, 1, 1, 3,n]*dqdx[3, 1, n] + Gv[2, 4, 1, 3,n]*dqdx[3, 4, n] )

    Fv[1, 3, n] = (Gv[3, 1, 1, 1,n]*dqdx[1, 1, n] + Gv[3, 3, 1, 1,n]*dqdx[1, 3, n]
                 + Gv[3, 1, 1, 2,n]*dqdx[2, 1, n] + Gv[3, 2, 1, 2,n]*dqdx[2, 2, n] )

    Fv[1, 4, n] = (Gv[4, 1, 1, 1,n]*dqdx[1, 1, n] + Gv[4, 4, 1, 1,n]*dqdx[1, 4, n]
                 + Gv[4, 1, 1, 3,n]*dqdx[3, 1, n] + Gv[4, 2, 1, 3,n]*dqdx[3, 2, n] )


    Fv[1, 5, n] = (Gv[5, 1, 1, 1,n]*dqdx[1, 1, n] + Gv[5, 2, 1, 1,n]*dqdx[1, 2, n]
                 + Gv[5, 3, 1, 1,n]*dqdx[1, 3, n] + Gv[5, 4, 1, 1,n]*dqdx[1, 4, n] + Gv[5, 5, 1, 1,n]*dqdx[1, 5, n]
                 + Gv[5, 1, 1, 2,n]*dqdx[2, 1, n] + Gv[5, 2, 1, 2,n]*dqdx[2, 2, n] + Gv[5, 3, 1, 2,n]*dqdx[2, 3, n]
                 + Gv[5, 1, 1, 3,n]*dqdx[3, 1, n] + Gv[5, 2, 1, 3,n]*dqdx[3, 2, n] + Gv[5, 4, 1, 3,n]*dqdx[3, 4, n] )

    Fv[2, 1, n] = 0.0

    Fv[2, 2, n] = (Gv[2, 1, 2, 1,n]*dqdx[1, 1, n] + Gv[2, 3, 2, 1,n]*dqdx[1, 3, n]
                 + Gv[2, 1, 2, 2,n]*dqdx[2, 1, n] + Gv[2, 2, 2, 2,n]*dqdx[2, 2, n] )

    Fv[2, 3, n] = (Gv[3, 1, 2, 1,n]*dqdx[1, 1, n] + Gv[3, 2, 2, 1,n]*dqdx[1, 2, n]
                 + Gv[3, 1, 2, 2,n]*dqdx[2, 1, n] + Gv[3, 3, 2, 2,n]*dqdx[2, 3, n] 
                 + Gv[3, 1, 2, 3,n]*dqdx[3, 1, n] + Gv[3, 4, 2, 3,n]*dqdx[3, 4, n] )

    Fv[2, 4, n] = (Gv[4, 1, 2, 2,n]*dqdx[2, 1, n] + Gv[4, 4, 2, 2,n]*dqdx[2, 4, n] 
                 + Gv[4, 1, 2, 3,n]*dqdx[3, 1, n] + Gv[4, 3, 2, 3,n]*dqdx[3, 3, n] )

    Fv[2, 5, n] = (Gv[5, 1, 2, 1,n]*dqdx[1, 1, n] + Gv[5, 2, 2, 1,n]*dqdx[1, 2, n] + Gv[5, 3, 2, 1,n]*dqdx[1, 3, n] 
                 + Gv[5, 1, 2, 2,n]*dqdx[2, 1, n] + Gv[5, 2, 2, 2,n]*dqdx[2, 2, n] 
                 + Gv[5, 3, 2, 2,n]*dqdx[2, 3, n] + Gv[5, 4, 2, 2,n]*dqdx[2, 4, n] + Gv[5, 5, 2, 2,n]*dqdx[2, 5, n]
                 + Gv[5, 1, 2, 3,n]*dqdx[3, 1, n] + Gv[5, 3, 2, 3,n]*dqdx[3, 3, n] + Gv[5, 4, 2, 3,n]*dqdx[3, 4, n] )

    Fv[3, 1, n] = 0.0

    Fv[3, 2, n] = (Gv[2, 1, 3, 1,n]*dqdx[1, 1, n] + Gv[2, 4, 3, 1,n]*dqdx[1, 4, n]
                 + Gv[2, 1, 3, 3,n]*dqdx[3, 1, n] + Gv[2, 2, 3, 3,n]*dqdx[3, 2, n] )

    Fv[3, 3, n] = (Gv[3, 1, 3, 2,n]*dqdx[2, 1, n] + Gv[3, 4, 3, 2,n]*dqdx[2, 4, n] 
                 + Gv[3, 1, 3, 3,n]*dqdx[3, 1, n] + Gv[3, 3, 3, 3,n]*dqdx[3, 3, n] )

    Fv[3, 4, n] = (Gv[4, 1, 3, 1,n]*dqdx[1, 1, n] + Gv[4, 2, 3, 1,n]*dqdx[1, 2, n] 
                 + Gv[4, 1, 3, 2,n]*dqdx[2, 1, n] + Gv[4, 3, 3, 2,n]*dqdx[2, 3, n] 
                 + Gv[4, 1, 3, 3,n]*dqdx[3, 1, n] + Gv[4, 4, 3, 3,n]*dqdx[3, 4, n] )

    Fv[3, 5, n] = (Gv[5, 1, 3, 1,n]*dqdx[1, 1, n] + Gv[5, 2, 3, 1,n]*dqdx[1, 2, n] + Gv[5, 4, 3, 1,n]*dqdx[1, 4, n] 
                 + Gv[5, 1, 3, 2,n]*dqdx[2, 1, n] + Gv[5, 3, 3, 2,n]*dqdx[2, 3, n] + Gv[5, 4, 3, 2,n]*dqdx[2, 4, n]
                 + Gv[5, 1, 3, 3,n]*dqdx[3, 1, n] + Gv[5, 2, 3, 3,n]*dqdx[3, 2, n] + Gv[5, 3, 3, 3,n]*dqdx[3, 3, n] 
                 + Gv[5, 4, 3, 3,n]*dqdx[3, 4, n] + Gv[5, 5, 3, 3,n]*dqdx[3, 5, n]  )
    # DEBUG ONLY
    # brutal loop
    # Fv0 = zeros(Tsol, 3, 5)
    # for d = 1 : dim
      # for d2 = 1 : dim
        # for row = 1 : numDofs
          # for col = 1 : numDofs
            # Fv0[d, row] += Gv[row, col, d, d2, n]*dqdx[d2, col, n]
          # end
        # end
      # end
    # end
    # diff = maximum(abs(real(Fv0 - ro_sview(Fv, :,:,n))))
    # if diff > 1.e-14
      # println(Fv0)
      # println(Fv)
    # end
    # DEBUG END
  end

  return nothing
end


@doc """

Compute the viscous flux. This is an element level function.

Input: 
  sbp   :
  q     : conservative variable at all element node
  dxidx :
  jac   :
  Fv    : viscous flux
Output
  jac = viscous flux jacobian at each node, dimension = (dim+2, dim+2, dim, dim, numNodes)
"""->
function calcFvis_elem(params::ParamType,
                       sbp::AbstractSBP,
                       q::AbstractArray{Tsol, 2},
                       dxidx::AbstractArray{Tmsh, 3},
                       jac::AbstractArray{Tmsh, 1},
                       Fv::AbstractArray{Tsol, 3}) where {Tsol, Tmsh}

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

  calcFvis(params, q, dqdx, Fv)

  return nothing
end

@doc """

Compute the viscous diffusion matrix. 2D version

Input: 
q = conservative variable at a node
Output
jac = viscous flux jacobian at each node, dimension = (dim+2, dim+2, dim, dim, numNodes)
"""->
function calcDiffusionTensor(params::ParamType{2, :conservative},
                             q::AbstractArray{Tsol, 2},
                             Gv::AbstractArray{Tsol, 5}) where Tsol
  @assert(size(q, 2) == size(Gv, 5))
  @assert(size(Gv, 4) == 2)
  @assert(size(Gv, 3) == 2)
  @assert(size(Gv, 2) == 2+2)
  @assert(size(Gv, 1) == 2+2)
  numNodes = size(q, 2)
  gamma = params.gamma
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

Compute the viscous diffusion matrix. 3D version TODO: check w/ JF

Input: 
q = conservative variable at a node
Output
jac = viscous flux jacobian at each node, dimension = (dim+2, dim+2, dim, dim, numNodes)
"""->
function calcDiffusionTensor(params::ParamType{3, :conservative},
                             q::AbstractArray{Tsol, 2},
                             Gv::AbstractArray{Tsol, 5}) where Tsol
  Tdim = 3
  @assert(size(q, 2) == size(Gv, 5))
  @assert(size(Gv, 4) == Tdim)
  @assert(size(Gv, 3) == Tdim)
  @assert(size(Gv, 2) == Tdim+2)
  @assert(size(Gv, 1) == Tdim+2)
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
  v1v3 = Array(Tsol, numNodes)
  v2v3 = Array(Tsol, numNodes)
  v    = Array(Tsol, 3, numNodes)
  rmuk = Array(Tsol, 2, numNodes)

  for n = 1 : numNodes
    v[1, n] = q[2,n]/q[1,n]
    v[2, n] = q[3,n]/q[1,n]
    v[3, n] = q[4,n]/q[1,n]
    v2[n]   = v[1,n]*v[1,n] + v[2,n]*v[2,n] + v[3,n]*v[3,n]
    v1v2[n] = v[1,n]*v[2,n]
    v2v3[n] = v[2,n]*v[3,n]
    v1v3[n] = v[1,n]*v[3,n]
    E[n]    = q[5,n]/q[1,n]
    T[n]    = gamma*gamma_1 * (E[n] - 0.5*v2[n] )
  end

  getMuK(T, rmuk)

  for n = 1 : numNodes
    coef = rmuk[1,n]/q[1,n]

    # Gv[1,1,1,1,n] = 0.0
    # Gv[1,2,1,1,n] = 0.0
    # Gv[1,3,1,1,n] = 0.0
    # Gv[1,4,1,1,n] = 0.0
    # Gv[1,5,1,1,n] = 0.0

    Gv[2,1,1,1,n] = -four3rd*v[1,n] * coef
    Gv[2,2,1,1,n] = four3rd  * coef
    # Gv[2,3,1,1,n] = 0.0
    # Gv[2,4,1,1,n] = 0.0
    # Gv[2,5,1,1,n] = 0.0

    Gv[3,1,1,1,n] = -v[2,n] * coef
    # Gv[3,2,1,1,n] = 0.0
    Gv[3,3,1,1,n] = coef
    # Gv[3,4,1,1,n] = 0.0
    # Gv[3,5,1,1,n] = 0.0

    Gv[4,1,1,1,n] = -v[3,n] * coef
    # Gv[4,2,1,1,n] = 0.0
    # Gv[4,3,1,1,n] = 0.0
    Gv[4,4,1,1,n] =  coef
    # Gv[4,5,1,1,n] = 0.0

    Gv[5,1,1,1,n] = -(one3rd*v[1,n]*v[1,n] + v2[n] + gamma_pr*(E[n] - v2[n])) * coef
    Gv[5,2,1,1,n] = (four3rd - gamma_pr)*v[1,n] * coef
    Gv[5,3,1,1,n] = (1.0 - gamma_pr)*v[2,n] * coef
    Gv[5,4,1,1,n] = (1.0 - gamma_pr)*v[3,n] * coef
    Gv[5,5,1,1,n] = gamma_pr * coef

    #
    # G12
    #
    # Gv[1,1,1,2,n] = 0.0
    # Gv[1,2,1,2,n] = 0.0
    # Gv[1,3,1,2,n] = 0.0
    # Gv[1,4,1,2,n] = 0.0
    # Gv[1,5,1,2,n] = 0.0

    Gv[2,1,1,2,n] = two3rd*v[2,n] * coef
    # Gv[2,2,1,2,n] = 0.0 
    Gv[2,3,1,2,n] = -two3rd * coef
    # Gv[2,4,1,2,n] = 0.0
    # Gv[2,5,1,2,n] = 0.0

    Gv[3,1,1,2,n] = -v[1,n] * coef
    Gv[3,2,1,2,n] = coef
    # Gv[3,3,1,2,n] = 0.0
    # Gv[3,4,1,2,n] = 0.0
    # Gv[3,5,1,2,n] = 0.0

    # Gv[4,1,1,2,n] = 0.0
    # Gv[4,2,1,2,n] = 0.0
    # Gv[4,3,1,2,n] = 0.0
    # Gv[4,4,1,2,n] = 0.0
    # Gv[4,5,1,2,n] = 0.0

    Gv[5,1,1,2,n] = -one3rd*v1v2[n] * coef
    Gv[5,2,1,2,n] = v[2,n] * coef
    Gv[5,3,1,2,n] = -two3rd*v[1,n] * coef 
    # Gv[5,4,1,2,n] = 0.0
    # Gv[5,5,1,2,n] = 0.0

    #
    # G13
    #
    # Gv[1,1,1,3,n] = 0.0
    # Gv[1,2,1,3,n] = 0.0
    # Gv[1,3,1,3,n] = 0.0
    # Gv[1,4,1,3,n] = 0.0
    # Gv[1,5,1,3,n] = 0.0

    Gv[2,1,1,3,n] = two3rd*v[3,n] * coef
    # Gv[2,2,1,3,n] = 0.0 
    # Gv[2,3,1,3,n] = 0.0
    Gv[2,4,1,3,n] = -two3rd * coef
    # Gv[2,5,1,3,n] = 0.0

    # Gv[3,1,1,3,n] = 0.0
    # Gv[3,2,1,3,n] = 0.0
    # Gv[3,3,1,3,n] = 0.0
    # Gv[3,4,1,3,n] = 0.0
    # Gv[3,5,1,3,n] = 0.0

    Gv[4,1,1,3,n] = -v[1,n] * coef
    Gv[4,2,1,3,n] = coef
    # Gv[4,3,1,3,n] = 0.0
    # Gv[4,4,1,3,n] = 0.0
    # Gv[4,5,1,3,n] = 0.0


    Gv[5,1,1,3,n] = -one3rd*v1v3[n] * coef
    Gv[5,2,1,3,n] = v[3,n] * coef
    # Gv[5,3,1,3,n] = 0.0
    Gv[5,4,1,3,n] = -two3rd*v[1,n]  * coef
    # Gv[5,5,1,3,n] = 0.0

    #
    # G21
    #
    # Gv[1,1,2,1,n] = 0.0
    # Gv[1,2,2,1,n] = 0.0
    # Gv[1,3,2,1,n] = 0.0
    # Gv[1,4,2,1,n] = 0.0
    # Gv[1,5,2,1,n] = 0.0

    Gv[2,1,2,1,n] = -v[2,n] * coef
    # Gv[2,2,2,1,n] = 0.0 
    Gv[2,3,2,1,n] = coef
    # Gv[2,4,2,1,n] = 0.0
    # Gv[2,5,2,1,n] = 0.0

    Gv[3,1,2,1,n] = two3rd*v[1,n] * coef
    Gv[3,2,2,1,n] = -two3rd * coef
    # Gv[3,3,2,1,n] = 0.0
    # Gv[3,4,2,1,n] = 0.0
    # Gv[3,5,2,1,n] = 0.0

    # Gv[4,1,2,1,n] = 0.0
    # Gv[4,2,2,1,n] = 0.0
    # Gv[4,3,2,1,n] = 0.0
    # Gv[4,4,2,1,n] = 0.0
    # Gv[4,5,2,1,n] = 0.0

    Gv[5,1,2,1,n] = -one3rd*v1v2[n] * coef
    Gv[5,2,2,1,n] = -two3rd*v[2,n] * coef
    Gv[5,3,2,1,n] = v[1,n]  * coef
    # Gv[4,4,2,1,n] = 0.0
    # Gv[5,5,2,1,n] = 0.0
    
    #
    # G22
    #
    # Gv[1,1,2,2,n] = 0.0
    # Gv[1,2,2,2,n] = 0.0
    # Gv[1,3,2,2,n] = 0.0
    # Gv[1,4,2,2,n] = 0.0
    # Gv[1,5,2,2,n] = 0.0

    Gv[2,1,2,2,n] = -v[1,n] * coef
    Gv[2,2,2,2,n] = coef
    # Gv[2,3,2,2,n] = 0.0
    # Gv[2,4,2,2,n] = 0.0
    # Gv[2,5,2,2,n] = 0.0

    Gv[3,1,2,2,n] = -four3rd*v[2,n] * coef
    # Gv[3,2,2,2,n] = 0.0
    Gv[3,3,2,2,n] = four3rd * coef
    # Gv[3,4,2,2,n] = 0.0
    # Gv[3,5,2,2,n] = 0.0

    Gv[4,1,2,2,n] = -v[1,n] * coef
    # Gv[4,2,2,2,n] = 0.0
    # Gv[4,3,2,2,n] = 0.0
    Gv[4,4,2,2,n] = coef
    # Gv[4,5,2,2,n] = 0.0

    Gv[5,1,2,2,n] = -(v2[n] + one3rd*v[2,n]*v[2,n] + gamma_pr*(E[n] - v2[n])) * coef
    Gv[5,2,2,2,n] = (1.0 - gamma_pr)*v[1,n] * coef
    Gv[5,3,2,2,n] = (four3rd - gamma_pr)*v[2,n] * coef
    Gv[5,4,2,2,n] = (1.0 - gamma_pr)*v[3,n] * coef
    Gv[5,5,2,2,n] = gamma_pr * coef

    #
    # G23
    #
    # Gv[1,1,2,3] = 0.0
    # Gv[1,2,2,3] = 0.0
    # Gv[1,3,2,3] = 0.0
    # Gv[1,4,2,3] = 0.0
    # Gv[1,5,2,3] = 0.0

    # Gv[2,1,2,3] = 0.0
    # Gv[2,2,2,3] = 0.0
    # Gv[2,3,2,3] = 0.0
    # Gv[2,4,2,3] = 0.0
    # Gv[2,5,2,3] = 0.0

    Gv[3,1,2,3] = two3rd * v[3,n] * coef
    # Gv[3,2,2,3] = 0.0
    # Gv[3,3,2,3] = 0.0
    Gv[3,4,2,3] = -two3rd * coef
    # Gv[3,5,2,3] = 0.0

    Gv[4,1,2,3] = v[2,n] * coef
    # Gv[4,2,2,3] = 0.0
    Gv[4,4,2,3] = coef
    # Gv[4,3,2,3] = 0.0
    # Gv[4,5,2,3] = 0.0

    Gv[5,1,2,3,n] = -one3rd*v2v3[n] * coef
    # Gv[5,2,2,3,n] = 0.0
    Gv[5,3,2,3,n] = v[3,n] * coef
    Gv[5,4,2,3,n] = -two3rd * v[1,n]  * coef
    # Gv[5,5,2,3,n] = 0.0

    #
    # G31
    #
    # Gv[1,1,3,1,n] = 0.0
    # Gv[1,2,3,1,n] = 0.0
    # Gv[1,3,3,1,n] = 0.0
    # Gv[1,4,3,1,n] = 0.0
    # Gv[1,5,3,1,n] = 0.0

    Gv[2,1,3,1,n] = -v[3,n] * coef
    # Gv[2,2,3,1,n] = 0.0 
    # Gv[2,3,3,1,n] = 0.0
    Gv[2,4,3,1,n] = coef
    # Gv[2,5,3,1,n] = 0.0

    # Gv[3,1,3,1,n] = 0.0
    # Gv[3,2,3,1,n] = 0.0
    # Gv[3,3,3,1,n] = 0.0
    # Gv[3,4,3,1,n] = 0.0
    # Gv[3,5,3,1,n] = 0.0

    Gv[4,1,3,1,n] = two3rd*v[1,n] * coef
    Gv[4,2,3,1,n] = -two3rd * coef
    # Gv[4,3,3,1,n] = 0.0
    # Gv[4,4,3,1,n] = 0.0
    # Gv[4,5,3,1,n] = 0.0


    Gv[5,1,3,1,n] = -one3rd*v1v3[n] * coef
    # Gv[5,2,3,1,n] = 0.0
    Gv[5,3,3,1,n] = -two3rd*v[3,n] * coef
    Gv[5,4,3,1,n] = v[2,n]  * coef
    # Gv[5,5,3,1,n] = 0.0
    #
    # G32
    #
    # Gv[1,1,3,2,n] = 0.0
    # Gv[1,2,3,2,n] = 0.0
    # Gv[1,3,3,2,n] = 0.0
    # Gv[1,4,3,2,n] = 0.0
    # Gv[1,5,3,2,n] = 0.0

    # Gv[2,1,3,2,n] = 0.0
    # Gv[2,2,3,2,n] = 0.0
    # Gv[2,3,3,2,n] = 0.0
    # Gv[2,4,3,2,n] = 0.0
    # Gv[2,5,3,2,n] = 0.0

    Gv[3,1,3,2,n] = -v[3,n] * coef
    # Gv[3,2,3,2,n] = 0.0
    # Gv[3,3,3,2,n] = 0.0
    Gv[3,4,3,2,n] = coef
    # Gv[3,5,3,2,n] = 0.0

    Gv[4,1,3,2,n] = two3rd*v[2,n]* coef
    # Gv[4,2,3,2,n] = 0.0
    Gv[4,3,3,2,n] = -two3rd* coef
    Gv[4,4,3,2,n] = 1.0* coef
    # Gv[4,5,3,2,n] = 0.0

    Gv[5,1,3,2,n] = -one3rd*v2v3[n]* coef
    # Gv[5,2,3,2,n] = 0.0
    Gv[5,3,3,2,n] = -two3rd * v[3,n]* coef
    Gv[5,4,3,2,n] = v[2,n]* coef 
    # Gv[5,5,3,2,n] = 0.0

    #
    # G33
    #
    # Gv[1,1,3,3,n] = 0.0
    # Gv[1,2,3,3,n] = 0.0
    # Gv[1,3,3,3,n] = 0.0
    # Gv[1,4,3,3,n] = 0.0
    # Gv[1,5,3,3,n] = 0.0

    Gv[2,1,3,3,n] = -v[1,n]* coef
    Gv[2,2,3,3,n] = coef
    # Gv[2,3,3,3,n] = 0.0
    # Gv[2,4,3,3,n] = 0.0
    # Gv[2,5,3,3,n] = 0.0

    Gv[3,1,3,3,n] = -v[2,n] * coef
    # Gv[3,2,3,3,n] = 0.0
    Gv[3,3,2,2,n] = coef
    # Gv[3,4,3,3,n] = 0.0
    # Gv[3,5,3,3,n] = 0.0

    Gv[4,1,2,3,n] = -four3rd * v[3,n]* coef
    # Gv[4,2,3,3,n] = 0.0
    # Gv[4,3,3,3,n] = 0.0
    Gv[4,4,3,3,n] = four3rd* coef
    # Gv[4,5,3,3,n] = 0.0

    Gv[5,1,3,3,n] = -(v2[n] + one3rd*v[3,n]*v[3,n] + gamma_pr*(E[n] - v2[n]))* coef
    Gv[5,2,3,3,n] = (1.0 - gamma_pr)*v[1,n]* coef
    Gv[5,3,3,3,n] = (1.0 - gamma_pr)*v[2,n]* coef
    Gv[5,4,3,3,n] = (four3rd - gamma_pr)*v[3,n] * coef
    Gv[5,5,3,3,n] = gamma_pr * coef

    # coef = rmuk[1,n]/q[1,n]

    # # TODO: hardcode this so that we don't loop over those zero values
    # for jDim = 1 : 3
      # for iDim = 1 : 3
        # for jDof = 1: 5
          # for iDof = 1 : 5
            # Gv[iDof, jDof, iDim, jDim ,n] *= coef 
          # end
        # end
      # end
    # end
  end
  return nothing
end
@doc """

Compute the viscous diffusion matrix for adiabatic wall. This matrix
is different from calcDiffusionTensor because it is changed to
satisfy ∇T = 0 condition. This is realized by setting gamma_pr = 0.0

Input: 
  q = conservative variable at a node
Output
  jac = viscous flux jacobian at each node, dimension = (dim+2, dim+2, dim, dim, numNodes)
"""->
function calcDiffusionTensorOnAdiabaticWall(params::ParamType{2, :conservative},
                                           q::AbstractArray{Tsol, 2},
                                           nrm::AbstractArray{Tmsh, 2},
                                           Gv::AbstractArray{Tsol, 5}) where {Tmsh, Tsol}

  Tdim = 2
  @assert(size(q, 2) == size(Gv, 5))
  @assert(size(Gv, 4) == Tdim)
  @assert(size(Gv, 3) == Tdim)
  @assert(size(Gv, 2) == Tdim+2)
  @assert(size(Gv, 1) == Tdim+2)
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
    #
    # cx and cy are use to impose n⋅∇T = 0.0
    # 
    # 
    cx = gamma_pr
    cy = 0.0
    if abs(cy) < 1.e-13
      cx = 0
      cy = 0
    else
      cy = -nrm[1,n] / nrm[2,n] * gamma_pr
    end
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

    Gv[4,1,1,1,n] = -(four3rd*v[1,n]*v[1,n] + v[2,n]*v[2,n] + cx*(E[n] - v2[n]))
    Gv[4,2,1,1,n] = (four3rd - cx)*v[1,n]
    Gv[4,3,1,1,n] = (1.0 - cx)*v[2,n]
    Gv[4,4,1,1,n] = cx

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

    Gv[4,1,2,1,n] = -one3rd*v1v2[n] - cy*(E[n] - v2[n])
    Gv[4,2,2,1,n] = -two3rd*v[2,n] - cy*v[1,n]
    Gv[4,3,2,1,n] = v[1,n] - cy*v[2,n]
    Gv[4,4,2,1,n] = cy 
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

    Gv[4,1,2,2,n] = -(v[1,n]*v[1,n] + four3rd*v[2,n]*v[2,n])
    Gv[4,2,2,2,n] = v[1,n]
    Gv[4,3,2,2,n] = four3rd*v[2,n]
    Gv[4,4,2,2,n] = 0.0

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

# NOTE: this is incorrect. states that grad T = 0 always by setting gamma_pr to 0, but what actually is
#       necessary is that n dot grad T = 0.
# TODO: extend the correct version of cDTOnAdiabaticWall (above) to 3D
function calcDiffusionTensor_adiabaticWall(params::ParamType{3, :conservative},
                                           q::AbstractArray{Tsol, 2},
                                           Gv::AbstractArray{Tsol, 5}) where Tsol

  Tdim = 3
  @assert(size(q, 2) == size(Gv, 5))
  @assert(size(Gv, 4) == Tdim)
  @assert(size(Gv, 3) == Tdim)
  @assert(size(Gv, 2) == Tdim+2)
  @assert(size(Gv, 1) == Tdim+2)
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
  v1v3 = Array(Tsol, numNodes)
  v2v3 = Array(Tsol, numNodes)
  v    = Array(Tsol, Tdim, numNodes)
  rmuk = Array(Tsol, 2, numNodes)

  for n = 1 : numNodes
    v[1, n] = q[2,n]/q[1,n]
    v[2, n] = q[3,n]/q[1,n]
    v[3, n] = q[4,n]/q[1,n]
    E[n]    = q[5,n]/q[1,n]
    v2[n]   = v[1,n]*v[1,n] + v[2,n]*v[2,n] + v[3,n]*v[3,n]
    v1v2[n] = v[1,n]*v[2,n]
    v1v3[n] = v[1,n]*v[3,n]
    v2v3[n] = v[2,n]*v[3,n]
    T[n]    = gamma*gamma_1 * (E[n] - 0.5*v2[n] )
  end

  getMuK(T, rmuk)

  for n = 1 : numNodes
    coef = rmuk[1,n]/q[1,n]

    # Gv[1,1,1,1,n] = 0.0
    # Gv[1,2,1,1,n] = 0.0
    # Gv[1,3,1,1,n] = 0.0
    # Gv[1,4,1,1,n] = 0.0
    # Gv[1,5,1,1,n] = 0.0

    Gv[2,1,1,1,n] = -four3rd*v[1,n] * coef
    Gv[2,2,1,1,n] = four3rd  * coef
    # Gv[2,3,1,1,n] = 0.0
    # Gv[2,4,1,1,n] = 0.0
    # Gv[2,5,1,1,n] = 0.0

    Gv[3,1,1,1,n] = -v[2,n] * coef
    # Gv[3,2,1,1,n] = 0.0
    Gv[3,3,1,1,n] = coef
    # Gv[3,4,1,1,n] = 0.0
    # Gv[3,5,1,1,n] = 0.0

    Gv[4,1,1,1,n] = -v[3,n] * coef
    # Gv[4,2,1,1,n] = 0.0
    # Gv[4,3,1,1,n] = 0.0
    Gv[4,4,1,1,n] =  coef
    # Gv[4,5,1,1,n] = 0.0

    Gv[5,1,1,1,n] = -(one3rd*v[1,n]*v[1,n] + v2[n] + gamma_pr*(E[n] - v2[n])) * coef
    Gv[5,2,1,1,n] = (four3rd - gamma_pr)*v[1,n] * coef
    Gv[5,3,1,1,n] = (1.0 - gamma_pr)*v[2,n] * coef
    Gv[5,4,1,1,n] = (1.0 - gamma_pr)*v[3,n] * coef
    Gv[5,5,1,1,n] = gamma_pr * coef

    #
    # G12
    #
    # Gv[1,1,1,2,n] = 0.0
    # Gv[1,2,1,2,n] = 0.0
    # Gv[1,3,1,2,n] = 0.0
    # Gv[1,4,1,2,n] = 0.0
    # Gv[1,5,1,2,n] = 0.0

    Gv[2,1,1,2,n] = two3rd*v[2,n] * coef
    # Gv[2,2,1,2,n] = 0.0 
    Gv[2,3,1,2,n] = -two3rd * coef
    # Gv[2,4,1,2,n] = 0.0
    # Gv[2,5,1,2,n] = 0.0

    Gv[3,1,1,2,n] = -v[1,n] * coef
    Gv[3,2,1,2,n] = coef
    # Gv[3,3,1,2,n] = 0.0
    # Gv[3,4,1,2,n] = 0.0
    # Gv[3,5,1,2,n] = 0.0

    # Gv[4,1,1,2,n] = 0.0
    # Gv[4,2,1,2,n] = 0.0
    # Gv[4,3,1,2,n] = 0.0
    # Gv[4,4,1,2,n] = 0.0
    # Gv[4,5,1,2,n] = 0.0

    Gv[5,1,1,2,n] = -one3rd*v1v2[n] * coef
    Gv[5,2,1,2,n] = v[2,n] * coef
    Gv[5,3,1,2,n] = -two3rd*v[1,n] * coef 
    # Gv[5,4,1,2,n] = 0.0
    # Gv[5,5,1,2,n] = 0.0

    #
    # G13
    #
    # Gv[1,1,1,3,n] = 0.0
    # Gv[1,2,1,3,n] = 0.0
    # Gv[1,3,1,3,n] = 0.0
    # Gv[1,4,1,3,n] = 0.0
    # Gv[1,5,1,3,n] = 0.0

    Gv[2,1,1,3,n] = two3rd*v[3,n] * coef
    # Gv[2,2,1,3,n] = 0.0 
    # Gv[2,3,1,3,n] = 0.0
    Gv[2,4,1,3,n] = -two3rd * coef
    # Gv[2,5,1,3,n] = 0.0

    # Gv[3,1,1,3,n] = 0.0
    # Gv[3,2,1,3,n] = 0.0
    # Gv[3,3,1,3,n] = 0.0
    # Gv[3,4,1,3,n] = 0.0
    # Gv[3,5,1,3,n] = 0.0

    Gv[4,1,1,3,n] = -v[1,n] * coef
    Gv[4,2,1,3,n] = coef
    # Gv[4,3,1,3,n] = 0.0
    # Gv[4,4,1,3,n] = 0.0
    # Gv[4,5,1,3,n] = 0.0


    Gv[5,1,1,3,n] = -one3rd*v1v3[n] * coef
    Gv[5,2,1,3,n] = v[3,n] * coef
    # Gv[5,3,1,3,n] = 0.0
    Gv[5,4,1,3,n] = -two3rd*v[1,n]  * coef
    # Gv[5,5,1,3,n] = 0.0

    #
    # G21
    #
    # Gv[1,1,2,1,n] = 0.0
    # Gv[1,2,2,1,n] = 0.0
    # Gv[1,3,2,1,n] = 0.0
    # Gv[1,4,2,1,n] = 0.0
    # Gv[1,5,2,1,n] = 0.0

    Gv[2,1,2,1,n] = -v[2,n] * coef
    # Gv[2,2,2,1,n] = 0.0 
    Gv[2,3,2,1,n] = coef
    # Gv[2,4,2,1,n] = 0.0
    # Gv[2,5,2,1,n] = 0.0

    Gv[3,1,2,1,n] = two3rd*v[1,n] * coef
    Gv[3,2,2,1,n] = -two3rd * coef
    # Gv[3,3,2,1,n] = 0.0
    # Gv[3,4,2,1,n] = 0.0
    # Gv[3,5,2,1,n] = 0.0

    # Gv[4,1,2,1,n] = 0.0
    # Gv[4,2,2,1,n] = 0.0
    # Gv[4,3,2,1,n] = 0.0
    # Gv[4,4,2,1,n] = 0.0
    # Gv[4,5,2,1,n] = 0.0

    Gv[5,1,2,1,n] = -one3rd*v1v2[n] * coef
    Gv[5,2,2,1,n] = -two3rd*v[2,n] * coef
    Gv[5,3,2,1,n] = v[1,n]  * coef
    # Gv[4,4,2,1,n] = 0.0
    # Gv[5,5,2,1,n] = 0.0
    
    #
    # G22
    #
    # Gv[1,1,2,2,n] = 0.0
    # Gv[1,2,2,2,n] = 0.0
    # Gv[1,3,2,2,n] = 0.0
    # Gv[1,4,2,2,n] = 0.0
    # Gv[1,5,2,2,n] = 0.0

    Gv[2,1,2,2,n] = -v[1,n] * coef
    Gv[2,2,2,2,n] = coef
    # Gv[2,3,2,2,n] = 0.0
    # Gv[2,4,2,2,n] = 0.0
    # Gv[2,5,2,2,n] = 0.0

    Gv[3,1,2,2,n] = -four3rd*v[2,n] * coef
    # Gv[3,2,2,2,n] = 0.0
    Gv[3,3,2,2,n] = four3rd * coef
    # Gv[3,4,2,2,n] = 0.0
    # Gv[3,5,2,2,n] = 0.0

    Gv[4,1,2,2,n] = -v[1,n] * coef
    # Gv[4,2,2,2,n] = 0.0
    # Gv[4,3,2,2,n] = 0.0
    Gv[4,4,2,2,n] = coef
    # Gv[4,5,2,2,n] = 0.0

    Gv[5,1,2,2,n] = -(v2[n] + one3rd*v[2,n]*v[2,n] + gamma_pr*(E[n] - v2[n])) * coef
    Gv[5,2,2,2,n] = (1.0 - gamma_pr)*v[1,n] * coef
    Gv[5,3,2,2,n] = (four3rd - gamma_pr)*v[2,n] * coef
    Gv[5,4,2,2,n] = (1.0 - gamma_pr)*v[3,n] * coef
    Gv[5,5,2,2,n] = gamma_pr * coef

    #
    # G23
    #
    # Gv[1,1,2,3] = 0.0
    # Gv[1,2,2,3] = 0.0
    # Gv[1,3,2,3] = 0.0
    # Gv[1,4,2,3] = 0.0
    # Gv[1,5,2,3] = 0.0

    # Gv[2,1,2,3] = 0.0
    # Gv[2,2,2,3] = 0.0
    # Gv[2,3,2,3] = 0.0
    # Gv[2,4,2,3] = 0.0
    # Gv[2,5,2,3] = 0.0

    Gv[3,1,2,3] = two3rd * v[3,n] * coef
    # Gv[3,2,2,3] = 0.0
    # Gv[3,3,2,3] = 0.0
    Gv[3,4,2,3] = -two3rd * coef
    # Gv[3,5,2,3] = 0.0

    Gv[4,1,2,3] = v[2,n] * coef
    # Gv[4,2,2,3] = 0.0
    Gv[4,4,2,3] = coef
    # Gv[4,3,2,3] = 0.0
    # Gv[4,5,2,3] = 0.0

    Gv[5,1,2,3,n] = -one3rd*v2v3[n] * coef
    # Gv[5,2,2,3,n] = 0.0
    Gv[5,3,2,3,n] = v[3,n] * coef
    Gv[5,4,2,3,n] = -two3rd * v[1,n]  * coef
    # Gv[5,5,2,3,n] = 0.0

    #
    # G31
    #
    # Gv[1,1,3,1,n] = 0.0
    # Gv[1,2,3,1,n] = 0.0
    # Gv[1,3,3,1,n] = 0.0
    # Gv[1,4,3,1,n] = 0.0
    # Gv[1,5,3,1,n] = 0.0

    Gv[2,1,3,1,n] = -v[3,n] * coef
    # Gv[2,2,3,1,n] = 0.0 
    # Gv[2,3,3,1,n] = 0.0
    Gv[2,4,3,1,n] = coef
    # Gv[2,5,3,1,n] = 0.0

    # Gv[3,1,3,1,n] = 0.0
    # Gv[3,2,3,1,n] = 0.0
    # Gv[3,3,3,1,n] = 0.0
    # Gv[3,4,3,1,n] = 0.0
    # Gv[3,5,3,1,n] = 0.0

    Gv[4,1,3,1,n] = two3rd*v[1,n] * coef
    Gv[4,2,3,1,n] = -two3rd * coef
    # Gv[4,3,3,1,n] = 0.0
    # Gv[4,4,3,1,n] = 0.0
    # Gv[4,5,3,1,n] = 0.0


    Gv[5,1,3,1,n] = -one3rd*v1v3[n] * coef
    # Gv[5,2,3,1,n] = 0.0
    Gv[5,3,3,1,n] = -two3rd*v[3,n] * coef
    Gv[5,4,3,1,n] = v[2,n]  * coef
    # Gv[5,5,3,1,n] = 0.0
    #
    # G32
    #
    # Gv[1,1,3,2,n] = 0.0
    # Gv[1,2,3,2,n] = 0.0
    # Gv[1,3,3,2,n] = 0.0
    # Gv[1,4,3,2,n] = 0.0
    # Gv[1,5,3,2,n] = 0.0

    # Gv[2,1,3,2,n] = 0.0
    # Gv[2,2,3,2,n] = 0.0
    # Gv[2,3,3,2,n] = 0.0
    # Gv[2,4,3,2,n] = 0.0
    # Gv[2,5,3,2,n] = 0.0

    Gv[3,1,3,2,n] = -v[3,n] * coef
    # Gv[3,2,3,2,n] = 0.0
    # Gv[3,3,3,2,n] = 0.0
    Gv[3,4,3,2,n] = coef
    # Gv[3,5,3,2,n] = 0.0

    Gv[4,1,3,2,n] = two3rd*v[2,n]* coef
    # Gv[4,2,3,2,n] = 0.0
    Gv[4,3,3,2,n] = -two3rd* coef
    Gv[4,4,3,2,n] = 1.0* coef
    # Gv[4,5,3,2,n] = 0.0

    Gv[5,1,3,2,n] = -one3rd*v2v3[n]* coef
    # Gv[5,2,3,2,n] = 0.0
    Gv[5,3,3,2,n] = -two3rd * v[3,n]* coef
    Gv[5,4,3,2,n] = v[2,n]* coef 
    # Gv[5,5,3,2,n] = 0.0

    #
    # G33
    #
    # Gv[1,1,3,3,n] = 0.0
    # Gv[1,2,3,3,n] = 0.0
    # Gv[1,3,3,3,n] = 0.0
    # Gv[1,4,3,3,n] = 0.0
    # Gv[1,5,3,3,n] = 0.0

    Gv[2,1,3,3,n] = -v[1,n]* coef
    Gv[2,2,3,3,n] = coef
    # Gv[2,3,3,3,n] = 0.0
    # Gv[2,4,3,3,n] = 0.0
    # Gv[2,5,3,3,n] = 0.0

    Gv[3,1,3,3,n] = -v[2,n] * coef
    # Gv[3,2,3,3,n] = 0.0
    Gv[3,3,2,2,n] = coef
    # Gv[3,4,3,3,n] = 0.0
    # Gv[3,5,3,3,n] = 0.0

    Gv[4,1,2,3,n] = -four3rd * v[3,n]* coef
    # Gv[4,2,3,3,n] = 0.0
    # Gv[4,3,3,3,n] = 0.0
    Gv[4,4,3,3,n] = four3rd* coef
    # Gv[4,5,3,3,n] = 0.0

    Gv[5,1,3,3,n] = -(v2[n] + one3rd*v[3,n]*v[3,n] + gamma_pr*(E[n] - v2[n]))* coef
    Gv[5,2,3,3,n] = (1.0 - gamma_pr)*v[1,n]* coef
    Gv[5,3,3,3,n] = (1.0 - gamma_pr)*v[2,n]* coef
    Gv[5,4,3,3,n] = (four3rd - gamma_pr)*v[3,n] * coef
    Gv[5,5,3,3,n] = gamma_pr * coef

    # coef = rmuk[1,n]/q[1,n]

    # # TODO: hardcode this so that we don't loop over those zero values
    # for jDim = 1 : 3
      # for iDim = 1 : 3
        # for jDof = 1: 5
          # for iDof = 1 : 5
            # Gv[iDof, jDof, iDim, jDim ,n] *= coef 
          # end
        # end
      # end
    # end
  end

  return nothing
end

abstract type AbstractBoundaryValueType end

@doc """

compute the adiabatic wall boundary value

Input:
q_in    :: conservative variables on boundary face
norm    :: outward unit normal of the boundary face
Output:
q_bnd    :: the boundary value

"""->
mutable struct AdiabaticWall <: AbstractBoundaryValueType
end
function (obj::AdiabaticWall)(
              q_in::AbstractArray{Tsol, 2},
              xy::AbstractArray{Tmsh, 2},
              norm::AbstractArray{Tmsh, 2},
              params::ParamType{Tdim, :conservative},
              q_bnd::AbstractArray{Tsol, 2}) where {Tsol, Tmsh, Tdim}
  @assert( size(q_in, 1) == size(q_bnd,  1))
  @assert( size(q_in, 2) == size(q_bnd,  2))
  @assert( size(q_in, 2) == size(norm,  2))

  dim = size(norm, 1)
  numNodes = size(q_in, 2)

  for n = 1 : numNodes
    q_bnd[1, n] = q_in[1, n]
    # noslip condition
    q_bnd[2:dim+1, n] = 0.0
    q_bnd[dim+2, n] = q_in[dim+2, n]
  end

  return nothing
end

mutable struct ExactChannel<: AbstractBoundaryValueType
end
function (obj::ExactChannel)(
                          q_in::AbstractArray{Tsol, 2},
                          xyz::AbstractArray{Tmsh, 2},
                          norm::AbstractArray{Tmsh, 2},
                          params::ParamType{3, :conservative},
                          qg::AbstractArray{Tsol, 2}) where {Tsol, Tmsh}

  gamma = params.gamma
  gamma_1 = params.gamma - 1
  sigma = 0.01
  aoa = params.aoa
  beta = params.sideslip_angle
  rhoInf = 1.0
  uInf = params.Ma * cos(beta) * cos(aoa)
  vInf = params.Ma * sin(beta) * -1
  wInf = params.Ma * cos(beta) * sin(aoa)
  TInf = 1.0
  numNodes = size(qg, 2)
  for n = 1 : numNodes
    x = xyz[1,n]
    y = xyz[2,n]
    z = xyz[3,n]

    rho = rhoInf * (1 + sigma*x*y*z)
    ux = sin(pi*x) + 1
    uy = sin(pi*y) + 1
    uz = sin(pi*z) + 1
    u  = (1 + sigma*ux * uy * uz )* uInf
    vx = sin(pi*x) + 1
    vy = sin(pi*y) + 1
    vz = sin(pi*z) + 1
    v  = (1 + sigma*vx * vy * vz )* vInf
    wx = sin(pi*x) + 1
    wy = sin(pi*y) + 1
    wz = sin(pi*z) + 1
    w  = (1 + sigma*wx * wy * wz) * wInf
    T  = TInf 

    if !params.isViscous
      u += 0.2 * uInf
      v += 0.2 * vInf
      w += 0.2 * wInf
    end

    qg[1,n] = rho
    qg[2,n] = rho*u
    qg[3,n] = rho*v
    qg[4,n] = rho*w
    qg[5,n] = T/(gamma * gamma_1) + 0.5 * (u*u + v*v + w*w)
    qg[5,n] *= rho
  end

  return nothing
end



function (obj::ExactChannel)(
                          q_in::AbstractArray{Tsol, 2},
                          xy::AbstractArray{Tmsh, 2},
                          norm::AbstractArray{Tmsh, 2},
                          params::ParamType{2, :conservative},
                          q_bnd::AbstractArray{Tsol, 2}) where {Tsol, Tmsh}
  @assert( size(q_in, 1) == size(q_bnd,  1))
  @assert( size(q_in, 2) == size(q_bnd,  2))
  @assert( size(q_in, 2) == size(norm,  2))
  dim = size(norm, 1)
  numNodes = size(q_in, 2)
  sigma = 0.1

  gamma = params.gamma
  gamma_1 = gamma - 1

  aoa = params.aoa
  qRef = zeros(Float64, dim+2)
  qRef[1] = 1.0
  qRef[2] = params.Ma*cos(aoa)
  qRef[3] = params.Ma*sin(aoa)
  qRef[4] = 1.0

  for n = 1 : numNodes
    x = xy[1, n]
    y = xy[2, n]
    rho = qRef[1] * (sigma*exp(sin(0.5*pi*(x+y))) +  1.0)
    ux  = (exp(x) * sin(pi*x) * sigma + 1) * qRef[2]
    uy  = exp(y) * sin(pi*y)
    u   = ux * uy
    vx  = (exp(x) * sin(pi*x) * sigma + 1) * qRef[3]
    vy  = exp(y) * sin(pi*y)
    v   = vx * vy
    T   = (1 + sigma*exp(0.1*x+0.1*y)) * qRef[4]
    # T   = qRef[4]
    if !params.isViscous
      u += 0.2 * qRef[2]
    end
    q_bnd[1, n] = rho 
    q_bnd[2, n] = rho * u
    q_bnd[3, n] = rho * v
    q_bnd[4, n] = T/(gamma*gamma_1) + 0.5 * (u*u + v*v)
    q_bnd[4, n] *= rho 
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
mutable struct Farfield <: AbstractBoundaryValueType
end

function (obj::Farfield)(
              q_in::AbstractArray{Tsol, 2},
              xy::AbstractArray{Tmsh, 2},
              norm::AbstractArray{Tmsh, 2},
              params::ParamType{2, :conservative},
              q_bnd::AbstractArray{Tsol, 2}) where {Tsol, Tmsh}
  @assert( size(q_in, 1) == size(q_bnd,  1))
  @assert( size(q_in, 2) == size(q_bnd,  2))
  @assert( size(q_in, 2) == size(norm,  2))

  dim      = size(norm, 1)
  numNodes = size(q_in, 2)
  lambda   = zeros(Float64, 3) 
  gamma    = params.gamma
  gamma_1  = params.gamma - 1
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
    u = q_in[2,n] / q_in[1,n]
    v = q_in[3,n] / q_in[1,n]
    vn = u*norm[1, n] + v*norm[2, n]
    v2 = u*u + v*v 
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
      q_bnd[1, n] = qInf[1]
      q_bnd[2, n] = qInf[2]
      q_bnd[3, n] = qInf[3]
      q_bnd[4, n] = p/gamma_1 + 0.5*MaInf*MaInf*qInf[1]
    else                        # subsonic outflow
      pInf = 1.0/gamma
      q_bnd[1, n] = q_in[1, n]
      q_bnd[2, n] = q_in[2, n]
      q_bnd[3, n] = q_in[3, n]
      q_bnd[4, n] = pInf/gamma_1 + 0.5*q_in[1,n] * v2
    end    
    
    # DEBUG ONLY
    # q_bnd[:,n] = qInf[:]
    # DEBUG END
  end

  return nothing
end

#
# TODO: combine this with the 2D version. The only difference
# between 2D and 3D is the freestream conditions.
#
function (obj::Farfield)(
              q_in::AbstractArray{Tsol, 2},
              xy::AbstractArray{Tmsh, 2},
              norm::AbstractArray{Tmsh, 2},
              params::ParamType{3, :conservative},
              q_bnd::AbstractArray{Tsol, 2}) where {Tsol, Tmsh}
  @assert( size(q_in, 1) == size(q_bnd,  1))
  @assert( size(q_in, 2) == size(q_bnd,  2))
  @assert( size(q_in, 2) == size(norm,  2))
  @assert( size(norm, 1) == 3 )

  dim      = size(norm, 1)
  numNodes = size(q_in, 2)
  lambda   = zeros(Float64, 3) 
  gamma    = params.gamma
  gamma_1  = params.gamma_1
  aoa      = params.aoa
  beta     = params.sideslip_angle
  gg_1     = (gamma*gamma_1)
  MaInf    = params.Ma    

  uInf = params.Ma * cos(beta) * cos(aoa)
  vInf = params.Ma * sin(beta) * -1
  wInf = params.Ma * cos(beta) * sin(aoa)
  # freestream values
  qInf = zeros(Float64, dim + 2)
  qInf[1] = 1.0
  qInf[2] = qInf[1]*uInf
  qInf[3] = qInf[1]*vInf
  qInf[4] = qInf[1]*wInf
  qInf[5] = qInf[1]*(1.0/gg_1 + 0.5*MaInf*MaInf)

  for n = 1 : numNodes
    #
    # contravariant velocity
    #
    vn = q_in[2, n]*norm[1, n] + q_in[3, n]*norm[2, n] + q_in[4,n] *norm[3,n]
    vn = vn/q_in[1,n]
    v2 = q_in[2,n]*q_in[2,n] + q_in[3,n]*q_in[3,n] + q_in[4,n]*q_in[4,n]
    v2 /= q_in[1,n]*q_in[1,n]
    T = gg_1*(q_in[5,n]/q_in[1,n] - 0.5*v2)
    a = sqrt(T)
    #
    # eigenvalues
    # 
    lambda[1] = real(vn)
    lambda[2] = real(vn + a)
    lambda[3] = real(vn - a)

    # depending on eigenvalue, set state boundary value to exterior or interior (depending on field)
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
    
    # debug only
    q_bnd[:,n] = qInf[:]        # for MMS in which q at boundary equals to qInf
    # end
  end

  return nothing
end

