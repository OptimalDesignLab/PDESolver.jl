# boundary condition functions

@doc """
Now actually we are integrating 
  ∫ G∇ϕ:[q] dΓ

Input: 
  mesh
  sbp
  eqn
  opts
Output:

"""->
function evalBoundaryIntegrals_vector(mesh::AbstractMesh{Tmsh},
                                      sbp::AbstractOperator,
                                      eqn::NSData{Tsol, Tres, Tdim},
                                      opts) where {Tmsh, Tsol, Tres, Tdim}

  sbpface = mesh.sbpface
  Dx = Array{Tmsh}((mesh.numNodesPerElement, mesh.numNodesPerElement, Tdim))
  R = sview(sbpface.interp[:,:])
  w = sview(sbpface.wface, :)
  res = sview(eqn.res, :,:,:)

  numNodes_elem = mesh.numNodesPerElement
  numNodes_face = mesh.numNodesPerFace
  stencilsize   = sbpface.stencilsize
  q_bnd = Array{Tsol}(mesh.numDofPerNode, mesh.numNodesPerFace)
  dq    = Array{Tsol}(mesh.numDofPerNode, mesh.numNodesPerFace)
  Gt = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode, Tdim, Tdim, mesh.numNodesPerFace)
  RDx = zeros(Tmsh, mesh.numNodesPerFace, mesh.numNodesPerElement, Tdim)
  GtRDx = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode, mesh.numNodesPerFace, mesh.numNodesPerElement)
  nrm     = Array{Tmsh}(Tdim, mesh.numNodesPerFace)
  nrm1 = Array{Tmsh}(Tdim, mesh.numNodesPerFace)
  area = Array{Tmsh}(mesh.numNodesPerFace)

  # loop over all the boundaries
  for bc = 1:mesh.numBC
    indx0 = mesh.bndry_offsets[bc]
    indx1 = mesh.bndry_offsets[bc+1] - 1

    for f = indx0:indx1
      bndry = mesh.bndryfaces[f]
      elem = bndry.element
      face = bndry.face
      p = sview(sbpface.perm, :, face)


      # compute RDx
      calcDx(mesh, sbp, elem, Dx)

      for i = 1 : length(RDx)
        RDx[i] = 0.0
      end

      for d =  1 : Tdim    
        for row = 1 : numNodes_face
          for col = 1 : numNodes_elem
            for s = 1 : stencilsize
              RDx[row, col, d] +=  R[s, row] * Dx[p[s], col, d]     
            end
          end
        end
      end

      vecflux = sview(eqn.vecflux_bndry, :,:,:,f)
      for i = 1 : numNodes_elem
        for j = 1 : numNodes_face
          for iDof = 2 : mesh.numDofPerNode
            # res[iDof, i, elem] +=  ( RDx[j, i, 1] * vecflux[1, iDof, j] 
            # + RDx[j, i, 2] * vecflux[2, iDof, j] ) * w[j]
            tmp = 0.0
            for iDim = 1 : Tdim
              tmp += RDx[j, i, iDim] * vecflux[iDim, iDof, j]
            end
            res[iDof, i, elem] += tmp * w[j]
          end
        end
      end
    end
  end

  return nothing
end  # end evalBoundaryIntegrals_vector


function calcViscousFlux_boundary(mesh::AbstractMesh{Tmsh},
                                  sbp::AbstractOperator,
                                  eqn::NSData{Tsol, Tres, Tdim},
                                  opts) where {Tmsh, Tsol, Tres, Tdim}
  # freestream info
  Ma = eqn.params.euler_params.Ma
  Re = eqn.params.Re
  gamma_1 = eqn.params.euler_params.gamma_1
  Pr = 0.72
  coef_nondim = Ma/Re

  p = opts["order"]
  sat_type = opts["SAT_type"]
  const_tii = (p + 1.0)*(p + Tdim)/Tdim
  sbpface = mesh.sbpface
  dq    = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)    
  dqn   = zeros(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerFace)    
  q_bnd = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)    
  pMat  = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode, mesh.numNodesPerFace)
  Gt    = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode, Tdim, Tdim, mesh.numNodesPerFace)
  Gt_bnd  = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode, Tdim, Tdim, mesh.numNodesPerFace)
  Fv_face = zeros(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerFace)
  Fv_bnd  = zeros(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerFace)
  vecflux = zeros(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerFace)

  nrm1 = Array{Tmsh}(Tdim, mesh.numNodesPerFace)
  area = Array{Tmsh}(mesh.numNodesPerFace)
  area_sum = sview(eqn.area_sum, :)

  # sigma = calcTraceInverseInequalityConst(sbp, sbpface)
  dqdx_elem = zeros(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerElement )
  dqdx_face = zeros(Tsol, Tdim, mesh.numDofPerNode, mesh.numNodesPerFace )
  for iBC = 1 : mesh.numBC
    indx0 = mesh.bndry_offsets[iBC]
    indx1 = mesh.bndry_offsets[iBC+1] - 1

    # specify boundary value function
    # TODO: Put it into a function 
    local bnd_functor::AbstractBoundaryValueType
    key_i = string("BC", iBC, "_name")
    val = opts[key_i]
    Gt_functor = calcDiffusionTensor
    if val == "FreeStreamBC"
      bnd_functor = Farfield()
    elseif val == "ExactChannelBC"
      bnd_functor = ExactChannel()
    elseif val == "nonslipBC"
      bnd_functor = AdiabaticWall()
      Gt_functor = calcDiffusionTensorOnAdiabaticWall
    elseif val == "noPenetrationBC"
      continue
    elseif val == "zeroPressGradientBC"
      bnd_functor = Farfield()
    else
      error("iBC = ", iBC, ", Only 'FreeStreamBC' and 'nonslipBC' available")
    end

    for f = indx0 : indx1

      flux  = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
      bndry = mesh.bndryfaces[f]
      elem = bndry.element
      face = bndry.face
      perm = sview(sbpface.perm, :, face)
      xy = sview(mesh.coords_bndry, :, :, f)

      # Compute geometric info on face
      nrm_xy = ro_sview(mesh.nrm_bndry, :, :, f)
      for n = 1 : mesh.numNodesPerFace
        area[n] = norm(ro_sview(nrm_xy, :, n)) 

        for i = 1 : Tdim
          nrm1[i,n] = nrm_xy[i,n] / area[n]
        end
      end

      # We need viscous flux and diffusion tensor on interfaces, and there
      # are different ways to compute them. For viscous flux:
      # 1) since we have Q on face nodes, if 𝛻Q is available on interface, then we are done.
      # 2) we can comoute viscous flux on volume nodes and then interpolate to interface node.
      # It's logically simple but computationally expensive.

      # Compute boundary viscous flux, F(q_b, ∇q) = G(q_b)∇q.
      # so we need viscousity tensor G, and derivative of q.
      q_face = sview(eqn.q_bndry, :, :, f)
      bnd_functor(q_face, xy, nrm1, eqn.params, q_bnd)

      # diffusion matrix used in penalty term should be computed from q_face rather than q_bnd
      if val == "nonslipBC"
        Gt_functor(eqn.params, q_bnd, nrm1, Gt)
      else
        Gt_functor(eqn.params, q_bnd, Gt)
      end
      q_elem = sview(eqn.q, :, :, elem)
      calcGradient(mesh, sbp, elem, q_elem, dqdx_elem)

      #
      # TODO: we can consider the first 2 dimension as a single dimension,
      # then we will not need slice here any more.
      #
      for d = 1 : Tdim
        q_x_node = sview(dqdx_elem, d, :, :)
        q_x_face = sview(dqdx_face, d, :, :)
        boundaryinterpolate(sbpface, bndry, q_x_node, q_x_face) 
      end

      calcFvis(eqn.params, Gt, dqdx_face, Fv_face)

      # compute penalty matrix
      cmptBPMat(mesh, sbp, eqn, opts, f, Gt, pMat)

      # Start to compute fluxes.  We have 3 terms on interfaces:
      # 1) -{Fv}⋅[ϕ]
      # 2) -{G^T ∇ϕ}⋅[q] 
      # 3) +δ{G}[q]:[ϕ]
      for n = 1 : mesh.numNodesPerFace
        for iDof = 1 : mesh.numDofPerNode
          dq[iDof, n] = q_face[iDof, n] - q_bnd[iDof, n]
          # dqn[1, iDof, n] = -dq[iDof, n]*nrm_xy[1, n]
          # dqn[2, iDof, n] = -dq[iDof, n]*nrm_xy[2, n]
        end
      end
      # DEBUG BEGIN
      # if maximum(abs.(real(dq))) > 1.e-11
        # println(real(dq))
      # end
      # DEBUG END


      # -----------------------------------------------
      # This part computes the contribution of
      # ∫ {G^T∇ϕ}:[q] dΓ = ∫ ∇ϕ⋅F dΓ , 
      # where 
      # [q] = (q+ - q-) ⊗ n, 
      # G = G(q_b) depends on boudanry value.
      # Then we can consider Δq⊗n as ∇q and F as viscous flux.
      # -----------------------------------------------

      # calcFvis(Gt, dqn, vecflux)
      fill!(vecflux, 0.0)
      for n = 1 : mesh.numNodesPerFace
        for iDof = 1 : mesh.numDofPerNode
          for iDim = 1 : Tdim
            for jDim = 1 : Tdim
              tmp = 0.0
              for jDof = 1 : mesh.numDofPerNode
                tmp += Gt[iDof, jDof, iDim, jDim, n]
              end
              vecflux[iDim, iDof, n] += tmp * nrm_xy[jDim, n]
            end
            vecflux[iDim,iDof,n] *=  dq[iDof,n]
          end
        end
      end

      for n = 1 : mesh.numNodesPerFace
        for iDof = 2 : mesh.numDofPerNode
          for iDim = 1 : Tdim
            flux[iDof, n] -= Fv_face[iDim, iDof, n]*nrm_xy[iDim,n] 
          end
        end
      end

      for n = 1 : mesh.numNodesPerFace
        for iDof = 2 : mesh.numDofPerNode
          for jDof = 1: mesh.numDofPerNode
            flux[iDof, n] +=  pMat[iDof, jDof, n]*dq[jDof, n]
          end
        end
      end

      # accumulate fluxes
      for n = 1 : mesh.numNodesPerFace
        for iDof = 2 : Tdim+2
          for iDim = 1 : Tdim
            eqn.vecflux_bndry[iDim, iDof, n, f] -=  vecflux[iDim, iDof, n]*coef_nondim
          end
          eqn.bndryflux[iDof, n, f] += flux[iDof, n]*coef_nondim
        end
      end
    end # loop over faces of one BC
  end # loop over BCs
  return nothing 
end


#------------------------------------------------------------------------------
# functions for computing qR

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
              params::ParamType{Tdim},
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
                          params::ParamType{3},
                          qg::AbstractArray{Tsol, 2}) where {Tsol, Tmsh}

  gamma = params.euler_params.gamma
  gamma_1 = params.euler_params.gamma - 1
  sigma = 0.01
  aoa = params.euler_params.aoa
  beta = params.sideslip_angle
  rhoInf = 1.0
  uInf = params.euler_params.Ma * cos(beta) * cos(aoa)
  vInf = params.euler_params.Ma * sin(beta) * -1
  wInf = params.euler_params.Ma * cos(beta) * sin(aoa)
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
                          params::ParamType{2},
                          q_bnd::AbstractArray{Tsol, 2}) where {Tsol, Tmsh}
  @assert( size(q_in, 1) == size(q_bnd,  1))
  @assert( size(q_in, 2) == size(q_bnd,  2))
  @assert( size(q_in, 2) == size(norm,  2))
  dim = size(norm, 1)
  numNodes = size(q_in, 2)
  sigma = 0.1

  gamma = params.euler_params.gamma
  gamma_1 = gamma - 1

  aoa = params.euler_params.aoa
  qRef = zeros(Float64, dim+2)
  qRef[1] = 1.0
  qRef[2] = params.euler_params.Ma*cos(aoa)
  qRef[3] = params.euler_params.Ma*sin(aoa)
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
              params::ParamType{2},
              q_bnd::AbstractArray{Tsol, 2}) where {Tsol, Tmsh}
  @assert( size(q_in, 1) == size(q_bnd,  1))
  @assert( size(q_in, 2) == size(q_bnd,  2))
  @assert( size(q_in, 2) == size(norm,  2))

  dim      = size(norm, 1)
  numNodes = size(q_in, 2)
  lambda   = zeros(Float64, 3) 
  gamma    = params.euler_params.gamma
  gamma_1  = params.euler_params.gamma - 1
  aoa      = params.euler_params.aoa
  gg_1     = (gamma*gamma_1)
  MaInf    = params.euler_params.Ma    

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
    lambda[1] = real(vn)  # TODO: taking the real part is not complex safe!
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
              params::ParamType{3},
              q_bnd::AbstractArray{Tsol, 2}) where {Tsol, Tmsh}
  @assert( size(q_in, 1) == size(q_bnd,  1))
  @assert( size(q_in, 2) == size(q_bnd,  2))
  @assert( size(q_in, 2) == size(norm,  2))
  @assert( size(norm, 1) == 3 )

  dim      = size(norm, 1)
  numNodes = size(q_in, 2)
  lambda   = zeros(Float64, 3) 
  gamma    = params.euler_params.gamma
  gamma_1  = params.euler_params.gamma_1
  aoa      = params.euler_params.aoa
  beta     = params.sideslip_angle
  gg_1     = (gamma*gamma_1)
  MaInf    = params.euler_params.Ma    

  uInf = params.euler_params.Ma * cos(beta) * cos(aoa)
  vInf = params.euler_params.Ma * sin(beta) * -1
  wInf = params.euler_params.Ma * cos(beta) * sin(aoa)
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



#------------------------------------------------------------------------------
# Boundary condition functors
mutable struct nonslipBC <: BCType
end
# low level function
function (obj::nonslipBC)(
              params::ParamType,
              q::AbstractArray{Tsol,1},  
              aux_vars::AbstractArray{Tres, 1},  
              x::AbstractArray{Tmsh,1}, 
              nrm_xy::AbstractArray{Tmsh,1}, 
              bndryflux::AbstractArray{Tres, 1}) where {Tmsh, Tsol, Tres}

  dim = length(nrm_xy)
	qg = params.qg
  # adiabatic wall
	qg[1] = q[1]
	qg[2:dim+1] = 0.0
	qg[dim+2] = q[dim+2]
  # isothermal wall
	# qg[1] = q[1]
	# rhoV2 = (q[2]*q[2] + q[3]*q[3])/q[1]
	# qg[2:dim+1] = 0.0
	# qg[dim+2] = q[4] - 0.5*rhoV2

	v_vals = params.v_vals
	convertFromNaturalToWorkingVars(params, qg, v_vals)
	# this is a problem: q is in conservative variables even if
	# params says we are using entropy variables
	calcEulerFlux(params, v_vals, aux_vars, nrm_xy, bndryflux)

	return nothing
end


mutable struct ExactChannelBC <: BCType
end
# low level function
function (obj::ExactChannelBC)(
              params::ParamType{3},
              q::AbstractArray{Tsol,1},  
              aux_vars::AbstractArray{Tres, 1},  
              xyz::AbstractArray{Tmsh,1}, 
              nrm_xy::AbstractArray{Tmsh,1}, 
              bndryflux::AbstractArray{Tres, 1}) where {Tmsh, Tsol, Tres}

  sigma = 0.01
  gamma = params.euler_params.gamma
  gamma_1 = params.euler_params.gamma - 1
  aoa = params.euler_params.aoa
  beta = params.sideslip_angle
  rhoInf = 1.0
  uInf = params.euler_params.Ma * cos(beta) * cos(aoa)
  vInf = params.euler_params.Ma * sin(beta) * -1
  wInf = params.euler_params.Ma * cos(beta) * sin(aoa)
  TInf = 1.0
  x = xyz[1]
  y = xyz[2]
  z = xyz[3]

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

  qg = Array{Tsol}(5)
	qg[1] = rho
	qg[2] = rho*u
	qg[3] = rho*v
	qg[4] = rho*w
  qg[5] = T/(gamma * gamma_1) + 0.5 * (u*u + v*v + w*w)
  qg[5] *= rho

	v_vals = params.v_vals
	convertFromNaturalToWorkingVars(params, qg, v_vals)
	# this is a problem: q is in conservative variables even if
	# params says we are using entropy variables
	# calcEulerFlux(params, v_vals, aux_vars, [nx2, ny2], bndryflux)
  RoeSolver(params, q, qg, aux_vars, nrm_xy, bndryflux)

	return nothing
end

# low level function
function (obj::ExactChannelBC)(
              params::ParamType{2},
              q::AbstractArray{Tsol,1},  
              aux_vars::AbstractArray{Tres, 1},  
              x::AbstractArray{Tmsh,1}, 
              nrm_xy::AbstractArray{Tmsh,1}, 
              bndryflux::AbstractArray{Tres, 1}) where {Tmsh, Tsol, Tres}

  # functor ExactChannel takes varibales on multiple nodes, so we need to reshape some variables
  xy = reshape(x, length(x), 1)
  norm = reshape(nrm_xy, length(nrm_xy), 1)
  q_in = reshape(q, length(q), 1)
  q_bnd = zeros(Tsol, 4, 1)
  bnd_functor = ExactChannel()
  bnd_functor(q_in, xy, norm, params, q_bnd)
  qg = reshape(q_bnd, length(q_bnd))

	v_vals = params.v_vals
	convertFromNaturalToWorkingVars(params, qg, v_vals)
  RoeSolver(params, q, qg, aux_vars, nrm_xy, bndryflux)

	return nothing
end

mutable struct zeroPressGradientBC <: BCType
end

# low level function
function (obj::zeroPressGradientBC)(
                                params::ParamType,
                                q::AbstractArray{Tsol,1},
                                aux_vars::AbstractArray{Tres, 1},
                                x::AbstractArray{Tmsh,1},
                                nrm_xy::AbstractArray{Tmsh,1},
                                bndryflux::AbstractArray{Tres, 1}) where {Tmsh, Tsol, Tres}


  dim = length(nrm_xy)

	gamma = params.euler_params.gamma
	gamma_1 = params.euler_params.gamma_1
	qg = params.qg
	dim = 2
  rhoV2 = (norm(view(q, 2:dim+1))) / q[1]
	# rhoV2 = (q[2]*q[2] + q[3]*q[3]) / q[1]
	pinf = 1./gamma
	qg[1:dim+1] = q[1:dim+1]
	qg[dim+2] = pinf/gamma_1 + 0.5*rhoV2

	v_vals = params.v_vals
	convertFromNaturalToWorkingVars(params, qg, v_vals)
	# this is a problem: q is in conservative variables even if
	# params says we are using entropy variables
	calcEulerFlux(params, v_vals, aux_vars, nrm_xy, bndryflux)

	return nothing
end

global const BCDict = Dict{String, BCType}(
  "nonslipBC" => nonslipBC(),
  "ExactChannelBC" => ExactChannelBC(),
  "zeroPressGradientBC" => zeroPressGradientBC(),
)

"""
### NavierStokesnMod.getBCFunctors

  This function uses the opts dictionary to populate mesh.bndry_funcs with
  the functors

    func(params::ParamType,
         q::AbstractArray{Tsol,1},
         aux_vars::AbstractArray{Tres, 1},  coords::AbstractArray{Tmsh,1},
         nrm_xy::AbstractArray{Tmsh,1},
         bndryflux::AbstractArray{Tres, 1},
         bndry::BoundaryNode=NullBoundaryNode)


  This is a high level function.
"""
# use this function to populate access the needed values in BCDict
function getBCFunctors(mesh::AbstractMesh, sbp::AbstractOperator, eqn::NSData, opts)

  # the boundary conditions should switch over to the same model as the 
  # other physics modules, but for now `calcViscousFlux_boundary` internally
  # picks the right functor
  # this will likely require each `eqn` object having its own `bndry_funcs`
  # array, rather than the mesh having a single one for all physics
  # In the meantime, `calcViscousFlux_boundary` is getting its functors
  # internally, and the invisicd functors are stored in the mesh.
  #
  # As of this writing, I believe the only boundary condition that works is
  # the FreeStreamBC, because it exists for both viscous and inviscid.
  # Eventually the inviscid part of every viscous boundary condition will
  # need to exist and both boundary conditions applied (related to the
  # comment above about each eqn object having its own BC functor array.
  error("getBCFunctors is not (currently) used for $PhysicsName")
  for i=1:mesh.numBC
    key_i = string("BC", i, "_name")
    val = opts[key_i]
    mesh.bndry_funcs[i] = BCDict[val]
  end

  return nothing
end # end function getBCFunctors


