# this file contains all functions related to performing edge stabilization
# edge stabilization is is executed from euler.jl


@doc """
### EulerEquationMod.stabscale

  This function calculates the edge stabilization scaling paramater at a 
  node and returns it

  Inputs:
    * u  : vector of conservative variables
    * dxidx : jacobian of xi wrt x coordinates at the node
    * nrm  : normal vector in xi space
    * mesh : AbstractMesh (only needed for order)
    * params : ParamType{2}

  This is a low level function
"""->
# this function is going to be deprecated soon
# low level function
function stabscale{T}(u::AbstractArray{T,1}, dxidx::AbstractArray{T,2}, nrm::AbstractArray{T,1}, mesh::AbstractMesh, params::ParamType{2})

#     println("==== entering stabscale ====")

    # grabbing conserved variables
    rho = u[1]
    vel_x = u[2]/rho
    vel_y = u[3]/rho
    Energy = u[4]

    # from JC's code below, eqn should still be in scope
    pressure = calcPressure(u, eqn.params)

    # solved eqn for e: E = rho*e + (1/2)*rho*u^2
    vel_squared = vel_x^2 + vel_y^2
    energy = Energy/rho - (1/2)*vel_squared

    # gamma stored in EulerData type
    gamma = eqn.gamma

#     println("pressure: ",pressure)
#     println("gamma: ",gamma)
#     println("rho: ",rho)
    # ideal gas law
    speed_sound = sqrt((gamma*pressure)/rho)

    # choice for edge stabilization constant: 
    #   refer to email from JH, 20150504:
    #   Anthony: there is little guidance in the literature regarding 
    #     gamma for the Euler equations.  I suggest starting with 
    #     gamma = 0.01.  If that fails (with a cfl of 1.0), then decrease 
    #     it by an order of magnitude at at time until RK is stable.  
    #     Once you find a value that works, try increasing it slowly.
    

     edge_stab_gamma = -0.01  # default
    #edge_stab_gamma = 0.0 
#     edge_stab_gamma = 0.00001

    # edge lengths component wise
    h_x = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
    h_y = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]

    # edge length
    h = sqrt(h_x^2 + h_y^2)

    # scaled velocity scalar
#     U = vel_x*(nrm[1]/h) + vel_y*(nrm[2]/h)
    U = vel_x*(h_x/h) + vel_y*(h_y/h)

    order = mesh.order
#     return (U + speed_sound)*edge_stab_gamma*h^2
    return (abs(U) + speed_sound)*edge_stab_gamma*h^(order+1)

  end


@doc """
### EulerEquationMod.edgestabilize!

Adds edge-stabilization (see Burman and Hansbo, doi:10.1016/j.cma.2003.12.032)
to a given residual.  Different methods are available depending on the rank of
`u`:

* For *scalar* fields, it is assumed that `u` is a rank-2 array, with the first
dimension for the local-node index, and the second dimension for the element
index.
* For *vector* fields, `u` is a rank-3 array, with the first dimension for the
index of the vector field, the second dimension for the local-node index, and
the third dimension for the element index.

Naturally, the number of entries in the dimension of `u` (and `res`)
corresponding to the nodes must be equal to the number of nodes in the SBP
operator `sbp`.

**Inputs**

* `sbp`: an SBP operator type
* `ifaces`: list of element interfaces stored as an array of `Interface`s
* `u`: the array of solution data
* `x`: Cartesian coordinates stored in (coord,node,element) format
* `dξdx`: scaled Jacobian of the mapping (as output from `mappingjacobian!`)
* `jac`: determinant of the Jacobian
* `α`: array of transformation terms (see below)
* `stabscale`: function to compute the edge-stabilization scaling (see below)

**In/Outs**

* `res`: where the result of the integration is stored

**Details**

The array `α` is used to compute the directional derivative normal to the faces.
For a 2-dimensional problem, it can be computed as follows:

      for k = 1:mesh.numelem
        for i = 1:sbp.numnodes
          for di1 = 1:2
            for di2 = 1:2
              α[di1,di2,i,k] = (dξdx[di1,1,i,k].*dξdx[di2,1,i,k] + 
                                dξdx[di1,2,i,k].*dξdx[di2,2,i,k])*jac[i,k]
            end
          end
        end
      end

The function `stabscale` has the signature

  function stabscale(u, dξdx, nrm)

where `u` is the solution at a node, `dξdx` is the (scaled) Jacobian at the same
node, and `nrm` is the normal to the edge in reference space.  `stabscale`
should return the necessary scaling to ensure the edge-stabilization has the
desired asymptotic convergence rate.

"""->
function edgestabilize!{T}(sbp::SBPOperator{T}, ifaces::Array{Interface},
                           u::AbstractArray{T,2}, x::AbstractArray{T,3},
                           dξdx::AbstractArray{T,4}, jac::AbstractArray{T,2},
                           α::AbstractArray{T,4}, stabscale::Function,
                           res::AbstractArray{T,2})
  @assert( sbp.numnodes == size(u,1) == size(res,1) == size(dξdx,3) == size(x,2) 
          == size(α,3) )
  @assert( size(dξdx,4) == size(α,4) == size(u,2) == size(res,2) == size(x,3) )
  @assert( length(u) == length(res) )
  dim = size(sbp.Q, 3)


  # JEH: temporary, until nbrnodeindex is part of sbp type
  nbrnodeindex = Array(sbp.numfacenodes:-1:1)

  Dn = zero(T)
  dirL = zeros(T, (dim))
  dirR = zeros(T, (dim))
  workvec = zeros(T, (dim))
  tmpL = zero(T)
  tmpR = zero(T)
  EDn = zeros(T, (sbp.numfacenodes) )
  for face in ifaces
    fill!(EDn, zero(T))
    for i = 1:sbp.numfacenodes
      # iL = element-local index for ith node on left element face
      # iR = element-local index for ith node on right element face
      iL = sbp.facenodes[i, face.faceL]::Int
      #iR = sbp.facenodes[getnbrnodeindex(sbp, face, i)::Int, face.faceR]::Int
      iR = sbp.facenodes[nbrnodeindex[i], face.faceR]::Int
      # apply the normal-derivative difference operator along the face
      smallmatvec!(view(α,:,:,iL,face.elementL), view(sbp.facenormal,:,face.faceL), dirL)
      Dn = zero(T)
      Dn = directionaldifferentiate!(sbp, dirL, view(u,:,face.elementL), iL)
      smallmatvec!(view(α,:,:,iR,face.elementR), view(sbp.facenormal,:,face.faceR), dirR)
      Dn += directionaldifferentiate!(sbp, dirR, view(u,:,face.elementR), iR)
      # get differential area element: need 1/ds for each Dn term (here and loop
      # below)to get unit normal, and then need ds for integration, so net
      # result is 1/ds
      ds = getdiffelementarea(view(sbp.facenormal,:,face.faceL),
                              view(dξdx,:,:,iL,face.elementL), workvec)::T
      # apply the scaling function
      Dn *= stabscale(u[iL,face.elementL], view(dξdx,:,:,iL,face.elementL),
                      view(sbp.facenormal,:,face.faceL))::T/ds # note that u[iL] = u[iR]
      # add the face-mass matrix contribution
      for j = 1:sbp.numfacenodes
        EDn[j] += sbp.wface[j,i]*Dn
      end
    end
    # here we use hand-coded reverse-mode to apply the transposed
    # normal-derivative difference operator
    for i = 1:sbp.numfacenodes
      iL = sbp.facenodes[i, face.faceL]::Int
      #iR = sbp.facenodes[getnbrnodeindex(sbp, face, i), face.faceR]::Int
      iR = sbp.facenodes[nbrnodeindex[i], face.faceR]::Int
      smallmatvec!(view(α,:,:,iL,face.elementL), view(sbp.facenormal,:,face.faceL), dirL)
      smallmatvec!(view(α,:,:,iR,face.elementR), view(sbp.facenormal,:,face.faceR), dirR)
      for di = 1:dim
        tmpL = dirL[di]*EDn[i]/sbp.w[iL]
        tmpR = dirR[di]*EDn[i]/sbp.w[iR]
        for j = 1:sbp.numnodes
          res[j,face.elementL] += sbp.Q[iL,j,di]*tmpL
          res[j,face.elementR] += sbp.Q[iR,j,di]*tmpR
        end
      end
    end
  end
end



function getdiffelementarea{T, T2, T3}(nrm::AbstractArray{T,1}, dxidx::AbstractArray{T2,2},
                                 workvec::AbstractArray{T3,1})
				   fill!(workvec, zero(T3))
				     for di1 = 1:size(nrm,1)
				           for di2 = 1:size(nrm,1)
					           workvec[di2] += nrm[di1]*dxidx[di1,di2]
						       end
						         end
							   return norm(workvec)
							 end



# for vector equations
function edgestabilize!{Tmsh,  Tsol, Tres}(sbp::SBPOperator, ifaces::Array{Interface},
                           u::AbstractArray{Tsol,3}, x::AbstractArray{Tmsh,3},
                           dξdx::AbstractArray{Tmsh,4}, jac::AbstractArray{Tmsh,2},
                           α::AbstractArray{Tmsh,4},stabscale::AbstractArray{Tres,2},
                           res::AbstractArray{Tres,3})

  @assert( sbp.numnodes == size(u,2) == size(res,2) == size(dξdx,3) == size(x,2) 
          == size(α,3) )
  @assert( size(dξdx,4) == size(α,4) == size(u,3) == size(res,3) == size(x,3) )
  @assert( length(u) == length(res) )
  dim = size(sbp.Q, 3)



  # JEH: temporary, until nbrnodeindex is part of sbp type
  nbrnodeindex = Array(sbp.numfacenodes:-1:1)

  Dn = zeros(Tsol, size(u,1))
  dirL = zeros(Tmsh, (dim))
  dirR = zeros(Tmsh, (dim))
  workvec = zeros(Tmsh, (dim))
  tmpL = zero(Dn)
  tmpR = zero(Dn)
  EDn = zeros(Tres, (size(u,1),sbp.numfacenodes) )
  for (facenum, face) in enumerate(ifaces)
#   for facenum = 1:length(ifaces)
#    face = ifaces[facenum]
    fill!(EDn, zero(Tres))
    for i = 1:sbp.numfacenodes
      # iL = element-local index for ith node on left element face
      # iR = element-local index for ith node on right element face
      iL = sbp.facenodes[i, face.faceL]
      #iR = sbp.facenodes[getnbrnodeindex(sbp, face, i), face.faceR]
      iR = sbp.facenodes[nbrnodeindex[i], face.faceR]
      # apply the normal-derivative difference operator along the face
      smallmatvec!(view(α,:,:,iL,face.elementL),view(sbp.facenormal,:,face.faceL), dirL)
      fill!(Dn, zero(Tres))
      directionaldifferentiate!(sbp, dirL, view(u,:,:,face.elementL), iL, Dn)
      smallmatvec!(view(α,:,:,iR,face.elementR), view(sbp.facenormal,:,face.faceR), dirR)
      directionaldifferentiate!(sbp, dirR, view(u,:,:,face.elementR), iR, Dn)
      # get differential area element: need 1/ds for each Dn term (here and loop
      # below) to get unit normals, and then need ds for integration, so net
      # result is 1/ds
      ds = getdiffelementarea(view(sbp.facenormal,:,face.faceL),
                              view(dξdx,:,:,iL,face.elementL), workvec)  # this assumes Tsbp is a lower type than the other
      # apply the scaling function
      scale = stabscale[i, facenum]
#      scale = stabscale(view(u,:,iL,face.elementL), view(dξdx,:,:,iL,face.elementL),
#                         view(sbp.facenormal,:,face.faceL))::T./ds # note that u[iL] = u[iR]
      for field = 1:size(u,1)
        Dn[field] *= scale
      end
      # add the face-mass matrix contribution
      for j = 1:sbp.numfacenodes
        for field = 1:size(u,1)
          EDn[field,j] += sbp.wface[j,i]*Dn[field]
        end
      end
    end
    for i = 1:sbp.numfacenodes
      iL = sbp.facenodes[i, face.faceL]
      #iR = sbp.facenodes[getnbrnodeindex(sbp, face, i), face.faceR]
      iR = sbp.facenodes[nbrnodeindex[i], face.faceR]
      smallmatvec!(view(α,:,:,iL,face.elementL), view(sbp.facenormal,:,face.faceL), dirL)
      smallmatvec!(view(α,:,:,iR,face.elementR), view(sbp.facenormal,:,face.faceR), dirR)
      # here we use hand-coded reverse-mode to apply the transposed
      # normal-derivative difference operator
      for di = 1:size(sbp.Q, 3)
        for field = 1:size(u,1)
          tmpL[field] = dirL[di]*EDn[field,i]/sbp.w[iL]
          tmpR[field] = dirR[di]*EDn[field,i]/sbp.w[iR]
        end
        for j = 1:sbp.numnodes
          for field = 1:size(u,1)
            res[field,j,face.elementL] += sbp.Q[iL,j,di]*tmpL[field]
            res[field,j,face.elementR] += sbp.Q[iR,j,di]*tmpR[field]
          end
        end
      end
    end
  end

end





@doc """
### EulerEquationMod.stabscale

  This function calculates the edge stabilization scalaing parameter at a node
  and returns it. Linear elements only.

   Inputs:
    * u : vector of conservative variables
    * dxidx : jacobian of xi wrt x coordinates at the node
    * nrm : normal vector in xi space
    * params : ParamType{2}

    This is a low level function
"""->
# low level function
function stabscale{Tmsh,  Tsol}(u::AbstractArray{Tsol,1}, dxidx::AbstractArray{Tmsh,2}, nrm::AbstractArray{Tmsh,1}, params::ParamType{2} )
# calculate stabscale for a single node

#     println("==== entering stabscale ====")

    # grabbing conserved variables
    rho = u[1]
    vel_x = u[2]/rho
    vel_y = u[3]/rho
    Energy = u[4]

    # from JC's code below, eqn should still be in scope
    pressure = calcPressure(u, params)

    # solved eqn for e: E = rho*e + (1/2)*rho*u^2
    vel_squared = vel_x^2 + vel_y^2
    energy = Energy/rho - (1/2)*vel_squared

    # gamma stored in EulerData type
    gamma = params.gamma

#     println("pressure: ",pressure)
#     println("gamma: ",gamma)
#     println("rho: ",rho)
    # ideal gas law
    speed_sound = sqrt((gamma*pressure)/rho)

    # choice for edge stabilization constant: 
    #   refer to email from JH, 20150504:
    #   Anthony: there is little guidance in the literature regarding 
    #     gamma for the Euler equations.  I suggest starting with 
    #     gamma = 0.01.  If that fails (with a cfl of 1.0), then decrease 
    #     it by an order of magnitude at at time until RK is stable.  
    #     Once you find a value that works, try increasing it slowly.

    edge_stab_gamma = params.edgestab_gamma  # default
    #edge_stab_gamma = 0.0 
#     edge_stab_gamma = 0.00001

    # edge lengths component wise
    h_x = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
    h_y = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]

    # edge length
    h = sqrt(h_x^2 + h_y^2)
#    println("h = ", h)

    # scaled velocity scalar
#     U = vel_x*(nrm[1]/h) + vel_y*(nrm[2]/h)
    U = vel_x*(h_x/h) + vel_y*(h_y/h)

#     return (U + speed_sound)*edge_stab_gamma*h^2
    return (abs(U) + speed_sound)*edge_stab_gamma*h^(2)

  end


@doc """
### EulerEquationMod.stabscale

  This function calculate the stabilization scaling parameter across the
  entire mesh by calling the low level method.

  This is a mid level function
"""->
# mid level function
function stabscale{Tmsh,  Tsol}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator, eqn::EulerData{Tsol})
# calculate stabscale for entire mesh

 nbrnodeindex = Array(sbp.numfacenodes:-1:1)

  for i=1:mesh.numInterfaces
    face_i = mesh.interfaces[i]
    for j=1:sbp.numfacenodes
      iL = sbp.facenodes[j, face_i.faceL]
      iR = sbp.facenodes[nbrnodeindex[j], face_i.faceR]
      q = view(eqn.q, :, iL, face_i.elementL)
      dxidx = view(mesh.dxidx, :, :, iL, face_i.elementL)
      nrm = view(sbp.facenormal, :, face_i.faceL)
      
      eqn.stabscale[j,i] = stabscale(q, dxidx, nrm, eqn.params)
    end
  end

  return nothing
end


@doc """
### EulerEquationMod.caclEdgeStabAlpha

  This function calculates the edge stabilization paramter alpha across the
  entire mesh.

  This is a mid level function.
"""
# used by EulerData Constructor - not that that matters for any reason
# mid level function
function calcEdgeStabAlpha{Tmsh,  Tsol, Tdim}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator, eqn::EulerData{Tsol, Tdim})
# calculate alpha, needed by edge stabilization


  numEl = mesh.numEl
  eqn.edgestab_alpha = Array(Tmsh,2,2,sbp.numnodes,numEl)
  dxidx = mesh.dxidx
  jac = mesh.jac

  # calculating alpha, required by edgestabilize!
  # this canbe compuated apriori
  for k = 1:numEl
    for i = 1:sbp.numnodes
      for di1 = 1:Tdim
        for di2 = 1:Tdim
          eqn.edgestab_alpha[di1,di2,i,k] = (dxidx[di1,1,i,k].*dxidx[di2,1,i,k] + dxidx[di1,2,i,k].*dxidx[di2,2,i,k])*jac[i,k]
        end
      end
    end
  end

  return nothing
end




