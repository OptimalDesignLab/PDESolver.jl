# file stabilization.jl
# this is the primary file for adding stabilization to the weak form.
# Stabilization functions may be put in separate files if they are included 
# here

include("SUPG.jl")


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
* `u`: the array of solution data, numDofPerNode x numNodesPerElement x numEl.
       for scalar fields (ie. numDofPerNode), this can be a 2D array.
* `x`: Cartesian coordinates stored in (coord,node,element) format
* `dxidx`: scaled Jacobian of the mapping (as output from `mappingjacobian!`)
* `jac`: determinant of the Jacobian array, numNodesPerElement x numEl
* `alpha`: array of transformation terms (see below)
* `stabscale`: numfaces x numNodes per face array of scaling parameter

**In/Outs**

* `res`: where the result of the integration is stored, same shape as u

**Details**

The array `alpha` is used to compute the directional derivative normal to the faces.
For a 2-dimensional problem, it can be computed as follows:

      for k = 1:mesh.numelem
        for i = 1:sbp.numnodes
          for di1 = 1:2
            for di2 = 1:2
              alpha[di1,di2,i,k] = (dxidx[di1,1,i,k].*dxidx[di2,1,i,k] + 
                                dxidx[di1,2,i,k].*dxidx[di2,2,i,k])*jac[i,k]
            end
          end
        end
      end

"""->
function edgestabilize!{T, Tres}(sbp::SBPOperator{T}, ifaces::Array{Interface},
                           u::AbstractArray{T,2}, x::AbstractArray{T,3},
                           dxidx::AbstractArray{T,4}, jac::AbstractArray{T,2},
                           alpha::AbstractArray{T,4}, 
                           stabscale::AbstractArray{Tres, 2},
                           res::AbstractArray{T,2})
# for scalar equations only!
  @assert( sbp.numnodes == size(u,1) == size(res,1) == size(dxidx,3) == size(x,2) 
          == size(alpha,3) )
  @assert( size(dxidx,4) == size(alpha,4) == size(u,2) == size(res,2) == size(x,3) )
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
  for (facenum, face) in enumerate(ifaces)
    fill!(EDn, zero(T))
    for i = 1:sbp.numfacenodes
      # iL = element-local index for ith node on left element face
      # iR = element-local index for ith node on right element face
      iL = sbp.facenodes[i, face.faceL]::Int
      #iR = sbp.facenodes[getnbrnodeindex(sbp, face, i)::Int, face.faceR]::Int
      iR = sbp.facenodes[nbrnodeindex[i], face.faceR]::Int
      # apply the normal-derivative difference operator along the face
      smallmatvec!(view(alpha,:,:,iL,face.elementL), 
                   view(sbp.facenormal,:,face.faceL), dirL)
      Dn = zero(T)
      Dn = directionaldifferentiate!(sbp, dirL, view(u,:,face.elementL), iL)
      smallmatvec!(view(alpha,:,:,iR,face.elementR), 
                   view(sbp.facenormal,:,face.faceR), dirR)
      Dn += directionaldifferentiate!(sbp, dirR, view(u,:,face.elementR), iR)
      # get differential area element: need 1/ds for each Dn term (here and loop
      # below)to get unit normal, and then need ds for integration, so net
      # result is 1/ds
      ds = calcDiffElementArea(view(sbp.facenormal,:,face.faceL),
                              view(dxidx,:,:,iL,face.elementL), workvec)::T
      # apply the scaling function

      scale = stabscale[i, facenum]
#      Dn *= stabscale(u[iL,face.elementL], view(dxidx,:,:,iL,face.elementL),
#                      view(sbp.facenormal,:,face.faceL))::T/ds # note that u[iL] = u[iR]
      Dn *= scale/ds
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
      smallmatvec!(view(alpha,:,:,iL,face.elementL), 
                   view(sbp.facenormal,:,face.faceL), dirL)
      smallmatvec!(view(alpha,:,:,iR,face.elementR), 
                   view(sbp.facenormal,:,face.faceR), dirR)
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


# for vector equations

function edgestabilize!{Tmsh,  Tsol, Tres}(sbp::SBPOperator, 
                        ifaces::Array{Interface}, u::AbstractArray{Tsol,3}, 
                        x::AbstractArray{Tmsh,3}, dxidx::AbstractArray{Tmsh,4},
                        jac::AbstractArray{Tmsh,2}, 
                        alpha::AbstractArray{Tmsh,4},
                        stabscale::AbstractArray{Tres,2},
                        res::AbstractArray{Tres,3})

  @assert( sbp.numnodes == size(u,2) == size(res,2) == size(dxidx,3) == size(x,2) 
          == size(alpha,3) )
  @assert( size(dxidx,4) == size(alpha,4) == size(u,3) == size(res,3) == size(x,3) )
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
    fill!(EDn, zero(Tres))
    for i = 1:sbp.numfacenodes
      # iL = element-local index for ith node on left element face
      # iR = element-local index for ith node on right element face
      iL = sbp.facenodes[i, face.faceL]
      #iR = sbp.facenodes[getnbrnodeindex(sbp, face, i), face.faceR]
      iR = sbp.facenodes[nbrnodeindex[i], face.faceR]

      # apply the normal-derivative difference operator along the face
      smallmatvec!(view(alpha,:,:,iL,face.elementL),
                   view(sbp.facenormal,:,face.faceL), dirL)

      fill!(Dn, zero(Tres))
      directionaldifferentiate!(sbp, dirL, view(u,:,:,face.elementL), iL, Dn)
      smallmatvec!(view(alpha,:,:,iR,face.elementR), 
                   view(sbp.facenormal,:,face.faceR), dirR)
      directionaldifferentiate!(sbp, dirR, view(u,:,:,face.elementR), iR, Dn)
      # get differential area element: need 1/ds for each Dn term (here and loop
      # below) to get unit normals, and then need ds for integration, so net
      # result is 1/ds
      ds = calcDiffElementArea(view(sbp.facenormal,:,face.faceL),
                              view(dxidx,:,:,iL,face.elementL), workvec)

      # apply the scaling function
      scale = stabscale[i, facenum]
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
      iR = sbp.facenodes[nbrnodeindex[i], face.faceR]
      smallmatvec!(view(alpha,:,:,iL,face.elementL), 
                   view(sbp.facenormal,:,face.faceL), dirL)
      smallmatvec!(view(alpha,:,:,iR,face.elementR), 
                   view(sbp.facenormal,:,face.faceR), dirR)

      # here we use hand-coded reverse-mode to apply the transposed
      # normal-derivative difference operator
      for di = 1:size(sbp.Q, 3)
        for field = 1:size(u,1)
          tmpL[field] = dirL[di]*EDn[field,i]/sbp.w[iL]
          tmpR[field] = dirR[di]*EDn[field,i]/sbp.w[iR]
        end
	elementL = face.elementL
	elementR = face.elementR
        for j = 1:sbp.numnodes
          for field = 1:size(u,1)
	    tmp2 =  sbp.Q[iL,j,di]*tmpL[field]
            res[field,j,face.elementL] += sbp.Q[iL,j,di]*tmpL[field]
	    tmp2 = sbp.Q[iR,j,di]*tmpR[field]
            res[field,j,face.elementR] += sbp.Q[iR,j,di]*tmpR[field]
          end
        end
      end
    end
  end

end


# for vector equations
# WIP: uses res_edge
#TODO: cleanup function signature
function edgestabilize!{Tmsh,  Tsol, Tres}(mesh, sbp::SBPOperator, eqn, 
                        ifaces::Array{Interface}, u::AbstractArray{Tsol,3}, 
                        x::AbstractArray{Tmsh,3}, dxidx::AbstractArray{Tmsh,4}, 
                        jac::AbstractArray{Tmsh,2}, 
                        alpha::AbstractArray{Tmsh,4},
                        stabscale::AbstractArray{Tres,2},
                        res::AbstractArray{Tres, 3}, res_edge::AbstractArray{Tres,4})

  @assert( sbp.numnodes == size(u,2) == size(res,2) == size(dxidx,3) == size(x,2) 
          == size(alpha,3) )
  @assert( size(dxidx,4) == size(alpha,4) == size(u,3) == size(res,3) == size(x,3) )
#  @assert( length(u) == length(res) )
# res is the residual array for element affecting themselves
# res_edge is the residual array of how elements affect each other
# this won't work for finite differences
# should there be a way to detect if this is a residual evaluation or a 
# differentiation and use the 2x faster version of this function?
  dim = size(sbp.Q, 3)

#  println("initially:")
#  println("res = \n", res)
#  println("res_edge = \n", res_edge)


  # JEH: temporary, until nbrnodeindex is part of sbp type
  nbrnodeindex = Array(sbp.numfacenodes:-1:1)

  Dn = zeros(Tsol, size(u,1))
  dirL = zeros(Tmsh, (dim))
  dirR = zeros(Tmsh, (dim))
  workvec = zeros(Tmsh, (dim))
  tmpL = zero(Dn)
  tmpR = zero(Dn)
  EDn = zeros(Tres, (size(u,1),sbp.numfacenodes) )

  # consider how elementR affects elementL and itself
  for (facenum, face) in enumerate(ifaces)
#   for facenum = 1:length(ifaces)
#    face = ifaces[facenum]
    fill!(EDn, zero(Tres))
    uL = real(view(u, :, :, face.elementL))
    uR = view(u, :, :, face.elementR)
#    println("elementL = ", face.elementL, ", elementR = ", face.elementR)
#    println("uL = ", uL, ", uR = ", uR) 
    for i = 1:sbp.numfacenodes  # consider making this its own function
      # iL = element-local index for ith node on left element face
      # iR = element-local index for ith node on right element face
      iL = sbp.facenodes[i, face.faceL]
      #iR = sbp.facenodes[getnbrnodeindex(sbp, face, i), face.faceR]
      iR = sbp.facenodes[nbrnodeindex[i], face.faceR]
      # apply the normal-derivative difference operator along the face
      smallmatvec!(view(alpha,:,:,iL,face.elementL),
                   view(sbp.facenormal,:,face.faceL), dirL)
      fill!(Dn, zero(Tres))
      directionaldifferentiate!(sbp, dirL, uL, iL, Dn)
      smallmatvec!(view(alpha,:,:,iR,face.elementR), 
                   view(sbp.facenormal,:,face.faceR), dirR)
      directionaldifferentiate!(sbp, dirR, uR, iR, Dn)
      # get differential area element: need 1/ds for each Dn term (here and loop
      # below) to get unit normals, and then need ds for integration, so net
      # result is 1/ds
      ds = calcDiffElementArea(view(sbp.facenormal,:,face.faceL),
                              view(dxidx,:,:,iL,face.elementL), workvec)  # this assumes Tsbp is a lower type than the other
      # apply the scaling function
      # calculate scale
      #=
      iL = sbp.facenodes[j, face_i.faceL]
      iR = sbp.facenodes[nbrnodeindex[j], face_i.faceR]
      q = view(eqn.q, :, iL, face_i.elementL)
      dxidx = view(mesh.dxidx, :, :, iL, face_i.elementL)
      nrm = view(sbp.facenormal, :, face_i.faceL)
      =#
      uL_node = real(view(eqn.q, :, iL, face.elementL)) # take real part here?
      dxidx_node = view(mesh.dxidx, :, :, iL, face.elementL)
      nrm_node = view(sbp.facenormal, :, face.faceL)
 
      scale = stabscale(uL_node, dxidx_node, nrm_node, eqn.params)
#      scale = stabscale[i, facenum]
#      scale = stabscale(view(u,:,iL,face.elementL), view(dxidx,:,:,iL,face.elementL),
#                         view(sbp.facenormal,:,face.faceL))::T./ds # note that u[iL] = u[iR]
#      println("before scaling, Dn = ", Dn)
#      println("scale = ", scale)
      for field = 1:size(u,1)
        Dn[field] *= scale
      end

#      println("after scaling Dn = ", Dn)
      # add the face-mass matrix contribution
      for j = 1:sbp.numfacenodes
        for field = 1:size(u,1)
          EDn[field,j] += sbp.wface[j,i]*Dn[field]
        end
      end
    end  # end loop over face nodes

#    println("EDn = ", EDn)

    for i = 1:sbp.numfacenodes
#      println("sbp.facenodes = ", sbp.facenodes)

      iL = sbp.facenodes[i, face.faceL]
      #iR = sbp.facenodes[getnbrnodeindex(sbp, face, i), face.faceR]
      iR = sbp.facenodes[nbrnodeindex[i], face.faceR]
      smallmatvec!(view(alpha,:,:,iL,face.elementL), 
                   view(sbp.facenormal,:,face.faceL), dirL)
      smallmatvec!(view(alpha,:,:,iR,face.elementR), 
                   view(sbp.facenormal,:,face.faceR), dirR)
#      println("dirL = ", dirL)
#      println("dirR = ", dirR)
      # here we use hand-coded reverse-mode to apply the transposed
      # normal-derivative difference operator
      for di = 1:size(sbp.Q, 3)
        for field = 1:size(u,1)
          tmpL[field] = dirL[di]*EDn[field,i]/sbp.w[iL]
          tmpR[field] = dirR[di]*EDn[field,i]/sbp.w[iR]
        end

#	println("tmpL = ", tmpL)
#	println("tmpR = ", tmpR)
	elementL = face.elementL
	elementR = face.elementR
	faceL = face.faceL
	faceR = face.faceR
        for j = 1:sbp.numnodes
          for field = 1:size(u,1)
	    # this is elementR affecting elementL
            tmp2 = sbp.Q[iL,j,di]*tmpL[field]
#            println("res_edge[$field, $j, $elementL, $faceL] += ", tmp2)

            res_edge[field,j,face.elementL, face.faceL] += sbp.Q[iL,j,di]*tmpL[field]
	    # this is elementR affecting itself
            tmp2 = sbp.Q[iR,j,di]*tmpR[field]
#            println("res[$field, $j, $elementR] += ", tmp2)


            res[field,j,face.elementR] += sbp.Q[iR,j,di]*tmpR[field]
	  end  # end loop over fields
	end  # end loop j=1:sbp.numnodes
      end  # end loop over directions di
    end  # end loop over i = 1:sbp.numfacenodes
  end  # end loop over interfaces

#  println("res = \n", res)
#  println("res_edge = \n", res_edge)
#  println("edgestabilization second calculation")
  # now consider how elementL affects elementR, and itself
  for (facenum, face) in enumerate(ifaces)
#   for facenum = 1:length(ifaces)
#    face = ifaces[facenum]
    fill!(EDn, zero(Tres))
    uL = view(u, :, :, face.elementL)
    uR = real(view(u, :, :, face.elementR))
#    println("elementL = ", face.elementL, ", elementR = ", face.elementR)
#    println("uL = ", uL, ", uR = ", uR) 
 
    for i = 1:sbp.numfacenodes  # consider making this its own function
      # iL = element-local index for ith node on left element face
      # iR = element-local index for ith node on right element face
      iL = sbp.facenodes[i, face.faceL]
      #iR = sbp.facenodes[getnbrnodeindex(sbp, face, i), face.faceR]
      iR = sbp.facenodes[nbrnodeindex[i], face.faceR]
      # apply the normal-derivative difference operator along the face
      smallmatvec!(view(alpha,:,:,iL,face.elementL),
                   view(sbp.facenormal,:,face.faceL), dirL)
      fill!(Dn, zero(Tres))
      directionaldifferentiate!(sbp, dirL, uL, iL, Dn)
      smallmatvec!(view(alpha,:,:,iR,face.elementR), 
                   view(sbp.facenormal,:,face.faceR), dirR)
      directionaldifferentiate!(sbp, dirR, uR, iR, Dn)
      # get differential area element: need 1/ds for each Dn term (here and loop
      # below) to get unit normals, and then need ds for integration, so net
      # result is 1/ds
      ds = calcDiffElementArea(view(sbp.facenormal,:,face.faceL),
                              view(dxidx,:,:,iL,face.elementL), workvec)
      # apply the scaling function
      uL_node = view(eqn.q, :, iL, face.elementL)
      dxidx_node = view(mesh.dxidx, :, :, iL, face.elementL)
      nrm_node = view(sbp.facenormal, :, face.faceL)
 
      scale = stabscale(uL_node, dxidx_node, nrm_node, eqn.params)

#      scale = stabscale[i, facenum]
#      scale = stabscale(view(u,:,iL,face.elementL), view(dxidx,:,:,iL,face.elementL),
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
    end  # end loop over face nodes


    for i = 1:sbp.numfacenodes
#      println("sbp.facenodes = ", sbp.facenodes)

      iL = sbp.facenodes[i, face.faceL]
      #iR = sbp.facenodes[getnbrnodeindex(sbp, face, i), face.faceR]
      iR = sbp.facenodes[nbrnodeindex[i], face.faceR]
      smallmatvec!(view(alpha,:,:,iL,face.elementL), 
                   view(sbp.facenormal,:,face.faceL), dirL)
      smallmatvec!(view(alpha,:,:,iR,face.elementR), 
                   view(sbp.facenormal,:,face.faceR), dirR)
      # here we use hand-coded reverse-mode to apply the transposed
      # normal-derivative difference operator
      for di = 1:size(sbp.Q, 3)
        for field = 1:size(u,1)
          tmpL[field] = dirL[di]*EDn[field,i]/sbp.w[iL]
          tmpR[field] = dirR[di]*EDn[field,i]/sbp.w[iR]
        end

	elementL = face.elementL
	elementR = face.elementR
	faceL = face.faceL
	faceR = face.faceR
        for j = 1:sbp.numnodes
          for field = 1:size(u,1)
	    # this is elementL affecting itself
            tmp2 = sbp.Q[iL,j,di]*tmpL[field]

#            println("res[$field, $j, $elementL] += ", tmp2)

            res[field,j,face.elementL] += sbp.Q[iL,j,di]*tmpL[field]
            tmp2 = sbp.Q[iR,j,di]*tmpR[field]
#            println("res_edge[$field, $j, $elementR, $faceR] += ", tmp2)
	    # this is elementL affecting elementR
            res_edge[field,j,face.elementR, face.faceR] += sbp.Q[iR,j,di]*tmpR[field]
	  end  # end loop over fields
	end  # end loop j=1:sbp.numnodes
      end  # end loop over directions di
    end  # end loop over i = 1:sbp.numfacenodes
  end  # end loop over interfaces

#  println("finished second edge stabilization calculation")
#  println("res = \n", res)
#  println("res_edge = \n", res_edge)

end  # end function





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

# u holds whatever the variable type that we are solving the equation in
# q_vals holds the conservative variables

#     println("==== entering stabscale ====")a
    q_vals = params.q_vals
    # convert to conservative variables if not already using them
    convertToConservative(params, u, q_vals)

    # grabbing conserved variables
    rho = q_vals[1]
    vel_x = q_vals[2]/rho
    vel_y = q_vals[3]/rho
    Energy = q_vals[4]

    # from JC's code below, eqn should still be in scope
    # calc pressure using the variables u and the params object
    # of type ParamType{2, :the_variable_type_of_u}
    pressure = calcPressure(params, u)

    # solved eqn for e: E = rho*e + (1/2)*rho*u^2
    vel_squared = vel_x^2 + vel_y^2
    energy = Energy/rho - (1/2)*vel_squared

    # gamma stored in EulerData type
    gamma = params.gamma

#     println("pressure: ",pressure)
#     println("gamma: ",gamma)
#     println("rho: ",rho)
    # ideal gas law
    speed_sound = calcSpeedofSound(params, u)
#    speed_sound = sqrt((gamma*pressure)/rho)

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
      #TODO: iL and iR need more descriptive names
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
### EulerEquationMod.calcEdgeStabAlpha

  This function calculates the edge stabilization paramter alpha across the
  entire mesh.

  This is a mid level function.
"""
# used by EulerData Constructor - not that that matters for any reason
# mid level function
function calcEdgeStabAlpha{Tmsh,  Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator, eqn::EulerData{Tsol, Tres, Tdim})
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



#####  Functions to for filtering ######################################

function applyFilter{Tmsh, Tsol}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator, eqn::AbstractSolutionData{Tsol}, arr, opts; trans=false)
# applies filter to array arr
# trans determine whether or not to transpose the filter matrix

  if trans
    filter_mat = eqn.params.filter_mat.'
  else
    filter_mat = eqn.params.filter_mat
  end

  q_filt = zeros(Tsol, mesh.numNodesPerElement, mesh.numDofPerNode)  # holds the filtered q variables
  len = mesh.numNodesPerElement*mesh.numDofPerNode
  for i=1:mesh.numEl
    q_vals = view(eqn.q, :, :, i)  # mesh.numDof x sbp.numnodes
    # apply filter matrix to q_vals transposed, so it is applied to
    # all the rho components, then all the x momentum, etc.
    smallmatmatT!(filter_mat, q_vals, q_filt)

    # copy values back into eqn.q, remembering to transpose q_filt
    for j=1:mesh.numNodesPerElement
      for k=1:mesh.numDofPerNode
	eqn.q[k, j, i] = q_filt[j, k]
      end
    end

  end  # end loop over elements

  return nothing
end



    
function calcFilter(sbp::SBPOperator, filter_name::ASCIIString, opts)
# calc the filter specified by filter_name


  filter_func = filter_dict[filter_name]
  filt_mat = filter_func(sbp, opts)

  # testing only
#  (m,n) = size(filt_mat)
#  filt_mat = eye(m)


  V = calcModalTransformationOp(sbp)
  println("filt_mat = \n", filt_mat)
  println("V = \n", V)
  println("cond(V) = ", cond(V))
  # calculate V*filter_mat*inv(V)
  # which converts from interpolating to modal, applies filter, then 
  # converts back to interpolating
  F_t = filt_mat.'
  V_t = V.'
  F_ret = (V_t\F_t).'
  F_ret = V*F_ret

  F_ret = V*filt_mat*inv(V)
  # for testing, return identity matrix
#  (m,n) = size(filt_mat)
#  F_ret = eye(m) 


  for i=1:length(F_ret)
    if abs(F_ret[i]) < 1e-15
      F_ret[i] = 0
    end
  end

  println("F_ret = \n", F_ret)
  return F_ret

end



function calcModalTransformationOp(sbp::SBPOperator)

  vtx = [-1. -1; 1 -1; -1 1]  # reference element
  x, y = SummationByParts.SymCubatures.calcnodes(sbp.cub, vtx)
  # loop over ortho polys up to degree d
  d = sbp.degree
  n = convert(Int, (d+1)*(d+2)/2)
  V = zeros(Float64, (sbp.numnodes, n))
  Vx = zeros(V)
  Vy = zeros(V)
  ptr = 0
  for r = 0:d
    for j = 0:r
      i = r-j
      V[:,ptr+1] = SummationByParts.OrthoPoly.proriolpoly(x, y, i, j)
      Vx[:,ptr+1], Vy[:,ptr+1] = SummationByParts.OrthoPoly.diffproriolpoly(x, y, i, j)
      ptr += 1
    end
  end


  Q, R = qr(V, thin=false)
  Qx, Rx = qr(Vx, thin=false)
  Qy, Ry = qr(Vy, thin=false)

  # make the V matrix square
  # for this to work, the QR decomposition *cannot* do column pivoting
  # if it does we will have to do a little more bookkepping
  V_full = zeros(sbp.numnodes, sbp.numnodes)
  V_full[:, 1:n] = V
  V_full[:, (n+1):end] = Q[:, (n+1):end]

  # make Vx, Vy square
  # also make them full rank by replacing linearly dependend columns with vectors from Q

  Vx_full = zeros(sbp.numnodes, sbp.numnodes)
  Vy_full = zeros(sbp.numnodes, sbp.numnodes)

  Vx_full[:, 1:n] = Vx
  Vy_full[:, 1:n] = Vy

  Vx_full[:, (n+1):end] = Qx[:, (n+1):end]
  Vy_full[:, (n+1):end] = Qy[:, (n+1):end]

  # also replace linearly dependent columns here?


  return V_full
end

function calcRaisedCosineFilter(sbp::SBPOperator, opts)
# calculates the 1D raised cosine filter
# from Spectral Methods for the Euler Equations: Part I - Fourier Methods and
# shock Capturing
# Hussaini, Koproiva, Salas, Zang

  filt = zeros(sbp.numnodes, sbp.numnodes)

  max_mode = getPascalLevel(sbp.numnodes)
  for i=1:sbp.numnodes
    mode_i = getPascalLevel(i)
    theta_i = (mode_i-1)/max_mode
    filt[i, i] = 0.5*(1 + cos(theta_i))
  end

  # get rid of nullspace component
  filt[sbp.numnodes, sbp.numnodes] = 0.0


  # increase the steepness of the filter
  for i=1:sbp.numnodes
    diff =  1 - filt[i, i]
    filt[i, i] -= 2*diff
  end

  for i=1:sbp.numnodes
    if filt[i,i] < 0
      filt[i,i] = 0
    end
  end

  
  return filt
end

function getPascalLevel(node::Integer)
# get the current polynomial order of some entry node in 
# Pascals triangle
# this assumes the tree is being traversed in order, from 1 to n

  level = 1
  for i=1:(node+1) # loop over level of triangles
            # looping all the way to i is excessive
    # get the maximum node in current level of triangle
#    println("Level = ", level)
    max_n = div(level*(1 + level), 2)
#    println("max_n = ", max_n)
    # break when we are at the right level
    if  node <= max_n
#      println("breaking")
      break
    end

    level += 1  # increment level
  end


  return level
end



function calcLowPassFilter(sbp::SBPOperator, opts)

  filt = zeros(sbp.numnodes, sbp.numnodes)
  for i=1:sbp.numnodes
    filt[i, i] = 1
  end

  filt[end, end] = 0.99

  
  return filt
end


global const filter_dict = Dict{ASCIIString, Function} (
"raisedCosineFilter" => calcRaisedCosineFilter,
"lowPassFilter" => calcLowPassFilter,
)


##### Artificial Dissipation Functions ######################################

function applyDissipation{Tmsh, Tsol, T}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator, eqn::AbstractSolutionData{Tsol}, arr::AbstractArray{T, 3}, opts)
# applies the artificial dissipation to the array arr
# arr must be mesh.numDofPerNode by sbp.numNodesPerElement by mesh.numEl
# trans determine whether or not to transpose the filter matrix


  q_filt = zeros(Tsol, mesh.numNodesPerElement, mesh.numDofPerNode)  # holds the filtered q variables
  len = mesh.numNodesPerElement*mesh.numDofPerNode
  for i=1:mesh.numEl
    filt_i = view(eqn.dissipation_mat, :, :, i)
    q_vals = view(eqn.q, :, :, i)  # mesh.numDof x sbp.numnodes
    # apply filter matrix to q_vals transposed, so it is applied to
    # all the rho components, then all the x momentum, etc.
    smallmatmatT!(filt_i, q_vals, q_filt)

    # update eqn.res, remembering to transpose q_filt
    for j=1:mesh.numNodesPerElement
      for k=1:mesh.numDofPerNode
	eqn.res[k, j, i] -= q_filt[j, k]
      end
    end

  end  # end loop over elements

  return nothing
end



function calcDissipationOperator{Tmsh, Tsol}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator, eqn::AbstractEulerData{Tsol}, dissipation_name::ASCIIString, opts)
# calculates and returns the artificial dissipation operator array

  epsilon = eqn.params.dissipation_const  # get the dissipation constant
  dissipation_func = dissipation_dict[dissipation_name]

  filt = getDissipationFilterOperator(sbp, dissipation_func)  # get the dissipation filter matrix for the reference element

  # store a sbp.numnodes square matrix for each element
  dissipation_mat = zeros(Tmsh, sbp.numnodes, sbp.numnodes, mesh.numEl)

    
  # filter matrix for a non-reference element
  filt_i = Array(Tmsh, sbp.numnodes, sbp.numnodes)
  hinv = inv(diagm(sbp.w))
  h = diagm(sbp.w)
  for i=1:mesh.numEl
    # modify filter matrix to be in real (non reference) space
#=    
    for col = 1:sbp.numnodes
      for row = 1:sbp.numnodes
	filt_i[row, col] = filt[row, col]/mesh.jac[row, i]
      end
    end
=#
    # scale the mass (digonal) mass matrix by jacobian determinent
    # then multiply by h = (1/jac)^(1/p)
#    h_jac = sbp.w./(mesh.jac[:, i].^2.0)

    h_jac = sbp.w./(mesh.jac[:, i].^1.5)
    # JC modification 11/4
#    h_jac = sbp.w./(mesh.jac[:, i])

    h_jac_inv = 1./h_jac

    # this is the proper artificial dissipation
    dissipation_mat[:, :, i] = epsilon*filt.'*diagm(h_jac)*filt

    # this is the used for preconditioning the iterative solver
#    dissipation_mat[:, :, i] = epsilon*filt.'*diagm(sbp.w)*filt
  end  # end loop over elements


  return dissipation_mat

end  # end function





function getDissipationFilterOperator{T}(sbp::TriSBP{T}, filter::Function)
# calculate the filter operator (including the conversion to and from
# the modal basis) used for artificial dissipation
# the filter function defines the filter kernel

  vtx = [-1. -1.; 1. -1.; -1. 1.]
  x, y = SummationByParts.SymCubatures.calcnodes(sbp.cub, vtx)
  # loop over ortho polys up to degree d
  d = sbp.degree
  eta_c = (d-1)/(d+1)
  s = d
  alpha = 36.0

  V = zeros(T, (sbp.numnodes, div((d+1)*(d+2), 2)) )
#  V = zeros(T, (sbp.numnodes, convert(Int, (d+1)*(d+2)/2)) )
  lambda = zeros(T, (sbp.numnodes) )
  ptr = 0
  for r = 0:d
    for j = 0:r
      i = r-j
      V[:,ptr+1] = SummationByParts.OrthoPoly.proriolpoly(x, y, i, j)
      lambda[ptr+1] = 1.0 - filter(r/(d+1), eta_c, alpha, s) 
      ptr += 1
    end
  end
  lambda[ptr+1:end] = 1.0
  #println("lambda = ",lambda)

  Z = nullspace(V')
  Vt = [V Z]
  F = Vt*diagm(lambda)*inv(Vt)
  return F
end


function damp1(eta, eta_c, alpha, s)
  if (eta <= eta_c)
    return 1.0
  else
    return exp(-alpha*((eta-eta_c)/(1-eta_c))^(2*s))
  end
end



global const dissipation_dict = Dict{ASCIIString, Function} (
"damp1" => damp1
)

