# file stabilization.jl
# this is the primary file for adding stabilization to the weak form.
# Stabilization functions may be put in separate files if they are included 
# here

include("GLS.jl")
include("filtering.jl")
include("artificial_dissipation.jl")

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
function edgestabilize!{T, Tres}(sbp::AbstractSBP{T}, ifaces::Array{Interface},
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

function edgestabilize!{Tmsh,  Tsol, Tres}(sbp::AbstractSBP, 
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
function edgestabilize!{Tmsh,  Tsol, Tres}(mesh, sbp::AbstractSBP, eqn, 
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
  and returns it.

   Inputs:
    * q : vector of conservative variables
    * dxidx : jacobian of xi wrt x coordinates at the node
    * nrm : normal vector in xi space
    * params : ParamType{2}

    This is a low level function
"""->
# low level function
function stabscale{Tmsh, Tsol}(q::AbstractArray{Tsol,1}, 
                   dxidx::AbstractArray{Tmsh,2}, nrm::AbstractArray{Tmsh,1}, 
                   params::ParamType{2} )
# calculate stabscale for a single node

#     println("==== entering stabscale ====")
   
    gamma = params.gamma
    edge_stab_gamma = params.edgestab_gamma

    # u holds whatever the variable type that we are solving the equation in
    # q_vals holds the conservative variables
    # convert to conservative variables if not already using them
    q_vals = params.q_vals
    convertToConservative(params, q, q_vals)

    # calc pressure using the variables u and the params object
    # of type ParamType{2, :the_variable_type_of_u}
    pressure = calcPressure(params, q)
    speed_sound = calcSpeedofSound(params, q)

    # extract conserved variables
    rho = q_vals[1]
    vel_x = q_vals[2]/rho
    vel_y = q_vals[3]/rho
    Energy = q_vals[4]

    vel_squared = vel_x^2 + vel_y^2
    energy = Energy/rho - (1/2)*vel_squared


    # edge lengths component wise
    h_x = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
    h_y = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]
    h = sqrt(h_x^2 + h_y^2)  # edge length

    U = vel_x*(h_x/h) + vel_y*(h_y/h)

    return (abs(U) + speed_sound)*edge_stab_gamma*h^(2)

  end


@doc """
### EulerEquationMod.stabscale

  This function calculate the stabilization scaling parameter across the
  entire mesh by calling the low level method.  This populates eqn.stabscale.

  This is a mid level function
"""->
# mid level function
function stabscale{Tmsh,  Tsol}(mesh::AbstractMesh{Tmsh}, sbp::AbstractSBP, 
                                eqn::EulerData{Tsol})
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
function calcEdgeStabAlpha{Tmsh,  Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh}, 
                           sbp::AbstractSBP, eqn::EulerData{Tsol, Tres, Tdim})
# calculate alpha, needed by edge stabilization

  numEl = mesh.numEl
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




