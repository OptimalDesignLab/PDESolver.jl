# function related to calculating areas and volumes

@doc """
### calcVolumeContribution!

Returns the contribution from the given faces to the volume of enclosed region.
The Gauss-Divergence theorem is applied to (1/3) ∫∫∫ ∇⋅(x,y,z) dΩ = vol.

**Inputs**

* `sbpface`: an SBP face operator type
* `xsbp`: SBP-face nodes in physical space; [coord, sbp node, face]
* `nrm`: scaled face-normal at the sbpface nodes; [component, sbp node, face]
* `vol`: volume contribution from boundaries other than this set of faces

**Returns**

* `vol`: the updated value for the volume

"""->
function calcVolumeContribution!(sbpface::AbstractFace{Tsbp},
                                 xsbp::AbstractArray{Tmsh,3},
                                 nrm::AbstractArray{Tmsh,3}) where {Tsbp,Tmsh}
  @assert( sbpface.numnodes == size(xsbp,2) == size(nrm,2) )
  @assert( size(xsbp,3) == size(nrm,3) )
  @assert( size(xsbp,1) == size(nrm,1) )
  vol = zero(Tmsh)
  for f = 1:size(nrm,3)
    for i = 1:sbpface.numnodes
      for di = 1:size(nrm,1)# 3
        vol += (1/3)*sbpface.wface[i]*xsbp[di,i,f]*nrm[di,i,f]
      end
    end
  end
  return vol           
end

"""
  This function computes the volume contribution of a set of faces, where the
  set of faces is given by which boundary condition they have applied to them.
"""
function calcVolumeContribution!(mesh::AbstractMesh{Tmsh}, eqn::AbstractSolutionData, bndry_nums::Array{Int, 1}) where Tmsh
# compute the volume contribution for the faces associated with a set of
# boundary conditions

  vol = zero(Tmsh)
  for i=1:length(bndry_nums)
    start_idx = mesh.bndry_offsets[ bndry_nums[i] ]
    end_idx = mesh.bndry_offsets[ bndry_nums[i] + 1 ] - 1
    face_range = start_idx:end_idx
    bndry_faces = sview(mesh.bndryfaces, face_range)
    
    # calculate normal
    nrm = sview(mesh.nrm_bndry, :, :, face_range)
#    nrm = computeNormal(mesh, eqn, bndry_faces, mesh.dxidx_bndry)
    coords = sview(mesh.coords_bndry, :, :, face_range)

    # compute volume contribution
    vol += calcVolumeContribution!(mesh.sbpface, coords, nrm)
  end

  return vol
end

function calcVolumeContribution_rev!(mesh::AbstractMesh{Tmsh}, eqn::AbstractSolutionData, bndry_nums::Array{Int, 1}, vol_bar::Tmsh) where Tmsh

  vol = zero(Tmsh)
  for i=1:length(bndry_nums)
    start_idx = mesh.bndry_offsets[ bndry_nums[i] ]
    end_idx = mesh.bndry_offsets[ bndry_nums[i] + 1 ] - 1
    face_range = start_idx:end_idx
    bndry_faces = sview(mesh.bndryfaces, face_range)
    coords = sview(mesh.coords_bndry, :, :, face_range)
    coords_bar = zeros(coords)  #TODO: make this a field of mesh?

    nrm = sview(mesh.nrm_bndry, :, :, face_range)
    nrm_bar = sview(mesh.nrm_bndry_bar, :, :, face_range)
#    nrm = computeNormal(mesh, eqn, bndry_faces)
#    nrm_bar = zeros(Tmsh, mesh.dim, mesh.numNodesPerFace, length(bndry_faces))

    # compute volume contribution
    calcVolumeContribution_rev!(mesh.sbpface, coords, coords_bar, nrm, nrm_bar, vol_bar)

    # calculate normal
#    computeNormal_rev(mesh, eqn,  bndry_faces, nrm_bar)

  end

  return nothing
end

"""
  Computes the normal vector at each surface node of the given list of faces.
  Straight-sided meshes only!

  Inputs:
    mesh: an AbstractMesh
    sbp: an SBP operator
    bndryfaces: list of boundary faces

  Outputs: nrm: mesh.dim x mesh.numNodesPerFace x length(bndryfaces) array 
                containing the normal vectors
"""
function computeNormal(mesh::AbstractMesh{Tmsh}, eqn::AbstractSolutionData,
                       bndryfaces::AbstractArray{Boundary, 1},
                       dxidx::AbstractArray{Tmsh, 4}) where Tmsh

  if typeof(mesh) <: AbstractDGMesh
    @assert mesh.coord_order == 1
  end

  nfaces = length(bndryfaces)
  nrm = zeros(Tmsh, mesh.dim, mesh.numNodesPerFace, nfaces)

  for i=1:nfaces
    bndry_i = bndryfaces[i]
    nrm_xi = sview(mesh.sbpface.normal, :, bndry_i.face)
    for j=1:mesh.numNodesPerFace
      dxidx_j = sview(dxidx, :, :, j, i)
      nrm_xy_j = sview(nrm, :, j, i)
      calcBCNormal(eqn.params, dxidx_j, nrm_xi, nrm_xy_j)
    end
  end

  return nrm
end

"""
  This function uses reverse mode to back propigate perturbations in nrm to 
  dxidx_bndry
"""
function computeNormal_rev(mesh::AbstractMesh, eqn::AbstractSolutionData,
                           bndryfaces::AbstractArray{Boundary, 1},
                           nrm_bar::Abstract3DArray{Tmsh},
                           dxidx::AbstractArray{Tmsh, 4},
                           dxidx_bar::AbstractArray{Tmsh, 4}) where Tmsh
  if typeof(mesh) <: AbstractDGMesh
    @assert mesh.coord_order == 1
  end

  nfaces = length(bndryfaces)
  for i=1:nfaces
    bndry_i = bndryfaces[i]
    nrm_xi = sview(mesh.sbpface.normal, :, bndry_i.face)
    for j=1:mesh.numNodesPerFace
      dxidx_j = sview(dxidx, :, :, j, i)
      dxidxbar_j = sview(dxidx_bar, :, :, j, i)
      nrmbar_xy_j = sview(nrm_bar, :, j, i)
      calcBCNormal_revm(eqn.params, dxidx_j, nrm_xi, nrmbar_xy_j, dxidxbar_j)
    end
  end

  return nothing
end

@doc """
### calcVolumeContribution_rev!

This is the reverse differentiated version of calcVolumeContribution!.  See docs
of calcVolumeContribution! for further details of the primal method.  This
function is differentiated with respect to the primal version's `xsbp` and `nrm`
variables.

**Inputs**

* `sbpface`: an SBP face operator type
* `xsbp`: SBP-face nodes in physical space; [coord, sbp node, face]
* `nrm`: scaled face-normal at the sbpface nodes; [component, sbp node, face]
* `vol_bar`: left multiplies d(vol)/d(nrm) and d(vol)/d(xsbp); it is a scalar

**InOuts**

* `xsbp_bar`: result of vector Jacobian product; [component, sbp node, face]
* `nrm_bar`: result of vector Jacobian product; [component, sbp node, face]

"""->
function calcVolumeContribution_rev!(sbpface::AbstractFace{Tsbp},
                                     xsbp::AbstractArray{Tmsh,3},
                                     xsbp_bar::AbstractArray{Tmsh,3},
                                     nrm::AbstractArray{Tmsh,3},
                                     nrm_bar::AbstractArray{Tmsh,3},
                                     vol_bar::Tmsh) where {Tsbp,Tmsh}
  @assert( sbpface.numnodes == size(xsbp,2) == size(xsbp_bar,2) == size(nrm,2)
           == size(nrm_bar,2) )
  @assert( size(xsbp,3) == size(xsbp_bar,3) == size(nrm,3) == size(nrm_bar,3) )
  @assert( size(xsbp,1) == size(xsbp,1) == size(nrm,1) == size(nrm_bar,1) )
  for f = 1:size(nrm,3)
    for i = 1:sbpface.numnodes
      for di = 1:size(nrm,1) # 3
        # vol += (1/3)*sbpface.wface[i]*xsbp[di,i,f]*nrm[di,i,f]
        fac = vol_bar*(1/3)*sbpface.wface[i]
        xsbp_bar[di,i,f] += fac*nrm[di,i,f]
        nrm_bar[di,i,f] += fac*xsbp[di,i,f]
      end
    end
  end
end

@doc """
### calcProjectedAreaContribution!

  Returns the contribution from the given faces to the projected area onto the
  plane of coordiante `di`

  Because this function uses the absolute value of the normal vector, it is not
  differentiable at zero.

**Inputs**

* `sbpface`: an SBP face operator type
* `di`: the index that indicates the projection plane (e.g. `di`=1 --> yz plane)
* `nrm`: scaled face-normal at the sbpface nodes; [component, sbp node, face]
* `projarea`: projected area of faces from other boundaries.

**Returns**

* `projarea`: the updated value for projected area

"""->
function calcProjectedAreaContribution!(sbpface::AbstractFace{Tsbp},
                                        nrm::AbstractArray{Tmsh,3},
                                        di::Int) where {Tsbp,Tmsh}
  @assert( sbpface.numnodes == size(nrm,2) )
  projarea = zero(Tmsh)
  for f = 1:size(nrm,3)
    for i = 1:sbpface.numnodes
      projarea += 0.5*sbpface.wface[i]*absvalue(nrm[di,i,f])
    end
  end
  return projarea
end

"""
  This function computes the volume contribution of a set of faces, where the
  set of faces is given by which boundary condition they have applied to them.
"""
function calcProjectedAreaContribution!(mesh::AbstractMesh{Tmsh}, eqn::AbstractSolutionData, bndry_nums::Array{Int, 1}, di::Int) where Tmsh

  proj_area = zero(Tmsh)
  for i=1:length(bndry_nums)
    start_idx = mesh.bndry_offsets[ bndry_nums[i] ]
    end_idx = mesh.bndry_offsets[ bndry_nums[i] + 1 ] - 1
    face_range = start_idx:end_idx
    bndry_faces = sview(mesh.bndryfaces, face_range)
    
    # calculate normal
    nrm = sview(mesh.nrm_bndry, :, :, face_range)
#    nrm = computeNormal(mesh, eqn, bndry_faces)
#    coords = sview(mesh.coords_bndry, :, :, face_range)

    proj_area += calcProjectedAreaContribution!(mesh.sbpface, nrm, di)
  end

  return proj_area
end



@doc """
### calcProjectedAreaContribution_rev!

This is the reverse differentiated version of calcProjectedAreaContribution!.
See docs of calcProjectedAreaContribution! for further details of the primal
method.  This function is differentiated with respect to the primal version's
`nrm` variable.

**Inputs**

* `sbpface`: an SBP face operator type
* `di`: the index that indicates the projection plane (e.g. `di`=1 --> yz plane)
* `nrm`: scaled face-normal at the sbpface nodes; [component, sbp node, face]
* `projarea_bar`: left multiplies d(projarea)/d(nrm); it is a scalar

**InOuts**

* `nrm_bar`: result of vector Jacobian product; [component, sbp node, face]

"""->
function calcProjectedAreaContribution_rev!(sbpface::AbstractFace{Tsbp}, nrm::AbstractArray{Tmsh,3},
    nrm_bar::AbstractArray{Tmsh,3}, di::Int, projarea_bar::Tmsh) where {Tsbp,Tmsh}
  @assert( sbpface.numnodes == size(nrm,2) == size(nrm_bar,2) )
  for f = 1:size(nrm,3)
    for i = 1:sbpface.numnodes
      # projarea += 0.5*sbpface.wface[i]*abs(nrm[di,i,f])
      #  v1 = nrm[di,i,f]
      #  v2 = abs(v1)
      #  v3 = 0.5*sbpface.wface[i]
      #  v4 = v2*v3

      nrm_bar[di,i,f] += (projarea_bar*0.5*sbpface.wface[i]*
                          absvalue_deriv(nrm[di,i,f]))
    end
  end
end

"""
  Back propigate projarea_bar to mesh.dxidx_bndry_bar
"""
function calcProjectedAreaContribution_rev!(mesh::AbstractMesh, eqn::AbstractSolutionData, bndry_nums::Array{Int, 1}, di::Int,  projarea_bar::Tmsh) where Tmsh

  
  for i=1:length(bndry_nums)
    start_idx = mesh.bndry_offsets[ bndry_nums[i] ]
    end_idx = mesh.bndry_offsets[ bndry_nums[i] + 1 ] - 1
    face_range = start_idx:end_idx
    bndry_faces = sview(mesh.bndryfaces, face_range)

    nrm = sview(mesh.nrm_bndry, :, :, face_range)
    nrm_bar = sview(mesh.nrm_bndry_bar, :, :, face_range)
#    nrm = computeNormal(mesh, eqn, bndry_faces)
#    nrm_bar = zero(nrm)
    calcProjectedAreaContribution_rev!(mesh.sbpface, nrm, nrm_bar, di, projarea_bar)

    # back propigate to dxidx_bndry_bar
#    computeNormal_rev(mesh, eqn, bndry_faces, nrm_bar)
  end

  return nothing
end
  

function crossProd(x::AbstractArray{T,1}, y::AbstractArray{T,1},
                   cp::AbstractArray{T,1}) where T
  cp[1] = x[2]*y[3] - y[2]*x[3]
  cp[2] = x[3]*y[1] - y[3]*x[1]
  cp[3] = x[1]*y[2] - y[1]*x[2]
end

function crossProd_rev(x::AbstractArray{T,1}, x_bar::AbstractArray{T,1},
                       y::AbstractArray{T,1}, y_bar::AbstractArray{T,1},
                       cp_bar::AbstractArray{T,1}) where T
  #cp[1] = x[2]*y[3] - y[2]*x[3]
  x_bar[2] +=  y[3]*cp_bar[1]
  y_bar[3] +=  x[2]*cp_bar[1]
  y_bar[2] += -x[3]*cp_bar[1]
  x_bar[3] += -y[2]*cp_bar[1]
  #cp[2] = x[3]*y[1] - y[3]*x[1]
  x_bar[3] +=  y[1]*cp_bar[2]
  y_bar[1] +=  x[3]*cp_bar[2]
  y_bar[3] += -x[1]*cp_bar[2]
  x_bar[1] += -y[3]*cp_bar[2]
  #cp[3] = x[1]*y[2] - y[1]*x[2]
  x_bar[1] +=  y[2]*cp_bar[3]
  y_bar[2] +=  x[1]*cp_bar[3]
  y_bar[1] += -x[2]*cp_bar[3]
  x_bar[2] += -y[1]*cp_bar[3]
end


