# functions for calculating boundary integrals

include("bc_solvers.jl")  # Roe solvers and related things

"""
  Type that identified a particular node on a particular face of
  a particular boundary.  It also contains the index of the face within
  the array of Boundaries with the same boundary condition

  **Fields**

   * element: the element the face is part of
   * face: the local face number
   * faceidx: the index within the array of Boundaries with this BC
   * node: the node on the face

"""
struct BoundaryNode
  element::UInt32
  face::UInt8
  faceidx::Int  # index within array of faces that have this boundary condition
  node::Int
end

function BoundaryNode(bndry::Boundary, faceidx, node)
  return BoundaryNode(bndry.element, bndry.face, faceidx, node)
end

"""
  Null boundary node (all fields zero).  This is useful as a default value
  for a functiona argument (if the argument is unused).
"""
global const NullBoundaryNode = BoundaryNode(Boundary(0, 0), 0, 0)

@doc """
### EulerEquationMod.getBCFluxes

  This function calls other functions to calculate the boundary fluxes, passing
  them pieces of the array needed.  This populates eqn.bndryflux.  It also
  calls writeBoundary() to do any requested output.  If the options dictionary
  specifies not to precompute the boundary flux, this function will do the
  integration as well and update `eqn.res`.

  This is a mid level function.

  
"""->
# this is a mid level function
function getBCFluxes(mesh::AbstractMesh, sbp::AbstractSBP, eqn::EulerData, opts)
  #get all the fluxes for all the boundary conditions and save them in eqn.bndryflux

  #println("mesh.bndry_funcs = ", mesh.bndry_funcs)
  for i=1:mesh.numBC
  #  println("computing flux for boundary condition ", i)
    functor_i = mesh.bndry_funcs[i]
    start_index = mesh.bndry_offsets[i]
    end_index = mesh.bndry_offsets[i+1]
    idx_range = start_index:(end_index - 1)
    bndry_facenums_i = sview(mesh.bndryfaces, start_index:(end_index - 1))

    if opts["use_staggered_grid"]

      calcBoundaryFlux_nopre(mesh, mesh.mesh2, sbp, mesh.sbp2, eqn, functor_i,
                             idx_range, bndry_facenums_i)
      
    elseif opts["precompute_boundary_flux"]
      bndryflux_i = sview(eqn.bndryflux, :, :, start_index:(end_index - 1))

      # call the function that calculates the flux for this boundary condition
      # passing the functor into another function avoid type instability
      calcBoundaryFlux(mesh, sbp, eqn, functor_i, idx_range, bndry_facenums_i, bndryflux_i)
    else
      calcBoundaryFlux_nopre(mesh, sbp, eqn, functor_i, idx_range, bndry_facenums_i)
    end
  end

  writeBoundary(mesh, sbp, eqn, opts)

  return nothing
end

@doc """
### EulerEquationMod.writeBoundary

  This function writes information about the boundary faces and fluxes to files.
  It is controlled by the input argument writeboundary, of type Bool.

  It generates the files:
    * boundaryfaces.dat : writes mesh.bndryfaces, an array with eltype Boundary
                          to a file, one element per line
    * boundaryflux.dat  : writes the element, local node number and boundary
                          flux to a line in a human readable format
    * boundaryflux2.dat : writes the real part ofmesh.bndryflux to space
                          delimited file

   This is a high level function.
"""->
function writeBoundary(mesh, sbp, eqn, opts)

  if !eqn.params.writeboundary
    return nothing
  end

    face_name = "boundaryfaces.dat"
    flux_name = "boundaryflux.dat"
    flux_dlm = "boundaryflux2.dat"

    rmfile(face_name)
    rmfile(flux_name)
    rmfile(flux_dlm)


  # write boundaryfaces.dat
  f = open(face_name, "a+")
  for i=1:length(mesh.bndryfaces)
    println(f, mesh.bndryfaces[i])
  end
  close(f)

  # write boundaryflux.dat
  f = open(flux_name, "a+")
  for i=1:mesh.numBoundaryFaces
    el = mesh.bndryfaces[i].element
    face = mesh.bndryfaces[i].face
    for j=1:mesh.numNodesPerFace  # TODO: should be mesh.numNodesPerFace?
      jb = mesh.facenodes[j, face]
      println(f, "el ", el, ", node_index ", jb, ", flux = ",
               real(eqn.bndryflux[:, j, i]))
    end
  end
  close(f)

  # write boundaryflux2.dat
  writedlm(flux_dlm, real(eqn.bndryflux))

  return nothing
end

@doc """
### EulerEquationMod.interpolateBoundary

  Interpolates the solution variables to the exterior boundary of the mesh
  and calculates any additional quantities at the boundary of the mesh.
  DG only

  Inputs:
    mesh: an AbstractDGMesh
    sbp
    eqn
    opts
    q : the 3D array of solution variables for all elements, numdofpernode x
        numNodesPerElement x numEl

  Inputs/Outputs:
    q_bndry: the array to be populated with the solution interpolated to
             the boundary, numdofpernode x numNodesPerFace x num boundary faces

    eqn.aux_vars_bndry is also populated

    Aliasing restrictions: none
"""->
function interpolateBoundary(mesh::AbstractDGMesh, sbp, eqn, opts, q::Abstract3DArray, q_bndry::Abstract3DArray)

  # interpolate solutions
  boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, eqn.q, eqn.q_bndry)

  # calculate aux_vars
  for i=1:mesh.numBoundaryFaces
    for j=1:mesh.numNodesPerFace
      q_vals = ro_sview(eqn.q_bndry, :, j, i)
      eqn.aux_vars_bndry[1, j, i] = calcPressure(eqn.params, q_vals)
    end
  end

end


@doc """
### EulerEquationMod.calcBoundaryFlux

  This function calculates the boundary flux for the portion of the boundary
  with a particular boundary condition.  The eqn.q are converted to
  conservative variables if needed

  Inputs:
  mesh : AbstractMesh
  sbp : AbstractSBP
  eqn : EulerEquation
  functor : a callable object that calculates the boundary flux at a node
  idx_range: the Range describing which Boundaries have the current BC
  bndry_facenums:  An array with elements of type Boundary that tell which
                   element faces have the boundary condition
  Outputs:
  bndryflux : the array to store the boundary flux, corresponds to
              bndry_facenums

  The functor must have the signature
  functor( q, aux_vars, x, nrm_xy, bndryflux_i, eqn.params)
  where q are the *conservative* variables.
  where all arguments (except params) are vectors of values at a node.

  params is the ParamType associated with the the EulerEquation object
  nrm = mesh.sbpface.normal[:, current_node]

  This is a mid level function.
"""->
# mid level function
function calcBoundaryFlux( mesh::AbstractCGMesh{Tmsh},
       sbp::AbstractSBP, eqn::EulerData{Tsol},
       functor::BCType, idx_range::UnitRange,
       bndry_facenums::AbstractArray{Boundary,1},
       bndryflux::AbstractArray{Tres, 3}) where {Tmsh,  Tsol, Tres}
# calculate the boundary flux for the boundary condition evaluated by the functor

#  println("enterted calcBoundaryFlux")


  nfaces = length(bndry_facenums)
  q2 = zeros(Tsol, mesh.numDofPerNode)
  nrm_xy = zeros(Tmsh, mesh.dim)
  for i=1:nfaces  # loop over faces with this BC
    bndry_i = bndry_facenums[i]
#    println("boundary ", i, "element = ", bndry_i.element, ", face = ", bndry_i.face)
#    println("interface ", i)
    for j = 1:mesh.numNodesPerFace
      k = mesh.facenodes[j, bndry_i.face]

      # get components
      q = ro_sview(eqn.q, :, k, bndry_i.element)
      # convert to conservative variables if needed
      convertToConservative(eqn.params, q, q2)
      aux_vars = ro_sview(eqn.aux_vars, :, k, bndry_i.element)
      x = ro_sview(mesh.coords, :, k, bndry_i.element)
      dxidx = ro_sview(mesh.dxidx, :, :, k, bndry_i.element)
      nrm = ro_sview(mesh.sbpface.normal, :, bndry_i.face)
      calcBCNormal(eqn.params, dxidx, nrm, nrm_xy)
      bndryflux_i = sview(bndryflux, :, j, i)

      functor(eqn.params, q2, aux_vars, x, nrm_xy, bndryflux_i)

    end

  end


  return nothing
end


# DG version
function calcBoundaryFlux( mesh::AbstractDGMesh{Tmsh},
       sbp::AbstractSBP, eqn::EulerData{Tsol},
       functor::BCType, idx_range::UnitRange,
       bndry_facenums::AbstractArray{Boundary,1},
       bndryflux::AbstractArray{Tres, 3}) where {Tmsh,  Tsol, Tres}
  # calculate the boundary flux for the boundary condition evaluated by the
  # functor

  nfaces = length(bndry_facenums)
  q2 = zeros(Tsol, mesh.numDofPerNode)
  params = eqn.params
  for i=1:nfaces  # loop over faces with this BC
    bndry_i = bndry_facenums[i]
    global_facenum = idx_range[i]

    for j = 1:mesh.numNodesPerFace
      # get components
      q = ro_sview(eqn.q_bndry, :, j, global_facenum)
      # convert to conservative variables if needed
      convertToConservative(eqn.params, q, q2)
      aux_vars = ro_sview(eqn.aux_vars_bndry, :, j, global_facenum)
      coords = ro_sview(mesh.coords_bndry, :, j, global_facenum)
      nrm_xy = ro_sview(mesh.nrm_bndry, :, j, global_facenum)
      bndryflux_i = sview(bndryflux, :, j, i)

      bndry_node = BoundaryNode(bndry_i, i, j)

      functor(params, q2, aux_vars, coords, nrm_xy, bndryflux_i, bndry_node)

    end
  end

  return nothing
end

"""
  Like calcBoundaryFlux, but performs the integration and updates res rather
  than storing the flux.
"""
function calcBoundaryFlux_nopre( mesh::AbstractDGMesh{Tmsh},
                          sbp::AbstractSBP, eqn::EulerData{Tsol, Tres},
                          functor::BCType, idx_range::UnitRange,
                          bndry_facenums::AbstractArray{Boundary,1}) where {Tmsh,  Tsol, Tres}
  # calculate the boundary flux for the boundary condition evaluated by the
  # functor

  nfaces = length(bndry_facenums)
  q2 = zeros(Tsol, mesh.numDofPerNode)
  params = eqn.params
  flux_face = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)
  for i=1:nfaces  # loop over faces with this BC
    bndry_i = bndry_facenums[i]
    global_facenum = idx_range[i]
    for j = 1:mesh.numNodesPerFace

      # get components
      #TODO: this doesn't work if precompute_q_bndr == false ?
      q = ro_sview(eqn.q_bndry, :, j, global_facenum)
      # convert to conservative variables if needed
      convertToConservative(eqn.params, q, q2)
      aux_vars = ro_sview(eqn.aux_vars_bndry, :, j, global_facenum)
      coords = ro_sview(mesh.coords_bndry, :, j, global_facenum)
      nrm_xy = ro_sview(mesh.nrm_bndry, :, j, global_facenum)
      bndryflux_i = sview(flux_face, :, j)

     
      bndry_node = BoundaryNode(bndry_i, i, j)

      functor(params, q2, aux_vars, coords, nrm_xy, bndryflux_i, bndry_node)
    end

    res_i = sview(eqn.res, :, :, bndry_i.element)
    boundaryFaceIntegrate!(mesh.sbpface, bndry_i.face, flux_face, res_i,
                           SummationByParts.Subtract())
  end

  return nothing
end


"""
  Staggered grid version
"""
function calcBoundaryFlux_nopre(mesh_s::AbstractDGMesh{Tmsh},
                          mesh_f::AbstractDGMesh{Tmsh},
                          sbp_s::AbstractSBP, sbp_f::AbstractSBP,
                          eqn::EulerData{Tsol, Tres},
                          functor::BCType, idx_range::UnitRange,
                          bndry_facenums::AbstractArray{Boundary,1}) where {Tmsh,  Tsol, Tres}
  # calculate the boundary flux for the boundary condition evaluated by the
  # functor

  nfaces = length(bndry_facenums)
  q_face = zeros(Tsol, mesh_f.numDofPerNode, mesh_f.numNodesPerFace)
  params = eqn.params
  flux_face = zeros(Tres, mesh_f.numDofPerNode, mesh_f.numNodesPerFace)
  res_f = zeros(Tres, mesh_f.numDofPerNode, mesh_f.numNodesPerElement)
  res_s = zeros(Tres, mesh_s.numDofPerNode, mesh_s.numNodesPerElement)
  aux_vars = zeros(Tsol, 1)
  for i=1:nfaces  # loop over faces with this BC
    bndry_i = bndry_facenums[i]
    global_facenum = idx_range[i]
    q_vol = ro_sview(eqn.q_flux, :, :, bndry_i.element)
    boundaryFaceInterpolate!(mesh_f.sbpface, bndry_i.face, q_vol, q_face)

    # interpolate to face
    for j = 1:mesh_f.numNodesPerFace

      # get components
      q_j = sview(q_face, :, j)

      # convert to conservative variables if needed
      aux_vars[1] = calcPressure(eqn.params, q_j)
      coords = ro_sview(mesh_f.coords_bndry, :, j, global_facenum)
      nrm_xy = ro_sview(mesh_f.nrm_bndry, :, j, global_facenum)
      bndryflux_i = sview(flux_face, :, j)

      functor(params, q_j, aux_vars, coords, nrm_xy, bndryflux_i)
    end  # end loop j

    # integrate
    fill!(res_f, 0.0)
    fill!(res_s, 0.0)
    boundaryFaceIntegrate!(mesh_f.sbpface, bndry_i.face, flux_face, res_f,
                           SummationByParts.Subtract())

    # interpolate back
    smallmatmat!(res_f, mesh_s.I_S2F, res_s)


    # accumulate into res
    @simd for j=1:mesh_s.numNodesPerElement
      @simd for k=1:mesh_s.numDofPerNode
        eqn.res[k, j, bndry_i.element] += res_s[k, j]
      end
    end

  end  # end loop i

  return nothing
end

#-----------------------------------------------------------------------------
# Boundary condition definitions
# TODO: move these to separate file

@doc """
### EulerEquationMod.isentropicVortexBC <: BCTypes

  This type and the associated call method define a functor to calculate
  the flux using the Roe Solver using the exact InsentropicVortex solution
  as boundary state.  See calcBoundaryFlux for the arguments all functors
  must support.

  This is a low level functor.

**Arguments**

*  `obj` : Object of type BCType used for multiple dispatch. Every new boundary
           condition needs to have its own type and entered in BCDict
*  `q`   : Solution variable
*  `aux_vars` : Auxiliary variables
*  `x`        : physical coordinates of the SBP node
*  `nrm_xy`      : sclaed face normal in physical space
*  `bndryflux` : Computed flux value at the boundary

"""->
mutable struct isentropicVortexBC <: BCType
end

function (obj::isentropicVortexBC)(params::ParamType,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  gamma = params.gamma
  gami = params.gamma_1

  # getting qg
  qg = params.qg
  calcIsentropicVortex(params, coords, qg) # Get the boundary value

  v_vals = params.q_vals
  convertFromNaturalToWorkingVars(params, q, v_vals)

  # Getting SAT terms
  specific_vol = 1.0/v_vals[1]
  u = v_vals[2]*specific_vol
  v = v_vals[3]*specific_vol
  phi = 0.5*(u*u + v*v)
  H = gamma*v_vals[4]*specific_vol - gami*phi # Total Enthalpy

#  dq = zeros(Tsol, 4)
  dq = v_vals - qg  #!!! this allocates a new vector dq every time
#  nrm2 = params.nrm2
#  calcBCNormal(params, dxidx, nrm, nrm2)
  sat = params.sat_vals
  roe_vars = params.roe_vars
  roe_vars[1] = u
  roe_vars[2] = v
  roe_vars[3] = H
  calcSAT(params, roe_vars, dq, nrm_xy, sat)

  euler_flux = zeros(Tsol, 4) # params.flux_vals1
  calcEulerFlux(params, v_vals, aux_vars, nrm_xy, euler_flux)

  sat_fac = 1.0 # Multiplier for SAT term
  for i=1:4
    bndryflux[i] = euler_flux[i] + sat_fac*sat[i]
  end

  return nothing
end

#=
# low level function
function (obj::isentropicVortexBC){Tmsh, Tsol, Tres}(params::ParamType,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, coords::AbstractArray{Tmsh,1},
               nrm::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode)

  qg = params.qg
  calcIsentropicVortex(params, coords, qg)
  RoeSolver(params, q, qg, aux_vars, nrm, bndryflux)

  return nothing

end # ends the function isentropicVortexBC
=#

@doc """
###EulerEquationMod.isentropicVortexBC_revm

Reverse mode for isentropicVortexBC.

**Inputs**

* `obj` : Type of the Boundary condition being evaluated. Its a subtype of
          BCType_revm
* `q`   : Solution variable
* `aux_vars` : Auxiliary variables
* `x`     : Node coordinates
* `nrm`   : scaled normal vector in x-y space
* `bndryflux_bar` : Input flux value seed that is used to compute the reverse
                    mode derivative.
* `params`        : equation object parameters

**Output**

* `nrm_bar` : Derivative of flux w.r.t the nrm

"""->

mutable struct isentropicVortexBC_revm <: BCType_revm
end

function (obj::isentropicVortexBC_revm)(params::ParamType2,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},  coords::AbstractArray{Tmsh,1},
              nrm::AbstractArray{Tmsh,1},
              nrm_bar::AbstractVector{Tmsh},
              bndryflux_bar::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  # Forward sweep
  gamma = params.gamma
  gami = params.gamma_1

  # getting qg
  qg = params.qg
  calcIsentropicVortex(params, coords, qg) # Get the boundary value
  v_vals = params.q_vals
  convertFromNaturalToWorkingVars(params, q, v_vals)

  # Getting SAT terms
  specific_vol = 1.0/v_vals[1]
  u = v_vals[2]*specific_vol
  v = v_vals[3]*specific_vol
  phi = 0.5*(u*u + v*v)
  H = gamma*v_vals[4]*specific_vol - gami*phi # Total Enthalpy

  dq = zeros(Tsol, 4)
  dq = v_vals - qg
#  nrm2 = params.nrm2
#  calcBCNormal(params, dxidx, nrm, nrm2)

  sat_fac = 1.0 # Multiplier for SAT term
  # for i=1:4
  #   bndryflux[i] = euler_flux[i] + sat_fac*sat[i]
  # end

  # Reverse sweep
  sat_bar = zeros(Tsol, 4)
  euler_flux_bar = zeros(Tsol, 4)
  for i = 1:4
    sat_bar[i] = sat_fac*bndryflux_bar[i]
    euler_flux_bar[i] = bndryflux_bar[i]
  end

#  nrm2_bar = zeros(Tmsh, 2)
  calcEulerFlux_revm(params, v_vals, aux_vars, nrm, euler_flux_bar, nrm_bar)
  calcSAT_revm(params, nrm, dq, [u,v], H, sat_bar, nrm_bar)
#  calcBCNormal_revm(params, dxidx, nrm, nrm2_bar, dxidx_bar)

  return nothing
end

@doc """
### EulerEquationMod.isentropicVortexBC_physical <: BCTypes

  This type and the associated call method define a functor to calculate
  the actual Euler flux  using the exact InsentropicVortex solution
  as boundary state.  See calcBoundaryFlux for the arguments all functors
  must support.

  This is a low level functor.
"""->
mutable struct isentropicVortexBC_physical <: BCType
end

function (obj::isentropicVortexBC_physical)(
              params::ParamType2,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  calcEulerFlux(params, q, aux_vars, nrm_xy, bndryflux)

  return nothing

end # end function isentropicVortexBC_physical


@doc """
### EulerEquationMod.noPenetrationBC <: BCTypes

  This functor uses the Roe solver to calculate the flux for a boundary
  state where the fluid velocity is projected into the wall.

  Works in 2D, untested in 3D.

  This is a low level functor

**Arguments**

*  `obj` : Object of type BCType used for multiple dispatch. Every new boundary
           condition needs to have its own type and entered in BCDict
*  `q`   : Solution variable
*  `aux_vars` : Auxiliary variables
*  `x`        : physical coordinates of the SBP node
*  `nrm_xy`      : scaled normal vector in x-y space
*  `bndryflux` : Computed flux value at the boundary

"""
mutable struct noPenetrationBC <: BCType
end

# low level function
function (obj::noPenetrationBC)(params::ParamType2,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},  coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}
# a clever optimizing compiler will clean this up
# there might be a way to do this with fewer flops using the tangent vector

  # calculate normal vector in xy space
  nx = nrm_xy[1]
  ny = nrm_xy[2]
#  nx2 = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
#  ny2 = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]
  fac = 1.0/(sqrt(nx*nx + ny*ny))
  # normalize normal vector
  nx *= fac
  ny *= fac

  # Get the normal momentum
  Unrm = nx*q[2] + ny*q[3]

  qg = params.qg
  for i=1:length(q)
    qg[i] = q[i]
  end

  # Subtract the normal component of the momentum from \xi & \eta components
  # of the momentum
  qg[2] -= nx*Unrm
  qg[3] -= ny*Unrm

  v_vals = params.v_vals
  convertFromNaturalToWorkingVars(params, qg, v_vals)
  # this is a problem: q is in conservative variables even if
  # params says we are using entropy variables


#  RoeSolver(params, q, qg, aux_vars, nrm_xy, bndryflux)
#  calcLFFlux(params, q, v_vals, aux_vars, nrm_xy, bndryflux)
  calcEulerFlux(params, v_vals, aux_vars, nrm_xy, bndryflux)

  return nothing
end


function (obj::noPenetrationBC)(params::ParamType3,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},  coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}
# a clever optimizing compiler will clean this up
# there might be a way to do this with fewer flops using the tangent vector


  # calculate normal vector in xy space
  nx = nrm_xy[1]
  ny = nrm_xy[2]
  nz = nrm_xy[3]
#  nx = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2] + dxidx[3,1]*nrm[3]
#  ny = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2] + dxidx[3,2]*nrm[3]
#  nz = dxidx[1,3]*nrm[1] + dxidx[2,3]*nrm[2] + dxidx[3,3]*nrm[3]
  fac = 1.0/(sqrt(nx*nx + ny*ny + nz*nz))
  # normalize normal vector
  nx = nx*fac
  ny = ny*fac
  nz = nz*fac

  # this is momentum, not velocity?
  Unrm = nx*q[2] + ny*q[3] + nz*q[4]

  qg = params.qg
  for i=1:length(q)
    qg[i] = q[i]
  end

  # calculate normal velocity
  qg[2] -= nx*Unrm
  qg[3] -= ny*Unrm
  qg[4] -= nz*Unrm

  # call Roe solver
  #RoeSolver(params, q, qg, aux_vars, nrm_xy, bndryflux)
#  nx2 = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2] + dxidx[3,1]*nrm[3]
#  ny2 = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2] + dxidx[3,2]*nrm[3]
#  nz2 = dxidx[1,3]*nrm[1] + dxidx[2,3]*nrm[2] + dxidx[3,3]*nrm[3]

  v_vals = params.v_vals
  convertFromNaturalToWorkingVars(params, qg, v_vals)
  # this is a problem: q is in conservative variables even if
  # params says we are using entropy variables
  calcEulerFlux(params, v_vals, aux_vars, nrm_xy, bndryflux)
#  calcLFFlux(params, q, v_vals, aux_vars, nrm_xy, bndryflux)

  return nothing
end

mutable struct noPenetrationESBC <: BCType
end

# low level function
function (obj::noPenetrationESBC)(params::ParamType2,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},  coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}
# a clever optimizing compiler will clean this up
# there might be a way to do this with fewer flops using the tangent vector


  # calculate normal vector in xy space
  nx = nrm_xy[1]
  ny = nrm_xy[2]
#  nx2 = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
#  ny2 = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]
  fac = 1.0/(sqrt(nx*nx + ny*ny))
  # normalize normal vector
  nx *= fac
  ny *= fac

  # Get the normal momentum
  Unrm = nx*q[2] + ny*q[3]

  # this is equivalent to:
  #   1. computing the normal and tangential components
  #   2. negating the normal component
  #   3. combining the negative normal and non-negated tangent component
  qg = params.qg
  qg[1] = q[1]
  qg[2] = -2*Unrm*nx + q[2]
  qg[3] = -2*Unrm*ny + q[3]
  qg[4] = q[4]

  calcLFFlux(params, q, qg, aux_vars,nrm_xy, bndryflux)

  return nothing
end

function (obj::noPenetrationESBC)(params::ParamType3,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},  coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}
# a clever optimizing compiler will clean this up
# there might be a way to do this with fewer flops using the tangent vector


  # calculate normal vector in xy space
  nx = nrm_xy[1]
  ny = nrm_xy[2]
  nz = nrm_xy[3]
#  nx = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2] + dxidx[3,1]*nrm[3]
#  ny = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2] + dxidx[3,2]*nrm[3]
#  nz = dxidx[1,3]*nrm[1] + dxidx[2,3]*nrm[2] + dxidx[3,3]*nrm[3]
  fac = 1.0/(sqrt(nx*nx + ny*ny + nz*nz))
  # normalize normal vector
  nx = nx*fac
  ny = ny*fac
  nz = nz*fac

  # this is momentum, not velocity?
  Unrm = nx*q[2] + ny*q[3] + nz*q[4]

  qg = params.qg
  qg[1] = q[1]
  qg[2] = -2*Unrm*nx + q[2]
  qg[3] = -2*Unrm*ny + q[3]
  qg[4] = -2*Unrm*nz + q[4]
  qg[5] = q[5]

  calcLFFlux(params, q, qg, aux_vars,nrm_xy, bndryflux)

  return nothing
end



@doc """
###EulerEquationMod.noPenetrationBC_revm

Reverse mode for noPenetrationBC.

**Input**

* `obj` : Type of the Boundary condition being evaluated. Its a subtype of
          BCType_revm
* `q`   : Solution variable
* `aux_vars` : Auxiliary variables
* `x`     : Node coordinates
* `dxidx` : Mapping jacobian matrix for the SBP node
* `nrm`   : sbpface normal vector
* `bndryflux_bar` : Input flux value seed that is used to compute the reverse
                    mode derivative.
* `params`        : equation object parameters

"""->

mutable struct noPenetrationBC_revm <: BCType_revm
end

function (obj::noPenetrationBC_revm)(params::ParamType2,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},  coords::AbstractArray{Tmsh,1},
              nrm::AbstractArray{Tmsh,1}, nrm_bar::AbstractVector{Tmsh},
              bndryflux_bar::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  # Forward sweep
  n1 = nrm[1]
  n2 = nrm[2]
#  n1 = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
#  n2 = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]
  fac = 1.0/(sqrt(n1*n1 + n2*n2))
  nx = n1*fac
  ny = n2*fac
  Unrm = nx*q[2] + ny*q[3]
  qg = params.qg
  for i=1:length(q)
    qg[i] = q[i]
  end

  # Subtract the normal component of the momentum from \xi & \eta components
  # of the momentum
  qg[2] = qg[2] - nx*Unrm
  qg[3] = qg[3] - ny*Unrm

  v_vals = params.v_vals
  convertFromNaturalToWorkingVars(params, qg, v_vals)

  # Reverse sweep
#  nrm2_bar = zeros(Tmsh, 2)
  qg_bar = zeros(Tsol, 4)
  calcEulerFlux_revm(params, v_vals, aux_vars, nrm, bndryflux_bar, nrm_bar)
  calcEulerFlux_revq(params, v_vals, aux_vars, nrm, bndryflux_bar, qg_bar)

  # TODO: reverse mode convertFromNaturalToWorkingVars(params, qg, v_vals)
  n1_bar = nrm_bar[1]
  n2_bar = nrm_bar[2]


  # qg[2] = qg[2] - nx*Unrm
  # qg[3] = qg[3] - ny*Unrm
  ny_bar = -qg_bar[3]*Unrm
  nx_bar = -qg_bar[2]*Unrm
  Unrm_bar = -qg_bar[3]*ny -qg_bar[2]*nx

  # Unrm = nx*q[2] + ny*q[3]
  nx_bar += Unrm_bar*q[2]
  ny_bar += Unrm_bar*q[3]

  # nx = n1*fac
  # ny = n2*fac
  n2_bar += ny_bar*fac
  n1_bar += nx_bar*fac
  fac_bar = ny_bar*n2 + nx_bar*n1

  # fac = 1.0/(sqrt(n1*n1 + n2*n2))
  n1_bar -= fac_bar*n1*((n1*n1 + n2*n2)^(-1.5))
  n2_bar -= fac_bar*n2*((n1*n1 + n2*n2)^(-1.5))

  # because n1_bar and n2_bar came out of nrm_bar earlier, this is an assigment
  nrm_bar[1] = n1_bar
  nrm_bar[2] = n2_bar

  #=
  # n1 = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
  dxidx_bar[1,1] += n1_bar*nrm[1]
  dxidx_bar[2,1] += n1_bar*nrm[2]
  # n2 = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]
  dxidx_bar[1,2] += n2_bar*nrm[1]
  dxidx_bar[2,2] += n2_bar*nrm[2]
  =#
  return nothing
end
#=
function (obj::noPenetrationBC_revm){Tmsh, Tsol, Tres}(params::ParamType3,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},  coords::AbstractArray{Tmsh,1},
              dxidx::AbstractArray{Tmsh,2}, dxidx_bar::AbstractArray{Tmsh, 2},
              nrm::AbstractArray{Tmsh,1}, bndryflux_bar::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode)

  # Forward sweep
  nx = zero(Tmsh)
  ny = zero(Tmsh)
  nz = zero(Tmsh)
  nx2 = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2] + dxidx[3,1]*nrm[3]
  ny2 = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2] + dxidx[3,2]*nrm[3]
  nz2 = dxidx[1,3]*nrm[1] + dxidx[2,3]*nrm[2] + dxidx[3,3]*nrm[3]
  fac = 1.0/(sqrt(nx2*nx2 + ny2*ny2 + nz2*nz2))
  # normalize normal vector
  nx = nx2*fac
  ny = ny2*fac
  nz = nz2*fac

  # this is momentum, not velocity?
  Unrm = nx*q[2] + ny*q[3] + nz*q[4]

  qg = params.qg
  for i=1:length(q)
    qg[i] = q[i]
  end

  # calculate normal velocity
  qg[2] -= nx*Unrm
  qg[3] -= ny*Unrm
  qg[4] -= nz*Unrm

  v_vals = params.v_vals
  convertFromNaturalToWorkingVars(params, qg, v_vals)

  # Reverse Sweep
  nrm2_bar = zeros(Tmsh, 3)
  v_vals_bar = zeros(Tsol, 5)
  calcEulerFlux_revm(params, v_vals, aux_vars, [nx2, ny2, nz2], bndryflux_bar, nrm2_bar)
  calcEulerFlux_revq(params, v_vals, aux_vars, [nx2, ny2, nz2], bndryflux_bar, v_vals_bar)

  nz2_bar = nrm2_bar[3]
  ny2_bar = nrm2_bar[2]
  nx2_bar = nrm2_bar[1]

  # TODO: reverse mode convertFromNaturalToWorkingVars(params, qg, v_vals)
  qg_bar = v_vals_bar

  # qg[4] = qg[4] - nz*Unrm
  nz_bar = -qg_bar[4]*Unrm
  Unrm_bar = -qg_bar[4]*nz
  # qg[3] -= ny*Unrm
  ny_bar = -qg_bar[3]*Unrm
  Unrm_bar -= qg_bar[3]*ny
  # qg[2] -= nx*Unrm
  nx_bar = -qg_bar[2]*Unrm
  Unrm_bar -= qg_bar[2]*nx

  # Unrm = nx*q[2] + ny*q[3] + nz*q[4]
  nx_bar += Unrm_bar*q[2]
  ny_bar += Unrm_bar*q[3]
  nz_bar += Unrm_bar*q[4]

  # nz = nz2*fac
  nz2_bar += nz_bar*fac
  fac_bar = nz_bar*nz2
  # ny = ny2*fac
  ny2_bar += ny_bar*fac
  fac_bar += ny_bar*ny2
  # nx = nx2*fac
  nx2_bar += nx_bar*fac
  fac_bar += nx_bar*nx2

  # fac = 1.0/(sqrt(nx2*nx2 + ny2*ny2 + nz2*nz2))
  nx2_bar -= fac_bar*nx2*((nx2*nx2 + ny2*ny2 + nz2*nz2)^(-1.5))
  ny2_bar -= fac_bar*ny2*((nx2*nx2 + ny2*ny2 + nz2*nz2)^(-1.5))
  nz2_bar -= fac_bar*nz2*((nx2*nx2 + ny2*ny2 + nz2*nz2)^(-1.5))

  # nz2 = dxidx[1,3]*nrm[1] + dxidx[2,3]*nrm[2] + dxidx[3,3]*nrm[3]
  dxidx_bar[1,3] += nz2_bar*nrm[1]
  dxidx_bar[2,3] += nz2_bar*nrm[2]
  dxidx_bar[3,3] += nz2_bar*nrm[3]
  # ny2 = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2] + dxidx[3,2]*nrm[3]
  dxidx_bar[1,2] += ny2_bar*nrm[1]
  dxidx_bar[2,2] += ny2_bar*nrm[2]
  dxidx_bar[3,2] += ny2_bar*nrm[3]
  # nx2 = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2] + dxidx[3,1]*nrm[3]
  dxidx_bar[1,1] += nx2_bar*nrm[1]
  dxidx_bar[2,1] += nx2_bar*nrm[2]
  dxidx_bar[3,1] += nx2_bar*nrm[3]

  return nothing
end # End noPenetrationBC_revm 3D
=#
@doc """
### EulerEquationMod.unsteadyVortexBC <: BCTypes

  This type and the associated call method define a functor to calculate
  the flux using the Roe Solver using the exact InsentropicVortex solution
  as boundary state.  See calcBoundaryFlux for the arguments all functors
  must support.

  This is a low level functor.

**Arguments**

*  `obj` : Object of type BCType used for multiple dispatch. Every new boundary
           condition needs to have its own type and entered in BCDict
*  `q`   : Solution variable
*  `aux_vars` : Auxiliary variables
*  `x`        : physical coordinates of the SBP node
*  `dxidx`    : Mapping jacobian matrix for the SBP node
*  `nrm_xy`      : SBP face normal
*  `bndryflux` : Computed flux value at the boundary

"""->
mutable struct unsteadyVortexBC <: BCType
end

# low level function
function (obj::unsteadyVortexBC)(params::ParamType,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},  coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}


#  println("entered isentropicOvrtexBC (low level)")
#  println("Tsol = ", Tsol)
  # getting qg
  qg = params.qg
  calcUnsteadyVortex(params, coords, qg)

  RoeSolver(params, q, qg, aux_vars, nrm_xy, bndryflux)

  return nothing

end # ends the function unsteadyVortex BC

mutable struct unsteadyVortex2BC <: BCType
end

# low level function
function (obj::unsteadyVortex2BC)(params::ParamType,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},  coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}


#  println("entered isentropicOvrtexBC (low level)")
#  println("Tsol = ", Tsol)
  # getting qg
  qg = params.qg
  calcUnsteadyVortex2(params, coords, qg)

  RoeSolver(params, q, qg, aux_vars, nrm_xy, bndryflux)

  return nothing

end # ends the function unsteadyVortex BC


@doc """
### EulerEquationMod.Rho1E2U1VW0BC <: BCTypes

  This functor uses the Roe solver to calculate the flux for a boundary
  state where the fluid density is 1, energy = 2, and 
    u = 1.0, v = 0.0, and w = 0.0 (if 3D).

  It should work for 2D and 3D meshes.

  This is a low level functor.

**Arguments**

*  `obj` : Object of type BCType used for multiple dispatch. Every new boundary
           condition needs to have its own type and entered in BCDict
*  `q`   : Solution variable
*  `aux_vars` : Auxiliary variables
*  `x`        : physical coordinates of the SBP node
*  `dxidx`    : Mapping jacobian matrix for the SBP node
*  `nrm`      : SBP face normal
*  `bndryflux` : Computed flux value at the boundary

"""
mutable struct Rho1E2U1VW0BC <: BCType
end

# low level function
function (obj::Rho1E2U1VW0BC)(params::ParamType,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},
              coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}


  #println("in Rho1E2BCU1V0W0")
  qg = params.qg

  calcRho1Energy2U1VW0(params, coords, qg)

  #println("qg = ", qg)
  # call Roe solver
  RoeSolver(params, q, qg, aux_vars, nrm_xy, bndryflux)

  return nothing

end

@doc """
### EulerEquationMod.Rho1E2BC <: BCTypes

  This functor uses the Roe solver to calculate the flux for a boundary
  state where the fluid density is 1, energy = 2, and u = v = 0.0

  This is a low level functor

**Arguments**

*  `obj` : Object of type BCType used for multiple dispatch. Every new boundary
           condition needs to have its own type and entered in BCDict
*  `q`   : Solution variable
*  `aux_vars` : Auxiliary variables
*  `x`        : physical coordinates of the SBP node
*  `dxidx`    : Mapping jacobian matrix for the SBP node
*  `nrm`      : SBP face normal
*  `bndryflux` : Computed flux value at the boundary

"""
mutable struct Rho1E2BC <: BCType
end

# low level function
function (obj::Rho1E2BC)(params::ParamType,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},
              coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}


  #println("in Rho1E2BC")
  qg = params.qg

  calcRho1Energy2(params, coords, qg)

  #println("qg = ", qg)
  # call Roe solver
  RoeSolver(params, q, qg, aux_vars, nrm_xy, bndryflux)

  return nothing

end

@doc """
### EulerEquationMod.Rho1E2U3BC <: BCTypes

  This functor uses the Roe solver to calculate the flux for a boundary
  state where the fluid density is 1, energy = 2, and x velocity = 0.5

  This is a low level functor

**Arguments**

*  `obj` : Object of type BCType used for multiple dispatch. Every new boundary
           condition needs to have its own type and entered in BCDict
*  `q`   : Solution variable
*  `aux_vars` : Auxiliary variables
*  `x`        : physical coordinates of the SBP node
*  `dxidx`    : Mapping jacobian matrix for the SBP node
*  `nrm`      : SBP face normal
*  `bndryflux` : Computed flux value at the boundary

"""
mutable struct Rho1E2U3BC <: BCType
end

# low level function
function (obj::Rho1E2U3BC)(params::ParamType,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},
              coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}


  #println("in Rho1E2U3Bc")
  qg = params.qg

  calcRho1Energy2U3(params, coords, qg)

  #println("qg = ", qg)
  # call Roe solver
  RoeSolver(params, q, qg, aux_vars, nrm_xy, bndryflux)

  return nothing

end

@doc """
### EulerEquationMod.FreeStreamBC <: BCTypes

  This functor uses the Roe solver to calculate the flux for a boundary
  state corresponding to the free stream velocity, using rho_free, Ma, aoa, and E_free

  Works in 2D and 3D

  This is a low level functor

**Arguments**

*  `obj` : Object of type BCType used for multiple dispatch. Every new boundary
           condition needs to have its own type and entered in BCDict
*  `q`   : Solution variable
*  `aux_vars` : Auxiliary variables
*  `x`        : physical coordinates of the SBP node
*  `dxidx`    : Mapping jacobian matrix for the SBP node
*  `nrm_xy`      : SBP face normal
*  `bndryflux` : Computed flux value at the boundary

"""
mutable struct FreeStreamBC <: BCType
end

function (obj::FreeStreamBC)(params::ParamType,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},  coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  qg = params.qg

  calcFreeStream(params, coords, qg)
  RoeSolver(params, q, qg, aux_vars, nrm_xy, bndryflux)

  return nothing
end

@doc """
###EulerEquationMod.FreeStreamBC_revm

Reverse mode for FreeStreamBC.

**Inputs**

* `obj` : Type of the Boundary condition being evaluated. Its a subtype of
          BCType_revm
* `q`   : Solution variable
* `aux_vars` : Auxiliary variables
* `x`     : Node coordinates
* `nrm_xy`   : scaled normal vector in x-y space
* `bndryflux_bar` : Input flux value seed that is used to compute the reverse
                    mode derivative.
* `params`        : equation object parameters

**Output**

* `nrm_bar` : Derivative of bndryflux_bar w.r.t the mapping jacobian

"""->

mutable struct FreeStreamBC_revm <: BCType_revm
end

function (obj::FreeStreamBC_revm)(params::ParamType,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},  coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1}, nrm_bar::AbstractVector{Tmsh},
              bndryflux_bar::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  # Forward sweep
  qg = params.qg
  calcFreeStream(params, coords, qg)

  # Reverse sweep
  RoeSolver_revm(params, q, qg, aux_vars, nrm_xy, bndryflux_bar, nrm_bar)

  return nothing
end

@doc """
### EulerEquationMod.FreeStreamBC_dAlpha <: BCTypes

  This functor uses the Roe solver to calculate the flux for a boundary
  state corresponding to the free stream velocity, using rho_free, Ma, aoa, and E_free

  This is a low level functor

**Arguments**

*  `obj` : Object of type BCType used for multiple dispatch. Every new boundary
           condition needs to have its own type and entered in BCDict
*  `q`   : Solution variable
*  `aux_vars` : Auxiliary variables
*  `x`        : physical coordinates of the SBP node
*  `nrm_xy`      : scaled normal vector in x-y space
*  `bndryflux` : Computed flux value at the boundary

"""->
mutable struct FreeStreamBC_dAlpha <: BCType
end

function (obj::FreeStreamBC_dAlpha)(params::ParamType2,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},  coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1}, nrm_bar::AbstractVector{Tmsh},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  qg = params.qg

  calcFreeStream_dAlpha(params, coords, qg)
  RoeSolver(params, q, qg, aux_vars, nrm, bndryflux)

  return nothing
end


@doc """
### EulerEquationMod.allOnesBC <: BCTypes

  This functor uses the Roe solver to calculate the flux for a boundary
  state where all the conservative variables have a value 1.0

  This is a low level functor
"""->

mutable struct allOnesBC <: BCType
end

function (obj::allOnesBC)(params::ParamType2,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  qg = zeros(Tsol, 4)
  calcOnes(params, coords, qg)

  RoeSolver(params, q, qg, aux_vars, nrm_xy, bndryflux)

  # println("bndryflux = ", bndryflux)
  return nothing
end # end function call

@doc """
### EulerEquationMod.allZerosBC <: BCTypes

  This functor uses the Roe solver to calculate the flux for a boundary
  state where all the conservative variables have a value 0.0

  This is a low level functor
"""->

mutable struct allZerosBC <: BCType
end

function (obj::allZerosBC)(params::ParamType2,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  qg = zeros(Tsol, 4)
  calcZeros(params, coords, qg)

  RoeSolver(params, q, qg, aux_vars, nrm_xy, bndryflux)

  # println("bndryflux = ", bndryflux)
  return nothing
end # end function

mutable struct ExpBC <: BCType
end

function (obj::ExpBC)(params::ParamType, q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  qg = params.qg
  calcExp(params, coords, qg)
  RoeSolver(params, q, qg, aux_vars, nrm_xy, bndryflux)

  # println("bndryflux = ", bndryflux)
  return nothing
end # end function

mutable struct ExpBC_revm <: BCType_revm
end

function (obj::ExpBC_revm)(params::ParamType,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},  coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1}, nrm_bar::AbstractVector{Tmsh},
              bndryflux_bar::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  # Forward Sweep
  qg = params.qg
  calcExp(params, coords, qg)

  # Reverse Sweep

  # RoeSolver(params, q, qg, aux_vars, nrm_xy, bndryflux)
  RoeSolver_revm(params, q, qg, aux_vars, nrm_xy, bndryflux_bar, nrm_bar)

  return nothing
end

mutable struct PeriodicMMSBC <: BCType
end

function (obj::PeriodicMMSBC)(params::ParamType,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}
# use the exact solution as the boundary condition for the PeriodicMMS
# solutions

  qg = params.qg
  calcPeriodicMMS(params, coords, qg)
  use_efix = 0
  RoeSolver(params, q, qg, aux_vars, nrm_xy, bndryflux, use_efix)

  # println("bndryflux = ", bndryflux)
  return nothing
end # end function

mutable struct ChannelMMSBC <: BCType
end

function (obj::ChannelMMSBC)(params::ParamType,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}
# use the exact solution as the boundary condition for the ChannelMMS
# solutions

  qg = params.qg
  calcChannelMMS(params, coords, qg)
  RoeSolver(params, q, qg, aux_vars, nrm_xy, bndryflux)

  return nothing
end # end function


mutable struct defaultBC <: BCType
end

function (obj::defaultBC)(params::ParamType,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  calcEulerFlux(params, q, aux_vars, nrm_xy, bndryflux)

  return nothing
end

mutable struct SubsonicInflowBC <: BCType
end

function (obj::SubsonicInflowBC)(params::ParamType2,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}

  #See NASA/TM-2011-217181: Inflow/Outflow Boundary Conditions with Application
  #                         to FUN3D by Carlson 
  # The derivation has some algebraic mistakes, but the approach is correct

  pt = 102010.0/params.p_free  # boundary stagnation pressure
  Tt = 288.6/params.T_free  # boundary stagnation temperature
  # i = interior quantity
  # b = boundary state

  # need normalized outward normal vector
  nrm_fac = 1/sqrt(nrm_xy[1]*nrm_xy[1] + nrm_xy[2]*nrm_xy[2])

  gamma = params.gamma
  gamma_1 = params.gamma_1

  pressi = calcPressure(params, q)
  # magnitude of velocity (negative sign because the normal is outward but
  # the velocity should be inward
  Ui = -nrm_fac*(q[2]*nrm_xy[1] + q[3]*nrm_xy[2])/q[1]
#  vi = q[3]/q[1]
  ai2 = gamma*pressi/q[1]  # speed of sound squared

  # stagnation enthalpy (specific)
  hti = ai2/gamma_1 + Ui*Ui

  # Riemann invarient for the characteristic exiting the domain
  Ri = Ui - 2*sqrt(ai2)/gamma_1

  # this step uses the adiabatic assumption + the Riemann invarient Rb to
  # form a quandratic equation for ab
  # pick the larger of the two roots

  a = 0.5 + 2/gamma_1
  b = 2*Ri
  c = 0.5*gamma_1*(Ri*Ri - hti)

  tmp1 = -0.5*b/a
  tmp2 = 0.5*sqrt(b*b - 4*a*c)/a
  ab1 = tmp1 + tmp2
  ab2 = tmp1 - tmp2

  ab = max(ab1, ab2)  # maximum root is the physically correct one

  # use Riemann invarient to find velocity magnitude on the boundary side
  Ub = Ri + 2*ab/gamma_1
  Mb = Ub/ab

  @assert Mb < 1.0

  operand = 1/(1 + 0.5*gamma_1*Mb*Mb)
  pb = pt*(operand)^(gamma/gamma_1)
  Tb = Tt*operand


  # convert back to conservative variables
  qg = params.qg
  rho1 = gamma*pb/(ab*ab)  # this is numerically equivalent to rho2 below,
                           # which is weuird because in this case Tb is never
                           # used
  rho2 = pb/(params.R_ND*Tb)
  qg[1] = rho2  # R is not nondimenstionalized
  qg[2] = -Ub*nrm_xy[1]*nrm_fac*qg[1]  # negative sign because the normal is neg
  qg[3] = -Ub*nrm_xy[2]*nrm_fac*qg[1]
  qg[4] = pb/gamma_1 + 0.5*qg[1]*Ub*Ub

#  R_computed = pb/(qg[1]*Tb)
  RoeSolver(params, q, qg, aux_vars, nrm_xy, bndryflux)


  return nothing
end

mutable struct SubsonicOutflowBC <: BCType
end

function(obj::SubsonicOutflowBC)(params::ParamType2,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}


  pb = 101300.0/params.p_free  # nondimensionalized pressure
  gamma = params.gamma
  gamma_1 = params.gamma_1

  pressi = calcPressure(params, q)
  # verify Mach number < 1
  ai2 = gamma*pressi/q[1]  # speed of sound squared
  # need normalized outward normal vector
  nrm_fac = 1/sqrt(nrm_xy[1]*nrm_xy[1] + nrm_xy[2]*nrm_xy[2])
  Un = (q[2]*nrm_xy[1] + q[3]*nrm_xy[2])*nrm_fac/q[1]

  @assert Un >= 0  # this should be outflow, not inflow
  @assert (Un*Un)/ai2 < 1

  qg = params.qg
  qg[1] = q[1]
  qg[2] = q[2]
  qg[3] = q[3]
  # compute energy from the specified pressure
  qg[4] = pb/gamma_1 + 0.5*(q[2]*q[2] + q[3]*q[3])/q[1]
#  qg[4] = pb/gamma_1 + 0.5*q[1]*(q[2]*q[2]

  RoeSolver(params, q, qg, aux_vars, nrm_xy, bndryflux)

  return nothing
end

mutable struct inviscidChannelFreeStreamBC <: BCType
end

# low level function
function(obj::inviscidChannelFreeStreamBC)(             
              params::ParamType,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},  coords::AbstractArray{Tmsh,1},
              nrm_xy::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1},
              bndry::BoundaryNode=NullBoundaryNode) where {Tmsh, Tsol, Tres}


#  println("entered isentropicOvrtexBC (low level)")
#  println("Tsol = ", Tsol)
  # getting qg
  qg = params.qg
  calcInvChannelFreeStream(params, coords, qg)
  RoeSolver(params, q, qg, aux_vars, nrm_xy, bndryflux)

  return nothing

end # ends the function unsteadyVortex BC

"""
   This is a special boundary condition used for imposing a numerical
   solution as a boundary condition.

   **Fields**

    * bc_vals: an array of size numDofPerNode x numNodesPerFace x numFaces with
               this boundary condition on it.  This array can be accessed
               using the `faceidx` field of [`BoundaryNode`](@ref).


  This BC is special because it has to store the boundary values it is imposing.
  Whenever this boundary condition is needed, the user should construct a
  new object, with the data inside it, and then register it just before
  constructing the `EulerData` object.
"""
mutable struct reanalysisBC{Tsol} <: BCType
  bc_vals::Array{Tsol, 3}
end

function reanalysisBC()
  bc_vals = Array(Float64, 0, 0, 0)
  return reanalysisBC{Float64}(bc_vals)
end

function (obj::reanalysisBC)(params::AbstractParamType{Tdim},
        q::AbstractArray{Tsol, 1},
        aux_vars::AbstractArray{Tres, 1},
        coords::AbstractVector{Tmsh},
        nrm_xy::AbstractVector{Tmsh},
        bndryflux::AbstractArray{Tres},
        bndry::EulerEquationMod.BoundaryNode) where {Tmsh, Tsol, Tres, Tdim}


  # get numerical data out of the array
  qg = sview(obj.bc_vals, :, bndry.node, bndry.faceidx)
  EulerEquationMod.RoeSolver(params, q, qg, aux_vars, nrm_xy, bndryflux)
  
  return nothing
end

include("bc_viscous.jl")

# every time a new boundary condition is created,
# add it to the dictionary
#const isentropicVortexBC_ = isentropicVortexBC()
#const noPenetrationBC_ = noPenetrationBC()

"""
  Maps boundary conditions names to the functor objects.
  Each functor should be callable with the signature
"""
global const BCDict = Dict{String, BCType}(
"isentropicVortexBC" => isentropicVortexBC(),
"noPenetrationBC" => noPenetrationBC(),
"noPenetrationESBC" => noPenetrationESBC(),
"Rho1E2BC" => Rho1E2BC(),
"Rho1E2U1VW0BC" => Rho1E2U1VW0BC(),
"Rho1E2U3BC" => Rho1E2U3BC(),
"isentropicVortexBC_physical" => isentropicVortexBC_physical(),
"FreeStreamBC" => FreeStreamBC(),
"allOnesBC" => allOnesBC(),
"unsteadyVortexBC" => unsteadyVortexBC(),
"unsteadyVortex2BC" => unsteadyVortex2BC(),
"ExpBC" => ExpBC(),
"PeriodicMMSBC" => PeriodicMMSBC(),
"ChannelMMSBC" => ChannelMMSBC(),
"subsonicInflowBC" => SubsonicInflowBC(),
"subsonicOutflowBC" => SubsonicOutflowBC(),
"inviscidChannelFreeStreamBC" => inviscidChannelFreeStreamBC(),
"reanalysisBC" => reanalysisBC(),
"defaultBC" => defaultBC(),
"nonslipBC" => nonslipBC(),
"ExactChannelBC" => ExactChannelBC(),
"zeroPressGradientBC" => zeroPressGradientBC(),
)

@doc """
### EulerEquationMod.getBCFunctors

  This function uses the opts dictionary to populate mesh.bndry_funcs with
  the the functors

    func(params::ParamType,
         q::AbstractArray{Tsol,1},
         aux_vars::AbstractArray{Tres, 1},  coords::AbstractArray{Tmsh,1},
         nrm_xy::AbstractArray{Tmsh,1},
         bndryflux::AbstractArray{Tres, 1},
         bndry::BoundaryNode=NullBoundaryNode)


  This is a high level function.
"""->
# use this function to populate access the needed values in BCDict
function getBCFunctors(mesh::AbstractMesh, sbp::AbstractSBP, eqn::EulerData, opts)
# populate the array mesh.bndry_funcs with the functors for the boundary condition types

#  println("Entered getBCFunctors")

  for i=1:mesh.numBC
    key_i = string("BC", i, "_name")
    val = opts[key_i]
#    println("BCDict[val] = ", BCDict[val])
    mesh.bndry_funcs[i] = BCDict[val]
  end

  return nothing
end # ENd function getBCFunctors

global const BCDict_revm = Dict{String, BCType_revm}(
"noPenetrationBC" => noPenetrationBC_revm(),
"FreeStreamBC" => FreeStreamBC_revm(),
"ExpBC" => ExpBC_revm(),
"isentropicVortexBC" => isentropicVortexBC_revm(),
)

function getBCFunctors_revm(mesh::AbstractMesh, sbp::AbstractSBP, eqn::EulerData, opts)

  for i = 1:mesh.numBC
    key_i = string("BC", i, "_name")
    val = opts[key_i]
#    println("BCDict_revm[$val] = ", BCDict_revm[val])
    mesh.bndry_funcs_revm[i] = BCDict_revm[val]
  end # End for i = 1:mesh.numBC

  return nothing
end # End function getBCFunctors_revm
