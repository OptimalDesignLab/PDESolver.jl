# functions for calculating boundary integrals

include("bc_solvers.jl")  # Roe solvers and related things


@doc """
### EulerEquationMod.getBCFluxes

  This function calls other functions to calculate the boundary fluxes, passing
  them pieces of the array needed.  This populates eqn.bndryflux.  It also
  calls writeBoundary() to do any requested output.

  This is a mid level function
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
    idx_range = start_index:end_index  # TODO: should this be start_index:(end_index - 1) ?
    bndry_facenums_i = sview(mesh.bndryfaces, start_index:(end_index - 1))
    bndryflux_i = sview(eqn.bndryflux, :, :, start_index:(end_index - 1))

    # call the function that calculates the flux for this boundary condition
    # passing the functor into another function avoid type instability
    calcBoundaryFlux(mesh, sbp, eqn, functor_i, idx_range, bndry_facenums_i, bndryflux_i)
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
    for j=1:sbp.numfacenodes  # TODO: should be mesh.numNodesPerFace?
      jb = sbp.facenodes[j, face]
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
      q_vals = sview(eqn.q_bndry, :, j, i)
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
  functor( q, aux_vars, x, dxidx, nrm, bndryflux_i, eqn.params)
  where q are the *conservative* variables.
  where all arguments (except params and nrm) are vectors of values at a node.

  params is the ParamType associated with the the EulerEquation object
  nrm = sbp.facenormal[:, current_node]

  This is a mid level function.
"""->
# mid level function
function calcBoundaryFlux{Tmsh,  Tsol, Tres}( mesh::AbstractCGMesh{Tmsh},
                          sbp::AbstractSBP, eqn::EulerData{Tsol},
                          functor::BCType, idx_range::UnitRange,
                          bndry_facenums::AbstractArray{Boundary,1},
                          bndryflux::AbstractArray{Tres, 3})
# calculate the boundary flux for the boundary condition evaluated by the functor

#  println("enterted calcBoundaryFlux")


  nfaces = length(bndry_facenums)
  q2 = zeros(Tsol, mesh.numDofPerNode)
  for i=1:nfaces  # loop over faces with this BC
    bndry_i = bndry_facenums[i]
#    println("boundary ", i, "element = ", bndry_i.element, ", face = ", bndry_i.face)
#    println("interface ", i)
    for j = 1:sbp.numfacenodes
      k = sbp.facenodes[j, bndry_i.face]

      # get components
      q = sview(eqn.q, :, k, bndry_i.element)
      # convert to conservative variables if needed
      convertToConservative(eqn.params, q, q2)
      aux_vars = sview(eqn.aux_vars, :, k, bndry_i.element)
      x = sview(mesh.coords, :, k, bndry_i.element)
      dxidx = sview(mesh.dxidx, :, :, k, bndry_i.element)
      nrm = sview(sbp.facenormal, :, bndry_i.face)
      bndryflux_i = sview(bndryflux, :, j, i)

      functor(q2, aux_vars, x, dxidx, nrm, bndryflux_i, eqn.params)

    end

  end


  return nothing
end

# DG version
function calcBoundaryFlux{Tmsh,  Tsol, Tres}( mesh::AbstractDGMesh{Tmsh},
                          sbp::AbstractSBP, eqn::EulerData{Tsol},
                          functor::BCType, idx_range::UnitRange,
                          bndry_facenums::AbstractArray{Boundary,1},
                          bndryflux::AbstractArray{Tres, 3})
  # calculate the boundary flux for the boundary condition evaluated by the
  # functor

  nfaces = length(bndry_facenums)
  q2 = zeros(Tsol, mesh.numDofPerNode)
  nrm = zeros(Tmsh, size(sbp.facenormal,1))
  for i=1:nfaces  # loop over faces with this BC
    bndry_i = bndry_facenums[i]
    global_facenum = idx_range[i]
    for j = 1:mesh.numNodesPerFace

      # get components
      q = sview(eqn.q_bndry, :, j, global_facenum)
      # convert to conservative variables if needed
      convertToConservative(eqn.params, q, q2)
      aux_vars = sview(eqn.aux_vars_bndry, :, j, global_facenum)
      x = sview(mesh.coords_bndry, :, j, global_facenum)
      dxidx = sview(mesh.dxidx_bndry, :, :, j, global_facenum)
      nrm[:] = sbp.facenormal[:,bndry_i.face]
      bndryflux_i = sview(bndryflux, :, j, i)

      # DEBUGGING: use analytical solution (avoid interpolation inexactness)
#      calcPeriodicMMS(x, eqn.params, q2)
#      println("after overwriting q with analytical solution, q_nodes = \n", q2)

#      println("coords = ", x)
      functor(q2, aux_vars, x, dxidx, nrm, bndryflux_i, eqn.params)
    end
  end

  return nothing
end


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
*  `dxidx`    : Mapping jacobian matrix for the SBP node
*  `nrm`      : SBP face normal
*  `bndryflux` : Computed flux value at the boundary

"""->
type isentropicVortexBC <: BCType
end

function call{Tmsh, Tsol, Tres}(obj::isentropicVortexBC,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, x::AbstractArray{Tmsh,1},
              dxidx::AbstractArray{Tmsh,2}, nrm::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1}, params::ParamType)

  gamma = params.gamma
  gami = params.gamma_1

  # getting qg
  qg = params.qg
  calcIsentropicVortex(x, params, qg) # Get the boundary value

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
  nrm2 = params.nrm2
  calcBCNormal(params, dxidx, nrm, nrm2)
  sat = params.sat_vals
  calcSAT(params, nrm2, dq, sat, [u, v], H)

  euler_flux = zeros(Tsol, 4) # params.flux_vals1
  calcEulerFlux(params, v_vals, aux_vars, nrm2, euler_flux)

  sat_fac = 1.0 # Multiplier for SAT term
  for i=1:4
    bndryflux[i] = euler_flux[i] + sat_fac*sat[i]
  end

  return nothing
end

#=
# low level function
function call{Tmsh, Tsol, Tres}(obj::isentropicVortexBC,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, x::AbstractArray{Tmsh,1},
              dxidx::AbstractArray{Tmsh,2}, nrm::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1}, params::ParamType)

  qg = params.qg
  calcIsentropicVortex(x, params, qg)
  RoeSolver(params, q, qg, aux_vars, dxidx, nrm, bndryflux)

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
* `dxidx` : Mapping jacobian matrix for the SBP node
* `nrm`   : sbpface normal vector
* `bndryflux_bar` : Input flux value seed that is used to compute the reverse
                    mode derivative.
* `params`        : equation object parameters

**Output**

* `dxidx_bar` : Derivative of bndryflux_bar w.r.t the mapping jacobian

"""->

type isentropicVortexBC_revm <: BCType_revm
end

function call{Tmsh, Tsol, Tres}(obj::isentropicVortexBC_revm, q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},  x::AbstractArray{Tmsh,1},
              dxidx::AbstractArray{Tmsh,2}, dxidx_bar::AbstractArray{Tmsh, 2},
              nrm::AbstractArray{Tmsh,1}, bndryflux_bar::AbstractArray{Tres, 1},
              params::ParamType{2})

  # Forward sweep
  gamma = params.gamma
  gami = params.gamma_1

  # getting qg
  qg = params.qg
  calcIsentropicVortex(x, params, qg) # Get the boundary value
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
  nrm2 = params.nrm2
  calcBCNormal(params, dxidx, nrm, nrm2)

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

  nrm2_bar = zeros(Tmsh, 2)
  calcEulerFlux_revm(params, v_vals, aux_vars, nrm2, euler_flux_bar, nrm2_bar)
  calcSAT_revm(params, nrm2, dq, [u,v], H, sat_bar, nrm2_bar)
  calcBCNormal_revm(params, dxidx, nrm, nrm2_bar, dxidx_bar)

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
type isentropicVortexBC_physical <: BCType
end

function call{Tmsh, Tsol, Tres}(obj::isentropicVortexBC_physical,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, x::AbstractArray{Tmsh,1},
              dxidx::AbstractArray{Tmsh,2}, nrm::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1}, params::ParamType{2})

  nx = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
  ny = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]

  calcEulerFlux(params, q, aux_vars, [nx, ny], bndryflux)

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
*  `dxidx`    : Mapping jacobian matrix for the SBP node
*  `nrm`      : SBP face normal
*  `bndryflux` : Computed flux value at the boundary

"""
type noPenetrationBC <: BCType
end

# low level function
function call{Tmsh, Tsol, Tres}(obj::noPenetrationBC, q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},  x::AbstractArray{Tmsh,1},
              dxidx::AbstractArray{Tmsh,2}, nrm::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1}, params::ParamType{2})
# a clever optimizing compiler will clean this up
# there might be a way to do this with fewer flops using the tangent vector


  # calculate normal vector in xy space
  nx = zero(Tmsh)
  ny = zero(Tmsh)
  nx2 = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
  ny2 = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]
  fac = 1.0/(sqrt(nx2*nx2 + ny2*ny2))
  # normalize normal vector
  nx = nx2*fac
  ny = ny2*fac

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


  calcEulerFlux(params, v_vals, aux_vars, [nx2, ny2], bndryflux)

  return nothing
end


function call{Tmsh, Tsol, Tres}(obj::noPenetrationBC, q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},  x::AbstractArray{Tmsh,1},
              dxidx::AbstractArray{Tmsh,2}, nrm::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1}, params::ParamType{3})
# a clever optimizing compiler will clean this up
# there might be a way to do this with fewer flops using the tangent vector


  # calculate normal vector in xy space
  nx = zero(Tmsh)
  ny = zero(Tmsh)
  nz = zero(Tmsh)
  tngt = Array(Tmsh, 2)  # tangent vector
  nx = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2] + dxidx[3,1]*nrm[3]
  ny = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2] + dxidx[3,2]*nrm[3]
  nz = dxidx[1,3]*nrm[1] + dxidx[2,3]*nrm[2] + dxidx[3,3]*nrm[3]
  fac = 1.0/(sqrt(nx*nx + ny*ny + nz*nz))
  # normalize normal vector
  nx *= fac
  ny *= fac
  nz *= fac

  # this is momentum, not velocity?
  Unrm = nx*q[2] + ny*q[3] + nz*q[4]

  qg = params.qg
  for i=1:length(q)
    qg[i] = q[i]
  end

  #qg = copy(q)

  # calculate normal velocity
  qg[2] -= nx*Unrm
  qg[3] -= ny*Unrm
  qg[4] -= nz*Unrm

  # call Roe solver
  #RoeSolver(params, q, qg, aux_vars, dxidx, nrm, bndryflux)
  nx2 = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2] + dxidx[3,1]*nrm[3]
  ny2 = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2] + dxidx[3,2]*nrm[3]
  nz2 = dxidx[1,3]*nrm[1] + dxidx[2,3]*nrm[2] + dxidx[3,3]*nrm[3]

  v_vals = params.v_vals
  convertFromNaturalToWorkingVars(params, qg, v_vals)
  # this is a problem: q is in conservative variables even if
  # params says we are using entropy variables
  calcEulerFlux(params, v_vals, aux_vars, [nx2, ny2, nz2], bndryflux)

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

type noPenetrationBC_revm <: BCType_revm
end

function call{Tmsh, Tsol, Tres}(obj::noPenetrationBC_revm, q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},  x::AbstractArray{Tmsh,1},
              dxidx::AbstractArray{Tmsh,2}, dxidx_bar::AbstractArray{Tmsh, 2},
              nrm::AbstractArray{Tmsh,1}, bndryflux_bar::AbstractArray{Tres, 1},
              params::ParamType{2})

  # Forward sweep
  n1 = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
  n2 = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]
  fac = 1.0/(sqrt(n1*n1 + n2*n2))
  nx = n1*fac
  ny = n2*fac
  Unrm = nx*q[2] + ny*q[3]

  # Subtract the normal component of the momentum from \xi & \eta components
  # of the momentum
  q[2] = q[2] - nx*Unrm
  q[3] = q[3] - ny*Unrm

  v_vals = params.v_vals
  convertFromNaturalToWorkingVars(params, q, v_vals)

  # Reverse sweep
  nrm2_bar = zeros(Tmsh, 2)
  q_bar = zeros(Tsol, 4)
  calcEulerFlux_revm(params, v_vals, aux_vars, [n1, n2], bndryflux_bar, nrm2_bar)
  calcEulerFlux_revq(params, v_vals, aux_vars, [n1, n2], bndryflux_bar, q_bar)

  # TODO: reverse mode convertFromNaturalToWorkingVars(params, qg, v_vals)
  n1_bar = nrm2_bar[1]
  n2_bar = nrm2_bar[2]


  # q[2] = q[2] - nx*Unrm
  # q[3] = q[3] - ny*Unrm
  ny_bar = -q_bar[3]*Unrm
  nx_bar = -q_bar[2]*Unrm
  Unrm_bar = -q_bar[3]*ny -q_bar[2]*nx

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

  # n1 = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
  dxidx_bar[1,1] += n1_bar*nrm[1]
  dxidx_bar[2,1] += n1_bar*nrm[2]
  # n2 = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]
  dxidx_bar[1,2] += n2_bar*nrm[1]
  dxidx_bar[2,2] += n2_bar*nrm[2]

  return nothing
end # End noPenetrationBC_revm 2D

function call{Tmsh, Tsol, Tres}(obj::noPenetrationBC_revm, q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},  x::AbstractArray{Tmsh,1},
              dxidx::AbstractArray{Tmsh,2}, dxidx_bar::AbstractArray{Tmsh, 2},
              nrm::AbstractArray{Tmsh,1}, bndryflux_bar::AbstractArray{Tres, 1},
              params::ParamType{3})

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
*  `nrm`      : SBP face normal
*  `bndryflux` : Computed flux value at the boundary

"""->
type unsteadyVortexBC <: BCType
end

# low level function
function call{Tmsh, Tsol, Tres}(obj::unsteadyVortexBC, q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},  x::AbstractArray{Tmsh,1},
              dxidx::AbstractArray{Tmsh,2}, nrm::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1}, params::ParamType{2})


#  println("entered isentropicOvrtexBC (low level)")
#  println("Tsol = ", Tsol)

  # getting qg
  qg = params.qg
  calcUnsteadyVortex(x, params, qg)

  RoeSolver(params, q, qg, aux_vars, dxidx, nrm, bndryflux)

  return nothing

end # ends the function unsteadyVortex BC

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
type Rho1E2U3BC <: BCType
end

# low level function
function call{Tmsh, Tsol, Tres}(obj::Rho1E2U3BC,
              q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},
              x::AbstractArray{Tmsh,1}, dxidx::AbstractArray{Tmsh,2},
              nrm::AbstractArray{Tmsh,1}, bndryflux::AbstractArray{Tres, 1},
              params::ParamType{2})



  #println("in Rho1E2U3Bc")
  qg = params.qg

  calcRho1Energy2U3(x, params, qg)

  #println("qg = ", qg)
  # call Roe solver
  RoeSolver(params, q, qg, aux_vars, dxidx, nrm, bndryflux)

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
*  `nrm`      : SBP face normal
*  `bndryflux` : Computed flux value at the boundary

"""
type FreeStreamBC <: BCType
end

function call{Tmsh, Tsol, Tres}(obj::FreeStreamBC, q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},  x::AbstractArray{Tmsh,1},
              dxidx::AbstractArray{Tmsh,2}, nrm::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1}, params::ParamType)

  qg = params.qg

  calcFreeStream(x, params, qg)
  RoeSolver(params, q, qg, aux_vars, dxidx, nrm, bndryflux)

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
* `dxidx` : Mapping jacobian matrix for the SBP node
* `nrm`   : sbpface normal vector
* `bndryflux_bar` : Input flux value seed that is used to compute the reverse
                    mode derivative.
* `params`        : equation object parameters

**Output**

* `dxidx_bar` : Derivative of bndryflux_bar w.r.t the mapping jacobian

"""->

type FreeStreamBC_revm <: BCType_revm
end

function call{Tmsh, Tsol, Tres}(obj::FreeStreamBC_revm, q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},  x::AbstractArray{Tmsh,1},
              dxidx::AbstractArray{Tmsh,2}, dxidx_bar::AbstractArray{Tmsh, 2},
              nrm::AbstractArray{Tmsh,1}, bndryflux_bar::AbstractArray{Tres, 1},
              params::ParamType)

  # Forward sweep
  qg = params.qg
  calcFreeStream(x, params, qg)

  # Reverse sweep
  RoeSolver_revm(params, q, qg, aux_vars, dxidx, nrm, bndryflux_bar, dxidx_bar)

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
*  `dxidx`    : Mapping jacobian matrix for the SBP node
*  `nrm`      : SBP face normal
*  `bndryflux` : Computed flux value at the boundary

"""->
type FreeStreamBC_dAlpha <: BCType
end

function call{Tmsh, Tsol, Tres}(obj::FreeStreamBC_dAlpha, q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},  x::AbstractArray{Tmsh,1},
              dxidx::AbstractArray{Tmsh,2}, nrm::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1}, params::ParamType{2})

  qg = params.qg

  calcFreeStream_dAlpha(x, params, qg)
  RoeSolver(params, q, qg, aux_vars, dxidx, nrm, bndryflux)

  return nothing
end


@doc """
### EulerEquationMod.allOnesBC <: BCTypes

  This functor uses the Roe solver to calculate the flux for a boundary
  state where all the conservative variables have a value 1.0

  This is a low level functor
"""->

type allOnesBC <: BCType
end

function call{Tmsh, Tsol, Tres}(obj::allOnesBC, q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, x::AbstractArray{Tmsh,1},
              dxidx::AbstractArray{Tmsh,2}, nrm::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1}, params::ParamType{2})

  qg = zeros(Tsol, 4)
  calcOnes(x, params, qg)

  RoeSolver(params, q, qg, aux_vars, dxidx, nrm, bndryflux)

  # println("bndryflux = ", bndryflux)
  return nothing
end # end function call

@doc """
### EulerEquationMod.allZerosBC <: BCTypes

  This functor uses the Roe solver to calculate the flux for a boundary
  state where all the conservative variables have a value 0.0

  This is a low level functor
"""->

type allZerosBC <: BCType
end

function call{Tmsh, Tsol, Tres}(obj::allZerosBC, q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, x::AbstractArray{Tmsh,1},
              dxidx::AbstractArray{Tmsh,2}, nrm::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1}, params::ParamType{2})

  qg = zeros(Tsol, 4)
  calcZeros(x, params, qg)

  RoeSolver(params, q, qg, aux_vars, dxidx, nrm, bndryflux)

  # println("bndryflux = ", bndryflux)
  return nothing
end # end function call

type ExpBC <: BCType
end

function call{Tmsh, Tsol, Tres}(obj::ExpBC, q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, coords::AbstractArray{Tmsh,1},
              dxidx::AbstractArray{Tmsh,2}, nrm::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1}, params::ParamType)

  qg = params.qg
  calcExp(coords, params, qg)
  RoeSolver(params, q, qg, aux_vars, dxidx, nrm, bndryflux)

  # println("bndryflux = ", bndryflux)
  return nothing
end # end function call

type ExpBC_revm <: BCType_revm
end

function call{Tmsh, Tsol, Tres}(obj::ExpBC_revm, q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},  x::AbstractArray{Tmsh,1},
              dxidx::AbstractArray{Tmsh,2}, dxidx_bar::AbstractArray{Tmsh, 2},
              nrm::AbstractArray{Tmsh,1}, bndryflux_bar::AbstractArray{Tres, 1},
              params::ParamType)

  # Forward Sweep
  qg = params.qg
  calcExp(x, params, qg)
  # RoeSolver(params, q, qg, aux_vars, dxidx, nrm, bndryflux)

  # Reverse Sweep
  RoeSolver_revm(params, q, qg, aux_vars, dxidx, nrm, bndryflux_bar, dxidx_bar)

  return nothing
end

type PeriodicMMSBC <: BCType
end

function call{Tmsh, Tsol, Tres}(obj::PeriodicMMSBC, q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, coords::AbstractArray{Tmsh,1},
              dxidx::AbstractArray{Tmsh,2}, nrm::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 1}, params::ParamType)
# use the exact solution as the boundary condition for the PeriodicMMS
# solutions

  qg = params.qg
  calcPeriodicMMS(coords, params, qg)
  use_efix = 0
  RoeSolver(params, q, qg, aux_vars, dxidx, nrm, bndryflux, use_efix)

  # println("bndryflux = ", bndryflux)
  return nothing
end # end function call



# every time a new boundary condition is created,
# add it to the dictionary
#const isentropicVortexBC_ = isentropicVortexBC()
#const noPenetrationBC_ = noPenetrationBC()
global const BCDict = Dict{ASCIIString, BCType}(
"isentropicVortexBC" => isentropicVortexBC(),
"noPenetrationBC" => noPenetrationBC(),
"Rho1E2U3BC" => Rho1E2U3BC(),
"isentropicVortexBC_physical" => isentropicVortexBC_physical(),
"FreeStreamBC" => FreeStreamBC(),
"FreeStreamBC_dAlpha" => FreeStreamBC_dAlpha(),
"allOnesBC" => allOnesBC(),
"unsteadyVortexBC" => unsteadyVortexBC(),
"ExpBC" => ExpBC(),
"PeriodicMMSBC" => PeriodicMMSBC(),
)

@doc """
### EulerEquationMod.getBCFunctors

  This function uses the opts dictionary to populate mesh.bndry_funcs with
  the the functors

  This is a high level function.
"""->
# use this function to populate access the needed values in BCDict
function getBCFunctors(mesh::AbstractMesh, sbp::AbstractSBP, eqn::EulerData, opts)
# populate the array mesh.bndry_funcs with the functors for the boundary condition types

#  println("Entered getBCFunctors")

  for i=1:mesh.numBC
    key_i = string("BC", i, "_name")
    val = opts[key_i]
    println("BCDict[val] = ", BCDict[val])
    mesh.bndry_funcs[i] = BCDict[val]
  end

  return nothing
end # ENd function getBCFunctors

global const BCDict_revm = Dict{ASCIIString, BCType_revm}(
"noPenetrationBC" => noPenetrationBC_revm(),
"FreeStreamBC" => FreeStreamBC_revm(),
"ExpBC" => ExpBC_revm(),
"isentropicVortexBC" => isentropicVortexBC_revm(),
)

function getBCFunctors_revm(mesh::AbstractMesh, sbp::AbstractSBP, eqn::EulerData, opts)

  for i = 1:mesh.numBC
    key_i = string("BC", i, "_name")
    val = opts[key_i]
    println("BCDict_revm[$val] = ", BCDict_revm[val])
    mesh.bndry_funcs_revm[i] = BCDict_revm[val]
  end # End for i = 1:mesh.numBC

  return nothing
end # End function getBCFunctors_revm
