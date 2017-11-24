# bc.jl

@doc """

  This function calculates the boundary flux for the portion of the boundary
  with a particular boundary condition.  The eqn.q are converted to 
  conservative variables if needed.  For the DG version, eqn.q_bndry must
  already be populated with the q variables interpolated to the boundary

  Inputs:

  mesh : AbstractMesh
  sbp : AbstractSBP
  eqn : AdvectionEquation
  functor : a callable object that calculates the boundary flux at a node
  idx_range: the Range describing which Boundaries have the current BC
  bndry_facenums:  An array with elements of type Boundary that tell which
                   element faces have the boundary condition
  Outputs:

  bndryflux : the array to store the boundary flux, corresponds to 
              bndry_facenums

  note that bndry_facenums and bndryflux must be only the portion of the 
  their parent arrays that correspond to the Boundaries that have the 
  current boundary condition applied.

  The functor must have the signature:
  functor( q, aux_vars, x, dxidx, nrm, bndryflux_i, eqn.params)
  where q are the *conservative* variables.
  where all arguments (except params and nrm) are vectors of values at a node.

  params is the ParamType associated with the the EulerEquation object
  nrm = mesh.sbpface.normal[:, current_node]

  This is a mid level function.
"""->
# mid level function
function calcBoundaryFlux{Tmsh,  Tsol, Tres}( mesh::AbstractCGMesh{Tmsh}, 
                          sbp::AbstractSBP, eqn::AdvectionData{Tsol}, 
                          functor::BCType, idx_range::UnitRange,
                          bndry_facenums::AbstractArray{Boundary,1}, 
                          bndryflux::AbstractArray{Tres, 3})
# calculate the boundary flux for the boundary condition evaluated by the functor

#  println("enterted calcBoundaryFlux CG")
#TODO: update this to use mesh.nrm_bndry

  t = eqn.t
  nfaces = length(bndry_facenums)
  nrm_scaled = zeros(Tmsh, mesh.dim)
  for i=1:nfaces  # loop over faces with this BC
    bndry_i = bndry_facenums[i]
    for j = 1:mesh.numNodesPerFace
      k = mesh.facenodes[j, bndry_i.face]

      # get components
      q = eqn.q[ 1, k, bndry_i.element]
      alpha_x = eqn.params.alpha_x
      alpha_y = eqn.params.alpha_y
      coords = ro_sview(mesh.coords, :, k, bndry_i.element)
      dxidx = ro_sview(mesh.dxidx, :, :, k, bndry_i.element)
      nrm = ro_sview(mesh.sbpface.normal, :, bndry_i.face)
      calcBCNormal(eqn.params, dxidx, nrm, nrm_scaled)
      bndryflux[1, j, i] = -functor(eqn.params, q, coords, nrm_scaled, t)
    end
  end

  return nothing
end


# DG version
function calcBoundaryFlux{Tmsh,  Tsol, Tres}( mesh::AbstractDGMesh{Tmsh}, 
                          sbp::AbstractSBP, eqn::AdvectionData{Tsol}, 
                          functor::BCType, idx_range::UnitRange,
                          bndry_facenums::AbstractArray{Boundary,1}, 
                          bndryflux::AbstractArray{Tres, 3})

  
# calculate the boundary flux for the boundary condition evaluated by the functor

#  println("entered calcBoundaryFlux DG")
  t = eqn.t
  println("t = ", t)
  nfaces = length(bndry_facenums)
  for i=1:nfaces  # loop over faces with this BC
    println("interface ", i)
    bndry_i = bndry_facenums[i]
    global_facenum = idx_range[i]
    for j = 1:mesh.numNodesPerFace
      println("face node ", j)

      # get components
      q = eqn.q_bndry[ 1, j, global_facenum]
      alpha_x = eqn.params.alpha_x
      alpha_y = eqn.params.alpha_y
      coords = ro_sview(mesh.coords_bndry, :, j, global_facenum)
      nrm_scaled = ro_sview(mesh.nrm_bndry, :, j, global_facenum)
      println("q = ", q)
      println("alpha_x = ", alpha_x)
      println("alpha_y = ", alpha_y)
      println("coords = ", coords)
      println("nrm_scaled = ", nrm_scaled)
      bndryflux[1, j, i] = -functor(eqn.params, q, coords, nrm_scaled, t)
      println("flux = ", bndryflux[1, j, i])
    end
  end

  return nothing
end

"""
  This function computes the boundary integrals (and should probably be renamed)
  without using eqn.q_bndry of eqn.bndryflux.  eqn.res is updated with
  the results.

  See calcBoundaryFlux for the meaning of the arguments
"""
function calcBoundaryFlux_nopre{Tmsh,  Tsol, Tres}( mesh::AbstractDGMesh{Tmsh}, 
                          sbp::AbstractSBP, eqn::AdvectionData{Tsol}, 
                          functor::BCType, idx_range::UnitRange,
                          bndry_facenums::AbstractArray{Boundary,1}, 
                          bndryflux::AbstractArray{Tres, 3})

  t = eqn.t
  nfaces = length(bndry_facenums)
  q_face = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
  flux_face = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)

  for i=1:nfaces  # loop over faces with this BC
    bndry_i = bndry_facenums[i]
    global_facenum = idx_range[i]

    # interpolate to face
    q_el = ro_sview(eqn.q, :, :, bndry_i.element)
    boundaryFaceInterpolate!(mesh.sbpface, bndry_i.face, q_el, q_face)

    for j = 1:mesh.numNodesPerFace
      # get components
      q = q_face[ 1, j]
      coords = ro_sview(mesh.coords_bndry, :, j, global_facenum)
      nrm_scaled = ro_sview(mesh.nrm_bndry, :, j, i)
      flux_face[1, j] = -functor(eqn.params, q, coords, nrm_scaled, t)
    end

    res_i = sview(eqn.res, :, :, bndry_i.element)
    boundaryFaceIntegrate!(mesh.sbpface, bndry_i.face, flux_face, res_i)
  end

  return nothing
end


@doc """
### AdvectionEquationMod.x5plusy5BC

Calculates q at the boundary which is equal to x^5 + y^5. It is a nodal 
level function.

**Inputs**

*  `u` : Advection variable (eqn.q)
*  `params`: the equation ParamType
*  `coords` : Nodal coordinates
*  `nrm_scaled`    : scaled face normal vector in x-y space
*  `t`:  current time value

**Outputs**
*  `bndryflux` : Flux at the boundary

*  None

"""->

type x5plusy5BC <: BCType
end

function call{Tmsh, Tsol}(obj::x5plusy5BC, params::ParamType2, u::Tsol, 
              coords::AbstractArray{Tmsh,1}, 
              nrm_scaled::AbstractArray{Tmsh,1}, t)

  u_bc = calc_x5plusy5(params, coords, t) # Calculate the actual analytic value of u at the bondary
  bndryflux = RoeSolver(params, u, u_bc, nrm_scaled)

  return bndryflux
end

type constantBC <: BCType
end

function call{Tmsh, Tsol}(obj::constantBC, params::ParamTypes, u::Tsol, 
              coords::AbstractArray{Tmsh,1}, 
              nrm_scaled::AbstractArray{Tmsh,1}, t)

  u_bc = 2
  bndryflux = RoeSolver(params, u, u_bc, nrm_scaled)
  return bndryflux
end



@doc """
### AdvectionEquationMod.exp_xplusyBC

Calculates q at the boundary which is equal to exp(x+y). It is a nodal 
level function.

**Inputs**

*  `u` : Advection variable (eqn.q)
*  `alpha_x` & `alpha_y` : velocities in the X & Y directions
*  `coords` : Nodal coordinates
*  `dxidx`  : Mapping Jacobian
*  `nrm`    : SBP face-normal vectors
*  `bndryflux` : Flux at the boundary

**Outputs**

*  None

"""->
type exp_xplusyBC <: BCType
end

function call{Tmsh, Tsol}(obj::exp_xplusyBC, params::ParamType2, u::Tsol, 
              coords::AbstractArray{Tmsh,1}, 
              nrm_scaled::AbstractArray{Tmsh,1}, t)

  u_bc = calc_exp_xplusy(params, coords, t)
  bndryflux = RoeSolver(params, u, u_bc, nrm_scaled)

  return bndryflux
end

@doc """
### AdvectionEquationMod.sinwave_BC

  Uses the Roe solver to calculate the boundary flux using calc_sinewave to
  get the boundary state
"""->
type sinwave_BC <: BCType
end

function call{Tmsh, Tsol}(obj::sinwave_BC, params::ParamType2, u::Tsol,
              coords::AbstractArray{Tmsh,1},
              nrm_scaled::AbstractArray{Tmsh,1}, t)

  u_bc = calc_sinwave(params, coords, t)
#  println("  u_bc = ", u_bc)
  bndryflux = RoeSolver(params, u, u_bc, nrm_scaled)

  return bndryflux
end

@doc """
### AdvectionEquationMod.sinwavey_BC

  Uses the Roe solver to calculate the boundary flux using calc_sinewavey to
  get the boundary state
"""->
type sinwavey_BC <: BCType
end

function call{Tmsh, Tsol}(obj::sinwavey_BC, params::ParamType2, u::Tsol,
              coords::AbstractArray{Tmsh,1},
              nrm_scaled::AbstractArray{Tmsh,1}, t)

  u_bc = calc_sinwavey(params, coords, t)
#  println("  u_bc = ", u_bc)
  bndryflux = RoeSolver(params, u, u_bc, nrm_scaled)

  return bndryflux
end

@doc """
### AdvectionEquationMod.sinwavey_pertBC

  Uses the Roe solver to calculate the boundary flux using calc_sinewave_pert to
  get the boundary state
"""->

type sinwavey_pertBC <: BCType
end

function call{Tmsh, Tsol}(obj::sinwavey_pertBC, params::ParamType2, u::Tsol,
              coords::AbstractArray{Tmsh,1},
              nrm_scaled::AbstractArray{Tmsh,1}, t)

  u_bc = calc_sinwavey_pert(params, coords, t)
#  println("  u_bc = ", u_bc)
  bndryflux = RoeSolver(params, u, u_bc, nrm_scaled)

  return bndryflux
end

@doc """
### AdvectionEquationMod.sinwave_ampl_BC

  Uses the Roe solver to calculate the boundary flux using calc_sinewave_ampl to
  get the boundary state
"""->
type sinwave_ampl_BC <: BCType
end

function call{Tmsh, Tsol}(obj::sinwave_ampl_BC, params::ParamType2, u::Tsol, 
              coords::AbstractArray{Tmsh,1}, nrm_scaled::AbstractArray{Tmsh,1}, t)

  u_bc = calc_sinwave_ampl(params, coords, t)
  bndryflux = RoeSolver(params, u, u_bc, nrm_scaled)

  return bndryflux
end

@doc """
### AdvectionEquationMod.mms1BC

  Uses the Roe solver to calculate the boundary flux using calc_mms1 to get
  the boundary state.
"""->
type mms1BC <: BCType
end

function call{Tmsh, Tsol}(obj::mms1BC, params::ParamType2, u::Tsol,
              coords::AbstractArray{Tmsh,1},
              nrm_scaled::AbstractArray{Tmsh,1}, t)

  u_bc = calc_mms1(params, coords, t)
#  println("  u_bc = ", u_bc)
  bndryflux = RoeSolver(params, u, u_bc, nrm_scaled)

  return bndryflux
end

@doc """
### AdvectionEquationMod.x4BC

  Uses the Roe solver to calculate the boundary flux using calc_x4 to
  get the boundary state.
"""->
type x4BC <: BCType
end

function call{Tmsh, Tsol}(obj::x4BC, params::ParamType2, u::Tsol,
              coords::AbstractArray{Tmsh,1},
              nrm_scaled::AbstractArray{Tmsh,1}, t)

  u_bc = calc_x4(params, coords, t)
#  println("  u_bc = ", u_bc)
  bndryflux = RoeSolver(params, u, u_bc, nrm_scaled)

  return bndryflux
end

@doc """
### AdvectionEquationMod.p0BC

  Uses the Roe solver to calculate the boundary flux using calc_p0 to
  get the boundary state
"""->

type p0BC <: BCType
end

function call{Tmsh, Tsol}(obj::p0BC, params::ParamTypes, u::Tsol,
              coords::AbstractArray{Tmsh,1},
              nrm_scaled::AbstractArray{Tmsh,1}, t)

  u_bc = calc_p0(params, coords, t)
#  println("  u_bc = ", u_bc)
  bndryflux = RoeSolver(params, u, u_bc, nrm_scaled)
  return bndryflux
end


@doc """
### AdvectionEquationMod.p1BC

  Uses the Roe solver to calculate the boundary flux using calc_p1 to
  get the boundary state
"""->

type p1BC <: BCType
end

function call{Tmsh, Tsol}(obj::p1BC, params::ParamTypes, u::Tsol,
              coords::AbstractArray{Tmsh,1},
              nrm_scaled::AbstractArray{Tmsh,1}, t)

  u_bc = calc_p1(params, coords, t)
#  println("  u_bc = ", u_bc)
  bndryflux = RoeSolver(params, u, u_bc, nrm_scaled)
  return bndryflux
end

@doc """
### AdvectionEquationMod.p2BC

  Uses the Roe solver to calculate the boundary flux using calc_p2 to
  get the boundary state
"""->
type p2BC <: BCType
end

function call{Tmsh, Tsol}(obj::p2BC, params::ParamTypes, u::Tsol, 
              coords::AbstractArray{Tmsh,1},
              nrm_scaled::AbstractArray{Tmsh,1}, t)

  u_bc = calc_p2(params, coords, t)
#  println("  u_bc = ", u_bc)
  bndryflux = RoeSolver(params, u, u_bc, nrm_scaled)

  return bndryflux
end

@doc """
### AdvectionEquationMod.p3BC

  Uses the Roe solver to calculate the boundary flux using calc_p3 to
  get the boundary state
"""->
type p3BC <: BCType
end

function call{Tmsh, Tsol}(obj::p3BC, params::ParamTypes, u::Tsol,
              coords::AbstractArray{Tmsh,1},
              nrm_scaled::AbstractArray{Tmsh,1}, t)

  u_bc = calc_p3(params, coords, t)
#  println("  u_bc = ", u_bc)
  bndryflux = RoeSolver(params, u, u_bc, nrm_scaled)

  return bndryflux
end


@doc """
### AdvectionEquationMod.p4BC

  Uses the Roe solver to calculate the boundary flux using calc_p4 to
  get the boundary state.
"""->
type p4BC <: BCType
end

function call{Tmsh, Tsol}(obj::p4BC, params::ParamTypes, u::Tsol,
              coords::AbstractArray{Tmsh,1},
              nrm_scaled::AbstractArray{Tmsh,1}, t)

  u_bc = calc_p4(params, coords, t)
#  println("  u_bc = ", u_bc)
  bndryflux = RoeSolver(params, u, u_bc, nrm_scaled)

  return bndryflux
end

@doc """
### AdvectionEquationMod.p5BC

  Uses the Roe solver to calculate the boundary flux using calc_p5 to
  get the boundary state.
"""->

type p5BC <: BCType
end

function call{Tmsh, Tsol}(obj::p5BC, params::ParamTypes, u::Tsol,
              coords::AbstractArray{Tmsh,1},
              nrm_scaled::AbstractArray{Tmsh,1}, t)
  # this is really slow: use Horner's rule!
  u_bc = calc_p5(params, coords, t)
#  println("calc_p5 @time printed above")
#  println("  u_bc = ", u_bc)
  bndryflux = RoeSolver(params, u, u_bc, nrm_scaled)
#  println("    RoeSolver @time printed above")

  return bndryflux
end

@doc """
### AdvectionEquationMod.exp5xplus4yplus2BC

Uses the Roe solver to calculate the boundary flux using calc_exp5xplus4yplus2 
to get the boundary state.

"""->

type exp5xplus4yplus2BC <: BCType
end

function call{Tmsh, Tsol}(obj::exp5xplus4yplus2BC, params::ParamType2, u::Tsol,
              coords::AbstractArray{Tmsh,1},
              nrm_scaled::AbstractArray{Tmsh,1}, t)

  u_bc = calc_exp5xplus4yplus2(params, coords, t)
  bndryflux = RoeSolver(params, u, u_bc, nrm_scaled)

  return bndryflux
end

@doc """
### AdvectionEquationMod.exp5xplusyBC

Uses the Roe solver to calculate the boundary flux using calc_exp5xplusy to get
the boundary state.

"""->

type exp5xplusyBC <:BCType
end

function call{Tmsh, Tsol}(obj::exp5xplusyBC, params::ParamType2, u::Tsol,
              coords::AbstractArray{Tmsh,1},
              nrm_scaled::AbstractArray{Tmsh,1}, t)

  u_bc = calc_exp5xplusy(params, coords, t)
  bndryflux = RoeSolver(params, u, u_bc, nrm_scaled)

  return bndryflux
end

@doc """
### AdvectionEquationMod.exp3xplusyBC

Uses the Roe solver to calculate the boundary flux using calc_exp3xplusy to get
the boundary state.

"""->

type exp3xplusyBC <:BCType
end

function call{Tmsh, Tsol}(obj::exp3xplusyBC, params::ParamType2, u::Tsol,
              coords::AbstractArray{Tmsh,1},
              nrm_scaled::AbstractArray{Tmsh,1}, t)

  u_bc = calc_exp3xplusy(params, coords, t)
  bndryflux = RoeSolver(params, u, u_bc, nrm_scaled)

  return bndryflux
end

@doc """
### AdvectionEquationMod.exp2xplus2yBC

Uses the Roe solver to calculate the boundary flux using calc_exp2xplus2y to get
the boundary state.

"""->

type exp2xplus2yBC <:BCType
end

function call{Tmsh, Tsol}(obj::exp2xplus2yBC, params::ParamType2,  u::Tsol,
              coords::AbstractArray{Tmsh,1},
              nrm_scaled::AbstractArray{Tmsh,1}, t)

  u_bc = calc_exp2xplus2y(params, coords, t)
  bndryflux = RoeSolver(params, u, u_bc, nrm_scaled)

  return bndryflux
end

@doc """
### AdvectionEquationMod.exp_xyBC

Uses the Roe solver to calculate the boundary flux using calc_exp_xy to get the
boundary state

"""->

type exp_xyBC <: BCType
end

function call{Tmsh, Tsol}(obj::exp_xyBC, params::ParamType2, u::Tsol,
              coords::AbstractArray{Tmsh,1},
              nrm_scaled::AbstractArray{Tmsh,1}, t)

  u_bc = calc_exp_xy(params, coords, t)
  bndryflux = RoeSolver(params, u, u_bc, nrm_scaled)

  return bndryflux
end

@doc """
### AdvectionEquationMod.xplusyBC

Uses Roe solver to calculate the boundary flux using calc_xplusy to get the 
boundary state

"""->

type xplusyBC <: BCType
end

function call{Tmsh, Tsol}(obj::xplusyBC, params::ParamType2, u::Tsol,
              coords::AbstractArray{Tmsh,1},
              nrm_scaled::AbstractArray{Tmsh,1}, t)

  u_bc = calc_xplusy(params, coords, t)
  bndryflux = RoeSolver(params, u, u_bc, nrm_scaled)

  return bndryflux
end


"""
  BC for unsteadymms
"""
type unsteadymmsBC <: BCType
end

function call{Tmsh, Tsol}(obj::unsteadymmsBC, params::ParamType, u::Tsol,
              coords::AbstractArray{Tmsh,1},
              nrm_scaled::AbstractArray{Tmsh,1}, t)

  u_bc = calc_unsteadymms(params, coords, t)
  bndryflux = RoeSolver(params, u, u_bc, nrm_scaled)

  return bndryflux
end

"""
  BC for unsteadypoly
"""
type unsteadypolyBC <: BCType
end

function call{Tmsh, Tsol}(obj::unsteadypolyBC, params::ParamType, u::Tsol,
              coords::AbstractArray{Tmsh,1},
              nrm_scaled::AbstractArray{Tmsh,1}, t)

  u_bc = calc_unsteadypoly(params, coords, t)
  bndryflux = RoeSolver(params, u, u_bc, nrm_scaled)

  return bndryflux
end

"""
  Default BC to calculate the boundary face integral (no numerical flux
  functions)
"""
type defaultBC <: BCType
end

function call{Tmsh, Tsol, Tdim}(obj::defaultBC, params::ParamType{Tdim},
              u::Tsol, coords::AbstractArray{Tmsh,1},
              nrm_scaled::AbstractArray{Tmsh,1}, t)

  alpha_nrm = params.alpha_x*nrm_scaled[1] + params.alpha_y*nrm_scaled[2]
  if Tdim == 3  # static analysis
    alpha_nrm += params.alpha_z*nrm_scaled[3]
  end

  return alpha_nrm*u
end


@doc """
### AdvectionEquationMod.BCDict

It stores all the possible boundary condition dictionary options. Whenever a 
new boundary condition is created, it should get added to BCDict.

"""->
global const BCDict = Dict{ASCIIString, BCType}(
"constantBC" => constantBC(),
"x5plusy5BC" => x5plusy5BC(),
"exp_xplusyBC" => exp_xplusyBC(),
"sinwaveBC" => sinwave_BC(),
"sinwaveyBC" => sinwavey_BC(),
"sinwavey_pertBC" => sinwavey_pertBC(),
"sinwaveamplBC" => sinwave_ampl_BC(),
"mms1BC" => mms1BC(),
"x4BC" => x4BC(),
"p0BC" => p0BC(),
"p1BC" => p1BC(),
"p2BC" => p2BC(),
"p3BC" => p3BC(),
"p4BC" => p4BC(),
"p5BC" => p5BC(),
"exp5xplus4yplus2BC" => exp5xplus4yplus2BC(),
"exp5xplusyBC" => exp5xplusyBC(),
"exp3xplusyBC" => exp3xplusyBC(),
"exp2xplus2yBC" => exp2xplus2yBC(),
"exp_xyBC" => exp_xyBC(),
"xplusyBC" => xplusyBC(),
"unsteadymmsBC" => unsteadymmsBC(),
"unsteadypolyBC" => unsteadypolyBC(),
"defaultBC" => defaultBC(),
)


@doc """
### AdvectionEquationMod.getBCFunctors

This function uses the opts dictionary to populate mesh.bndry_funcs with
the the functors

This is a high level function.

**Inputs**

*  `mesh` : Abstract mesh type
*  `sbp`  : Summation-by-parts operator
*  `eqn`  : Advection equation object
*  `opts` : Input dictionary options

**Outputs**

*  None

"""->
# use this function to populate access the needed values in BCDict
function getBCFunctors{Tmsh, Tsol, Tdim}(mesh::AbstractMesh{Tmsh}, sbp::AbstractSBP, 
                       eqn::AdvectionData{Tsol, Tdim}, opts)
# populate the array mesh.bndry_funcs with the functors for the boundary condition types

  println("Entered getBCFunctors")

  for i=1:mesh.numBC
    key_i = string("BC", i, "_name")
    val = opts[key_i]
    println("BCDict[val] = ", BCDict[val])
    mesh.bndry_funcs[i] = BCDict[val]
  end

  return nothing
end
