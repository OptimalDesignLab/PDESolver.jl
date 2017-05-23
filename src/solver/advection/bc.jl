# bc.jl

@doc """
### AdvectionEquationMod.calcBoundaryFlux

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
  
  t = eqn.t
  nfaces = length(bndry_facenums)
  for i=1:nfaces  # loop over faces with this BC
    bndry_i = bndry_facenums[i]
    for j = 1:sbp.numfacenodes
      k = sbp.facenodes[j, bndry_i.face]

      # get components
      q = eqn.q[ 1, k, bndry_i.element]
      alpha_x = eqn.params.alpha_x
      alpha_y = eqn.params.alpha_y
      coords = sview(mesh.coords, :, k, bndry_i.element)
      dxidx = sview(mesh.dxidx, :, :, k, bndry_i.element)
      nrm = sview(mesh.sbpface.normal, :, bndry_i.face)
      bndryflux[1, j, i] = -functor(q, eqn.params, coords, dxidx, nrm, t)
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
  nfaces = length(bndry_facenums)
  for i=1:nfaces  # loop over faces with this BC
    bndry_i = bndry_facenums[i]
    global_facenum = idx_range[i]
    for j = 1:mesh.numNodesPerFace

      # get components
      q = eqn.q_bndry[ 1, j, global_facenum]
      alpha_x = eqn.params.alpha_x
      alpha_y = eqn.params.alpha_y
      coords = sview(mesh.coords_bndry, :, j, global_facenum)
      dxidx = sview(mesh.dxidx_bndry, :, :, j, global_facenum)
      nrm = sview(mesh.sbpface.normal, :, bndry_i.face)
      bndryflux[1, j, i] = -functor(q, eqn.params, coords, dxidx, nrm, t)
    end
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
*  `dxidx`  : Mapping Jacobian
*  `nrm`    : SBP face-normal vectors
*  `bndryflux` : Flux at the boundary

**Outputs**

*  None

"""->

type x5plusy5BC <: BCType
end

function call{Tmsh, Tsol}(obj::x5plusy5BC, u::Tsol, 
              params::ParamType2, coords::AbstractArray{Tmsh,1}, 
              dxidx::AbstractArray{Tmsh,2}, nrm::AbstractArray{Tmsh,1}, t)

  u_bc = calc_x5plusy5(coords, params, t) # Calculate the actual analytic value of u at the bondary
  bndryflux = RoeSolver(u, u_bc, params, nrm, dxidx)

  return bndryflux
end

type constantBC <: BCType
end

function call{Tmsh, Tsol}(obj::constantBC, u::Tsol, 
              params::ParamTypes, coords::AbstractArray{Tmsh,1}, 
              dxidx::AbstractArray{Tmsh,2}, nrm::AbstractArray{Tmsh,1}, t)

  u_bc = 2
  bndryflux = RoeSolver(u, u_bc, params, nrm, dxidx)
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

function call{Tmsh, Tsol}(obj::exp_xplusyBC, u::Tsol, 
              params::ParamType2, coords::AbstractArray{Tmsh,1}, 
              dxidx::AbstractArray{Tmsh,2}, nrm::AbstractArray{Tmsh,1}, t)

  u_bc = calc_exp_xplusy(coords, params, t)
  bndryflux = RoeSolver(u, u_bc, params, nrm, dxidx)

  return bndryflux
end

@doc """
### AdvectionEquationMod.sinwave_BC

  Uses the Roe solver to calculate the boundary flux using calc_sinewave to
  get the boundary state
"""->
type sinwave_BC <: BCType
end

function call{Tmsh, Tsol}(obj::sinwave_BC, u::Tsol, params::ParamType2,
              coords::AbstractArray{Tmsh,1}, dxidx::AbstractArray{Tmsh, 2},
              nrm::AbstractArray{Tmsh,1}, t)

  u_bc = calc_sinwave(coords, params, t)
#  println("  u_bc = ", u_bc)
  bndryflux = RoeSolver(u, u_bc, params, nrm, dxidx)

  return bndryflux
end

@doc """
### AdvectionEquationMod.sinwavey_BC

  Uses the Roe solver to calculate the boundary flux using calc_sinewavey to
  get the boundary state
"""->
type sinwavey_BC <: BCType
end

function call{Tmsh, Tsol}(obj::sinwavey_BC, u::Tsol, params::ParamType2,
              coords::AbstractArray{Tmsh,1}, dxidx::AbstractArray{Tmsh, 2},
              nrm::AbstractArray{Tmsh,1}, t)

  u_bc = calc_sinwavey(coords, params, t)
#  println("  u_bc = ", u_bc)
  bndryflux = RoeSolver(u, u_bc, params, nrm, dxidx)

  return bndryflux
end

@doc """
### AdvectionEquationMod.sinwavey_pertBC

  Uses the Roe solver to calculate the boundary flux using calc_sinewave_pert to
  get the boundary state
"""->

type sinwavey_pertBC <: BCType
end

function call{Tmsh, Tsol}(obj::sinwavey_pertBC, u::Tsol, params::ParamType2,
              coords::AbstractArray{Tmsh,1}, dxidx::AbstractArray{Tmsh, 2},
              nrm::AbstractArray{Tmsh,1}, t)

  u_bc = calc_sinwavey_pert(coords, params, t)
#  println("  u_bc = ", u_bc)
  bndryflux = RoeSolver(u, u_bc, params, nrm, dxidx)

  return bndryflux
end

@doc """
### AdvectionEquationMod.mms1BC

  Uses the Roe solver to calculate the boundary flux using calc_mms1 to get
  the boundary state.
"""->
type mms1BC <: BCType
end

function call{Tmsh, Tsol}(obj::mms1BC, u::Tsol, params::ParamType2,
              coords::AbstractArray{Tmsh,1}, dxidx::AbstractArray{Tmsh, 2},
              nrm::AbstractArray{Tmsh,1}, t)

  u_bc = calc_mms1(coords, params, t)
#  println("  u_bc = ", u_bc)
  bndryflux = RoeSolver(u, u_bc, params, nrm, dxidx)

  return bndryflux
end

@doc """
### AdvectionEquationMod.x4BC

  Uses the Roe solver to calculate the boundary flux using calc_x4 to
  get the boundary state.
"""->
type x4BC <: BCType
end

function call{Tmsh, Tsol}(obj::x4BC, u::Tsol, params::ParamType2,
              coords::AbstractArray{Tmsh,1}, dxidx::AbstractArray{Tmsh, 2},
              nrm::AbstractArray{Tmsh,1}, t)

  u_bc = calc_x4(coords, params, t)
#  println("  u_bc = ", u_bc)
  bndryflux = RoeSolver(u, u_bc, params, nrm, dxidx)

  return bndryflux
end

@doc """
### AdvectionEquationMod.p1BC

  Uses the Roe solver to calculate the boundary flux using calc_p1 to
  get the boundary state
"""->

type p1BC <: BCType
end

function call{Tmsh, Tsol}(obj::p1BC, u::Tsol, params::ParamTypes,
              coords::AbstractArray{Tmsh,1}, dxidx::AbstractArray{Tmsh, 2},
              nrm::AbstractArray{Tmsh,1}, t)

  u_bc = calc_p1(coords, params, t)
#  println("  u_bc = ", u_bc)
  bndryflux = RoeSolver(u, u_bc, params, nrm, dxidx)
  return bndryflux
end

@doc """
### AdvectionEquationMod.p2BC

  Uses the Roe solver to calculate the boundary flux using calc_p2 to
  get the boundary state
"""->
type p2BC <: BCType
end

function call{Tmsh, Tsol}(obj::p2BC, u::Tsol, params::ParamTypes,
              coords::AbstractArray{Tmsh,1}, dxidx::AbstractArray{Tmsh, 2},
              nrm::AbstractArray{Tmsh,1}, t)

  u_bc = calc_p2(coords, params, t)
#  println("  u_bc = ", u_bc)
  bndryflux = RoeSolver(u, u_bc, params, nrm, dxidx)

  return bndryflux
end

@doc """
### AdvectionEquationMod.p3BC

  Uses the Roe solver to calculate the boundary flux using calc_p3 to
  get the boundary state
"""->
type p3BC <: BCType
end

function call{Tmsh, Tsol}(obj::p3BC, u::Tsol, params::ParamTypes,
              coords::AbstractArray{Tmsh,1}, dxidx::AbstractArray{Tmsh, 2},
              nrm::AbstractArray{Tmsh,1}, t)

  u_bc = calc_p3(coords, params, t)
#  println("  u_bc = ", u_bc)
  bndryflux = RoeSolver(u, u_bc, params, nrm, dxidx)

  return bndryflux
end


@doc """
### AdvectionEquationMod.p4BC

  Uses the Roe solver to calculate the boundary flux using calc_p4 to
  get the boundary state.
"""->
type p4BC <: BCType
end

function call{Tmsh, Tsol}(obj::p4BC, u::Tsol, params::ParamTypes,
              coords::AbstractArray{Tmsh,1}, dxidx::AbstractArray{Tmsh, 2},
              nrm::AbstractArray{Tmsh,1}, t)

  u_bc = calc_p4(coords, params, t)
#  println("  u_bc = ", u_bc)
  bndryflux = RoeSolver(u, u_bc, params, nrm, dxidx)

  return bndryflux
end

@doc """
### AdvectionEquationMod.p5BC

  Uses the Roe solver to calculate the boundary flux using calc_p5 to
  get the boundary state.
"""->

type p5BC <: BCType
end

function call{Tmsh, Tsol}(obj::p5BC, u::Tsol, params::ParamTypes,
              coords::AbstractArray{Tmsh,1}, dxidx::AbstractArray{Tmsh, 2},
              nrm::AbstractArray{Tmsh,1}, t)
  # this is really slow: use Horner's rule!
  u_bc = calc_p5(coords, params, t)
#  println("calc_p5 @time printed above")
#  println("  u_bc = ", u_bc)
  bndryflux = RoeSolver(u, u_bc, params, nrm, dxidx)
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

function call{Tmsh, Tsol}(obj::exp5xplus4yplus2BC, u::Tsol, params::ParamType2,
              coords::AbstractArray{Tmsh,1}, dxidx::AbstractArray{Tmsh, 2},
              nrm::AbstractArray{Tmsh,1}, t)

  u_bc = calc_exp5xplus4yplus2(coords, params, t)
  bndryflux = RoeSolver(u, u_bc, params, nrm, dxidx)

  return bndryflux
end

@doc """
### AdvectionEquationMod.exp5xplusyBC

Uses the Roe solver to calculate the boundary flux using calc_exp5xplusy to get
the boundary state.

"""->

type exp5xplusyBC <:BCType
end

function call{Tmsh, Tsol}(obj::exp5xplusyBC, u::Tsol, params::ParamType2,
              coords::AbstractArray{Tmsh,1}, dxidx::AbstractArray{Tmsh, 2},
              nrm::AbstractArray{Tmsh,1}, t)

  u_bc = calc_exp5xplusy(coords, params, t)
  bndryflux = RoeSolver(u, u_bc, params, nrm, dxidx)

  return bndryflux
end

@doc """
### AdvectionEquationMod.exp3xplusyBC

Uses the Roe solver to calculate the boundary flux using calc_exp3xplusy to get
the boundary state.

"""->

type exp3xplusyBC <:BCType
end

function call{Tmsh, Tsol}(obj::exp3xplusyBC, u::Tsol, params::ParamType2,
              coords::AbstractArray{Tmsh,1}, dxidx::AbstractArray{Tmsh, 2},
              nrm::AbstractArray{Tmsh,1}, t)

  u_bc = calc_exp3xplusy(coords, params, t)
  bndryflux = RoeSolver(u, u_bc, params, nrm, dxidx)

  return bndryflux
end

@doc """
### AdvectionEquationMod.exp2xplus2yBC

Uses the Roe solver to calculate the boundary flux using calc_exp2xplus2y to get
the boundary state.

"""->

type exp2xplus2yBC <:BCType
end

function call{Tmsh, Tsol}(obj::exp2xplus2yBC, u::Tsol, params::ParamType2,
              coords::AbstractArray{Tmsh,1}, dxidx::AbstractArray{Tmsh, 2},
              nrm::AbstractArray{Tmsh,1}, t)

  u_bc = calc_exp2xplus2y(coords, params, t)
  bndryflux = RoeSolver(u, u_bc, params, nrm, dxidx)

  return bndryflux
end

@doc """
### AdvectionEquationMod.exp_xyBC

Uses the Roe solver to calculate the boundary flux using calc_exp_xy to get the
boundary state

"""->

type exp_xyBC <: BCType
end

function call{Tmsh, Tsol}(obj::exp_xyBC, u::Tsol, params::ParamType2,
              coords::AbstractArray{Tmsh,1}, dxidx::AbstractArray{Tmsh, 2},
              nrm::AbstractArray{Tmsh,1}, t)

  u_bc = calc_exp_xy(coords, params, t)
  bndryflux = RoeSolver(u, u_bc, params, nrm, dxidx)

  return bndryflux
end

@doc """
### AdvectionEquationMod.xplusyBC

Uses Roe solver to calculate the boundary flux using calc_xplusy to get the 
boundary state

"""->

type xplusyBC <: BCType
end

function call{Tmsh, Tsol}(obj::xplusyBC, u::Tsol, params::ParamType2,
              coords::AbstractArray{Tmsh,1}, dxidx::AbstractArray{Tmsh, 2},
              nrm::AbstractArray{Tmsh,1}, t)

  u_bc = calc_xplusy(coords, params, t)
  bndryflux = RoeSolver(u, u_bc, params, nrm, dxidx)

  return bndryflux
end


"""
  BC for unsteadymms
"""
type unsteadymmsBC <: BCType
end

function call{Tmsh, Tsol}(obj::unsteadymmsBC, u::Tsol, params::ParamType,
              coords::AbstractArray{Tmsh,1}, dxidx::AbstractArray{Tmsh, 2},
              nrm::AbstractArray{Tmsh,1}, t)

  u_bc = calc_unsteadymms(coords, params, t)
  bndryflux = RoeSolver(u, u_bc, params, nrm, dxidx)

  return bndryflux
end

"""
  BC for unsteadypoly
"""
type unsteadypolyBC <: BCType
end

function call{Tmsh, Tsol}(obj::unsteadypolyBC, u::Tsol, params::ParamType,
              coords::AbstractArray{Tmsh,1}, dxidx::AbstractArray{Tmsh, 2},
              nrm::AbstractArray{Tmsh,1}, t)

  u_bc = calc_unsteadypoly(coords, params, t)
  bndryflux = RoeSolver(u, u_bc, params, nrm, dxidx)

  return bndryflux
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
"mms1BC" => mms1BC(),
"x4BC" => x4BC(),
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
