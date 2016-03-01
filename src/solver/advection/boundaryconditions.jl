# boundaryconditions.jl

@doc """
### AdvectionEquationMod.calcBoundaryFlux

  This function calculates the boundary flux for the portion of the boundary
  with a particular boundary condition.  The eqn.q are converted to 
  conservative variables if needed

  Inputs:
  mesh : AbstractMesh
  sbp : SBPOperator
  eqn : AdvectionEquation
  functor : a callable object that calculates the boundary flux at a node
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
function calcBoundaryFlux{Tmsh,  Tsol, Tres}( mesh::AbstractMesh{Tmsh}, 
                          sbp::SBPOperator, eqn::AdvectionData{Tsol}, 
                          functor::BCType, 
                          bndry_facenums::AbstractArray{Boundary,1}, 
                          bndryflux::AbstractArray{Tres, 3})
# calculate the boundary flux for the boundary condition evaluated by the functor

#  println("enterted calcBoundaryFlux")
  
  t = eqn.t
  nfaces = length(bndry_facenums)
  for i=1:nfaces  # loop over faces with this BC
    bndry_i = bndry_facenums[i]
#    println("element = ", bndry_i.element, ", face = ", bndry_i.face)
#    println("face ", i)

#      println("evaluating boundary flux for element ", bndry_i.element, " face ", bndry_i.face)
    for j = 1:sbp.numfacenodes
#      println("  node ", j)
      k = sbp.facenodes[j, bndry_i.face]

      # get components
      q = eqn.q[ 1, k, bndry_i.element]
      alpha_x = eqn.alpha_x[1, k, bndry_i.element]
      alpha_y = eqn.alpha_y[ 1, k, bndry_i.element]
      # flux_parametric = view(eqn.flux_parametric, :, k, bndry_i.element, :)
      # aux_vars = view(eqn.aux_vars, :, k, bndry_i.element)
      coords = view(mesh.coords, :, k, bndry_i.element)
      dxidx = view(mesh.dxidx, :, :, k, bndry_i.element)
      nrm = view(sbp.facenormal, :, bndry_i.face)
      #println("eqn.bndryflux = ", eqn.bndryflux)
      # functor(u, flux_parametric, aux_vars, coords, dxidx, nrm, bndryflux_i, eqn.params)
#      println("  q = ", q)
#      println("  alpha_x = ", alpha_x)
#      println("  alpha_Y = ", alpha_y)
#      println("  coords = ", coords)
#      println("  dxidx = ", dxidx)
#      println("  nrm = ", nrm)
#      println("  t = ", t)
      bndryflux[1, j, i] = -functor(q, alpha_x, alpha_y, coords, dxidx, nrm, t)
#      println("  bndryflux = ", bndryflux[1, j, i])
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
*  `alpha_x` & `alpha_y` : velocities in the X & Y directions
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
              alpha_x, alpha_y, coords::AbstractArray{Tmsh,1}, 
              dxidx::AbstractArray{Tmsh,2}, nrm::AbstractArray{Tmsh,1}, t)

  u_bc = calc_x5plusy5(coords) # Calculate the actual analytic value of u at the bondary
  bndryflux = RoeSolver(u, u_bc, alpha_x, alpha_y, nrm, dxidx)

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
              alpha_x, alpha_y, coords::AbstractArray{Tmsh,1}, 
              dxidx::AbstractArray{Tmsh,2}, nrm::AbstractArray{Tmsh,1}, t)

  u_bc = calc_exp_xplusy(coords, alpha_x, alpha_y, t)
  bndryflux = RoeSolver(u, u_bc, alpha_x, alpha_y, nrm, dxidx)

  return bndryflux
end

@doc """
### AdvectionEquationMod.sinwave_BC

  Uses the Roe solver to calculate the boundary flux using calc_sinewave to
  get the boundary state
"""->
type sinwave_BC <: BCType
end

function call{Tmsh, Tsol}(obj::sinwave_BC, u::Tsol, alpha_x, alpha_y,
              coords::AbstractArray{Tmsh,1}, dxidx::AbstractArray{Tmsh, 2},
              nrm::AbstractArray{Tmsh,1}, t)

  u_bc = calc_sinwave(coords, alpha_x, alpha_y, t)
#  println("  u_bc = ", u_bc)
  bndryflux = RoeSolver(u, u_bc, alpha_x, alpha_y, nrm, dxidx)

  return bndryflux
end

@doc """
### AdvectionEquationMod.sinwavey_BC

  Uses the Roe solver to calculate the boundary flux using calc_sinewavey to
  get the boundary state
"""->
type sinwavey_BC <: BCType
end

function call{Tmsh, Tsol}(obj::sinwavey_BC, u::Tsol, alpha_x, alpha_y,
              coords::AbstractArray{Tmsh,1}, dxidx::AbstractArray{Tmsh, 2},
              nrm::AbstractArray{Tmsh,1}, t)

  u_bc = calc_sinwavey(coords, alpha_x, alpha_y, t)
#  println("  u_bc = ", u_bc)
  bndryflux = RoeSolver(u, u_bc, alpha_x, alpha_y, nrm, dxidx)

  return bndryflux
end

@doc """
### AdvectionEquationMod.sinwavey_pertBC

  Uses the Roe solver to calculate the boundary flux using calc_sinewave_pert to
  get the boundary state
"""->

type sinwavey_pertBC <: BCType
end

function call{Tmsh, Tsol}(obj::sinwavey_pertBC, u::Tsol, alpha_x, alpha_y,
              coords::AbstractArray{Tmsh,1}, dxidx::AbstractArray{Tmsh, 2},
              nrm::AbstractArray{Tmsh,1}, t)

  u_bc = calc_sinwavey_pert(coords, alpha_x, alpha_y, t)
#  println("  u_bc = ", u_bc)
  bndryflux = RoeSolver(u, u_bc, alpha_x, alpha_y, nrm, dxidx)

  return bndryflux
end

@doc """
### AdvectionEquationMod.mms1BC

  Uses the Roe solver to calculate the boundary flux using calc_mms1 to get
  the boundary state.
"""->
type mms1BC <: BCType
end

function call{Tmsh, Tsol}(obj::mms1BC, u::Tsol, alpha_x, alpha_y,
              coords::AbstractArray{Tmsh,1}, dxidx::AbstractArray{Tmsh, 2},
              nrm::AbstractArray{Tmsh,1}, t)

  u_bc = calc_mms1(coords, alpha_x, alpha_y, t)
#  println("  u_bc = ", u_bc)
  bndryflux = RoeSolver(u, u_bc, alpha_x, alpha_y, nrm, dxidx)

  return bndryflux
end

@doc """
### AdvectionEquationMod.x4BC

  Uses the Roe solver to calculate the boundary flux using calc_x4 to
  get the boundary state.
"""->
type x4BC <: BCType
end

function call{Tmsh, Tsol}(obj::x4BC, u::Tsol, alpha_x, alpha_y,
              coords::AbstractArray{Tmsh,1}, dxidx::AbstractArray{Tmsh, 2},
              nrm::AbstractArray{Tmsh,1}, t)

  u_bc = calc_x4(coords, alpha_x, alpha_y, t)
#  println("  u_bc = ", u_bc)
  bndryflux = RoeSolver(u, u_bc, alpha_x, alpha_y, nrm, dxidx)

  return bndryflux
end

@doc """
### AdvectionEquationMod.p1BC

  Uses the Roe solver to calculate the boundary flux using calc_p1 to
  get the boundary state
"""->

type p1BC <: BCType
end

function call{Tmsh, Tsol}(obj::p1BC, u::Tsol, alpha_x, alpha_y,
              coords::AbstractArray{Tmsh,1}, dxidx::AbstractArray{Tmsh, 2},
              nrm::AbstractArray{Tmsh,1}, t)

  u_bc = calc_p1(coords, alpha_x, alpha_y, t)
#  println("  u_bc = ", u_bc)
  bndryflux = RoeSolver(u, u_bc, alpha_x, alpha_y, nrm, dxidx)

  return bndryflux
end

@doc """
### AdvectionEquationMod.p2BC

  Uses the Roe solver to calculate the boundary flux using calc_p2 to
  get the boundary state
"""->
type p2BC <: BCType
end

function call{Tmsh, Tsol}(obj::p2BC, u::Tsol, alpha_x, alpha_y,
              coords::AbstractArray{Tmsh,1}, dxidx::AbstractArray{Tmsh, 2},
              nrm::AbstractArray{Tmsh,1}, t)

  u_bc = calc_p2(coords, alpha_x, alpha_y, t)
#  println("  u_bc = ", u_bc)
  bndryflux = RoeSolver(u, u_bc, alpha_x, alpha_y, nrm, dxidx)

  return bndryflux
end

@doc """
### AdvectionEquationMod.p3BC

  Uses the Roe solver to calculate the boundary flux using calc_p3 to
  get the boundary state
"""->
type p3BC <: BCType
end

function call{Tmsh, Tsol}(obj::p3BC, u::Tsol, alpha_x, alpha_y,
              coords::AbstractArray{Tmsh,1}, dxidx::AbstractArray{Tmsh, 2},
              nrm::AbstractArray{Tmsh,1}, t)

  u_bc = calc_p3(coords, alpha_x, alpha_y, t)
#  println("  u_bc = ", u_bc)
  bndryflux = RoeSolver(u, u_bc, alpha_x, alpha_y, nrm, dxidx)

  return bndryflux
end


@doc """
### AdvectionEquationMod.p4BC

  Uses the Roe solver to calculate the boundary flux using calc_p4 to
  get the boundary state.
"""->
type p4BC <: BCType
end

function call{Tmsh, Tsol}(obj::p4BC, u::Tsol, alpha_x, alpha_y,
              coords::AbstractArray{Tmsh,1}, dxidx::AbstractArray{Tmsh, 2},
              nrm::AbstractArray{Tmsh,1}, t)

  u_bc = calc_p4(coords, alpha_x, alpha_y, t)
#  println("  u_bc = ", u_bc)
  bndryflux = RoeSolver(u, u_bc, alpha_x, alpha_y, nrm, dxidx)

  return bndryflux
end

@doc """
### AdvectionEquationMod.p5BC

  Uses the Roe solver to calculate the boundary flux using calc_p5 to
  get the boundary state.
"""->

type p5BC <: BCType
end

function call{Tmsh, Tsol}(obj::p5BC, u::Tsol, alpha_x, alpha_y,
              coords::AbstractArray{Tmsh,1}, dxidx::AbstractArray{Tmsh, 2},
              nrm::AbstractArray{Tmsh,1}, t)

  u_bc = calc_p5(coords, alpha_x, alpha_y, t)
#  println("  u_bc = ", u_bc)
  bndryflux = RoeSolver(u, u_bc, alpha_x, alpha_y, nrm, dxidx)

  return bndryflux
end


@doc """
### AdvectionEquationMod.BCDict

It stores all the possible boundary condition dictionary options. Whenever a 
new boundary condition is created, it should get added to BCDict.

"""->
global const BCDict = Dict{ASCIIString, BCType}(
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
function getBCFunctors{Tmsh, Tsol, Tdim}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator, 
                       eqn::AdvectionData{Tsol, Tdim}, opts)
# populate the array mesh.bndry_funcs with the functors for the boundary condition types

#  println("Entered getBCFunctors")

  for i=1:mesh.numBC
    key_i = string("BC", i, "_name")
    val = opts[key_i]
    println("BCDict[val] = ", BCDict[val])
    mesh.bndry_funcs[i] = BCDict[val]
  end

  return nothing
end
