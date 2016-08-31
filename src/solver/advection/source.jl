# source.jl
# this file defines the source term functors for the advection equation,
# similar to the boundary condition functors

@doc """
### AdvectionEquationMod.SRC0

  This is the zero source term.  This is the default of source term
  is specified
"""->
type SRC0 <: SRCType
end

function call(obj::SRC0, coords::AbstractVector, params::ParamType2, t)
  return 0
end

@doc """
### AdvectionEquationMod.SRC1

  This source term returns 1 everywhere.
"""->
type SRC1 <: SRCType
end

function call(obj::SRC1, coords::AbstractVector, params::ParamType2, t)
  return 1
end

@doc """
### AdvectionEquationMod.SRC1

  This source term that returns: f = x
"""->
type SRCx <: SRCType
end

function call(obj::SRCx, coords::AbstractVector, params::ParamType2, t)
  return coords[1]
end

@doc """
### AdvectionEquationMod.SRC1

  This source term that returns: the derivative of mms1
"""->
type SRCmms1 <: SRCType
end

function call(obj::SRCmms1, coords::AbstractVector, params::ParamType2, t)
  return params.alpha_x*calc_mms1dx(coords, params, t)
end

@doc """
### AdvectionEquationMod.SRCx4

  This source term that returns: the source term for a manufactured solution
  using a 4th order polynomial
"""->
type SRCx4 <: SRCType
end

function call(obj::SRCx4, coords::AbstractVector, params::ParamType2, t)
  return params.alpha_x*calc_x4der(coords, params, t)
end

@doc """
### AdvectionEquationMod.SRCp1

  This source term that returns: the source term for a manufactured solution
  using a 1st order polynomial
"""->
type SRCp1 <: SRCType
end

function call(obj::SRCp1, coords::AbstractVector, params::ParamType2, t)
  return params.alpha_x*calc_p1dx(coords, params, t) + params.alpha_y*calc_p1dy(coords, params, t)

end

function call(obj::SRCp1, coords::AbstractVector, params::ParamType3, t)
  return params.alpha_x*calc_p1dx(coords, params, t) + params.alpha_y*calc_p1dy(coords, params, t) + params.alpha_z*calc_p1dz(coords, params, t)
end


@doc """
### AdvectionEquationMod.SRCp3

  This source term that returns: the source term for a manufactured solution
  using a 2nd order polynomial
"""->
type SRCp2 <: SRCType
end

function call(obj::SRCp2, coords::AbstractVector, params::ParamType2, t)
  return params.alpha_x*calc_p2dx(coords, params, t) + params.alpha_y*calc_p2dy(coords, params, t)

end

function call(obj::SRCp2, coords::AbstractVector, params::ParamType3, t)
  return params.alpha_x*calc_p2dx(coords, params, t) + params.alpha_y*calc_p2dy(coords, params, t) + params.alpha_z*calc_p2dz(coords, params, t)

end


@doc """
### AdvectionEquationMod.SRCp3

  This source term that returns: the source term for a manufactured solution
  using a 3rd order polynomial
"""->
type SRCp3 <: SRCType
end

function call(obj::SRCp3, coords::AbstractVector, params::ParamType2, t)
  return params.alpha_x*calc_p3dx(coords, params, t) + params.alpha_y*calc_p3dy(coords, params, t)

end

function call(obj::SRCp3, coords::AbstractVector, params::ParamType3, t)
  return params.alpha_x*calc_p3dx(coords, params, t) + params.alpha_y*calc_p3dy(coords, params, t) + params.alpha_z*calc_p3dz(coords, params, t)

end


@doc """
### AdvectionEquationMod.SRCp4

  This source term that returns: the source term for a manufactured solution
  using a 4th order polynomial
"""->
type SRCp4 <: SRCType
end

function call(obj::SRCp4, coords::AbstractVector, params::ParamType2, t)
  return params.alpha_x*calc_p4dx(coords, params, t) + params.alpha_y*calc_p4dy(coords, params, t)
end

function call(obj::SRCp4, coords::AbstractVector, params::ParamType3, t)
  return params.alpha_x*calc_p4dx(coords, params, t) + params.alpha_y*calc_p4dy(coords, params, t) + params.alpha_z*calc_p4dz(coords, params, t)
end


@doc """
### AdvectionEquationMod.SRCp5

  This source term that returns: the source term for a manufactured solution
  using a 5th order polynomial
"""->
type SRCp5 <: SRCType
end

function call(obj::SRCp5, coords::AbstractVector, params::ParamType2, t)
  return params.alpha_x*calc_p5dx(coords, params, t) + params.alpha_y*calc_p5dy(coords, params, t)
end

function call(obj::SRCp5, coords::AbstractVector, params::ParamType3, t)
  return params.alpha_x*calc_p5dx(coords, params, t) + params.alpha_y*calc_p5dy(coords, params, t) + params.alpha_z*calc_p5dz(coords, params, t)
end


@doc """
### AdvectionEquationMod.SRCexp_xplusy

  This is a source term that returns a source term for e^(x+y)
"""->
type SRCexp_xplusy <: SRCType
end

function call(obj::SRCexp_xplusy, coords::AbstractVector, params::ParamType2, t)
  u = calc_exp_xplusy(coords, params, t)
  return params.alpha_x*u + params.alpha_y*u 
end

@doc """
### AdvectionEquationMod.SRCx5plusy5

  This is a source term that returns a source term for e^(x+y)
"""->
type SRCx5plusy5 <: SRCType
end

function call(obj::SRCx5plusy5, coords::AbstractVector, params::ParamType2, t)
  # u = calc_x5plusy5(coords, params, t)
  x = coords[1]
  y = coords[2]

  return params.alpha_x*5*(x^4) + params.alpha_y*5*(y^4)
end

@doc """
### AdvectionEquationMod.SRCexp5xplus4yplus2

Calculates the source term for q = exp(5*x + 4*y +2)

"""->
type SRCexp5xplus4yplus2 <: SRCType
end

function call(obj::SRCexp5xplus4yplus2, coords::AbstractVector, params::ParamType2, t)

  u = calc_exp5xplus4yplus2(coords, params, t)
  return params.alpha_x*5*u + params.alpha_y*4*u
end

@doc """
### AdvectionEquationMod.SRCexp5xplusy

Calculates the source term for q = exp(5*x + y)
"""->
type SRCexp5xplusy <: SRCType
end

function call(obj::SRCexp5xplusy, coords::AbstractVector, params::ParamType2, t)
  u = calc_exp5xplusy(coords, params, t)
  return params.alpha_x*5*u + params.alpha_y*u
end

@doc """
### AdvectionEquationMod.SRCexp3xplusy

Calculates the source term for q = exp(3*x + y)
"""->
type SRCexp3xplusy <: SRCType
end

function call(obj::SRCexp3xplusy, coords::AbstractVector, params::ParamType2, t)
  u = calc_exp3xplusy(coords, params, t)
  return params.alpha_x*3*u + params.alpha_y*u
end

@doc """
### AdvectionEquationMod.SRCexp2xplus2y

Calculates the source term for q = exp(2*x + 2*y)
"""->
type SRCexp2xplus2y <: SRCType
end

function call(obj::SRCexp2xplus2y, coords::AbstractVector, params::ParamType2, t)
  u = calc_exp2xplus2y(coords, params, t)
  return params.alpha_x*2*u + params.alpha_y*2*u
end

@doc """
### AdvectionEquationMod.SRCexp_xy

Calculates the source term for q = exp(x*y)
"""->

type SRCexp_xy <: SRCType
end

function call(obj::SRCexp_xy, coords::AbstractVector, params::ParamType2, t)
  x = coords[1]
  y = coords[2]
  u = calc_exp_xy(coords, params, t)
  return params.alpha_x*y*u + params.alpha_y*x*u
end

@doc """
### AdvectionEquationMod.SRCxplusy

calculates the source term for q = x + y

"""->
type SRCxplusy <: SRCType
end

function call(obj::SRCxplusy, coords::AbstractVector, params::ParamType2, t)

  return params.alpha_x + params.alpha_y
end


@doc """
### AdvectionEquationMod.SRCDict

  It stores all the possible boundary condition dictionary options. Whenever a 
  new boundary condition is created, it should get added to BCDict.

  All functors must have the signature:

  src_func(coords, params::ParamType, t)

  where coords is the vector of length 2 containing the x and y coordinates
  of the node, params.alpha_x and params.alpha_y are the advection velocities in the x an y
  directions, and t is the current time
"""->
global const SRCDict = Dict{ASCIIString, SRCType}(
"SRC0" => SRC0(),
"SRC1" => SRC1(),
"SRCx" => SRCx(),
"SRCmms1" => SRCmms1(),
"SRCx4" => SRCx4(),
"SRCp1" => SRCp1(),
"SRCp2" => SRCp2(),
"SRCp3" => SRCp3(),
"SRCp4" => SRCp4(),
"SRCp5" => SRCp5(),
"SRCexp_xplusy" => SRCexp_xplusy(),
"SRCx5plusy5" => SRCx5plusy5(),
"SRCexp5xplus4yplus2" => SRCexp5xplus4yplus2(),
"SRCexp5xplusy" => SRCexp5xplusy(),
"SRCexp3xplusy" => SRCexp3xplusy(),
"SRCexp2xplus2y" => SRCexp2xplus2y(),
"SRCexp_xy" => SRCexp_xy(),
"SRCxplusy" => SRCxplusy(),
)

@doc """
### AdvectionEquationMod.getSRCFunctors

  This function gets the functor specified by opts["SRCname"] and stores
  it to the equation object.  Currently one 1 source functor is allowed.

"""->
function getSRCFunctors(mesh::AbstractMesh, sbp::AbstractSBP, 
                        eqn::AdvectionData, opts)

  # currently we only allow 1 source functor
  eqn.src_func = SRCDict[opts["SRCname"]]
  println("using source term functor ", eqn.src_func)
  return nothing
end
