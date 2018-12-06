# source.jl
# this file defines the source term functors for the advection equation,
# similar to the boundary condition functors

@doc """
### AdvectionEquationMod.SRC0

  This is the zero source term.  This is the default of source term
  is specified
"""->
mutable struct SRC0 <: SRCType
end

function (obj::SRC0)(params::ParamTypes, coords::AbstractVector, t)
  return 0
end

@doc """
### AdvectionEquationMod.SRC1

  This source term returns 1 everywhere.
"""->
mutable struct SRC1 <: SRCType
end

function (obj::SRC1)(params::ParamType2, coords::AbstractVector, t)
  return 1
end

@doc """
### AdvectionEquationMod.SRC2

  This source term returns 1 everywhere.
"""->
mutable struct SRC2 <: SRCType
end

function (obj::SRC2)(params::ParamType2, coords::AbstractVector, t)
  return 2
end

@doc """
### AdvectionEquationMod.SRCx

  This source term that returns: f = x
"""->
mutable struct SRCx <: SRCType
end

function (obj::SRCx)(params::ParamType2, coords::AbstractVector, t)
  return coords[1]
end

@doc """
### AdvectionEquationMod.SRC1

  This source term that returns: the derivative of mms1
"""->
mutable struct SRCmms1 <: SRCType
end

function (obj::SRCmms1)(params::ParamType2, coords::AbstractVector, t)
  return params.alpha_x*calc_mms1dx(params, coords, t)
end

@doc """
### AdvectionEquationMod.SRCx4

  This source term that returns: the source term for a manufactured solution
  using a 4th order polynomial
"""->
mutable struct SRCx4 <: SRCType
end

function (obj::SRCx4)(params::ParamType2, coords::AbstractVector, t)
  return params.alpha_x*calc_x4der(params, coords, t)
end

@doc """
### AdvectionEquationMod.SRCp1

  This source term that returns: the source term for a manufactured solution
  using a 1st order polynomial
"""->
mutable struct SRCp1 <: SRCType
end

function (obj::SRCp1)(params::ParamType2, coords::AbstractVector, t)
  return params.alpha_x*calc_p1dx(params, coords, t) + params.alpha_y*calc_p1dy(params, coords, t)

end

function (obj::SRCp1)(params::ParamType3, coords::AbstractVector, t)
  return params.alpha_x*calc_p1dx(params, coords, t) + params.alpha_y*calc_p1dy(params, coords, t) + params.alpha_z*calc_p1dz(params, coords, t)
end


@doc """
### AdvectionEquationMod.SRCp3

  This source term that returns: the source term for a manufactured solution
  using a 2nd order polynomial
"""->
mutable struct SRCp2 <: SRCType
end

function (obj::SRCp2)(params::ParamType2, coords::AbstractVector, t)
  return params.alpha_x*calc_p2dx(params, coords, t) + params.alpha_y*calc_p2dy(params, coords, t)

end

function (obj::SRCp2)(params::ParamType3, coords::AbstractVector, t)
  return params.alpha_x*calc_p2dx(params, coords, t) + params.alpha_y*calc_p2dy(params, coords, t) + params.alpha_z*calc_p2dz(params, coords, t)

end


@doc """
### AdvectionEquationMod.SRCp3

  This source term that returns: the source term for a manufactured solution
  using a 3rd order polynomial
"""->
mutable struct SRCp3 <: SRCType
end

function (obj::SRCp3)(params::ParamType2, coords::AbstractVector, t)
  return params.alpha_x*calc_p3dx(params, coords, t) + params.alpha_y*calc_p3dy(params, coords, t)

end

function (obj::SRCp3)(params::ParamType3, coords::AbstractVector, t)
  return params.alpha_x*calc_p3dx(params, coords, t) + params.alpha_y*calc_p3dy(params, coords, t) + params.alpha_z*calc_p3dz(params, coords, t)

end


@doc """
### AdvectionEquationMod.SRCp4

  This source term that returns: the source term for a manufactured solution
  using a 4th order polynomial
"""->
mutable struct SRCp4 <: SRCType
end

function (obj::SRCp4)(params::ParamType2, coords::AbstractVector, t)
  return params.alpha_x*calc_p4dx(params, coords, t) + params.alpha_y*calc_p4dy(params, coords, t)
end

function (obj::SRCp4)(params::ParamType3, coords::AbstractVector, t)
  return params.alpha_x*calc_p4dx(params, coords, t) + params.alpha_y*calc_p4dy(params, coords, t) + params.alpha_z*calc_p4dz(params, coords, t)
end


@doc """
### AdvectionEquationMod.SRCp5

  This source term that returns: the source term for a manufactured solution
  using a 5th order polynomial
"""->
mutable struct SRCp5 <: SRCType
end

function (obj::SRCp5)(params::ParamType2, coords::AbstractVector, t)
  return params.alpha_x*calc_p5dx(params, coords, t) + params.alpha_y*calc_p5dy(params, coords, t)
end

function (obj::SRCp5)(params::ParamType3, coords::AbstractVector, t)
  return params.alpha_x*calc_p5dx(params, coords, t) + params.alpha_y*calc_p5dy(params, coords, t) + params.alpha_z*calc_p5dz(params, coords, t)
end


@doc """
### AdvectionEquationMod.SRCexp_xplusy

  This is a source term that returns a source term for e^(x+y)
"""->
mutable struct SRCexp_xplusy <: SRCType
end

function (obj::SRCexp_xplusy)(params::ParamType2, coords::AbstractVector, t)
  u = calc_exp_xplusy(params, coords, t)
  return params.alpha_x*u + params.alpha_y*u 
end

@doc """
### AdvectionEquationMod.SRCx5plusy5

  This is a source term that returns a source term for e^(x+y)
"""->
mutable struct SRCx5plusy5 <: SRCType
end

function (obj::SRCx5plusy5)(params::ParamType2, coords::AbstractVector, t)
  # u = calc_x5plusy5(params, coords, t)
  x = coords[1]
  y = coords[2]

  return params.alpha_x*5*(x^4) + params.alpha_y*5*(y^4)
end

@doc """
### AdvectionEquationMod.SRCexp5xplus4yplus2

Calculates the source term for q = exp(5*x + 4*y +2)

"""->
mutable struct SRCexp5xplus4yplus2 <: SRCType
end

function (obj::SRCexp5xplus4yplus2)(params::ParamType2, coords::AbstractVector, t)

  u = calc_exp5xplus4yplus2(params, coords, t)
  return params.alpha_x*5*u + params.alpha_y*4*u
end

@doc """
### AdvectionEquationMod.SRCexp5xplusy

Calculates the source term for q = exp(5*x + y)
"""->
mutable struct SRCexp5xplusy <: SRCType
end

function (obj::SRCexp5xplusy)(params::ParamType2, coords::AbstractVector, t)
  u = calc_exp5xplusy(params, coords, t)
  return params.alpha_x*5*u + params.alpha_y*u
end

@doc """
### AdvectionEquationMod.SRCexp3xplusy

Calculates the source term for q = exp(3*x + y)
"""->
mutable struct SRCexp3xplusy <: SRCType
end

function (obj::SRCexp3xplusy)(params::ParamType2, coords::AbstractVector, t)
  u = calc_exp3xplusy(params, coords, t)
  return params.alpha_x*3*u + params.alpha_y*u
end

@doc """
### AdvectionEquationMod.SRCexp2xplus2y

Calculates the source term for q = exp(2*x + 2*y)
"""->
mutable struct SRCexp2xplus2y <: SRCType
end

function (obj::SRCexp2xplus2y)(params::ParamType2, coords::AbstractVector, t)
  u = calc_exp2xplus2y(params, coords, t)
  return params.alpha_x*2*u + params.alpha_y*2*u
end

@doc """
### AdvectionEquationMod.SRCexp_xy

Calculates the source term for q = exp(x*y)
"""->

mutable struct SRCexp_xy <: SRCType
end

function (obj::SRCexp_xy)(params::ParamType2, coords::AbstractVector, t)
  x = coords[1]
  y = coords[2]
  u = calc_exp_xy(params, coords, t)
  return params.alpha_x*y*u + params.alpha_y*x*u
end

@doc """
### AdvectionEquationMod.SRCxplusy

calculates the source term for q = x + y

"""->
mutable struct SRCxplusy <: SRCType
end

function (obj::SRCxplusy)(params::ParamType2, coords::AbstractVector, t)

  return params.alpha_x + params.alpha_y
end

"""
  Source term for unsteady mms
"""
mutable struct SRCunsteadymms <: SRCType
end

function (obj::SRCunsteadymms)(params::ParamType2, coords::AbstractVector, t)

  alpha_x = params.alpha_x
  alpha_y = params.alpha_y

  dudt = calc_unsteadymmsdt(params, coords, t)
  dudx = calc_unsteadymmsdx(params, coords, t)
  dudy = calc_unsteadymmsdy(params, coords, t)
  return  dudt + alpha_x*dudx + alpha_y*dudy
end

function (obj::SRCunsteadymms)(params::ParamType3, coords::AbstractVector, t)

  alpha_x = params.alpha_x
  alpha_y = params.alpha_y
  alpha_z = params.alpha_z

  dudt = calc_unsteadymmsdt(params, coords, t)
  dudx = calc_unsteadymmsdx(params, coords, t)
  dudy = calc_unsteadymmsdy(params, coords, t)
  dudz = calc_unsteadymmsdz(params, coords, t)
  return  dudt + alpha_x*dudx + alpha_y*dudy + alpha_z*dudz
end

"""
  Source term for unsteady poly
"""
mutable struct SRCunsteadypoly <: SRCType
end

function (obj::SRCunsteadypoly)(params::ParamType2, coords::AbstractVector, t)

  alpha_x = params.alpha_x
  alpha_y = params.alpha_y

  dudt = calc_unsteadypolydt(params, coords, t)
  dudx = calc_unsteadypolydx(params, coords, t)
  dudy = calc_unsteadypolydy(params, coords, t)
 

  return dudt + alpha_x*dudx + alpha_y*dudy
end

@doc """
### AdvectionEquationMod.SRCDict

  It stores all the possible source term dictionary options. Whenever a 
  new source is created, it should get added to SRCDict.

  All functors must have the signature:

  src_func(params, coords::ParamType, t)

  where coords is the vector of length 2 containing the x and y coordinates
  of the node, params.alpha_x and params.alpha_y are the advection velocities in the x an y
  directions, and t is the current time
"""->
global const SRCDict = Dict{String, SRCType}(
"SRC0" => SRC0(),
"SRC1" => SRC1(),
"SRC2" => SRC2(),
"SRCx" => SRCx(),
"SRCmms1" => SRCmms1(),
"SRCx4" => SRCx4(),
"SRCp0" => SRC0(),  # this *should* avoid calculating teh source term?
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
"SRCunsteadymms" => SRCunsteadymms(),
"SRCunsteadypoly" => SRCunsteadypoly(),
)

@doc """

  This function gets the functor specified by opts["SRCname"] and stores
  it to the equation object.  Currently one 1 source functor is allowed.

"""->
function getSRCFunctors(mesh::AbstractMesh, sbp::AbstractOperator, 
                        eqn::AdvectionData, opts)

  # currently we only allow 1 source functor
  eqn.src_func = SRCDict[opts["SRCname"]]
#  println("using source term functor ", eqn.src_func)
  return nothing
end
