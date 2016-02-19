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

function call(obj::SRC0, coords::AbstractVector, alpha_x, alpha_y, t)
  return 0
end

@doc """
### AdvectionEquationMod.SRC1

  This source term returns 1 everywhere.
"""->
type SRC1 <: SRCType
end

function call(obj::SRC1, coords::AbstractVector, alpha_x, alpha_y, t)
  return 1
end

@doc """
### AdvectionEquationMod.SRC1

  This source term that returns: f = x
"""->
type SRCx <: SRCType
end

function call(obj::SRCx, coords::AbstractVector, alpha_x, alpha_y, t)
  return coords[1]
end

@doc """
### AdvectionEquationMod.SRC1

  This source term that returns: the derivative of mms1
"""->
type SRCmms1 <: SRCType
end

function call(obj::SRCmms1, coords::AbstractVector, alpha_x, alpha_y, t)
  return alpha_x*calc_mms1dx(coords, alpha_x, alpha_y, t)
end

@doc """
### AdvectionEquationMod.SRCx4

  This source term that returns: the source term for a manufactured solution
  using a 4th order polynomial
"""->
type SRCx4 <: SRCType
end

function call(obj::SRCx4, coords::AbstractVector, alpha_x, alpha_y, t)
  return alpha_x*calc_x4der(coords, alpha_x, alpha_y, t)
end

@doc """
### AdvectionEquationMod.SRCp1

  This source term that returns: the source term for a manufactured solution
  using a 1st order polynomial
"""->
type SRCp1 <: SRCType
end

function call(obj::SRCp1, coords::AbstractVector, alpha_x, alpha_y, t)
  return alpha_x*calc_p1dx(coords, alpha_x, alpha_y, t) + alpha_y*calc_p1dy(coords, alpha_x, alpha_y, t)

end

@doc """
### AdvectionEquationMod.SRCp3

  This source term that returns: the source term for a manufactured solution
  using a 2nd order polynomial
"""->
type SRCp2 <: SRCType
end

function call(obj::SRCp2, coords::AbstractVector, alpha_x, alpha_y, t)
  return alpha_x*calc_p2dx(coords, alpha_x, alpha_y, t) + alpha_y*calc_p2dy(coords, alpha_x, alpha_y, t)

end

@doc """
### AdvectionEquationMod.SRCp3

  This source term that returns: the source term for a manufactured solution
  using a 3rd order polynomial
"""->
type SRCp3 <: SRCType
end

function call(obj::SRCp3, coords::AbstractVector, alpha_x, alpha_y, t)
  return alpha_x*calc_p3dx(coords, alpha_x, alpha_y, t) + alpha_y*calc_p3dy(coords, alpha_x, alpha_y, t)

end

@doc """
### AdvectionEquationMod.SRCp4

  This source term that returns: the source term for a manufactured solution
  using a 4th order polynomial
"""->
type SRCp4 <: SRCType
end

function call(obj::SRCp4, coords::AbstractVector, alpha_x, alpha_y, t)
  return alpha_x*calc_p4dx(coords, alpha_x, alpha_y, t) + alpha_y*calc_p4dy(coords, alpha_x, alpha_y, t)
end

@doc """
### AdvectionEquationMod.SRCp5

  This source term that returns: the source term for a manufactured solution
  using a 5th order polynomial
"""->
type SRCp5 <: SRCType
end

function call(obj::SRCp5, coords::AbstractVector, alpha_x, alpha_y, t)
  return alpha_x*calc_p5dx(coords, alpha_x, alpha_y, t) + alpha_y*calc_p5dy(coords, alpha_x, alpha_y, t)
end


@doc """
### AdvectionEquationMod.SRCDict

  It stores all the possible boundary condition dictionary options. Whenever a 
  new boundary condition is created, it should get added to BCDict.

  All functors must have the signature:

  src_func(coords, alpha_x, alpha_y, t)

  where coords is the vector of length 2 containing the x and y coordinates
  of the node, alpha_x and alpha_y are the advection velocities in the x an y
  directions, and t is the current time
"""->
global const SRCDict = Dict{ASCIIString, SRCType} (
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
)

@doc """
### AdvectionEquationMod.getSRCFunctors

  This function gets the functor specified by opts["SRCname"] and stores
  it to the equation object.  Currently one 1 source functor is allowed.

"""->
function getSRCFunctors(mesh::AbstractMesh, sbp::SBPOperator, 
                        eqn::AdvectionData, opts)

  # currently we only allow 1 source functor
  eqn.src_func = SRCDict[opts["SRCname"]]
  println("using source term functor ", eqn.src_func)
  return nothing
end
