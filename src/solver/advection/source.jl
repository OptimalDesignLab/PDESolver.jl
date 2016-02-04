# source.jl
# this file defines the source term functors for the advection equation,
# similar to the boundary condition functors


type SRC0 <: SRCType
end

function call(obj::SRC0, coords::AbstractVector, alpha_x, alpha_y, t)
  return 0
end

type SRC1 <: SRCType
end

function call(obj::SRC1, coords::AbstractVector, alpha_x, alpha_y, t)
  return 1
end

type SRCx <: SRCType
end

function call(obj::SRCx, coords::AbstractVector, alpha_x, alpha_y, t)
  return coords[1]
end

type SRCmms1 <: SRCType
end

function call(obj::SRCmms1, coords::AbstractVector, alpha_x, alpha_y, t)
  return alpha_x*calc_mms1dx(coords, alpha_x, alpha_y, t)
end



@doc """
### AdvectionEquationMod.SRCDict

It stores all the possible boundary condition dictionary options. Whenever a 
new boundary condition is created, it should get added to BCDict.

"""->
global const SRCDict = Dict{ASCIIString, SRCType} (
"SRC0" => SRC0(),
"SRC1" => SRC1(),
"SRCx" => SRCx(),
"SRCmms1" => SRCmms1(),
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
