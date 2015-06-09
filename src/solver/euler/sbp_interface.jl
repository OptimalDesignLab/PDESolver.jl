# this provides an interface for any SBP function that take a function as an argument
# this interfaces uses scoping to allow the underlying function to take
# more arguments than SBP allows

using SummationByParts
#using Equation
#include("abstract_types.jl")
using PdePumiInterface
using PDESolverCommon


import SummationByParts.boundaryintegrate2!
import SummationByParts.edgestabilize!

# any function passed into this interface should take all the same arguments
# required by SBP, plus an AbstractMesh and and AbstractEquation as the last two


function boundaryintegrate2!{T}(sbp::SBPOperator{T}, bndryfaces::Array{Boundary},
                               u::AbstractArray{T,3}, x::AbstractArray{T,3},
                               dxidx::AbstractArray{T,4}, bndryflux::Function,
                               res::AbstractArray{T,3}, mesh::AbstractMesh, eqn::AbstractEquation)

  # define a nested function
  function bndryflux_wrapper{T}(q::AbstractArray{T,1}, x::AbstractArray{T,1}, dxidx::AbstractArray{T,2}, nrm::AbstractArray{T,1}, flux::AbstractArray{T,1})
    # call the boundary flux function with all the arguments
    bndryflux(q, x, dxidx, nrm, flux, mesh, eqn)



    return nothing

  end
  # call the SBP boundaryintegrate! with bndryflux_wrapper
  boundaryintegrate!(sbp, bndryfaces, u, x, dxidx, bndryflux_wrapper, res)

  return nothing
end



function edgestabilize!{T}(sbp::SBPOperator{T}, ifaces::Array{Interface},
                           u::AbstractArray{T,3}, x::AbstractArray{T,3},
                           dxidx::AbstractArray{T,4}, jac::AbstractArray{T,2},
                           alpha::AbstractArray{T,4}, stabscale::Function,
                           res::AbstractArray{T,3}, mesh::AbstractMesh, eqn::AbstractEquation)

  function stabscale_wrapper{T}(u::AbstractArray{T,1}, dxidx::AbstractArray{T,2}, nrm::AbstractArray{T,1})


    return stabscale(u, dxidx, nrm, mesh, eqn)
  end

  # call the SBP edgestabilize with stabscale_wrapper
  edgestabilize!(sbp, ifaces, u, x, dxidx, jac, alpha, stabscale_wrapper, res)

  return nothing
end
