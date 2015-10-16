module AdvectionEquationMod

using ArrayViews
using PDESolverCommon
using SummationByParts
using PdePumiInterface
using ForwardDiff
export AdvectionData, AdvectionData_, getMass

# include("advectionFunctions.jl")
include("getMass.jl")

abstract AbstractAdvectionData{Tsol} <: AbstractSolutionData{Tsol}
abstract AdvectionData{Tsol, Tdim} <: AbstractAdvectionData{Tsol}

@doc """
### AdvectionEquationMod.AdvectionData_

  This type is an implimentation of the abstract AdvectionData.  It is
  paramterized by the residual type Tres and the mesh type Tmsh
  because it stores some arrays of those types.  Tres is the 'maximum' type of
  Tsol and Tmsh, where Tsol is the type of the conservative variables.

"""->

type AdvectionData_{Tsol, Tres, Tdim, Tmsh} <: AdvectionData{Tsol, Tdim}

  # params::ParamType{Tdim}
  res_type::DataType  # type of res
  u::Array{Tsol, 3}
  res::Array{Tsol, 3}
  aux_vars::Array{Tres, 3}  # storage for auxiliary variables 
  F_xi::Array{Tsol,4}  # flux in xi direction
  res::Array{Tres, 3}  # result of computation
  SL::Array{Tres, 1}  # result of computation in vector form
  SL0::Array{Tres,1}  # initial condition in vector form
  bndryflux::Array{Tsol, 2}  # boundary flux
  M::Array{Float64, 2}  # mass matrix
  disassembleSolution::Function # function SL0 -> eqn.q
  assembleSolution::Function  # function : eqn.res -> SL

  function AdvectionData_(mesh::PumiMesh2, sbp::SBPOperator, opts)
    println("\nConstruction AdvectionData object")
    println("  Tsol = ", Tsol)
    println("  Tres = ", Tres)
    println("  Tdim = ", Tdim)
    println("  Tmsh = ", Tmsh)
    eqn = new()  # incomplete initilization
    eqn.res_type = Tres
    eqn.disassembleSolution = disassembleSolution
    eqn.assembleSolution = assembleSolution
    eqn.M = getMass(sbp, mesh)
    eqn.u = zeros(Tsol, 1, sbp.numnodes, mesh.numEl)
    eqn.res = zeros(Tsol, 1, sbp.numnodes, mesh.numEl)
    eqn.SL = zeros(Tres, mesh.numDof)
    eqn.SL0 = zeros(Tres, mesh.numDof)
    eqn.bndryflux = zeros(Tsol, 1, sbp.numfacenodes, mesh.numBoundaryEdges)

    return eqn
  end # ends the constructer AdvectionData_

end # End type AdvectionData_

@doc """
### AdvectionEquationMod.assembleSolution

  This function takes the 2D array of variables in arr and 
  reassmbles is into the vector SL.  Note that
  This is a reduction operation and requires eqn.SL to be zerod before
  calling this function.

  This is a mid level function, and does the right thing regardless of
  equation dimension

"""->

function assembleSolution{Tmsh, Tsol, Tres}(mesh::AbstractMesh{Tmsh}, 
                                            sbp::SBPOperator, 
                                            eqn::AdvectionData{Tsol}, opts, arr,
                                            SL::AbstractArray{Tres,1})

  for i=1:mesh.numEl  # loop over elements
    for j=1:mesh.numNodesPerElement
      dofnum_k = mesh.dofs[1, j, i]
      SL[dofnum_k] += arr[1,j,i]
    end
  end

  return nothing
end # end function assembleSolution

@doc """
### AdvectionEquationMod.disassembleSolution

  This takes eqn.SL0 (the initial state), and disassembles it into eqn.q, the
  3 dimensional array of conservative variables.  This function uses mesh.dofs
  to speed the process.

  This is a mid level function, and does the right thing regardless of equation
  dimension.
"""->

function disassembleSolution{Tmsh, Tsol, Tdim}(mesh::AbstractMesh{Tmsh}, sbp,
                                               eqn::AdvectionData{Tsol, Tdim},
                                               opts, 
                                               array::AbstractArray{Tsol, 1})
  # disassemble SL0 into eqn.
  for i=1:mesh.numEl  # loop over elements
    for j = 1:mesh.numNodesPerElement
      dofnum_k = mesh.dofs[1, j, i]
      eqn.u[1, j, i] = array[dofnum_k]
    end
  end

  # writeQ(mesh, sbp, eqn, opts)

  return nothing
end

end # end module

