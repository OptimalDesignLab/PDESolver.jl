module AdvectionEquationMod

using ArrayViews
using PDESolverCommon
using SummationByParts
using PdePumiInterface
using ForwardDiff
export AbstractEulerData, EulerData, EulerData_

abstract AdvectionData {Tsol, Tdim} <: AbstractAdvectionData{Tsol}

@doc """
### AdvectionEquationMod.AdvectionData_

  This type is an implimentation of the abstract AdvectionData.  It is
  paramterized by the residual type Tres and the mesh type Tmsh
  because it stores some arrays of those types.  Tres is the 'maximum' type of
  Tsol and Tmsh, where Tsol is the type of the conservative variables.

"""->

type AdvectionData_{Tsol, Tres, Tdim, Tmsh} <: AdvectionData{Tsol, Tdim}

  params::ParamType{Tdim}
  res_type::DataType  # type of res
  u::Array{Tsol,2}  
  aux_vars::Array{Tres, 3}  # storage for auxiliary variables 
  F_xi::Array{Tsol,4}  # flux in xi direction
  res::Array{Tres, 3}  # result of computation
  SL::Array{Tres, 1}  # result of computation in vector form
  SL0::Array{Tres,1}  # initial condition in vector form
  bndryflux::Array{Tsol, 3}  # boundary flux
  M::Array{Float64, 1}  # mass matrix
  disassembleSolution::Function # function SL0 -> eqn.q
  assembleSolution::Function  # function : eqn.res -> SL

  function AdvectionData_(mesh::PumiMesh2, sbp::SBPOperator, opts)

  end # ends the constructer AdvectionData_



end # End type AdvectionData_