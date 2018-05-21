# some functions used by all physics modules, particularly durin ginitialization

module SolverCommon

using PDESolver
using Utils
using ODLCommonTools
using LinearSolvers
using NonlinearSolvers
using Input
using MPI
using SummationByParts
using PdePumiInterface

# from common.jl
export createMeshAndOperator, loadRestartState, call_nlsolver, getDataTypes

include("common.jl")


end # end module
