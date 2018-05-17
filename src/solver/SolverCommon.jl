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

# from initialization.jl
export createMeshAndOperator, loadRestartState, call_nlsolver

include("common.jl")


end # end module
