# module to hold optimization interface

module OptimizationInterface

using PDESolver  # evalFunctional, calcFunctionalDeriv
using SummationByParts
using ODLCommonTools
using Utils
using LinearSolvers
using NonlinearSolvers  # need jacobian calculation function  # unneded?
using PdePumiInterface  # need write vtk files
using MPI

# steady_adjoint.jl
export calcAdjoint

# error_estimate.jl
export doHAdaptation, solveAdaptive, calcErrorEstimate, AdaptOpts

include("steady_adjoint.jl")
include("error_estimate.jl")
include("target_els.jl")

end  # end module
