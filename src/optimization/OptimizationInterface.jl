# module to hold optimization interface

module OptimizationInterface

using PDESolver  # evalFunctional, calcFunctionalDeriv
using SummationByParts
using ODLCommonTools
using Utils
using LinearSolvers
using NonlinearSolvers  # need jacobian calculation function

# steady_adjoint.jl
export calcAdjoint

include("steady_adjoint.jl")

end  # end module
