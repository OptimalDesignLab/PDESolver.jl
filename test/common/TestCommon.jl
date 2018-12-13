module TestCommon
# tests that can be run by any physics module

using Base.Test

using ODLCommonTools
using PDESolver
using Utils
using OptimizationInterface
using PdePumiInterface


# fd_tests.jl
export testMeshDerivatives, testResidualDerivatives, testDJDx

include("fd_tests.jl")



end  # end module
