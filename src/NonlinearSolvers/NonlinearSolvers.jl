module NonlinearSolvers

#pde_pumi_interface_path = joinpath(Pkg.dir("PumiInterface"), "src")
#push!(LOAD_PATH, pde_pumi_interface_path)
using ArrayViews
import ArrayViews.view
using PDESolver
using PdePumiInterface  # needed to write vtk files
using ODLCommonTools
using Utils
import ODLCommonTools.sview
import MPI
using PETSc2  #TODO: see if this can be removed
using SummationByParts
using Jacobian
using LinearSolvers

import LinearSolvers: calcPC, applyPC, applyPCTranspose, calcLinearOperator,
                      applyLinearOperator, applyLinearOperatorTranspose,
                      needParallelData

include(joinpath(Pkg.dir("PDESolver"), "src/solver/debug.jl"))  # debug macro
include("rk4.jl")
include("lserk.jl")
export lserk54

include("jac_recalc.jl")
include("newton.jl")
export getNewtonPCandLO
include("preconditioning_types.jl")

include("crank_nicolson.jl")
include("preconditioning.jl")
include("globalization.jl")
include("predictor_corrector.jl")
include("p_homotopy.jl")
#include("newton_fd_old.jl")

#-----------------------------------------------------------
# NON WORKING AND EXPERIMENTAL!!!
# unsteady adjoint crank nicolson code:
include("crank_nicolson_uadj/crank_nicolson_uadj.jl")
include("crank_nicolson_uadj/crank_nicolson_uadj_jacandrhs.jl")
include("crank_nicolson_uadj/crank_nicolson_uadj_objective.jl")
include("crank_nicolson_uadj/crank_nicolson_uadj_checkpointstraight.jl")
#-----------------------------------------------------------

# predictor_corrector.jl
export predictorCorrectorHomotopy

# p_homotopy.jl
export pHomotopy

end
