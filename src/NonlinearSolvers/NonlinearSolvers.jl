module NonlinearSolvers

pde_pumi_interface_path = joinpath(Pkg.dir("PumiInterface"), "src")
push!(LOAD_PATH, pde_pumi_interface_path)
using ArrayViews
import ArrayViews.view
using PDESolver
using PdePumiInterface  # needed to write vtk files
using ODLCommonTools
using Utils
import ODLCommonTools.sview
import MPI
using PETSc2
using Utils
using SummationByParts
using LinearSolvers

import LinearSolvers: calcPC, applyPC, applyPCTranspose, calcLinearOperator,
                      applyLinearOperator, applyLinearOperatorTranspose,
                      needParallelData

include(joinpath(Pkg.dir("PDESolver"), "src/solver/debug.jl"))  # debug macro
include("rk4.jl")
include("lserk.jl")
  export lserk54
include("lserk_ds.jl")
  export lserk54_ds

include("explicit_euler.jl")
  export explicit_euler

include("drag_timeaverage.jl")
  export calcDragTimeAverage, calcFinalIter, calcQuadWeight

include("jac_recalc.jl")
include("preconditioning_types.jl")
include("newton.jl")
export getNewtonPCandLO

include("crank_nicolson.jl")
include("preconditioning.jl")
include("globalization.jl")
include("predictor_corrector.jl")
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

# jacobian.jl
export AssembleElementData, assembleElement, assembleInterface,
       assembleSharedFace, assembleBoundary, NullAssembleElementData

end
