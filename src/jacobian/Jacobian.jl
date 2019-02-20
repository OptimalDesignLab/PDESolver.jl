module Jacobian

using PDESolver
using ODLCommonTools
using SummationByParts
using PETSc2
using Utils

# residual_evaluation.jl
export physicsRhs, assembleResidual

# sparse.jl
export getBlockSparsityCounts, INVISCID, VISCOUS, COLORING, VISCOUSTIGHT

# jacobian.jl - main functions
export physicsJac

# jacobian.jl - explicit Jacobian calculation
export _AssembleElementData, assembleElement, assembleInterface,
       assembleInterfaceVisc, assembleSharedFace, assembleBoundary,
       assembleBoundaryFull, assembleInterfaceFull, assembleSharedFaceFull,
       NullAssembleElementData, calcJacCol

# jacobian_diag.jl
export DiagJac, diagMatVec, AssembleDiagJacData

# experimental stuff
export filterDiagJac, removeUnstableModes!, findStablePerturbation!

include("residual_evaluation.jl")
include("sparse.jl")
include("jacobian.jl")
include("jacobian_diag.jl")

end  # module
