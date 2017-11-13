module NonlinearSolvers

pde_pumi_interface_path = joinpath(Pkg.dir("PumiInterface"), "src")
push!(LOAD_PATH, pde_pumi_interface_path)
using ArrayViews
using Debug
using PdePumiInterface  # needed to write vtk files
using ODLCommonTools
using Utils
import ODLCommonTools.sview
import MPI
using PETSc
using Utils
using SummationByParts

include(joinpath(Pkg.dir("PDESolver"), "src/solver/debug.jl"))  # debug macro
include("rk4.jl")
include("lserk.jl")
  export lserk54

include("crank_nicolson.jl")
include("preconditioning_types.jl")
include("newton.jl")
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

end
