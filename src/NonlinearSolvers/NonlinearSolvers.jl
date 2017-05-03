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
include("newton.jl")
include("globalization.jl")
#include("newton_fd_old.jl")

end
