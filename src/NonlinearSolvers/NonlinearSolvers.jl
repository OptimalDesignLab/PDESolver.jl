module NonlinearSolvers

pde_pumi_interface_path = joinpath(Pkg.dir("PumiInterface"), "src")
push!(LOAD_PATH, pde_pumi_interface_path)
using ArrayViews
using Debug
using PdePumiInterface  # needed to write vtk files
using ODLCommonTools
import MPI
using PETSc
using Utils

include("rk4.jl")
include("newton.jl")
include("globalization.jl")
#include("newton_fd_old.jl")

end
