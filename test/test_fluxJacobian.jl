# test if the flux jacobians are correct

push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/solver/euler"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/NonlinearSolvers"))

#include("complexify.jl")   # TODO: include location needs to be reconsidered

using ODLCommonTools
using PdePumiInterface  # common mesh interface - pumi
using SummationByParts  # SBP operators
using EulerEquationMod
using ForwardDiff
using NonlinearSolvers   # non-linear solvers
using ArrayViews

include(joinpath(Pkg.dir("PDESolver"),"src/solver/euler/output.jl"))  # printing results to files
include(joinpath(Pkg.dir("PDESolver"), "src/input/read_input.jl"))

#function runtest(flag::Int)
println("ARGS = ", ARGS)
println("size(ARGS) = ", size(ARGS))
opts = read_input(ARGS[1])
#opts = read_input("input_vals_channel2.jl")

# flag determines whether to calculate u, dR/du, or dR/dx (1, 2, or 3)
flag = opts["run_type"]

# timestepping parameters
delta_t = opts["delta_t"]   # delta_t: timestep for RK
t_max = opts["t_max"]       # t_max: maximum time for RK
order = opts["order"]       # order of accuracy

# types of the mesh, SBP, Equation objects
if flag == 1 || flag == 8  # normal run
  Tmsh = Float64
  Tsbp = Float64
  Tsol = Float64
  Tres = Float64
elseif flag == 2  # calculate dR/du
  Tmsh = Float64
  Tsbp = Float64
  Tsol = Float64
  Tres = Float64
elseif flag == 3  # calcualte dR/dx using forward mode
  Tmsh = Dual{Float64}
  Tsbp = Float64
  Tsol = Dual{Float64}
  Tres = Dual{Float64}
elseif flag == 4  # use Newton method using finite difference
  Tmsh = Float64
  Tsbp = Float64
  Tsol = Float64
  Tres = Float64
elseif flag == 5  # use complex step dR/du
  Tmsh = Float64
  Tsbp = Float64
  Tsol = Complex128
  Tres = Complex128
elseif flag == 6 || flag == 7  # evaluate residual error and print to paraview
  Tmsh = Float64
  Tsbp = Float64
  Tsol = Complex128
  Tres = Complex128
end

# create SBP object
println("\nConstructing SBP Operator")
sbp = TriSBP{Tsbp}(degree=order)  # create linear sbp operator

dmg_name = opts["dmg_name"]
smb_name = opts["smb_name"]

# create linear mesh with 4 dof per node
mesh = PumiMesh2{Tmsh}(dmg_name, smb_name, order, sbp, arg_dict; dofpernode=1)
eqn = EulerData_{Tsol, Tres, 2, Tmsh}(mesh, sbp, opts)