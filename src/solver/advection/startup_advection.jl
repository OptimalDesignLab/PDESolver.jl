# Startup file for advection equation

push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/solver/euler"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/nl_solvers"))

using PDESolverCommon
using PumiInterface # pumi interface
using PdePumiInterface  # common mesh interface - pumi
using SummationByParts  # SBP operators
using EulerEquationMod
using ForwardDiff
using nl_solvers   # non-linear solvers
using ArrayViews

include(joinpath(Pkg.dir("PDESolver"),"src/solver/euler/output.jl"))  # printing results to files
include(joinpath(Pkg.dir("PDESolver"), "src/input/read_input.jl"))

function getResType(Tmsh::DataType, Tsbp::DataType, Tsol::DataType )
  # figure out what type eqn.res needs to be, taking into account
  # algorithmic differentiation of the different arguments
  # to support reverse mode, will need to consider the type of Tres as an input
  
  if Tsol <: DualNumbers.Dual # differentiating wrt eqn.q
    Tres = Tsol
  elseif  Tmsh <: DualNumbers.Dual  # differentiating wrt mesh.coords
    Tres = Tmsh
    Tsol = Tmsh  # derivative of coordinates will end up in eqn.q
  else  # no algorithmic differntiation
    Tres = Tsol
  end

  return Tres
end

#function runtest(flag::Int)
println("ARGS = ", ARGS)
println("size(ARGS) = ", size(ARGS))
opts = read_input(ARGS[1])
#opts = read_input("input_vals_channel2.jl")

# flag determines whether to calculate u, dR/du, or dR/dx (1, 2, or 3)
flag = opts["run_type"]

alpha_x = 1.0
alpha_y = 0.0

# timestepping parameters
delta_t = opts["delta_t"]
t_max = opts["t_max"]

order = opts["order"]  # order of accuracy

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

sbp = TriSBP{Tsbp}(degree=order)  # create linear sbp operator

# create mesh with 1 dofpernode
dmg_name = opts["dmg_name"]
smb_name = opts["smb_name"]
mesh = PumiMesh2{Tmsh}(dmg_name, smb_name, order, sbp, arg_dict ; dofpernode=1)

u_i = zeros(mesh.numDof)  # create u vector (current timestep)
u_i_1 = zeros(mesh.numDof) # u at next timestep

for i = 1:mesh.numEl
  for j = 1:sbp.numnodes
    x_j = coords[1,j,i]
    u_j = sin(x_j)
    dofnum = getGlobalNodeNumber(mesh, i, j)
    u_i[dofnum] = u_j
  end # end for j = 1:sbp.numnodes
end # end for i=1:mesh.numEl

mass_matrix = getMass(sbp, mesh)
println("mass_matrix = ", mass_matrix)
global int_advec = 1

u