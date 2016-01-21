# Startup file for 1 dof advection equation

push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/solver/advection"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/NonlinearSolvers"))

using ODLCommonTools
using PumiInterface # pumi interface
using PdePumiInterface  # common mesh interface - pumi
using SummationByParts  # SBP operators
using AdvectionEquationMod # Advection equation module
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

# create linear mesh with 4 dof per node
mesh = PumiMesh2{Tmsh}(dmg_name, smb_name, order, sbp, opts; dofpernode=4, 
                       coloring_distance=opts["coloring_distance"])
if opts["jac_type"] == 3 || opts["jac_type"] == 4
  pmesh = PumiMesh2Preconditioning(mesh, sbp, opts; coloring_distance=opts["coloring_distance_prec"])
else
  pmesh = mesh
end

# Create advection equation object
eqn = AdvectionData_{Tsol, Tres, 2, Tmsh}(mesh, sbp, opts)

u_vec = eqn.u_vec

# Initialize the advection equation
init(mesh, sbp, eqn, opts)

# Populate with initial conditions
println("\nEvaluating initial condition")
ICfunc_name = opts["IC_name"]
ICfunc = ICDict[ICfunc_name]
println("ICfunc = ", ICfunc)
ICfunc(mesh, sbp, eqn, opts, u_vec) 
println("finished initializing u")

if opts["calc_error"]
  println("\ncalculating error of file ", opts["calc_error_infname"], 
          " compared to initial condition")
  vals = readdlm(opts["calc_error_infname"])
  @assert length(vals) == mesh.numDof

  err_vec = vals - eqn.u_vec
  err = calcNorm(eqn, err_vec)
  outname = opts["calc_error_outfname"]
  println("printint err = ", err, " to file ", outname)
  f = open(outname, "w")
  println(f, err)
  close(f)
end

if opts["calc_trunc_error"]  # calculate truncation error
  println("\nCalculating residual for truncation error")
  tmp = calcResidual(mesh, sbp, eqn, opts, evalAdvection)

  f = open("error_trunc.dat", "w")
  println(f, tmp)
  close(f)
end

if opts["perturb_ic"]
  println("\nPerturbing initial condition")
  perturb_mag = opts["perturb_mag"]
  for i=1:mesh.numDof
    u_vec[i] += perturb_mag*rand()
  end
end

res_vec_exact = deepcopy(u_vec)

rmfile("IC.dat")
writedlm("IC.dat", real(u_vec))
saveSolutionToMesh(mesh, u_vec)

global int_advec = 1

#------------------------------------------------------------------------------
#=
# Calculate the recommended delta t
CFLMax = 1      # Maximum Recommended CFL Value
const alpha_x = 1.0 # advection velocity in x direction
const alpha_y = 0.0 # advection velocity in y direction
Dt = zeros(mesh.numNodesPerElement,mesh.numEl) # Array of all possible delta t
for i = 1:mesh.numEl
  for j = 1:mesh.numNodesPerElement
    h = 1/sqrt(mesh.jac[j,i])
    V = sqrt((alpha_x.^2) + (alpha_y.^2))
    Dt[j,i] = CFLMax*h/(V)
  end # end for j = 1:mesh.numNodesPerElement
end   # end for i = mesh.numEl
RecommendedDT = minimum(Dt)
println("Recommended delta t = ", RecommendedDT) =#
#------------------------------------------------------------------------------

t = 0.0
# evalAdvection(mesh, sbp, eqn, opts, t)

# Now put it into the RK4 solver
rk4(evalAdvection, delta_t, t_max, mesh, sbp, eqn, opts, res_tol=opts["res_abstol"])