# Name: startup.jl
# Description: startup script for solving an equation

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
mesh = PumiMesh2{Tmsh}(dmg_name, smb_name, order, sbp, arg_dict; dofpernode=4)
# TODO: input argument for dofpernode

# create euler equation
eqn = EulerData_{Tsol, Tres, 2, Tmsh}(mesh, sbp, opts)

# TODO: needs comment
init(mesh, sbp, eqn, opts)

res_vec = eqn.res_vec         # solution at previous timestep
q_vec = eqn.q_vec       # solution at current timestep

# calculate residual of some other function for res_reltol0
# TODO: add a boolean options here?
Relfunc_name = opts["Relfunc_name"]
if haskey(ICDict, Relfunc_name)
  println("\ncalculating residual for relative residual tolerance")
  Relfunc = ICDict[Relfunc_name]
  println("Relfunc = ", Relfunc)
  Relfunc(mesh, sbp, eqn, opts, q_vec)

#  println("eqn.q_vec = ", eqn.q_vec)
  res_real = zeros(mesh.numDof)
  tmp = calcResidual(mesh, sbp, eqn, opts, evalEuler, res_real)
#  println("res_real = \n", res_real)
#  println("eqn.res_vec = ", eqn.res_vec)
#  println("res_real = ", res_real)
  opts["res_reltol0"] = tmp
  println("res_reltol0 = ", tmp)
  
#  writedlm("relfunc_res.dat", eqn.res)
#  writedlm("relfunc_resvec.dat", res_real)
  saveSolutionToMesh(mesh, res_real)
  writeVisFiles(mesh, "solution_relfunc")
end

# populate u0 with initial condition
println("\nEvaluating initial condition")
ICfunc_name = opts["IC_name"]
ICfunc = ICDict[ICfunc_name]
println("ICfunc = ", ICfunc)
ICfunc(mesh, sbp, eqn, opts, q_vec)

# TODO: cleanup 20151009 start

if opts["calc_error"]
  println("\ncalculating error of file ", opts["calc_error_infname"], 
          " compared to initial condition")
  vals = readdlm(opts["calc_error_infname"])
  @assert length(vals) == mesh.numDof

  err_vec = vals - eqn.q_vec
  err = calcNorm(eqn, err_vec)
  outname = opts["calc_error_outfname"]
  println("printint err = ", err, " to file ", outname)
  f = open(outname, "w")
  println(f, err)
  close(f)
end

if opts["calc_trunc_error"]  # calculate truncation error
  println("\nCalculating residual for truncation error")
  res_real = zeros(mesh.numDof)
  tmp = calcResidual(mesh, sbp, eqn, opts, evalEuler, res_real)

  f = open("error_trunc.dat", "w")
  println(f, tmp)
  close(f)
end

if opts["perturb_ic"]
  println("\nPerturbing initial condition")
  perturb_mag = opts["perturb_mag"]
  for i=1:mesh.numDof
    q_vec[i] += perturb_mag*rand()
  end
end

res_vec_exact = deepcopy(q_vec)

rmfile("IC.dat")
writedlm("IC.dat", real(q_vec))
saveSolutionToMesh(mesh, q_vec)

writeVisFiles(mesh, "solution_ic")

# initialize some variables in nl_solvers module
initializeTempVariables(mesh)

#------------------------------------------------------------------------------
include("checkEigenValues.jl")
# include("artificialViscosity.jl")

# Calculate the recommended delta t
res_0 = zeros(eqn.res_vec)
res_0_norm = calcResidual(mesh, sbp, eqn, opts, evalEuler, res_0)
CFLMax = 1      # Maximum Recommended CFL Value
Dt = zeros(mesh.numNodesPerElement,mesh.numEl) # Array of all possible delta t

for i = 1:mesh.numEl
  for j = 1:mesh.numNodesPerElement
    h = 1/sqrt(mesh.jac[j,i])
    velocities = zeros(2) # Nodal velocities
    velocities[1] = eqn.q[2,j,i]/eqn.q[1,j,i]
    velocities[2] = eqn.q[3,j,i]/eqn.q[1,j,i] 
    vmax = norm(velocities)
    q = view(eqn.q,:,j,i)
    T = (q[4] - 0.5*(q[2]*q[2] + q[3]*q[3])/q[1])*(1/(q[1]*eqn.params.cv))
    c = sqrt(eqn.params.gamma*eqn.params.R*T) # Speed of sound
    Dt[j,i] = CFLMax*h/(vmax + c)
  end
end
RecommendedDT = minimum(Dt)
println("Recommended delta t = ", RecommendedDT)

#elementEigenValues(mesh, sbp, eqn)

#------------------------------------------------------------------------------

# call timestepper
if opts["solve"]
  
  if flag == 1 # normal run
   rk4(evalEuler, delta_t, t_max, mesh, sbp, eqn, opts, 
       res_tol=opts["res_abstol"])

   println("finish rk4")
   printSolution("rk4_solution.dat", eqn.res_vec)
  # println("rk4 @time printed above")
  elseif flag == 2 # forward diff dR/du

    # define nested function
    function dRdu_rk4_wrapper(u_vals::AbstractVector, res_vec::AbstractVector)
      eqn.q_vec = u_vals
      eqn.q_vec = res_vec
      rk4(evalEuler, delta_t, t_max, mesh, sbp, eqn)
      return nothing
    end

    # use ForwardDiff package to generate function that calculate jacobian
    calcdRdu! = forwarddiff_jacobian!(dRdu_rk4_wrapper, Float64, 
                fadtype=:dual, n = mesh.numDof, m = mesh.numDof)

    jac = zeros(Float64, mesh.numDof, mesh.numDof)  # array to be populated
    calcdRdu!(eqn.q_vec, jac)

  elseif flag == 3 # calculate dRdx

    # dRdx here

  elseif flag == 4 || flag == 5
    @time newton(evalEuler, mesh, sbp, eqn, opts, itermax=opts["itermax"], 
                 step_tol=opts["step_tol"], res_abstol=opts["res_abstol"], 
                 res_reltol=opts["res_reltol"], res_reltol0=opts["res_reltol0"])

    println("total solution time printed above")
    printSolution("newton_solution.dat", eqn.res_vec)

  elseif flag == 6
    newton_check(evalEuler, mesh, sbp, eqn, opts)
    vals = abs(real(eqn.res_vec))  # remove unneded imaginary part
    saveSolutionToMesh(mesh, vals)
    writeVisFiles(mesh, "solution_error")
    printBoundaryEdgeNums(mesh)
    printSolution(mesh, vals)

  elseif flag == 7
    jac_col = newton_check(evalEuler, mesh, sbp, eqn, opts, 1)
    writedlm("solution.dat", jac_col)

  elseif flag == 8
    jac_col = newton_check_fd(evalEuler, mesh, sbp, eqn, opts, 1)
    writedlm("solution.dat", jac_col)

  end       # end of if/elseif blocks checking flag

  if opts["write_finalsolution"]
    writedlm("solution_final.dat", real(eqn.q_vec))
  end

  if opts["write_finalresidual"]
    writedlm("residual_final.dat", real(eqn.res_vec))
  end


##### Do postprocessing ######
println("\nDoing postprocessing")
  if flag == 1
      res_vec_diff = res_vec - res_vec_exact
      step = q_vec - res_vec_exact
      step_norm = norm(step)/mesh.numDof
      println("step_norm = ", step_norm)
      res_vec_norm = calcNorm(eqn, res_vec)
      #res_vec_side_by_side = [res_vec_exact  res_vec]

      #=
      println("\n\n\n")
      println("res_vec_diff: \n")
      for i=1:size(res_vec_diff)[1]
	println(res_vec_diff[i,:])
      end
      println("res_vec_side_by_side: \n")
      for i=1:size(res_vec_side_by_side)[1]
	println(i, " ", res_vec_side_by_side[i,:])
      end
      =#
      println("res_vec_norm: \n",res_vec_norm,"\n")

  end

      saveSolutionToMesh(mesh, real(eqn.q_vec))
      printSolution(mesh, real(eqn.q_vec))
      printCoordinates(mesh)
      writeVisFiles(mesh, "solution_done")

end  # end if (opts[solve])

#runtest(1)
