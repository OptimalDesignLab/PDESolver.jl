# Name: startup.jl
# Description: startup script for solving an equation

push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/solver/euler"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/NonlinearSolvers"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Debugging"))
#include("complexify.jl")   # TODO: include location needs to be reconsidered

using ODLCommonTools
using PdePumiInterface  # common mesh interface - pumi
using SummationByParts  # SBP operators
using EulerEquationMod
using ForwardDiff
using NonlinearSolvers   # non-linear solvers
using ArrayViews
#using Debugging   # some debugging utils.
include(joinpath(Pkg.dir("PDESolver"),"src/solver/euler/output.jl"))  # printing results to files
include(joinpath(Pkg.dir("PDESolver"), "src/input/read_input.jl"))

#function runtest()
println("ARGS = ", ARGS)
println("size(ARGS) = ", size(ARGS))
opts = read_input(ARGS[1])
#opts = read_input("input_vals_channel2.jl")

# flag determines whether to calculate u, dR/du, or dR/dx (1, 2, or 3)
flag = opts["run_type"]

# timestepping parameters
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

# record these choices in the dictionary
opts["Tmsh"] = Tmsh
opts["Tsbp"] = Tsbp
opts["Tsol"] = Tsol
opts["Tres"] = Tres

# create SBP object
println("\nConstructing SBP Operator")
sbp = TriSBP{Tsbp}(degree=order)  # create linear sbp operator

dmg_name = opts["dmg_name"]
smb_name = opts["smb_name"]
Tdim = opts["dimensions"]

# create linear mesh with 4 dof per node
mesh = PumiMesh2{Tmsh}(dmg_name, smb_name, order, sbp, opts; dofpernode=4, coloring_distance=opts["coloring_distance"])

if opts["jac_type"] == 3 || opts["jac_type"] == 4
  pmesh = PumiMesh2Preconditioning(mesh, sbp, opts; coloring_distance=opts["coloring_distance_prec"])
else
  pmesh = mesh
end

# TODO: input argument for dofpernode

# create euler equation
var_type = opts["variable_type"]
eqn = EulerData_{Tsol, Tres, 2, Tmsh, var_type}(mesh, sbp, opts)

# initialize physics module and populate any fields in mesh and eqn that
# depend on the physics module
init(mesh, sbp, eqn, opts, pmesh)

#delta_t = opts["delta_t"]   # delta_t: timestep for RK


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
 
  if var_type == :entropy
    println("converting to entropy variables")
    for i=1:mesh.numDofPerNode:mesh.numDof
      q_view = view(q_vec, i:(i+mesh.numDofPerNode-1))
      convertFromNaturalToWorkingVars(eqn.params, q_view, q_view)
    end
  end
#  println("eqn.q_vec = ", eqn.q_vec)
  tmp = calcResidual(mesh, sbp, eqn, opts, evalEuler)
  res_real = real(eqn.res_vec)
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

if var_type == :entropy
  for i=1:mesh.numDofPerNode:mesh.numDof
    q_view = view(q_vec, i:(i+mesh.numDofPerNode-1))
    convertFromNaturalToWorkingVars(eqn.params, q_view, q_view)
  end
end

# TODO: cleanup 20151009 start

if opts["calc_error"]
  println("\ncalculating error of file ", opts["calc_error_infname"], 
          " compared to initial condition")
  vals = readdlm(opts["calc_error_infname"])
  @assert length(vals) == mesh.numDof

  err_vec = vals - eqn.q_vec
  err = calcNorm(eqn, err_vec)

  # calculate avg mesh size
  jac_3d = reshape(mesh.jac, 1, mesh.numNodesPerElement, mesh.numEl)
  jac_vec = zeros(Tmsh, mesh.numNodes)
  EulerEquationMod.assembleArray(mesh, sbp, eqn, opts, jac_3d, jac_vec)
  # scale by the minimum distance between nodes on a reference element
  # this is a bit of an assumption, because for distorted elements this
  # might not be entirely accurate
  println("mesh.min_node_distance = ", mesh.min_node_dist)
  h_avg = sum(1./sqrt(jac_vec))/length(jac_vec)
#  println("h_avg = ", h_avg)
  h_avg *= mesh.min_node_dist
#  println("h_avg = ", h_avg)

  outname = opts["calc_error_outfname"]
  println("printint err = ", err, " to file ", outname)
  f = open(outname, "w")
  println(f, err, " ", h_avg)
  close(f)
end

if opts["calc_trunc_error"]  # calculate truncation error
  println("\nCalculating residual for truncation error")
  tmp = calcResidual(mesh, sbp, eqn, opts, evalEuler)

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

wave_speed = EulerEquationMod.calcMaxWaveSpeed(mesh, sbp, eqn, opts)
println("max wave speed = ", wave_speed)
delta_t = opts["CFL"]*opts["mesh_size"]/wave_speed
println("for a CFL of ", opts["CFL"], " delta_t = ", delta_t)
opts["delta_t"] = delta_t

#DEBUGGING
if opts["test_GLS2"]
  calcResidual(mesh, sbp, eqn, opts, evalEuler)
end

#------------------------------------------------------------------------------
#=
include("checkEigenValues.jl")
# include("artificialViscosity.jl")
# include("SUPG.jl")

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

#=
FluxJacobian(mesh, sbp, eqn) # Calculate the euler flux jacobian  
tau = zeros(Tsol, mesh.numNodesPerElement, mesh.numEl) # Stabilization term
calcStabilizationTerm(mesh, sbp, eqn, tau) =#

#elementEigenValues(mesh, sbp, eqn)
#SUPG2(mesh, sbp, eqn)
# residualComparison(mesh, sbp, eqn, opts)

=#

#res_0_norm = calcResidual(mesh, sbp, eqn, opts, evalEuler)

#------------------------------------------------------------------------------

# call timestepper
if opts["solve"]
  
  if flag == 1 # normal run
   @time rk4(evalEuler, delta_t, t_max, mesh, sbp, eqn, opts, res_tol=opts["res_abstol"], real_time=opts["real_time"])
#   @time rk4(evalEuler, delta_t, t_max, eqn.q_vec, eqn.res_vec, 
#              (mesh, sbp, eqn), opts, majorIterationCallback=eqn.majorIterationCallback, 
#              res_tol=opts["res_abstol"], real_time=opts["real_time"])

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
    @time newton(evalEuler, mesh, sbp, eqn, opts, pmesh, itermax=opts["itermax"], 
                 step_tol=opts["step_tol"], res_abstol=opts["res_abstol"], 
                 res_reltol=opts["res_reltol"], res_reltol0=opts["res_reltol0"])

    printSolution("newton_solution.dat", eqn.res_vec)

  elseif flag == 6
    @time newton_check(evalEuler, mesh, sbp, eqn, opts)
    vals = abs(real(eqn.res_vec))  # remove unneded imaginary part
    saveSolutionToMesh(mesh, vals)
    writeVisFiles(mesh, "solution_error")
    printBoundaryEdgeNums(mesh)
    printSolution(mesh, vals)

  elseif flag == 7
    @time jac_col = newton_check(evalEuler, mesh, sbp, eqn, opts, 1)
    writedlm("solution.dat", jac_col)

  elseif flag == 8
    @time jac_col = newton_check_fd(evalEuler, mesh, sbp, eqn, opts, 1)
    writedlm("solution.dat", jac_col)

  end       # end of if/elseif blocks checking flag

  println("total solution time printed above")
  # evaluate residual at final q value
  eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
  evalEuler( mesh, sbp, eqn, opts, eqn.params.t)

  eqn.res_vec[:] = 0.0
  eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)


  if opts["write_finalsolution"]
    println("writing final solution")
    writedlm("solution_final.dat", real(eqn.q_vec))
  end

  if opts["write_finalresidual"]
    writedlm("residual_final.dat", real(eqn.res_vec))
  end


  ##### Do postprocessing ######
  println("\nDoing postprocessing")

  if opts["do_postproc"]
    exfname = opts["exact_soln_func"]
    if haskey(ICDict, exfname)
      exfunc = ICDict[exfname]
      q_exact = zeros(Tsol, mesh.numDof)
      exfunc(mesh, sbp, eqn, opts, q_exact)
#    if var_type == :entropy
#      println("converting to entropy variables")
#      for i=1:mesh.numDofPerNode:mesh.numDof
#        q_view = view(q_vec, i:(i+mesh.numDofPerNode-1))
#        convertFromNaturalToWorkingVars(eqn.params, q_view, q_view)
#      end
#    end

      q_diff = eqn.q_vec - q_exact
      diff_norm = calcNorm(eqn, q_diff)
      discrete_norm = norm(q_diff/length(q_diff))

      println("solution error norm = ", diff_norm)
      println("solution discrete L2 norm = ", discrete_norm)

      # print to file
      outname = opts["calc_error_outfname"]
      f = open(outname, "w")
      println(f, mesh.numEl, " ", diff_norm, " ", discrete_norm)
      close(f)

    end
  end

  saveSolutionToMesh(mesh, real(eqn.q_vec))
  printSolution(mesh, real(eqn.q_vec))
  printCoordinates(mesh)
  writeVisFiles(mesh, "solution_done")

end  # end if (opts[solve])

#  return mesh, sbp, eqn, opts
#end  # end function
#runtest(1)
