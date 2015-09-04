# startup script for solving an equation

# push!(LOAD_PATH, "/users/creanj/julialib_fork/PUMI.jl/src")
# push!(LOAD_PATH, "../../equation")
# push!(LOAD_PATH, "../../common")
# push!(LOAD_PATH, "/users/creanj/.julia/v0.4/PDESolver/src/solver/euler")
#push!(LOAD_PATH, "../../../../PUMI")



push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/solver/euler"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/nl_solvers"))
#include("complexify.jl")
using PDESolverCommon
using PumiInterface # pumi interface
using PdePumiInterface  # common mesh interface - pumi
using SummationByParts  # SBP operators
using EulerEquationMod
using ForwardDiff
using nl_solvers   # non-linear solvers

#include(joinpath(Pkg.dir("PDESolver"),"src/nl_solvers/rk4.jl"))  # timestepping

#include(joinpath(Pkg.dir("PDESolver"), "src/nl_solvers/newton_fd.jl"))  # timestepping
include(joinpath(Pkg.dir("PDESolver"),"src/solver/euler/output.jl"))  # printing results to files
include(joinpath(Pkg.dir("PDESolver"), "src/input/read_input.jl"))
#include(joinpath(Pkg.dir("PDESolver"), "src/tools/misc.jl"))  # assorted utilities

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

flag = opts["run_type"]
# flag determines whether to calculate u, dR/du, or dR/dx (1, 2, or 3)
# timestepping parameters
#delta_t = 0.005
delta_t = opts["delta_t"]
#t_max = 0.025
#t_max = 5.00
t_max = opts["t_max"]
#t_max = 1.0
#order = 1  # order of accuracy
order = opts["order"]


#set_bigfloat_precision(80)  # use 128 bit floats
# types of the mesh, SBP, Equation objects
if flag == 1 || flag == 8  # normal run
  Tmsh = Float64
  Tsbp = Float64
  Tsol = Float64
  Tres = Float64
elseif flag == 2  # calculate dR/du
  Tmsh = Float64
  Tsbp = Float64
  Tsol = Dual{Float64}
  Tres = Dual{Float64}
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

#  Tsol = BigFloat
#  Tres = BigFloat
elseif flag == 5  # use complex step dR/du

  Tmsh = Float64
  Tsbp = Float64
  Tsol = Complex128
  Tres = Complex128
#=
  Tmsh = Float64
  Tsbp = Float64
  Tsol = Complex{BigFloat}
  Tres = Complex{BigFloat}
=#
elseif flag == 6 || flag == 7  # evaluate residual error and print to paraview
  Tmsh = Float64
  Tsbp = Float64
  Tsol = Complex128
  Tres = Complex128

end


# create operator
sbp = TriSBP{Tsbp}(degree=order)  # create linear sbp operator

# create mesh
#dmg_name = ".null"
#smb_name = "../../mesh_files/quarter_vortex10l.smb"

#dmg_name = "../../mesh_files/vortex.dmg"
#smb_name = "../../mesh_files/vortex.smb"

dmg_name = opts["dmg_name"]
smb_name = opts["smb_name"]

#smb_name = "../../mesh_files/quarter_vortex33l.smb"
#smb_name = "../../mesh_files/quarter_vortex1000l.smb"
#smb_name = "../../mesh_files/quarter_vortex8l.smb"
#smb_name = "../../mesh_files/tri30l.smb"
mesh = PumiMesh2{Tmsh}(dmg_name, smb_name, order, sbp, arg_dict ; dofpernode=4)  #create linear mesh with 4 dof per node





# create euler equation
eqn = EulerData_{Tsol, Tres, 2, Tmsh}(mesh, sbp, opts)
#eqn = EulerEquation{Tsol}(mesh, sbp, Float64)

println("size(q) = ", size(eqn.q))

#SL0 = zeros(Tsol, mesh.numDof)  # solution at previous timestep
#SL = zeros(Tsol, mesh.numDof) # solution at current timestep
SL = eqn.SL
SL0 = eqn.SL0

# get BC functors
getBCFunctors(mesh, sbp, eqn, opts)

# calculate residual of some other function for res_reltol0
Relfunc_name = opts["Relfunc_name"]
if haskey(ICDict, Relfunc_name)
  Relfunc = ICDict[Relfunc_name]
  println("Relfunc = ", Relfunc)
  Relfunc(mesh, sbp, eqn, opts, SL0)

  res_real = zeros(mesh.numDof)
  println("calculating residual for relative residual tolerance")
  tmp = calcResidual(mesh, sbp, eqn, opts, evalEuler, res_real)

  opts["res_reltol0"] = tmp
  println("res_reltol0 = ", tmp)
  print("\n")
end



# populate u0 with initial condition
ICfunc_name = opts["IC_name"]
ICfunc = ICDict[ICfunc_name]
println("ICfunc = ", ICfunc)
ICfunc(mesh, sbp, eqn, opts, SL0)
# ICZero(mesh, sbp, eqn, SL0)
# ICLinear(mesh, sbp, eqn, SL0)
# ICIsentropicVortex(mesh, sbp, eqn, SL0)
#ICRho1E2(mesh, sbp, eqn, SL0)
#ICRho1E2U3(mesh, sbp, eqn, SL0)

#ICVortex(mesh, sbp, eqn, SL0)
#ICIsentropicVortex(mesh, sbp, eqn, SL0)

if opts["calc_error"]
  println("calculating error of file ", opts["calc_error_infname"], " compared to initial condition")
  vals = readdlm(opts["calc_error_infname"])
  @assert length(vals) == mesh.numDof

  err_vec = vals - eqn.SL0
  err = calcNorm(eqn, err_vec)
#  err = norm(vals - eqn.SL0)/mesh.numDof
  outname = opts["calc_error_outfname"]
  println("printint err = ", err, " to file ", outname)
  f = open(outname, "w")
  println(f, err)
  close(f)
  print("\n")
end

if opts["calc_trunc_error"]  # calculate truncation error

  res_real = zeros(mesh.numDof)
  println("calculating residual for truncation error")
  tmp = calcResidual(mesh, sbp, eqn, opts, evalEuler, res_real)

  tmp = 0
  # calculate a norm
  for i=1:mesh.numDof
    tmp = res_real[i]*eqn.Minv[i]*res_real[i]
  end

  tmp = sqrt(tmp)

  f = open("error_trunc.dat", "w")
  println(f, tmp)
  close(f)
end




if opts["perturb_ic"]
  perturb_mag = opts["perturb_mag"]
  for i=1:mesh.numDof
    SL0[i] += perturb_mag*rand()
  end
end


SL_exact = deepcopy(SL0)

rmfile("IC.dat")
writedlm("IC.dat", real(SL0))
saveSolutionToMesh(mesh, SL0)

writeVisFiles(mesh, "solution_ic")

# initialize some variables in nl_solvers module
initializeTempVariables(mesh)

# call timestepper
if opts["solve"]
  
  if flag == 1 # normal run
   rk4(evalEuler, delta_t, t_max, mesh, sbp, eqn, opts, res_tol=opts["res_tol"])
   println("finish rk4")
   printSolution("rk4_solution.dat", eqn.SL)
  # println("rk4 @time printed above")
  elseif flag == 2 # forward diff dR/du

    # define nested function
    function dRdu_rk4_wrapper(u_vals::AbstractVector, SL::AbstractVector)
      eqn.SL0 = u_vals
      eqn.SL0 = SL
      rk4(evalEuler, delta_t, t_max, mesh, sbp, eqn)
      return nothing
    end

    # use ForwardDiff package to generate function that calculate jacobian
    calcdRdu! = forwarddiff_jacobian!(drDu_rk4_wrapper, Float64, fadtype=:dual, n = mesh.numDof, m = mesh.numDof)

    jac = zeros(Float64, mesh.numDof, mesh.numDof)  # array to be populated
    calcdRdu!(eqn.SL0, jac)

  elseif flag == 3 # calculate dRdx

    # dRdx here

  elseif flag == 4 || flag == 5
    @time newton(evalEuler, mesh, sbp, eqn, opts, itermax=opts["itermax"], step_tol=opts["step_tol"], res_abstol=opts["res_abstol"], res_reltol=opts["res_reltol"], res_reltol0=opts["res_reltol0"])
    println("total solution time printed above")
    printSolution("newton_solution.dat", eqn.SL)
#=
  elseif flag == 5
  #  newton_complex(evalEuler, mesh, sbp, eqn, opts, itermax=200, step_tol=1e-6, res_tol=1e-8)

    newton_complex(evalEuler, mesh, sbp, eqn, opts, itermax=opts["itermax"], step_tol=opts["step_tol"], res_tol=opts["res_tol"])
=#
  elseif flag == 6
    newton_check(evalEuler, mesh, sbp, eqn, opts)

    vals = abs(real(eqn.SL))  # remove unneded imaginary part
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

  end



  if opts["write_finalsolution"]
    writedlm("solution_final.dat", real(eqn.SL0))
  end

  if opts["write_finalresidual"]
    writedlm("residual_final.dat", real(eqn.SL))
  end



  if flag == 1


      SL_diff = SL - SL_exact
      step = SL0 - SL_exact
      step_norm = norm(step)/mesh.numDof
      println("step_norm = ", step_norm)
      SL_norm = norm(SL)/mesh.numDof
      #SL_side_by_side = [SL_exact  SL]

      #=
      println("\n\n\n")
      println("SL_diff: \n")
      for i=1:size(SL_diff)[1]
	println(SL_diff[i,:])
      end
      println("SL_side_by_side: \n")
      for i=1:size(SL_side_by_side)[1]
	println(i, " ", SL_side_by_side[i,:])
      end
      =#
      println("SL_norm: \n",SL_norm,"\n")


  end

      saveSolutionToMesh(mesh, real(eqn.SL0))
      printSolution(mesh, real(eqn.SL0))
      printCoordinates(mesh)
      writeVisFiles(mesh, "solution_done")

end  # end if (opts[solve])
  #end


#runtest(1)
