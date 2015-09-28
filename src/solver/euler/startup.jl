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
using ArrayViews

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

#------------------------------------------------------------------------------
include("checkEigenValues.jl")


# Calculate the recommended delta t

res_0 = zeros(eqn.SL)

res_0_norm = calcResidual(mesh, sbp, eqn, opts, evalEuler, res_0)
# println("eqn.SL \n", eqn.SL)

CFLMax = 1 # Maximum Recommended CFL Value
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
pert = 1e-6
#println(eqn.SL0)

numDofPerElement = mesh.numDofPerNode*mesh.numNodesPerElement
elem_eigen_vals = zeros(numDofPerElement)
elem_count = 0 # Count the number of elements with negative eigen values

for i = 1:mesh.numEl
  elem_jac = zeros(numDofPerElement,numDofPerElement)
  elemSL = zeros(numDofPerElement) # Element wise residual
  elemSL0 = zeros(elemSL)  # 1D array containing element wise conservative variables
  elemRes0 = zeros(elemSL) # Collects the original residual
  entry_orig = zero(eltype(eqn.SL0)) # stores one element of original SL) at a time
  orig_SL = copy(eqn.SL)

  orig_SL0 = copy(eqn.SL)
  for j = 1:mesh.numNodesPerElement # Populating elemRes0 & elemSL0
    for k = 1:mesh.numDofPerNode
      elemRes0[(j-1)*mesh.numDofPerNode + k] = eqn.res[k,j,i]
      elemSL0[(j-1)*mesh.numDofPerNode + k] = eqn.q[k,j,i]
    end
  end

  # println("elemRes0 - eqn.SL = \n", elemRes0 - eqn.SL)
  # println("elemSL0 - eqn.SL0 = \n", elemSL0 - eqn.SL0)

  for j = 1:numDofPerElement # Perturb all dofs in the element  
    if j == 1
      entry_orig = elemSL0[j]
      elemSL0[j] += pert  # Perturb the jth element of q
    else
      elemSL0[j-1] = entry_orig
      entry_orig = elemSL0[j]
      elemSL0[j] += pert
    end
    #println("elemSL0 - eqn.SL0 = \n", elemSL0 - eqn.SL0)
    # Disassembling SL0 to q at the element level
    for l = 1:mesh.numNodesPerElement
      for k = 1:mesh.numDofPerNode
        eqn.q[k,l,i] = elemSL0[(l-1)*mesh.numDofPerNode + k]
      end
    end
    evalEuler(mesh, sbp, eqn, opts) # Evaluate function with perturbed q
    
    # Populate elemSL with perturbed residuals
    for l = 1:mesh.numNodesPerElement
      for k = 1:mesh.numDofPerNode
        elemSL[(l-1)*mesh.numDofPerNode + k] = eqn.res[k,l,i]
      end
    end
    
    fill!(eqn.SL, 0.0)
    eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.res,  eqn.SL)
    # println("elemSL - eqn.SL = \n", elemSL - eqn.SL)
    # println("orig_SL - eqn.SL = \n", orig_SL - eqn.SL)

    nl_solvers.calcJacRow(unsafe_view(elem_jac, :, j), elemRes0, elemSL, pert)
  end  # End for j = = 1:numDofPerElement
  
  # println("elemSL = \n", elemSL)
  eqn.q[mesh.numDofPerNode,mesh.numNodesPerElement,i] = entry_orig
  
  elem_eigen_vals = eigvals(elem_jac)
  checkEigenValues(elem_eigen_vals, i, elem_count)
  
  #= println("eigen values of element ", i)
  counter = 0
  for kappa = 1:length(elem_eigen_vals)
    if elem_eigen_vals[kappa] < 0 && abs(elem_eigen_vals[kappa]) > 1e-12
      counter += 1
    end
  end =#

  #println(elem_eigen_vals, '\n', '\n')
  #println("number of negative element eigen values = ", counter)
  
end  

println("Number of elements with negative eigen values = ", elem_count)

res_0_norm = calcResidual(mesh, sbp, eqn, opts, evalEuler, res_0)
jac = zeros(Float64, mesh.numDof, mesh.numDof)
res_0 = copy(eqn.SL)
nl_solvers.calcJacFD(mesh, sbp, eqn, opts, evalEuler, res_0, pert, jac)
println("Jacobian successfully computed") 


eigen_values = eigvals(jac)
counter = 0
for i = 1:length(eigen_values)
  if eigen_values[i] < 0 && abs(eigen_values[i]) > 1e-12
    counter += 1
  end
end
#println(eigen_values)
println("number of negative eigen values for assembled residual = ", counter)
println("Number of eigen values of assembled residual = ", length(eigen_values))
printSolution("eigen_values.dat", eigen_values)  


#------------------------------------------------------------------------------

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
    calcdRdu! = forwarddiff_jacobian!(dRdu_rk4_wrapper, Float64, fadtype=:dual, n = mesh.numDof, m = mesh.numDof)

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
