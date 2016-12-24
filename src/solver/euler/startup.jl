# Name: startup.jl
# Description: startup script for solving an equation

push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/solver/euler"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/NonlinearSolvers"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Debugging"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Utils"))
#include("complexify.jl")   # TODO: include location needs to be reconsidered

using ODLCommonTools
using PdePumiInterface  # common mesh interface - pumi
using SummationByParts  # SBP operators
using EulerEquationMod
using ForwardDiff
using NonlinearSolvers   # non-linear solvers
using ArrayViews
using Utils
import ODLCommonTools.sview
using MPI

#using Debugging   # some debugging utils.
include(joinpath(Pkg.dir("PDESolver"),"src/solver/euler/output.jl"))  # printing results to files
include(joinpath(Pkg.dir("PDESolver"), "src/input/read_input.jl"))


if !MPI.Initialized()
  MPI.Init()
end

#function runtest()
println("ARGS = ", ARGS)
println("size(ARGS) = ", size(ARGS))
opts = read_input(ARGS[1])
#opts = read_input("input_vals_channel2.jl")

# flag determines whether to calculate u, dR/du, or dR/dx (1, 2, or 3)
t_max = opts["t_max"]
flag = opts["run_type"]

# timestepping parameters
order = opts["order"]       # order of accuracy
Tdim = opts["dimensions"]
dofpernode = Tdim + 2

sbp, mesh, pmesh, Tsol, Tres, Tmsh, mesh_time = createMeshAndOperator(opts, dofpernode)

myrank = mesh.myrank
# TODO: input argument for dofpernode

# create euler equation
var_type = opts["variable_type"]
eqn = EulerData_{Tsol, Tres, Tdim, Tmsh, var_type}(mesh, sbp, opts)

# initialize physics module and populate any fields in mesh and eqn that
# depend on the physics module
init(mesh, sbp, eqn, opts, pmesh)

#delta_t = opts["delta_t"]   # delta_t: timestep for RK


res_vec = eqn.res_vec 
q_vec = eqn.q_vec       # solution at current timestep

# calculate residual of some other function for res_reltol0
# TODO: add a boolean options here?
Relfunc_name = opts["Relfunc_name"]
if haskey(ICDict, Relfunc_name)
  @mpi_master println("\ncalculating residual for relative residual tolerance")
  Relfunc = ICDict[Relfunc_name]
  @mpi_master println("Relfunc = ", Relfunc)
  Relfunc(mesh, sbp, eqn, opts, q_vec)
 
  if var_type == :entropy
    @mpi_master println("converting to entropy variables")
    for i=1:mesh.numDofPerNode:mesh.numDof
      q_view = sview(q_vec, i:(i+mesh.numDofPerNode-1))
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
@mpi_master println("\nEvaluating initial condition")
ICfunc_name = opts["IC_name"]
ICfunc = ICDict[ICfunc_name]
@mpi_master println("ICfunc = ", ICfunc)
ICfunc(mesh, sbp, eqn, opts, q_vec)

if var_type == :entropy
  for i=1:mesh.numDofPerNode:mesh.numDof
    q_view = sview(q_vec, i:(i+mesh.numDofPerNode-1))
    convertFromNaturalToWorkingVars(eqn.params, q_view, q_view)
  end
end

# TODO: cleanup 20151009 start

if opts["calc_error"]
  @mpi_master println("\ncalculating error of file ", opts["calc_error_infname"], 
          " compared to initial condition")
  vals = readdlm(opts["calc_error_infname"])
  @assert length(vals) == mesh.numDof

  err_vec = abs(vals - eqn.q_vec)
  err = calcNorm(eqn, err_vec)

  # calculate avg mesh size
  h_avg = calcMeshH(mesh, sbp, eqn, opts)

  @mpi_master begin
    outname = opts["calc_error_outfname"]
    println("printint err = ", err, " to file ", outname)
    f = open(outname, "w")
    println(f, err, " ", h_avg)
    close(f)
  end

  # write visualization
  saveSolutionToMesh(mesh, vec(err_vec))
  writeVisFiles(mesh, "solution_error")
end

if opts["calc_trunc_error"]  # calculate truncation error
  @mpi_master println("\nCalculating residual for truncation error")
  tmp = calcResidual(mesh, sbp, eqn, opts, evalEuler)

  @mpi_master begin
    f = open("error_trunc.dat", "w")
    println(f, tmp)
    close(f)
  end
end

if opts["perturb_ic"]
  @mpi_master println("\nPerturbing initial condition")
  perturb_mag = opts["perturb_mag"]
  for i=1:mesh.numDof
    q_vec[i] += perturb_mag*rand()
  end
end

res_vec_exact = deepcopy(q_vec)

rmfile("IC_$myrank.dat")
writedlm("IC_$myrank.dat", real(q_vec))
saveSolutionToMesh(mesh, q_vec)

writeVisFiles(mesh, "solution_ic")
writedlm("solution_ic.dat", real(eqn.q_vec))
if opts["calc_dt"]
  wave_speed = EulerEquationMod.calcMaxWaveSpeed(mesh, sbp, eqn, opts)
  @mpi_master println("max wave speed = ", wave_speed)
  @mpi_master println("min element size = ", mesh.min_el_size)
  delta_t = opts["CFL"]*mesh.min_el_size/wave_speed
  println("for a CFL of ", opts["CFL"], " delta_t = ", delta_t)
  opts["delta_t"] = delta_t
end


#------------------------------------------------------------------------------
#=
geometric_edge_number = 4
eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
force = EulerEquationMod.calcNumericalForce(mesh, sbp, eqn, opts, geometric_edge_number)
println("\nNumerical force on geometric edge ", geometric_edge_number, 
        " = ", norm(force,2))
println("force in the X-direction = ", force[1])
println("force in the Y-direction = ", force[2])
=#
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
    q = sview(eqn.q,:,j,i)
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
   @time rk4(evalEuler, opts["delta_t"], t_max, mesh, sbp, eqn, opts, res_tol=opts["res_abstol"], real_time=opts["real_time"])
#   @time rk4(evalEuler, delta_t, t_max, eqn.q_vec, eqn.res_vec, 
#              (mesh, sbp, eqn), opts, majorIterationCallback=eqn.majorIterationCallback, 
#              res_tol=opts["res_abstol"], real_time=opts["real_time"])

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

  elseif flag == 20

#     @time crank_nicolson(evalEuler, opts["delta_t"], t_max, mesh, sbp, 
#                          eqn, opts, res_tol=opts["res_abstol"], real_time=opts["real_time"])
    @time crank_nicolson(evalEuler, opts["delta_t"], t_max, mesh, sbp, eqn, 
                         opts, opts["res_abstol"], opts["real_time"])

  end       # end of if/elseif blocks checking flag

  println("total solution time printed above")
  # evaluate residual at final q value
  eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
  evalEuler( mesh, sbp, eqn, opts, eqn.params.t)

  eqn.res_vec[:] = 0.0
  eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)


  if opts["write_finalsolution"]
    @mpi_master println("writing final solution")
    writedlm("solution_final.dat", real(eqn.q_vec))
  end

  if opts["write_finalresidual"]
    writedlm("residual_final_$myrank.dat", real(eqn.res_vec))
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
#        q_view = sview(q_vec, i:(i+mesh.numDofPerNode-1))
#        convertFromNaturalToWorkingVars(eqn.params, q_view, q_view)
#      end
#    end

      myrank = mesh.myrank
      q_diff = eqn.q_vec - q_exact
      saveSolutionToMesh(mesh, abs(real(q_diff)))
      writeVisFiles(mesh, "solution_error")


      diff_norm = calcNorm(eqn, q_diff)
#      diff_norm = MPI.Allreduce(diff_norm, MPI.SUM, mesh.comm)
#      diff_norm = sqrt(diff_norm)


      @mpi_master println("solution error norm = ", diff_norm)
      h_avg = calcMeshH(mesh, sbp, eqn, opts)

      # print to file
      @mpi_master begin
        outname = opts["calc_error_outfname"]
        f = open(outname, "w")
        println(f, diff_norm, " ", h_avg)
        close(f)
      end

      #---- Calculate functional on a boundary  -----#
      if opts["calc_functional"]
        
        eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
        if mesh.isDG
          boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, eqn.q, eqn.q_bndry)
        end

        # Calculate functional over edges
        num_functionals = opts["num_functionals"]
        for j = 1:num_functionals
          # Geometric edge at which the functional needs to be integrated
          key_j = string("geom_edges_functional", j)
          functional_edges = opts[key_j]
          functional_name = getFunctionalName(opts, j)

          functional_val = zero(Tsol)
          functional_val = calcBndryFunctional(mesh, sbp, eqn, opts, 
                           functional_name, functional_edges)

          println("\nNumerical functional value on geometric edges ", 
                  functional_edges, " = ", functional_val)

          analytical_functional_val = opts["analytical_functional_val"]
          println("analytical_functional_val = ", analytical_functional_val)

          absolute_functional_error = norm((functional_val - 
                                           analytical_functional_val), 2)
          relative_functional_error = absolute_functional_error/
                                      norm(analytical_functional_val, 2)

          mesh_metric = 1/sqrt(mesh.numEl/2)  # TODO: Find a suitable mesh metric

          # write functional error to file
          outname = string(opts["functional_error_outfname"], j, ".dat")
          println("printed relative functional error = ", 
                  relative_functional_error, " to file ", outname, '\n')
          f = open(outname, "w")
          println(f, relative_functional_error, " ", mesh_metric)
          close(f)
        end  # End for i = 1:num_functionals
      end    # End if opts["calc_functional"]


      #----- Calculate Adjoint Vector For A Functional -----#
      if opts["calc_adjoint"]
        eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
        if mesh.isDG
          boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, eqn.q, eqn.q_bndry)
        end

        # TODO: Presently adjoint computation only for 1 functional. Figure out
        # API based on future use.
        j = 1
        key = string("geom_edges_functional", j)
        functional_edges = opts[key]
        functional_number = j
        functional_name = getFunctionalName(opts, j)
        
        adjoint_vec = zeros(Tsol, mesh.numDof)
        calcAdjoint(mesh, sbp, eqn, opts, functional_name, functional_number, adjoint_vec)


        # Write adjoint vector to file and mesh
        file_object = open("adjoint_vector.dat", "w")
        for iter = 1:length(adjoint_vec)
          println(file_object, real(adjoint_vec[iter]))
        end
        close(file_object)
        saveSolutionToMesh(mesh, real(adjoint_vec))
        writeVisFiles(mesh, "adjoint_field")

      end  # End if opts["calc_adjoint"]


    end
  end

  saveSolutionToMesh(mesh, real(eqn.q_vec))
  printSolution(mesh, real(eqn.q_vec))
  printCoordinates(mesh)
  writeVisFiles(mesh, "solution_done")
  writedlm("solution_done.dat", real(eqn.q_vec))

end  # end if (opts[solve])

fname = "timing_breakdown_$myrank"
write_timings(eqn.params.time, fname)


#  return mesh, sbp, eqn, opts
#end  # end function
#runtest(1)
