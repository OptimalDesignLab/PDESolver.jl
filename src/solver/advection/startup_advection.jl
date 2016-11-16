# Startup file for 1 dof advection equation

push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/solver/advection"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/NonlinearSolvers"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Debugging"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Utils"))
using ODLCommonTools
using PdePumiInterface     # common mesh interface - pumi
using SummationByParts     # SBP operators
using AdvectionEquationMod # Advection equation module
using ForwardDiff
using NonlinearSolvers     # non-linear solvers
using ArrayViews
using Utils
using MPI
include(joinpath(Pkg.dir("PDESolver"),"src/solver/advection/output.jl"))  # printing results to files
include(joinpath(Pkg.dir("PDESolver"), "src/input/read_input.jl"))


if !MPI.Initialized()
  MPI.Init()
end

#function runtest(flag::Int)
println("ARGS = ", ARGS)
println("size(ARGS) = ", size(ARGS))
opts = read_input(ARGS[1])

# timestepping parameters
delta_t = opts["delta_t"]
t_max = opts["t_max"]
dim = opts["dimensions"]
# flag determines whether to calculate u, dR/du, or dR/dx (1, 2, or 3)
flag = opts["run_type"]

dofpernode = 1

sbp, mesh, pmesh, Tsol, Tres, Tmsh, mesh_time = createMeshAndOperator(opts, dofpernode)

#mesh.numColors = 4
#mesh.maxColors = 4

if opts["write_timing"]
  MPI.Barrier(mesh.comm)
  if mesh.myrank == 0
    f = open("timing.dat", "a+")
    println(f, mesh_time)
    close(f)
  end
end

myrank = mesh.myrank

println("\ntypeof(mesh) = ", typeof(mesh))
println("is subtype of DG mesh = ", typeof(mesh) <: AbstractDGMesh)
println("mesh.isDG = ", mesh.isDG)

# Create advection equation object
Tdim = dim
eqn = AdvectionData_{Tsol, Tres, Tdim, Tmsh}(mesh, sbp, opts)

q_vec = eqn.q_vec


# Initialize the advection equation
init(mesh, sbp, eqn, opts)

# Populate with initial conditions
println("\nEvaluating initial condition")
ICfunc_name = opts["IC_name"]
ICfunc = ICDict[ICfunc_name]
println("ICfunc = ", ICfunc)
ICfunc(mesh, sbp, eqn, opts, q_vec) 
println("finished initializing q")

writedlm("solution_ic.dat", eqn.q_vec)

if opts["calc_error"]
  println("\ncalculating error of file ", opts["calc_error_infname"], 
          " compared to initial condition")
  vals = readdlm(opts["calc_error_infname"])
  @assert length(vals) == mesh.numDof

  err_vec = vals - eqn.q_vec
  err_local = calcNorm(eqn, err_vec)
  if myrank == 0
    outname = opts["calc_error_outfname"]
    println("printed err = ", err, " to file ", outname)
    f = open(outname, "w")
    println(f, err)
    close(f)
  end
end

if opts["calc_trunc_error"]  # calculate truncation error
  println("\nCalculating residual for truncation error")
  tmp = calcResidual(mesh, sbp, eqn, opts, evalAdvection)
  if myrank == 0
    f = open("error_trunc.dat", "w")
    println(f, tmp)
    close(f)
  end
end

if opts["calc_havg"]
  h_avg = calcMeshH(mesh, sbp, eqn, opts)
  if myrank == 0  # TODO: make this globally accurate
    rmfile("havg.dat")
    f = open("havg.dat", "w")
    println(f, h_avg)
    close(f)
  end
end



if opts["perturb_ic"]
  println("\nPerturbing initial condition")
  perturb_mag = opts["perturb_mag"]
  for i=1:mesh.numDof
    q_vec[i] += perturb_mag*rand()
  end
end

res_vec_exact = deepcopy(q_vec)

#rmfile("IC_$myrank.dat")
#writedlm("IC_$myrank.dat", real(q_vec))
saveSolutionToMesh(mesh, q_vec)
writeVisFiles(mesh, "solution_ic")
global int_advec = 1

if opts["calc_dt"]
  alpha_net = sqrt(eqn.params.alpha_x^2 + eqn.params.alpha_y^2)
  opts["delta_t"] = opts["CFL"]*mesh.min_el_size/alpha_net
end


if opts["test_GLS2"]
  calcResidual(mesh, sbp, eqn, opts, evalAdvection)
end
println("mesh.min_node_dist = ", mesh.min_node_dist)
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
MPI.Barrier( mesh.comm)
# evalAdvection(mesh, sbp, eqn, opts, t)
if opts["solve"]
  
  solve_time = @elapsed if flag == 1 # normal run
    # RK4 solver
    @time rk4(evalAdvection, delta_t, t_max, mesh, sbp, eqn, opts, 
              res_tol=opts["res_abstol"], real_time=opts["real_time"])
    println("finish rk4")
#    printSolution("rk4_solution.dat", eqn.res_vec)
  
  elseif flag == 2 # forward diff dR/du
    
    # define nested function
    function dRdu_rk4_wrapper(u_vals::AbstractVector, res_vec::AbstractVector)
      eqn.q_vec = u_vals
      eqn.q_vec = res_vec
      rk4(evalAdvection, delta_t, t_max, mesh, sbp, eqn)
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
    @time newton(evalAdvection, mesh, sbp, eqn, opts, pmesh, itermax=opts["itermax"], 
                 step_tol=opts["step_tol"], res_abstol=opts["res_abstol"], 
                 res_reltol=opts["res_reltol"], res_reltol0=opts["res_reltol0"])

    printSolution("newton_solution.dat", eqn.res_vec)

  elseif flag == 9
    # to non-pde rk4 run
    function pre_func(mesh, sbp, eqn,  opts)
      println("pre_func was called")
      return nothing
    end

    function post_func(mesh, sbp, eqn, opts)
#      for i=1:mesh.numDof
#        eqn.res_vec[i] *= eqn.Minv[i]
#      end
      nrm = norm(eqn.res_vec)
      println("post_func returning residual norm = ", nrm)
      return nrm
    end

    rk4(evalAdvection, delta_t, t_max, eqn.q_vec, eqn.res_vec, pre_func, post_func, (mesh, sbp, eqn), opts, majorIterationCallback=eqn.majorIterationCallback, real_time=opts["real_time"])

  elseif flag == 10
    function test_pre_func(mesh, sbp, eqn, opts)
      
      eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
    end

    function test_post_func(mesh, sbp, eqn, opts, calc_norm=true)
      return calcNorm(eqn, eqn.res_vec)
    end


    rk4(evalAdvection, delta_t, t_max, eqn.q_vec, eqn.res_vec, test_pre_func,
        test_post_func, (mesh, sbp, eqn), opts, 
        majorIterationCallback=eqn.majorIterationCallback, real_time=opts["real_time"])

  elseif flag == 20

    @time crank_nicolson(evalAdvection, opts["delta_t"], t_max, mesh, sbp, eqn, 
                         opts, opts["res_abstol"], opts["real_time"])

#   else
#     throw(ErrorException("No flag specified: no solve will take place"))
#     return nothing

  end       # end of if/elseif blocks checking flag

  println("total solution time printed above")

  if opts["write_timing"]
    MPI.Barrier(mesh.comm)
    if mesh.myrank == 0
      f = open("timing.dat", "a+")
      println(f, solve_time)
      close(f)
    end
  end

  # evaluate residual at final q value
  need_res = false
  if need_res
    eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
    evalAdvection(mesh, sbp, eqn, opts, eqn.t)

    eqn.res_vec[:] = 0.0
    eqn.assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
  end

  if opts["write_finalsolution"]
    println("writing final solution")
    writedlm("solution_final_$myrank.dat", real(eqn.q_vec))
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
      q_diff = eqn.q_vec - q_exact
      diff_norm = calcNorm(eqn, q_diff)
      discrete_norm = norm(q_diff/length(q_diff))
      discrete_norm = MPI.Allreduce(discrete_norm*discrete_norm, MPI.SUM, mesh.comm)
      discrete_norm = sqrt(discrete_norm)

      if myrank == 0
        println("solution error norm = ", diff_norm)
        println("solution discrete L2 norm = ", discrete_norm)
      end

      sol_norm = calcNorm(eqn, eqn.q_vec)
      exact_norm = calcNorm(eqn, q_exact)
      @mpi_master println("numerical solution norm = ", sol_norm)
      @mpi_master println("exact solution norm = ", exact_norm)

      # calculate the average mesh size
      h_avg = calcMeshH(mesh, sbp, eqn, opts)
      if myrank == 0
#        println("mesh.min_node_distance = ", mesh.min_node_dist)
        # print to file
        outname = opts["calc_error_outfname"]
        println("printed err = ", diff_norm, " to file ", outname)
        f = open(outname, "w")
        println(f, diff_norm, " ", h_avg)
        close(f)
      end

      #----  Calculate functional on a boundary  -----#
      
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
          functional_val = calcBndryfunctional(mesh, sbp, eqn, opts, 
                           functional_name, functional_edges)
          println("\nNumerical functional value on geometric edges ", 
                  functional_edges, " = ", functional_val)
          
          analytical_functional_val = opts["analytical_functional_val"]
          println("analytical_functional_val = ", analytical_functional_val)
          
          absolute_functional_error = norm((functional_val - 
                                           analytical_functional_val), 2)
          relative_functional_error = absolute_functional_error/
                                      norm(analytical_functional_val, 2)
          
          mesh_metric = 1/sqrt(mesh.numEl/2)  # Only for a square domain with triangular elements
          
          # write functional error to file
          outname = string(opts["functional_error_outfname"], j, ".dat")
          println("printed relative functional error = ", 
                  relative_functional_error, " to file ", outname, '\n')
          f = open(outname, "w")
          println(f, relative_functional_error, " ", mesh_metric)
          close(f)
        end  # End for j = 1:num_functional
      end    # End if opts["calc_functional"]


      #-----  Calculate adjoint vector for a functional  -----#

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
      end  # end if opts["calc_adjoint"]

      


#=
      outname = opts["calc_error_outfname"]
      f = open(outname, "w")
      println(f, mesh.numEl, " ", diff_norm, " ", discrete_norm)
      close(f)
=#
    end
  end

  myrank = mesh.myrank
#  f = open("profile_$myrank.dat", "a+")
#  Profile.print(f, format=:flat, C=true)
#  close(f)

  saveSolutionToMesh(mesh, real(eqn.q_vec))
#  printSolution(mesh, real(eqn.q_vec))
#  printCoordinates(mesh)
  writeVisFiles(mesh, "solution_done")

  # write timings
  params = eqn.params
  myrank = mesh.myrank
#  timings = [params.t_volume, params.t_face, params.t_source, params.t_sharedface, params.t_bndry, params.t_send, params.t_wait, params.t_allreduce, params.t_jacobian, params.t_solve, params.t_barrier, params.t_barrier2, params.t_barrier3]
#  writedlm("timing_breakdown_$myrank.dat", vcat(timings, params.t_barriers))
  fname = "timing_breakdown_$myrank"
  write_timings(params.time, fname)

  MPI.Barrier(mesh.comm)
  if opts["finalize_mpi"]
    MPI.Finalize()
  end
end  # end if (opts[solve])
