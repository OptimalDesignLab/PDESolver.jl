# Startup file for 1 dof advection equation

push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/solver/advection"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/NonlinearSolvers"))

push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Utils"))
using ODLCommonTools
using PdePumiInterface     # common mesh interface - pumi
using SummationByParts     # SBP operators
using AdvectionEquationMod # Advection equation module
using ForwardDiff
using NonlinearSolvers     # non-linear solvers
using ArrayViews
using Utils

include(joinpath(Pkg.dir("PDESolver"),"src/solver/advection/output.jl"))  # printing results to files
include(joinpath(Pkg.dir("PDESolver"), "src/input/read_input.jl"))


#function runtest(flag::Int)
println("ARGS = ", ARGS)
println("size(ARGS) = ", size(ARGS))
opts = read_input(ARGS[1])

# flag determines whether to calculate u, dR/du, or dR/dx (1, 2, or 3)
flag = opts["run_type"]

# timestepping parameters
delta_t = opts["delta_t"]
t_max = opts["t_max"]

order = opts["order"]  # order of accuracy

if flag == 1 || flag == 8  || flag == 9 || flag == 10  # normal run
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

if opts["use_DG"]
  println("\nConstructing SBP Operator")
  # create DG SBP operator with internal nodes only
  sbp = TriSBP{Tsbp}(degree=order, reorder=false, internal=true)
  ref_verts = [0. 1 0; 0 0 1]
  interp_op = SummationByParts.buildinterpolation(sbp, ref_verts)
  sbpface = TriFace{Float64}(order, sbp.cub, ref_verts.')

  # create linear mesh with 4 dof per node

  println("constructing DG mesh")
  mesh = PumiMeshDG2{Tmsh}(dmg_name, smb_name, order, sbp, opts, interp_op, sbpface; 
                   dofpernode=1, coloring_distance=opts["coloring_distance"])
  if opts["jac_type"] == 3 || opts["jac_type"] == 4
    pmesh = PumiMeshDG2Preconditioning(mesh, sbp, opts; 
                   coloring_distance=opts["coloring_distance_prec"])
  else
    pmesh = mesh
  end

else  # continuous Galerkin
  # create SBP object
  println("\nConstructing SBP Operator")
  sbp = TriSBP{Tsbp}(degree=order)  # create linear sbp operator
  # create linear mesh with 4 dof per node

  println("constructing CG mesh")
  mesh = PumiMesh2{Tmsh}(dmg_name, smb_name, order, sbp, opts; dofpernode=1, coloring_distance=opts["coloring_distance"])

  if opts["jac_type"] == 3 || opts["jac_type"] == 4
    pmesh = PumiMesh2Preconditioning(mesh, sbp, opts; coloring_distance=opts["coloring_distance_prec"])
  else
    pmesh = mesh
  end
end

println("\ntypeof(mesh) = ", typeof(mesh))
println("is subtype of DG mesh = ", typeof(mesh) <: AbstractDGMesh)
println("mesh.isDG = ", mesh.isDG)

# Create advection equation object
Tdim = 2
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

if opts["calc_error"]
  println("\ncalculating error of file ", opts["calc_error_infname"], 
          " compared to initial condition")
  vals = readdlm(opts["calc_error_infname"])
  @assert length(vals) == mesh.numDof

  err_vec = vals - eqn.q_vec
  err = calcNorm(eqn, err_vec)
  outname = opts["calc_error_outfname"]
  println("printed err = ", err, " to file ", outname)
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

if opts["calc_havg"]
  # calculate the average mesh size
  jac_3d = reshape(mesh.jac, 1, mesh.numNodesPerElement, mesh.numEl)
  jac_vec = zeros(Tmsh, mesh.numNodes)
  assembleArray(mesh, sbp, eqn, opts, jac_3d, jac_vec)
  # scale by the minimum distance between nodes on a reference element
  # this is a bit of an assumption, because for distorted elements this
  # might not be entirely accurate
  println("mesh.min_node_distance = ", mesh.min_node_dist)
  h_avg = sum(1./sqrt(jac_vec))/length(jac_vec)
#  println("h_avg = ", h_avg)
  h_avg *= mesh.min_node_dist
#  println("h_avg = ", h_avg)

  rmfile("havg.dat")
  f = open("havg.dat", "w")
  println(f, h_avg)
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
global int_advec = 1

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

# evalAdvection(mesh, sbp, eqn, opts, t)
if opts["solve"]
  
  if flag == 1 # normal run
    # RK4 solver
    @time rk4(evalAdvection, delta_t, t_max, mesh, sbp, eqn, opts, 
              res_tol=opts["res_abstol"], real_time=opts["real_time"])
    println("finish rk4")
    printSolution("rk4_solution.dat", eqn.res_vec)
  
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

  elseif flag == 6
    @time newton_check(evalAdvection, mesh, sbp, eqn, opts)
    vals = abs(real(eqn.res_vec))  # remove unneded imaginary part
    saveSolutionToMesh(mesh, vals)
    writeVisFiles(mesh, "solution_error")
    printBoundaryEdgeNums(mesh)
    printSolution(mesh, vals)

  elseif flag == 7
    @time jac_col = newton_check(evalAdvection, mesh, sbp, eqn, opts, 1)
    writedlm("solution.dat", jac_col)

  elseif flag == 8
    @time jac_col = newton_check_fd(evalAdvection, mesh, sbp, eqn, opts, 1)
    writedlm("solution.dat", jac_col)

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

    function test_post_func(mesh, sbp, eqn, opts)
      return calcNorm(eqn, eqn.res_vec)
    end


    rk4(evalAdvection, delta_t, t_max, eqn.q_vec, eqn.res_vec, test_pre_func, test_post_func, (mesh, sbp, eqn), opts, majorIterationCallback=eqn.majorIterationCallback, real_time=opts["real_time"])
  end       # end of if/elseif blocks checking flag

  println("total solution time printed above")
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
      q_diff = eqn.q_vec - q_exact
#      diff_norm = norm(q_diff, Inf)
      diff_norm = calcNorm(eqn, q_diff)
      discrete_norm = norm(q_diff/length(q_diff))

      println("solution error norm = ", diff_norm)
      println("solution discrete L2 norm = ", discrete_norm)

#      sol_norm = norm(eqn.q_vec, Inf)
#      exact_norm = norm(q_exact)

      sol_norm = calcNorm(eqn, eqn.q_vec)
      exact_norm = calcNorm(eqn, q_exact)
      println("numerical solution norm = ", sol_norm)
      println("exact solution norm = ", exact_norm)

      # calculate the average mesh size
      jac_3d = reshape(mesh.jac, 1, mesh.numNodesPerElement, mesh.numEl)
      jac_vec = zeros(Tmsh, mesh.numNodes)
      assembleArray(mesh, sbp, eqn, opts, jac_3d, jac_vec)
      # scale by the minimum distance between nodes on a reference element
      # this is a bit of an assumption, because for distorted elements this
      # might not be entirely accurate
      println("mesh.min_node_distance = ", mesh.min_node_dist)
      h_avg = sum(1./sqrt(jac_vec))/length(jac_vec)
    #  println("h_avg = ", h_avg)
      h_avg *= mesh.min_node_dist
    #  println("h_avg = ", h_avg)


      #----  Calculate functional on a boundary  -----
      
      # Create the adjoint vector
      
        #file_name = "./jacobian1.dat"
        #res_jac = readdlm(file_name) # get the residual Jacobian
        #@assert size(res_jac) == (mesh.numDof, mesh.numDof)
      #=
      newton_data = NonlinearSolvers.NewtonData(mesh, sbp, eqn, opts)
      res_jac = zeros(Tres, mesh.numDof, mesh.numDof)
      pert = complex(0, opts["epsilon"])
      NonlinearSolvers.calcJacobianComplex(newton_data, mesh, sbp, eqn, opts, evalAdvection, pert, res_jac)
      
      func_deriv = ones(Tsol, mesh.numDof)
      adjoint_vec = -res_jac\func_deriv
      
      file_object = open("adjoint_vector.dat", "w")
      for iter = 1:length(adjoint_vec)
        println(file_object, adjoint_vec[iter])
      end
      close(file_object)
      =#

      if opts["calc_force"]
        if mesh.isDG
          boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, eqn.q, eqn.q_bndry)
        end
        eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
        
        # calculate functional over edges
        functional_edges = [2,3] # geometric edge at which the functional needs to be integrated
        functional_val = 0.0
        for k = 1:length(functional_edges)
          geometric_edge_number = functional_edges[k]  
          functional_val += AdvectionEquationMod.calcBndryfunctional(mesh, sbp, eqn,
                           opts, geometric_edge_number)
        end
        println("\nNumerical functional value on geometric edges ", 
                functional_edges, " = ", functional_val)
        analytical_functional_val = (exp(6+6) - exp(6+2))
        # analytical_functional_val = 3^6 + (3^6)/6 - 3^5 - 1/6
        println("analytical_functional_val = ", analytical_functional_val)
        absolute_functional_error = norm((functional_val - analytical_functional_val), 2)
        relative_functional_error = absolute_functional_error/norm(analytical_functional_val, 2)
        
        # write force error to file
        outname = opts["force_error_outfname"]
        println("printed relative functional error = ", 
                relative_functional_error, " to file ", outname, '\n')
        f = open(outname, "w")
        println(f, relative_functional_error, " ", h_avg)
        close(f)
      end

      # print to file
      outname = opts["calc_error_outfname"]
      println("printed err = ", diff_norm, " to file ", outname)
      f = open(outname, "w")
      println(f, diff_norm, " ", h_avg)
      close(f)


#=
      outname = opts["calc_error_outfname"]
      f = open(outname, "w")
      println(f, mesh.numEl, " ", diff_norm, " ", discrete_norm)
      close(f)
=#
    end
  end

  saveSolutionToMesh(mesh, real(eqn.q_vec))
  printSolution(mesh, real(eqn.q_vec))
  printCoordinates(mesh)
  writeVisFiles(mesh, "solution_done")

end  # end if (opts[solve])
