# Description: startup function for solving an equation

"""
This function invokes the solver for the Euler equations, using the
specified input file

Inputs:
input_file: a string containing the path to an input file, or just a file
name.  If it is just a file name, it is taken to be in the
users pwd.

Outputs:
mesh: an AbstractMesh object
sbp: the SBP operator used in solving the equation
eqn: the AbstractSolutionData object used to solve the equation
opts: the options dictonary

"""
function run_euler(input_file::AbstractString)

  mesh, sbp, eqn, opts, pmesh = createObjects(input_file)
  solve_euler(mesh, sbp, eqn, opts, pmesh)

  return mesh, sbp, eqn, opts
end


"""
This function creates and initializes the mesh, sbp, eqn, and opts objects

Inputs:
file_name: input file name

Outputs:
mesh: an AbstractMesh.  The concrete type is determined by the options
dictionary
sbp: an AbstractSBP.  The concrete type is determined by the options
dictionary
eqn: an EulerData object
opts: the options dictionary
pmesh: mesh used for preconditioning, can be same object as mesh
"""
function createObjects(input_file::AbstractString)

  if !MPI.Initialized()
    MPI.Init()
  end

  opts = read_input(input_file)  # read input file and get default values
  checkOptions(opts)  # physics specific options checking
  #opts = read_input("input_vals_channel2.jl")

  Tdim = opts["dimensions"]
  dofpernode = Tdim + 2

  sbp, mesh, pmesh, Tsol, Tres, Tmsh, mesh_time = createMeshAndOperator(opts, dofpernode)

  # TODO: input argument for dofpernode

  # create euler equation
  var_type = opts["variable_type"]
  eqn = EulerData_{Tsol, Tres, Tdim, Tmsh, var_type}(mesh, sbp, opts)

  # initialize physics module and populate any fields in mesh and eqn that
  # depend on the physics module
  init(mesh, sbp, eqn, opts, pmesh)

  return mesh, sbp, eqn, opts, pmesh
end

"""
Given fully initialized mesh, sbp, eqn, opts, this function solves
the Euler equations.  The 4 object should be obtained from createObjects().


Specifically, it applies an initial condition and invokes a nonlinear
solver according to the options dictionary.

Inputs:
mesh: an AbstractMesh
sbp: an AbstractSBP
eqn: an AbstractEulerData
opts: the options dictionary.  This must be the options dictionary returned
by createObjects().  Changing values in the options dictionary after
calling createObjects() results in undefined behavior.
pmesh: mesh used for preconditioning, can be same object as mesh.
default value of mesh

"""
function solve_euler(mesh::AbstractMesh, sbp, eqn::AbstractEulerData, opts, pmesh=mesh)
  #delta_t = opts["delta_t"]   # delta_t: timestep for RK

  myrank = mesh.myrank
  # flag determines whether to calculate u, dR/du, or dR/dx (1, 2, or 3)
  t_max = opts["t_max"]
  flag = opts["run_type"]
  var_type = opts["variable_type"]

  # timestepping parameters
  order = opts["order"]       # order of accuracy


  # zero some arrays out in case the user forgot

  fill!(eqn.res, 0.0)
  fill!(eqn.res_vec, 0.0)
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
    tmp = calcResidual(mesh, sbp, eqn, opts, evalResidual)
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

  if opts["calc_error"]
    @mpi_master println("\ncalculating error of file ",
                        opts["calc_error_infname"],
                        " compared to initial condition")

    # read in this processors portion of the solution
    vals = readdlm(get_parallel_fname(opts["calc_error_infname"], myrank))
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
    tmp = calcResidual(mesh, sbp, eqn, opts, evalResidual)

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

  call_nlsolver(mesh, sbp, eqn, opts, pmesh)
  postproc(mesh, sbp, eqn, opts)

  cleanup(mesh, sbp, eqn, opts)

  MPI.Barrier(mesh.comm)
  if opts["finalize_mpi"]
    MPI.Finalize()
  end

  return mesh, sbp, eqn, opts
end  # end function
#runtest(1)

"""
This function does post processing, if requested by the input options.
Typical post processing includes calculation of errors, norms of important
quantities, writing of files. etc.

Inputs:
mesh
sbp
eqn
opts
"""
function postproc(mesh, sbp, eqn, opts)

  ##### Do postprocessing ######
  println("\nDoing postprocessing")

  Tsol = eltype(eqn.q)
  Tres = eltype(eqn.res)
  Tmsh = eltype(mesh.dxidx)
  myrank = mesh.myrank

  if haskey(opts, "Functional") && haskey(opts, "exactFunctional")
    exact_func = opts["exactFunctional"]
    func_name = opts["Functional"]
    functional = VolumeFunctionalDict[func_name]
    func_val = Array(Tsol, mesh.numDofPerNode)
    functional(mesh, sbp, eqn, opts, func_val)
    func_error = real(func_val[1]) - exact_func
    func_error = abs(func_error)
    fname = "functional.dat"
    f = open(fname, "w")
    println(f, func_error)
    println("functional error = ", func_error)
  end

  if haskey(opts, "exact_soln_func") && opts["solve"]
    exfname = opts["exact_soln_func"]
    if haskey(ICDict, exfname)
      println("calculating error...")
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
      # q_diff = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)
      # q_diff_vec = reshape(q_diff, mesh.numDof) 
      q_diff_vec = eqn.q_vec - q_exact
      # eqn.disassembleSolution(mesh, sbp, eqn, opts, q_diff, q_diff_vec)
      saveSolutionToMesh(mesh, abs(real(q_diff_vec)))
      writeVisFiles(mesh, "solution_error")

      diff_norm = calcNorm(eqn, q_diff_vec)
      #      diff_norm = MPI.Allreduce(diff_norm, MPI.SUM, mesh.comm)
      #      diff_norm = sqrt(diff_norm)
      q_error = zeros(Float64, mesh.numDofPerNode)
      for el = 1 : mesh.numEl
        for n = 1 : mesh.numNodesPerElement
          for dof = 1 : mesh.numDofPerNode
            dofnum = mesh.dofs[dof, n, el]
            q_error_j = real(q_diff_vec[dofnum])

            q_error[dof] += q_error_j  * q_error_j* sbp.w[n]/mesh.jac[n, el]
          end
        end
      end
      @mpi_master println("solution error norm = ", diff_norm)
      for dof = 1 : mesh.numDofPerNode
        q_error[dof] = sqrt(q_error[dof])
        @mpi_master println("solt_error[", dof, "] = ", q_error[dof])
      end
      # TODO: make this mesh.min_el_size?
      h_avg = calcMeshH(mesh, sbp, eqn, opts)

      # print to file
      @mpi_master begin
        outname = opts["calc_error_outfname"]
        f = open(outname, "w")
        println(f, diff_norm, " ", h_avg)
        close(f)
      end

    end
  end

  # Calculate functionals on a boundary
  if opts["calc_functional"]

    eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
    if mesh.isDG
      boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, eqn.q, eqn.q_bndry)
    end

    # Calculate functional over edges
    # num_functionals = opts["num_functionals"]
    # for j = 1:num_functionals
    # Geometric edge at which the functional needs to be integrated
    key_j = string("geom_edges_functional")
    functional_edges = opts[key_j]
    # functional_name = getFunctionalName(opts, j)
    # println(functional_name)

    force = zeros(Tsol,2)
    # functional_val = calcBndryFunctional(mesh, sbp, eqn, opts, 
    # functional_name, functional_edges)
    force = calcBndryFunctional(mesh, sbp, eqn, opts, functional_edges)

    println("\nNumerical functional value on geometric edges ", 
            functional_edges, " = ", real(force))

    # if haskey(opts, "analytical_functional_val")
      # analytical_functional_val = opts["analytical_functional_val"]
      # println("analytical_functional_val = ", analytical_functional_val)

      # functional_val = force[1]
      # absolute_functional_error = norm((functional_val - analytical_functional_val), 2)
      # relative_functional_error = absolute_functional_error/ norm(analytical_functional_val, 2)

      # mesh_metric = 1/sqrt(mesh.numEl/2)  # TODO: Find a suitable mesh metric

      # # write functional error to file
      # outname = string(opts["functional_error_outfname"], j, ".dat")
      # println("printed relative functional error = ", 
              # relative_functional_error, " to file ", outname, '\n')
      # f = open(outname, "w")
      # println(f, relative_functional_error, " ", mesh_metric)
      # close(f)
    # end
    # end  # End for i = 1:num_functionals
  end    # End if opts["calc_functional"]

  #----- Calculate Adjoint Vector For A Functional -----#
  if opts["calc_adjoint"]
    println("\nAjoint solving...\n")
    eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
    if mesh.isDG
      boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, eqn.q, eqn.q_bndry)
    end

    # TODO: Presently adjoint computation only for 1 functional. Figure out
    # API based on future use.
    key = string("geom_edges_functional")
    functional_edges = opts[key]
    # functional_number = j
    # functional_name = getFunctionalName(opts, j)

    adjoint_vec = zeros(Tsol, mesh.numDof, 2)
    # calcAdjoint(mesh, sbp, eqn, opts, functional_number, adjoint_vec)
    calcAdjoint(mesh, sbp, eqn, opts, adjoint_vec)
    saveSolutionToMesh(mesh, sview(real(adjoint_vec[:,1])))
    writeVisFiles(mesh, "adjoint_Cl")
    saveSolutionToMesh(mesh, sview(real(adjoint_vec[:,2])))
    writeVisFiles(mesh, "adjoint_Cd")

  end  # End if opts["calc_adjoint"]

  return nothing
end


