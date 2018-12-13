# Description: startup function for solving an equation

import PDESolver.solvePDE

#=
"""
This function invokes the solver for the Elliptic equations, using the
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
function run_elliptic(input_file::AbstractString)

  mesh, sbp, eqn, opts, pmesh = createObjects(input_file)
  solve_elliptic(mesh, sbp, eqn, opts, pmesh)

  return mesh, sbp, eqn, opts
end
=#

"""
This function creates and initializes the mesh, sbp, eqn, and opts objects

Inputs:
    opts: options dictionary

Outputs:
mesh: an AbstractMesh.  The concrete type is determined by the options
dictionary
sbp: an AbstractOperator.  The concrete type is determined by the options
dictionary
eqn: an EllipticData object
opts: the options dictionary
pmesh: mesh used for preconditioning, can be same object as mesh
"""
function createObjects(opts::Dict)

  Tdim = opts["dimensions"]
  dofpernode = 1 

  sbp, mesh, pmesh, Tsol, Tres, Tmsh, mesh_time = createMeshAndOperator(opts, dofpernode)

  # TODO: input argument for dofpernode

  # create elliptic equation
  # var_type = opts["variable_type"]
  # eqn = EllipticData_{Tsol, Tres, Tdim, Tmsh, var_type}(mesh, sbp, opts)
  eqn = EllipticData_{Tsol, Tres, Tdim, Tmsh}(mesh, sbp, opts)

  # initialize physics module and populate any fields in mesh and eqn that
  # depend on the physics module
  init(mesh, sbp, eqn, opts, pmesh)

  return mesh, sbp, eqn, opts, pmesh
end

"""
Given fully initialized mesh, sbp, eqn, opts, this function solves
the Elliptic equations.  The 4 object should be obtained from createObjects().


Specifically, it applies an initial condition and invokes a nonlinear
solver according to the options dictionary.

Inputs:
mesh: an AbstractMesh
sbp: an AbstractOperator
eqn: an AbstractEllipticData
opts: the options dictionary.  This must be the options dictionary returned
by createObjects().  Changing values in the options dictionary after
calling createObjects() results in undefined behavior.
pmesh: mesh used for preconditioning, can be same object as mesh.
default value of mesh

"""
function solvePDE(mesh::AbstractMesh, sbp::AbstractOperator, eqn::AbstractEllipticData, opts::Dict, pmesh::AbstractMesh=mesh)
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

    tmp = calcResidual(mesh, sbp, eqn, opts, evalResidual)
    res_real = real(eqn.res_vec)
    opts["res_reltol0"] = tmp
    println("res_reltol0 = ", tmp)

    saveSolutionToMesh(mesh, res_real)
    writeVisFiles(mesh, "solution_relfunc")
  end

  # populate u0 with initial condition
  @mpi_master println("\nEvaluating initial condition")
  ICfunc_name = opts["IC_name"]
  ICfunc = ICDict[ICfunc_name]
  @mpi_master println("ICfunc = ", ICfunc)
  ICfunc(mesh, sbp, eqn, opts, q_vec)

  res_vec_exact = deepcopy(q_vec)

  #writeSolutionFiles(mesh, sbp, eqn, opts, "IC")
  saveSolutionToMesh(mesh, q_vec)

  writeVisFiles(mesh, "solution_ic")

  call_nlsolver(mesh, sbp, eqn, opts, pmesh)
  postproc(mesh, sbp, eqn, opts)

  MPI.Barrier(mesh.comm)

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
  myrank = mesh.myrank

  @mpi_master println("\nDoing postprocessing")

  if haskey(opts, "exactSolution")
    l2norm, lInfnorm = calcErrorL2Norm(mesh, sbp, eqn, opts)
    @mpi_master begin
      println("L2Norm = ", l2norm)
      println("LinfNorm = ", lInfnorm)
      fname = "l2norm.dat"
      f = open(fname, "w")
      println(f, l2norm)
      close(f)
    end
  end

  if haskey(opts, "Functional")
    if haskey(opts, "exactFunctional")
      exactFunctional = opts["exactFunctional"]
    end
    functional_value = Array{typeof(eqn.q[1,1,1])}(mesh.numDofPerNode)
    eqn.functional(mesh, sbp, eqn, opts, functional_value)

    @mpi_master begin
      println("functional = ", abs(real(functional_value[1]) - exactFunctional))
      fname = "functional.dat"
      f = open(fname, "w")
      println(f, abs(real(functional_value[1]) - exactFunctional))
      close(f)
    end
  end

  return nothing
end


