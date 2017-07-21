"""
  Create a SBP operator and a mesh.  This is used by all physics modules
  to create the right type of operator and mesh based on the input options.
  It is type unstable, but that is ok.

  Inputs:
    opts: options dictonary
    dofpernode: number of degrees of freedom on each node

  Outputs
    sbp : an AbstractSBP
    mesh : an AbstractMesh
    pmesh : an AbstractMesh, used for preconditioning, may be same object as
            mesh
    Tsol : DataType that should be used for eqn.q
    Tres : DataType that should be used for eqn.res
    Tmsh : DataType of mesh.dxidx and friends
    mesh_time : time in seconds for creation of mesh (Float64)
"""
function createMeshAndOperator(opts, dofpernode)

  # create mesh with 1 dofpernode
  dmg_name = opts["dmg_name"]
  smb_name = opts["smb_name"]
  order = opts["order"]  # order of accuracy
  dim = opts["dimensions"]

  Tmsh, Tsbp, Tsol, Tres = getDataTypes(opts)

  opts["Tsol"] = Tsol
  opts["Tres"] = Tres
  opts["Tsbp"] = Tsbp
  opts["Tmsh"] = Tmsh

  sbp, sbpface, shape_type, topo = createSBPOperator(opts, Tsbp)
 
  mesh_time = @elapsed mesh, pmesh = createMesh(opts, sbp, sbpface, shape_type,
                                                topo, Tmsh, dofpernode)
  return sbp, mesh, pmesh, Tsol, Tres, Tmsh, mesh_time
end

"""
  This function determines the datatypes of the elements of the arrays of the
  mesh quantities, sbp operator, solution variables and residual.

  If the datatypes cannot be determined, an error is thrown.

  Inputs:
    opts: the options dictionary

  Outputs
    Tmsh
    Tsbp
    Tsol
    Tres
"""
function getDataTypes(opts::Dict)

  flag = opts["run_type"]
  if haskey(opts, "jac_method")
    jac_method = opts["jac_method"]
  end

  if flag == 1 || flag == 8  || flag == 9 || flag == 10 || flag == 30  # normal run
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
  elseif flag == 5 || flag == 40 || flag == 41
    if jac_method == 1 # use Newton method using finite difference  (former flag 4)
      # println("========== utils/initialization: flag 5, jac_method 1")
      Tmsh = Float64
      Tsbp = Float64
      Tsol = Float64
      Tres = Float64
    elseif jac_method == 2 # use complex step dR/du
      # println("========== utils/initialization: flag 5, jac_method 2")
      Tmsh = Float64
      Tsbp = Float64
      Tsol = Complex128
      Tres = Complex128
    else
      throw(ErrorException("Illegal or no jac_method specified for steady Newton initialization."))
    end
  elseif flag == 6 || flag == 7  # evaluate residual error and print to paraview
    Tmsh = Float64
    Tsbp = Float64
    Tsol = Complex128
    Tres = Complex128
  elseif flag == 11 # Same as Flag 5 but Tmsh is complex
    if jac_method == 1 # use Newton method using finite difference  (former flag 4)
      # println("========== utils/initialization: flag 11, jac_method 1")
      Tmsh = Float64
      Tsbp = Float64
      Tsol = Float64
      Tres = Float64
    elseif jac_method == 2 # use complex step dR/du
      # println("========== utils/initialization: flag 11, jac_method 2")
      Tmsh = Complex128
      Tsbp = Float64
      Tsol = Complex128
      Tres = Complex128
    else
      throw(ErrorException("Illegal or no jac_method specified for steady Newton initialization."))
    end
  elseif flag == 20 # jac_method needs to be symbol
    if jac_method == 1 # Crank-Nicolson, FD Jac
      Tmsh = Float64
      Tsbp = Float64
      Tsol = Float64
      Tres = Float64
    elseif jac_method == 2 # Crank-Nicolson, CS Jac
      Tmsh = Float64
      Tsbp = Float64
      Tsol = Complex128
      Tres = Complex128
    else
      throw(ErrorException("Illegal or no jac_method specified for CN initialization."))
    end
  else
    throw(ErrorException("Unrecognized run_type: $flag"))
  end

  return Tmsh, Tsbp, Tsol, Tres
end

"""
  This function constructs the SBP operator and the associated SBP face
  operator, as specified by the options dictionary.  It also determines
  the shape_type that PumiInterface uses to describe the SBP operator to
  Pumi.

  Inputs:
    opts: the options dictionary
    Tsbp: the DataType specifying the Tsbp passed to the SBP operator
          constructor

  Outputs:
    sbp: the SBP operator
    sbpface: the SBP face operator
    shape_type: an integer passed to the mesh constructor to describe the
                operator
    topo: in the 3D DG case, an ElementTopology describing the SBP reference
          element, otherwise the integer 0.
"""
function createSBPOperator(opts::Dict, Tsbp::DataType)
  # construct SBP operator and figure out shape_type needed by Pumi
  order = opts["order"]  # order of accuracy
  dim = opts["dimensions"]

  println("\nConstructing SBP Operator")
  topo = 0  # generally not needed, so return a default value
  if opts["use_DG"]
    if opts["operator_type"] == "SBPOmega"
      if dim == 2
        sbp = getTriSBPOmega(degree=order, Tsbp=Tsbp)
      else
        sbp = getTetSBPOmega(degree=order, Tsbp=Tsbp)
      end
      shape_type = 2
    elseif opts["operator_type"] == "SBPGamma"
      if dim == 2
        sbp = getTriSBPGamma(degree=order, Tsbp=Tsbp)
      else
        sbp = getTetSBPGamma(degree=order, Tsbp=Tsbp)
      end
      shape_type = 3
    elseif opts["operator_type"] == "SBPDiagonalE"
      if dim == 2
        sbp = getTriSBPWithDiagE(degree=order, Tsbp=Tsbp)
      else
        throw(ArgumentError("3 dimensional SBPDiagonalE no supported"))
      end
      shape_type = 4
    elseif opts["operator_type"] == "SBPDiagonalE2"  # no vert nodes
      if dim == 2
        sbp = getTriSBPWithDiagE(degree=order, Tsbp=Tsbp, vertices=false)
      else
        sbp = getTetSBPWithDiagE(degree=order, Tsbp=Tsbp)
      end
      shape_type = 5
    else
      op_type = opts["operator_type"]
      throw(ArgumentError("unrecognized operator type $op_type for DG mesh"))
    end
  else  # CG mesh
    # the CG varient of SBP gamma is the only option
    if opts["operator_type"] != "SBPGamma"
      op_type = opts["operator_type"]
      throw(ArgumentError("invalid operator type $op_type for CG"))
    end
    sbp = getTriSBPGamma(degree=order, Tsbp=Tsbp)
    shape_type = 1
  end
 
  println("\nConstructing SBP Face Operator")
  if opts["use_DG"]
    if dim == 2
      # TODO: use sbp.vtx instead
      ref_verts = [-1. 1 -1; -1 -1 1]
      if opts["operator_type"] == "SBPDiagonalE"
        sbpface = getTriFaceForDiagE(order, sbp.cub, ref_verts.')
      elseif opts["operator_type"] == "SBPDiagonalE2"
        println("getting TriFaceForDiagE2")
        sbpface = getTriFaceForDiagE(order, sbp.cub, ref_verts.', vertices=false)
      else
        sbpface = TriFace{Float64}(order, sbp.cub, ref_verts.')
      end
    else  # dim == 3
      ref_verts = sbp.vtx
      face_verts = SummationByParts.SymCubatures.getfacevertexindices(sbp.cub)
      topo = ElementTopology{3}(face_verts)

      if opts["operator_type"] == "SBPDiagonalE2"
        sbpface = getTetFaceForDiagE(order, sbp.cub, ref_verts)
      else
        sbpface = TetFace{Float64}(order, sbp.cub, ref_verts)
      end
    end  # end if dim == 2
  else   # CG
    if dim == 2
      sbpface = TriFace{Float64}(order, sbp.cub, sbp.vtx)
    else
      throw(ErrorException("3D CG not supported"))
    end
  end
 
  return sbp, sbpface, shape_type, topo
end

"""
  This function creates the mesh object and, optionally, a second mesh
  used for preconditioning

  Inputs:
    opts: the options dictionary
    sbp: an SBP operator
    sbpface: an SBP face operator
    topo: an ElementTopology describing the SBP reference element.  Only
          needed for 3D DG, otherwise can be any value
    Tmsh: the DataType of the elements of the mesh arrays (dxidx, jac, etc.)
    dofpernode: number of degrees of freedom on every node

  All arguments except opts are typically provided by 
  [`createSBPOperator`](@ref) and [`getDataTypes`](@ref)
"""
function createMesh(opts::Dict, sbp::AbstractSBP, sbpface, shape_type, topo,
                    Tmsh, dofpernode)

  dmg_name = opts["dmg_name"]
  smb_name = opts["smb_name"]
  order = opts["order"]  # order of accuracy
  dim = opts["dimensions"]


  if opts["use_DG"]
    println("constructing DG mesh")
    if dim == 2
      mesh = PumiMeshDG2{Tmsh}(dmg_name, smb_name, order, sbp, opts, sbpface; 
                               dofpernode=dofpernode, 
                               coloring_distance=opts["coloring_distance"],
                               shape_type=shape_type)
    else
      mesh = PumiMeshDG3{Tmsh}(dmg_name, smb_name, order, sbp, opts, sbpface,
                               topo; dofpernode=dofpernode,
                               coloring_distance=opts["coloring_distance"],
                               shape_type=shape_type)
    end

    # create preconditioning mesh
    if (opts["jac_type"] == 3 || opts["jac_type"] == 4) && opts["use_jac_precond"]
      @assert dim == 2
      pmesh = PumiMeshDG2Preconditioning(mesh, sbp, opts;
                     coloring_distance=opts["coloring_distance_prec"])
    else
      pmesh = mesh
    end

  else  # continuous Galerkin
    if dim == 2

      println("constructing CG mesh")
      mesh = PumiMesh2{Tmsh}(dmg_name, smb_name, order, sbp, opts, sbpface;
                             dofpernode=dofpernode,
                             coloring_distance=opts["coloring_distance"],
                             shape_type=shape_type)

      if opts["jac_type"] == 3 || opts["jac_type"] == 4
        pmesh = PumiMesh2Preconditioning(mesh, sbp, opts;
                               coloring_distance=opts["coloring_distance_prec"])
      else
        pmesh = mesh
      end

    else  # dim == 3
      throw(ErrorException("3D continuous Galerkin not supported"))
    end
  end  # end if DG


  return mesh, pmesh
end

"""
  This function takes in the 4 principle object, fully initialized, and calls
  a nonlinear solver on them, according to the options in the dictionary.
  The evalResidual function is passed to the nonlinear solver

  Inputs:
    mesh: a mesh object
    sbp: an SBP operator
    eqn: an equation object
    opts: options dictionary, used to determine which nonlinear solver to call
    pmesh: mesh used for calculating preconditioning jacobian in Newton's
           method, default to using mesh if not specified

  Outputs:
    none

  Aliasing restrictions: none (specificaly, mesh and pmesh *can* be the same
                         object)
"""
function call_nlsolver(mesh::AbstractMesh, sbp::AbstractSBP,
                       eqn::AbstractSolutionData, opts::Dict,
                       pmesh::AbstractMesh=mesh)
  flag = opts["run_type"]::Int
  t = 0.0  # for steady methods, t = 0.0 always, for unsteady, the time
           # stepper returns a new t value
  if opts["solve"]

    solve_time = @elapsed if flag == 1 # normal run
      # RK4 solver
      delta_t = opts["delta_t"]
      t_max = opts["t_max"]

      @time t = rk4(evalResidual, delta_t, t_max, mesh, sbp, eqn, opts, 
                res_tol=opts["res_abstol"], real_time=opts["real_time"])
      println("finish rk4")
  #    printSolution("rk4_solution.dat", eqn.res_vec)

    elseif flag == 2 # forward diff dR/du
    #=
      # define nested function
      function dRdu_rk4_wrapper(u_vals::AbstractVector, res_vec::AbstractVector)
        eqn.q_vec = u_vals
        eqn.q_vec = res_vec
        rk4(evalResidual, delta_t, t_max, mesh, sbp, eqn)
        return nothing
      end

      # use ForwardDiff package to generate function that calculate jacobian
      calcdRdu! = forwarddiff_jacobian!(dRdu_rk4_wrapper, Float64,
                  fadtype=:dual, n = mesh.numDof, m = mesh.numDof)

      jac = zeros(Float64, mesh.numDof, mesh.numDof)  # array to be populated
      calcdRdu!(eqn.q_vec, jac)
    =#
    elseif flag == 3 # calculate dRdx

      # dRdx here

    elseif flag == 4 || flag == 5 || flag == 11
      @time newton(evalResidual, mesh, sbp, eqn, opts, pmesh, itermax=opts["itermax"],
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
      delta_t = opts["delta_t"]
      t_max = opts["t_max"]

      t = rk4(evalResidual, delta_t, t_max, eqn.q_vec, eqn.res_vec, pre_func, post_func, (mesh, sbp, eqn), opts, majorIterationCallback=eqn.majorIterationCallback, real_time=opts["real_time"])

    elseif flag == 10
      function test_pre_func(mesh, sbp, eqn, opts)

        eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
      end

      function test_post_func(mesh, sbp, eqn, opts, calc_norm=true)
        return calcNorm(eqn, eqn.res_vec)
      end
      delta_t = opts["delta_t"]
      t_max = opts["t_max"]


      t = rk4(evalResidual, delta_t, t_max, eqn.q_vec, eqn.res_vec, test_pre_func,
          test_post_func, (mesh, sbp, eqn), opts, 
          majorIterationCallback=eqn.majorIterationCallback, real_time=opts["real_time"])


    elseif flag == 20

      @time t = crank_nicolson(evalResidual, opts["delta_t"], opts["t_max"],
                               mesh, sbp, eqn, opts, opts["res_abstol"],
                               opts["real_time"])

    elseif flag == 30  # lserk54

      t = lserk54(evalResidual, opts["delta_t"], opts["t_max"], eqn.q_vec, eqn.res_vec, (mesh, sbp, eqn), opts, eqn.params.time, majorIterationCallback=eqn.majorIterationCallback, res_tol=opts["res_abstol"], real_time=opts["real_time"])


    elseif flag == 40  # predictor-corrector newton

      predictorCorrectorHomotopy(evalResidual, evalHomotopy, mesh, sbp, eqn, opts, pmesh=pmesh)

    elseif flag == 41  # special mode: use regular Newton to solve homotopy

     @time newton(evalHomotopy, mesh, sbp, eqn, opts, pmesh, itermax=opts["itermax"],
                   step_tol=opts["step_tol"], res_abstol=opts["res_abstol"],
                   res_reltol=opts["res_reltol"], res_reltol0=opts["res_reltol0"])



    else
       throw(ErrorException("No flag specified: no solve will take place"))

    end       # end of if/elseif blocks checking flag

    println("total solution time printed above")
    params = eqn.params
    myrank = mesh.myrank

    if opts["write_timing"]
      MPI.Barrier(mesh.comm)
      fname = "timing_breakdown_$myrank"
      write_timings(params.time, fname)
    end

    # evaluate residual at final q value
    need_res = true
    if need_res
      eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
      # this will make sure the t value is stored into the equation object
      # this is important for calculating error norms later, to make sure
      # they exact solution is calculated at the right time
      # TODO: better way to update final time
      evalResidual(mesh, sbp, eqn, opts, t)

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

    myrank = mesh.myrank
  #  f = open("profile_$myrank.dat", "a+")
  #  Profile.print(f, format=:flat, C=true)
  #  close(f)

    saveSolutionToMesh(mesh, real(eqn.q_vec))
    writeVisFiles(mesh, "solution_done")

  end  # end if (opts[solve])

  return nothing
end
