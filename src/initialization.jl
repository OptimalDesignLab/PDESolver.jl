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
  # flag determines whether to calculate u, dR/du, or dR/dx (1, 2, or 3)
  flag = opts["run_type"]
  if haskey(opts, "jac_method")
    jac_method = opts["jac_method"]
  end
  dim = opts["dimensions"]

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
  elseif flag == 5  
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
    throw(ErrorException("Illegal flag or jac_method combination specified in input."))
  end
  # If the user specifies a flag other than the ones within the above if checks, 
  #   then an error will be thrown now because Tsol is not defined
  opts["Tsol"] = Tsol
  opts["Tres"] = Tres
  opts["Tsbp"] = Tsbp
  opts["Tmsh"] = Tmsh

  # figure out reorder, internal args for SBP, shape_type for Pumi
  # should shape_type live here or be encapsulated in Pumi?
  if opts["use_DG"]
    if opts["operator_type"] == "SBPOmega"
      reorder = false
      internal = true
      shape_type = 2
    elseif opts["operator_type"] == "SBPGamma"
      reorder = false
      internal = false
      shape_type = 3
    else
      op_type = opts["operator_type"]
      throw(ArgumentError("unrecognized operator type $op_type for DG mesh"))
    end
  else  # CG mesh
    if opts["operator_type"] != "SBPGamma"
      op_type = opts["operator_type"]
      throw(ArgumentError("invalid operator type $op_type for CG"))
    end
    # the CG varient of SBP gamma is the only option
    reorder = true
    internal = false
    shape_type = 1
  end


  mesh_time = @elapsed if opts["use_DG"]
    println("\nConstructing SBP Operator")
    # create DG SBP operator with internal nodes only
    if dim == 2
      sbp = TriSBP{Float64}(degree=order, reorder=reorder, internal=internal)
      # TODO: use sbp.vtx instead
      ref_verts = [-1. 1 -1; -1 -1 1]
      interp_op = SummationByParts.buildinterpolation(sbp, ref_verts)
      sbpface = TriFace{Float64}(order, sbp.cub, ref_verts.')
    else 
      sbp = TetSBP{Float64}(degree=order, reorder=reorder, internal=internal)
      ref_verts = sbp.vtx
      interp_op = SummationByParts.buildinterpolation(sbp, ref_verts.')
      face_verts = SummationByParts.SymCubatures.getfacevertexindices(sbp.cub)
      topo = ElementTopology{3}(face_verts)
      sbpface = TetFace{Float64}(order, sbp.cub, ref_verts)
    end

    # create mesh with 4 dof per node

    println("constructing DG mesh")
    if dim == 2

      mesh = PumiMeshDG2{Tmsh}(dmg_name, smb_name, order, sbp, opts, interp_op, sbpface; dofpernode=dofpernode, coloring_distance=opts["coloring_distance"], shape_type=shape_type)
    else
      mesh = PumiMeshDG3{Tmsh}(dmg_name, smb_name, order, sbp, opts, interp_op, sbpface, topo; dofpernode=dofpernode, coloring_distance=opts["coloring_distance"], shape_type=shape_type)
    end
    if (opts["jac_type"] == 3 || opts["jac_type"] == 4) && opts["use_jac_precond"]
      @assert dim == 2
      pmesh = PumiMeshDG2Preconditioning(mesh, sbp, opts; 
                     coloring_distance=opts["coloring_distance_prec"])
    else
      pmesh = mesh
    end

  else  # continuous Galerkin
    # create SBP object
    println("\nConstructing SBP Operator")
    sbp = TriSBP{Float64}(degree=order, reorder=reorder, internal=internal)  # create linear sbp operator
    sbpface = TriFace{Float64}(order, sbp.cub, sbp.vtx)
    # create linear mesh with 4 dof per node

    println("constructing CG mesh")
    mesh = PumiMesh2{Tmsh}(dmg_name, smb_name, order, sbp, opts, sbpface; dofpernode=dofpernode, coloring_distance=opts["coloring_distance"], shape_type=shape_type)

    if opts["jac_type"] == 3 || opts["jac_type"] == 4
      pmesh = PumiMesh2Preconditioning(mesh, sbp, opts; coloring_distance=opts["coloring_distance_prec"])
    else
      pmesh = mesh
    end
  end

  return sbp, mesh, pmesh, Tsol, Tres, Tmsh, mesh_time
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
  if opts["solve"]
    
    solve_time = @elapsed if flag == 1 # normal run
      # RK4 solver
      delta_t = opts["delta_t"]
      t_max = opts["t_max"]
      @time rk4(evalResidual, delta_t, t_max, mesh, sbp, eqn, opts, 
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

    elseif flag == 4 || flag == 5
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

      rk4(evalResidual, delta_t, t_max, eqn.q_vec, eqn.res_vec, pre_func, post_func, (mesh, sbp, eqn), opts, majorIterationCallback=eqn.majorIterationCallback, real_time=opts["real_time"])

    elseif flag == 10
      function test_pre_func(mesh, sbp, eqn, opts)
        
        eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
      end

      function test_post_func(mesh, sbp, eqn, opts, calc_norm=true)
        return calcNorm(eqn, eqn.res_vec)
      end
      delta_t = opts["delta_t"]
      t_max = opts["t_max"]

      rk4(evalResidual, delta_t, t_max, eqn.q_vec, eqn.res_vec, test_pre_func,
          test_post_func, (mesh, sbp, eqn), opts, 
          majorIterationCallback=eqn.majorIterationCallback, real_time=opts["real_time"])

    elseif flag == 20

      @time t = crank_nicolson(evalResidual, opts["delta_t"], opts["t_max"], 
                               mesh, sbp, eqn, opts, opts["res_abstol"], 
                               opts["real_time"])

      eqn.t = t

  #   else
  #     throw(ErrorException("No flag specified: no solve will take place"))
  #     return nothing

    end       # end of if/elseif blocks checking flag

    println("total solution time printed above")
    params = eqn.params
    myrank = mesh.myrank

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
      evalResidual(mesh, sbp, eqn, opts, eqn.t)

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
  #  printSolution(mesh, real(eqn.q_vec))
  #  printCoordinates(mesh)
    writeVisFiles(mesh, "solution_done")

    # write timings
  #  timings = [params.t_volume, params.t_face, params.t_source, params.t_sharedface, params.t_bndry, params.t_send, params.t_wait, params.t_allreduce, params.t_jacobian, params.t_solve, params.t_barrier, params.t_barrier2, params.t_barrier3]
  #  writedlm("timing_breakdown_$myrank.dat", vcat(timings, params.t_barriers))
    fname = "timing_breakdown_$myrank"
    write_timings(params.time, fname)

  end  # end if (opts[solve])

  return nothing
end
