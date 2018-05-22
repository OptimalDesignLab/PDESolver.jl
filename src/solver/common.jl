"""
  Create a SBP operator and a mesh.  This is used by all physics modules
  to create the right type of operator and mesh based on the input options.
  It is type unstable, but that is ok.

  If the options dictionary specifies a second SBP operator type, a second
  mesh and SBP operator will be created and stored in the `mesh2` and `sbp2`

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

  # if there is a second mesh/sbp pair, construct it first
  # this is important because Pumi is going to finalize the first mesh
  # when the second mesh is created

  mesh_time = 0.0
  if opts["use_staggered_grid"]
    println("constructing flux grid")

    sbp2, sbpface, shape_type, topo = createSBPOperator(opts, Tsbp, 2)

    mesh_time = @elapsed mesh2, pmesh2 = createMesh(opts, sbp2, sbpface, 
                                                  shape_type, topo, Tmsh,
                                                  dofpernode, 2)
    if !(mesh2 === pmesh2)
      throw(ErrorException("preconditioning mesh not supported with staggered grids"))
    end

  end

  println("constructing solution grid")
  sbp, sbpface, shape_type, topo = createSBPOperator(opts, Tsbp)
 
  mesh_time += @elapsed mesh, pmesh = createMesh(opts, sbp, sbpface, shape_type,
                                                topo, Tmsh, dofpernode)

  # store the second mesh and SBP operator inside the first mesh
  if opts["use_staggered_grid"]
    mesh.mesh2 = mesh2
    mesh.sbp2 = sbp2

    # build the interpolation operators
    I_S2F, I_F2S = buildinterpolation(sbp, sbp2)
    mesh.I_S2F = I_S2F
    mesh.I_S2FT = I_S2F.'
    mesh.I_F2S = I_F2S
    mesh.I_F2ST = I_F2S.'
  end

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

  **Options Keys**

   * run_type
   * jac_method
   * force_solution_complex
   * force_mesh_complex
"""
function getDataTypes(opts::Dict)

  flag = opts["run_type"]
  if haskey(opts, "jac_method")  #TODO: this should be handled in read_input()
    jac_method = opts["jac_method"]
  end

  calc_jac_explicit = opts["calc_jac_explicit"]

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

      if calc_jac_explicit
        Tsol = Float64
        Tres = Float64
      else
        Tsol = Complex128
        Tres = Complex128
      end
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
      if calc_jac_explicit
        Tsol = Float64
      else
        Tsol = Complex128
      end
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
      if calc_jac_explicit
        Tsol = Float64
        Tres = Float64
      else
        Tsol = Complex128
        Tres = Complex128
      end
    else
      throw(ErrorException("Illegal or no jac_method specified for CN initialization."))
    end
  elseif flag == 660 # jac_method needs to be symbol
    if jac_method == 1 # Crank-Nicolson, unsteady adjoint, FD Jac
      Tmsh = Float64
      Tsbp = Float64
      Tsol = Float64
      Tres = Float64
    elseif jac_method == 2 # Crank-Nicolson, unsteady adjoint, CS Jac
      Tmsh = Float64
      Tsbp = Float64
      Tsol = Complex128
      Tres = Complex128
    else
      throw(ErrorException("Illegal or no jac_method specified for CN uadj initialization."))
    end
  else
    throw(ErrorException("Unrecognized run_type: $flag"))
  end

  if opts["force_solution_complex"]
    Tsol = Complex128
    Tres = Complex128
  end

  if opts["force_mesh_complex"]
    Tmsh = Complex128
  end

  return Tmsh, Tsbp, Tsol, Tres
end

"""
  This function constructs the SBP operator and the associated SBP face
  operator, as specified by the options dictionary.  It also determines
  the shape_type that PumiInterface uses to describe the SBP operator to
  Pumi.

  **Inputs**

   * opts: the options dictionary
   * Tsbp: the DataType specifying the Tsbp passed to the SBP operator
          constructor
   * suffix: this suffix is added to all keys accessed in the options dictionary.
            Usually the suffix is either the empty string or an integer.  This
            provides a convenient way for the input file to specify several
            different SBP operator and have this operator construct them.
            Default value is the empty string.

  **Outputs**

   * sbp: the SBP operator
   * sbpface: the SBP face operator
   * shape_type: an integer passed to the mesh constructor to describe the
                 operator
   * topo: in the 3D DG case, an ElementTopology describing the SBP reference
           element, otherwise the integer 0.

  **DG Operator Names**

   * SBPOmega: nodes on the interior of the element only
   * SBPGamma: nodes on the faces of the element and the interior (similar
               to Lagrange finite elements)
   * SBPDiagonalE: operator with diagonal E matrix, with nodes on vertices,
                   faces, and interior (similar to Lagrange FE)
   * SBPDiagonalE2: operator with diagonal E matrix, with nodes on faces
                    (but not vertices)
   * SBPOmega2: attempt at optimized SBP Omega-type operators, probably not
                working
   * SBPOmega3: SBP Omega-type operator with degree 2p cubature rule for
                all degree operators (p=1 and 2 are the same as SBPOmega),
                unlike `SBPOmega2`, not optimized

  **CG Operator Names**

   * SBPGamma: see above, this operator can be used to CG as well

"""
function createSBPOperator(opts::Dict, Tsbp::DataType, suffix="")
  # construct SBP operator and figure out shape_type needed by Pumi
  order = opts["order$suffix"]  # order of accuracy
  dim = opts["dimensions"]

  println("\nConstructing SBP Operator $suffix")
  topo = 0  # generally not needed, so return a default value
  if opts["use_DG"]
    if opts["operator_type$suffix"] == "SBPOmega"
      if dim == 2
        sbp = getTriSBPOmega0(degree=order, Tsbp=Tsbp)
      else
        sbp = getTetSBPOmega(degree=order, Tsbp=Tsbp)
      end
      shape_type = 2
    elseif opts["operator_type$suffix"] == "SBPGamma"
      if dim == 2
        sbp = getTriSBPGamma(degree=order, Tsbp=Tsbp)
      else
        sbp = getTetSBPGamma(degree=order, Tsbp=Tsbp)
      end
      shape_type = 3
    elseif opts["operator_type$suffix"] == "SBPDiagonalE"
      if dim == 2
        sbp = getTriSBPDiagE(degree=order, Tsbp=Tsbp)
      else
        sbp = getTetSBPDiagE(degree=order, Tsbp=Tsbp)
#        throw(ArgumentError("3 dimensional SBPDiagonalE no supported"))
      end
      shape_type = 4
    elseif opts["operator_type$suffix"] == "SBPDiagonalE2"  # no vert nodes
      if dim == 2
        sbp = getTriSBPDiagE(degree=order, Tsbp=Tsbp, vertices=false)
      else
        sbp = getTetSBPDiagE(degree=order, Tsbp=Tsbp)
      end
      shape_type = 5
    elseif opts["operator_type$suffix"] == "SBPOmega2"
      shape_type = 6
      if dim == 2
        sbp = getTriSBPOmega2(degree=order, Tsbp=Tsbp)
      else
        throw(ArgumentError("3D SBPOmega2 not supported"))
      end
    elseif opts["operator_type$suffix"] == "SBPOmega3"
      shape_type = 7
      if dim == 2
        sbp = getTriSBPOmega(degree=order, Tsbp=Tsbp)
      else
        throw(ArgumentError("3D SBPOmega3 not supported"))
      end
    else
      op_type = opts["operator_type$suffix"]
      throw(ArgumentError("unrecognized operator type $op_type for DG mesh"))
    end
  else  # CG mesh
    # the CG varient of SBP gamma is the only option
    if opts["operator_type$suffix"] != "SBPGamma"
      op_type = opts["operator_type$suffix"]
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
      if opts["operator_type$suffix"] == "SBPDiagonalE"
        sbpface = getTriFaceForDiagE(order, sbp.cub, ref_verts.')
      elseif opts["operator_type$suffix"] == "SBPDiagonalE2"
        println("getting TriFaceForDiagE2")
        sbpface = getTriFaceForDiagE(order, sbp.cub, ref_verts.', vertices=false)
      else
        sbpface = TriFace{Tsbp}(order, sbp.cub, ref_verts.')
      end
    else  # dim == 3
      ref_verts = sbp.vtx
      face_verts = SummationByParts.SymCubatures.getfacevertexindices(sbp.cub)
      edge_verts = [1 2 1 1 2 3;  # TODO: SBP should provide this
                    2 3 3 4 4 4]
      topo2 = ElementTopology2()   #TODO: make this the correct one for SBP
      topo = ElementTopology{3}(face_verts, edge_verts, topo2=topo2)

      if opts["operator_type$suffix"] == "SBPDiagonalE"
        sbpface = getTetFaceForDiagE(order, sbp.cub, ref_verts)
      else
        sbpface = TetFace{Tsbp}(order, sbp.cub, ref_verts)
      end
    end  # end if dim == 2
  else   # CG
    if dim == 2
      sbpface = TriFace{Tsbp}(order, sbp.cub, sbp.vtx)
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
    suffix: suffix added to options dictionary keys that describe the SBP
            operator.  See [`createSBPOperator`](@ref)

  All arguments except opts are typically provided by 
  [`createSBPOperator`](@ref) and [`getDataTypes`](@ref)
"""
function createMesh(opts::Dict, sbp::AbstractSBP, sbpface, shape_type, topo,
                    Tmsh, dofpernode, suffix="")

  dmg_name = opts["dmg_name"]
  smb_name = opts["smb_name"]
  # the mesh constructor looks at the order option, so temporarily make
  # order = order$suffix
  order_orig = opts["order"]
  opts["order"] = opts["order$suffix"]
  order = opts["order$suffix"]  # order of accuracy
  dim = opts["dimensions"]


  if opts["use_DG"]
    println("constructing DG mesh")
    if dim == 2
      mesh = PumiMeshDG2(Tmsh, sbp, opts, sbpface, dofpernode=dofpernode,
                                                   shape_type=shape_type)
#      mesh = PumiMeshDG2{Tmsh, typeof(sbpface)}(dmg_name, smb_name, order, sbp, opts, sbpface; 
#                               dofpernode=dofpernode, 
#                               coloring_distance=opts["coloring_distance"],
#                               shape_type=shape_type)
    else
      mesh = PumiMeshDG3(Tmsh, sbp, opts, sbpface, topo, dofpernode=dofpernode,
                                                   shape_type=shape_type)
#      mesh = PumiMeshDG3{Tmsh, typeof(sbpface)}(dmg_name, smb_name, order, sbp, opts, sbpface,
#                               topo; dofpernode=dofpernode,
#                               coloring_distance=opts["coloring_distance"],
#                               shape_type=shape_type)
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
      mesh = PumiMesh2{Tmsh, typeof(sbpface)}(dmg_name, smb_name, order, sbp, opts, sbpface;
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


  # reset order
  opts["order"] = order_orig

  return mesh, pmesh
end

"""
  This function is used by all physics modules to load the most recently
  saved state when restarting.

  **Inputs**

   * mesh: the mesh
   * sbp: AbstractSBP
   * eqn: AbstractSolutionData, eqn.q_vec is overwritten with the saved state
   * opts: options dictionary

   The keys described in the [`Checkpointer`](@ref Utils.Checkpointer)
   documentation are used to determine the most recent complete checkpoint.

   Implementation notes:
     currently pmesh isn't used for anything because checkpointing does not
     support mesh adaptation.  When this changes, this function will have to
     be updated.
"""
function loadRestartState(mesh::AbstractMesh, sbp::AbstractSBP,
                        eqn::AbstractSolutionData, opts::Dict,
                        pmesh::AbstractMesh=mesh)

  chkpointer = Checkpointer(opts, mesh.myrank)
  loadLastCheckpoint(chkpointer, mesh, sbp, eqn, opts)

  return nothing
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

    t_nlsolve = @elapsed if flag == 1 # normal run
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
      @time newton(evalResidual, mesh, sbp, eqn, opts, pmesh)

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

        array1DTo3D(mesh, sbp, eqn, opts, eqn.q_vec, eqn.q)
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

     @time newton(evalHomotopy, mesh, sbp, eqn, opts, pmesh)

    elseif flag == 660    # Unsteady adjoint crank nicolson code. DOES NOT PRODUCE CORRECT RESULTS. See Anthony.
      # error("Unsteady adjoint Crank-Nicolson code called.\nThis code does run, but incorrect numerical results are obtained.\nTo run this, you must comment out this error message in initialization.jl.\n\n")

      if opts["adjoint_revolve"]
        error("adjoint_revolve not fully implemented yet.")
      end

      if opts["adjoint_straight"]

        if opts["uadj_global"]

          println(" GLOBAL: forming WWW, ZZZ")
          # dof_global = mesh.numDof*t_steps
          # blksz = 3   # testing 44
          blksz = mesh.numDof
          t_steps = 4
          dof_global = blksz*t_steps
          WWW = rand(dof_global)
          ZZZ = rand(dof_global)
          println(" GLOBAL: forming dRdu")
          dRdu_global_fwd = zeros(dof_global, dof_global)
          dRdu_global_rev = zeros(dof_global, dof_global)
          # PM stands for piecemeal. intended to do it step by step during the adj calc to test bookkeeping
          dRdu_global_rev_PM = zeros(dof_global, dof_global)
          println(" GLOBAL: size(dRdu_global_fwd): ", size(dRdu_global_fwd))
          writedlm("global_www.dat", WWW)
          writedlm("global_zzz.dat", ZZZ)
          writedlm("global_dRdu_fwd_initial.dat", dRdu_global_fwd)
          writedlm("global_dRdu_rev_initial.dat", dRdu_global_rev)
        else
          WWW = zeros(1,1)
          ZZZ = zeros(1,1)
          dRdu_global_fwd = zeros(1,1)
          dRdu_global_rev = zeros(1,1)
          dRdu_global_rev_PM = zeros(1,1)
        end      # end of if opts["uadj_global"]

        # forward sweep
        # @time t = crank_nicolson(evalResidual, opts["delta_t"], opts["t_max"],
                                 # mesh, sbp, eqn, opts, opts["res_abstol"], store_u_to_disk=true)
        println(" Calling CN, forward sweep.")
        @time t = crank_nicolson_uadj(evalResidual, opts["delta_t"], opts["t_max"],
                                 mesh, sbp, eqn, opts,
                                 WWW, ZZZ, dRdu_global_fwd, dRdu_global_rev, dRdu_global_rev_PM,
                                 opts["res_abstol"],
                                 store_u_to_disk=true)

        # reverse sweep
        # @time t = crank_nicolson(evalResidual, opts["delta_t"], opts["t_max"],
                                 # mesh, sbp, eqn, opts, opts["res_abstol"], neg_time=true)
        println(" Calling CN, reverse sweep.")
        @time t = crank_nicolson_uadj(evalResidual, opts["delta_t"], opts["t_max"],
                                 mesh, sbp, eqn, opts,
                                 WWW, ZZZ, dRdu_global_fwd, dRdu_global_rev, dRdu_global_rev_PM,
                                 opts["res_abstol"],
                                 neg_time=true)

        if opts["uadj_global"]
          dRdu_global_fwd_WWW = dRdu_global_fwd*WWW
          fwd_check_number = dot(dRdu_global_fwd_WWW, ZZZ)
          filename = "global_dRdu_check_fwd.dat"
          writedlm(filename, fwd_check_number)

          dRdu_global_rev_ZZZ = dRdu_global_rev*ZZZ
          rev_check_number = dot(dRdu_global_rev_ZZZ, WWW)
          filename = "global_dRdu_check_rev.dat"
          writedlm(filename, rev_check_number)
        end     # end of if opts["uadj_global"]

      else
        @time t = crank_nicolson_uadj(evalResidual, opts["delta_t"], opts["t_max"],
                                 mesh, sbp, eqn, opts, opts["res_abstol"])
      end      # end of if opts["adjoint_straight"]

      eqn.t = t

    else
       throw(ErrorException("No flag specified: no solve will take place"))

    end       # end of if/elseif blocks checking flag

    println("total solution time printed above")
    params = eqn.params
    params.time.t_nlsolve += t_nlsolve
    myrank = mesh.myrank

    if opts["write_timing"]
      MPI.Barrier(mesh.comm)
      fname = "timing_breakdown_$myrank"
      write_timings(params.time, fname)
    end

    # evaluate residual at final q value
    need_res = true
    if need_res
      array1DTo3D(mesh, sbp, eqn, opts, eqn.q_vec, eqn.q)
      # this will make sure the t value is stored into the equation object
      # this is important for calculating error norms later, to make sure
      # they exact solution is calculated at the right time
      # TODO: better way to update final time
      evalResidual(mesh, sbp, eqn, opts, t)

      eqn.res_vec[:] = 0.0
      array3DTo1D(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
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

