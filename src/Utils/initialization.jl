"""
### Utils.createMeshAndOperator

  Create a SBP operator and a mesh.

  Inputs:
    opts: options dictonary
    dofpernode: number of degrees of freedom on each node

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

