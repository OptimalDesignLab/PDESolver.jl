function make_input_mms(opts0, degree; dg=false, operator="SBPOmega")
  opts = copy(opts0)
  opts["order"] = degree
  opts["IC_name"] = "ICp$degree"
  opts["BC1_name"] = string("p", degree, "BC")
  opts["SRCname"] = "SRCp$degree"
  opts["operator_type"] = operator

  if dg
    opts["use_DG"] = true
    opts["Flux_name"] = "LFFlux"
  else
    opts["use_DG"] = false
  end

  read_input(opts)  # get default values

#  fname = "input_vals_mms$degree.jl"
#  rmfile(fname)
#  f = open(fname, "w")
#  print(f, "arg_dict = ")  # TODO: change this to use make_input
#  println(f, opts)
#  ARGS[1] = fname
#  close(f)
  return opts
end

"""
  Test that the residual is zero for polynomials up to a specified degree

  **Inputs**

   * opts: options dictionary
   * maxdegree: maximum degree polynomials to test
"""
function test_mms_poly(opts::Dict, maxdegree::Integer)

  mesh, sbp, eqn, opts = createObjects(opts)

  for d=0:maxdegree
    println("testing degree ", d, " polynomials")

    # get IC
    AdvectionEquationMod.ICDict["ICp$d"](mesh, sbp, eqn, opts, eqn.q_vec)
    disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)

    # set BCs
    opts["BC1_name"] = string("p", d, "BC")
    AdvectionEquationMod.getBCFunctors(mesh, sbp, eqn, opts)

    # set source
    opts["SRCname"] = string("SRCp", d)
    AdvectionEquationMod.getSRCFunctors(mesh, sbp, eqn, opts)

    # compute residual
    fill!(eqn.res, 0.0)
    AdvectionEquationMod.evalResidual(mesh, sbp, eqn, opts)
    assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
    @fact calcNorm(eqn, eqn.res_vec) --> roughly(0.0, atol=1e-12)
#    @fact eqn.res_vec --> roughly(zeros(mesh.numDof), atol=1e-12)
  end

  return nothing
end


"""
  Test manufactured solutions, including cases that should be exact, for both
  gamma and omega operators.
"""
function test_mms()
  facts("----- Testing using manufactured polynomials -----") do

    # test 2d DG
    println("  ----- Testing 2d -----")
    opts0 = read_input_file("input_vals_mms.jl")
    operator_names = ["SBPOmega", "SBPGamma", "SBPDiagonalE", "SBPDiagonalE2"]
    max_degree =     [         4,          4,              4,              2 ]

    for i=1:length(operator_names)
      for degree=1:max_degree[i]
        println(  "  ----- testing degree $degree $(operator_names[i])")
        opts = make_input_mms(opts0, degree, dg=true, operator=operator_names[i])
        test_mms_poly(opts, degree)
      end
    end

    println("  ----- Testing 3d -----")
    opts0 = read_input_file("input_vals_3d.jl")
    operator_names = ["SBPOmega", "SBPGamma", "SBPDiagonalE"]
    max_degree =     [         4,          4,              2]

    for i=1:length(operator_names)
      for degree=1:max_degree[i]
        println(  "  ----- testing degree $degree $(operator_names[i])")
        opts = make_input_mms(opts0, degree, dg=true, operator=operator_names[i])
        test_mms_poly(opts, degree)
      end
    end


#=

    println("  -----testing degree 1 polynomial -----")
    resize!(ARGS, 1)
    #TODO: uncomment when SBP is fixed
    
    ARGS[1] = "input_vals_mms.jl"
    mesh, sbp, eqn, opts = solvePDE(ARGS[1])
    #=
    fill!(eqn.res, 0.0)
    disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
    println("eqn.q_vec = ", eqn.q_vec)
    AdvectionEquationMod.evalResidual(mesh, sbp, eqn, opts)
    assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
    @fact eqn.res_vec --> roughly(zeros(mesh.numDof), atol=1e-12)
    =#

    ARGS[1] = "input_vals_mms1.jl"
    fname = make_input_mms(opts, 1, dg=true)
    mesh, sbp, eqn, opts = solvePDE(ARGS[1])
    fill!(eqn.res, 0.0)
    disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
    AdvectionEquationMod.evalResidual(mesh, sbp, eqn, opts)
    assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
    @fact eqn.res_vec --> roughly(zeros(mesh.numDof), atol=1e-12)

    fname = make_input_mms(opts, 1, dg=true, operator="SBPOmega")
    mesh, sbp, eqn, opts = solvePDE(ARGS[1])
    fill!(eqn.res, 0.0)
    disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
    AdvectionEquationMod.evalResidual(mesh, sbp, eqn, opts)
    assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
    @fact eqn.res_vec --> roughly(zeros(mesh.numDof), atol=1e-12)

    println("  -----testing degree 2 polynomial -----")
    #=
    fname = make_input_mms(opts, 2)
    ARGS[1] = fname
    mesh, sbp, eqn, opts = solvePDE(ARGS[1])
    fill!(eqn.res, 0.0)
    disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
    AdvectionEquationMod.evalResidual(mesh, sbp, eqn, opts)
    assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
    @fact eqn.res_vec --> roughly(zeros(mesh.numDof), atol=1e-12)
    =#
    fname = make_input_mms(opts, 2, dg=true)
    mesh, sbp, eqn, opts = solvePDE(ARGS[1])
    fill!(eqn.res, 0.0)
    disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
    AdvectionEquationMod.evalResidual(mesh, sbp, eqn, opts)
    assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
    @fact eqn.res_vec --> roughly(zeros(mesh.numDof), atol=1e-11)



    println("  -----testing degree 3 polynomial -----")
    #=
    fname = make_input_mms(opts, 3)
    ARGS[1] = fname
    mesh, sbp, eqn, opts = solvePDE(ARGS[1])
    fill!(eqn.res, 0.0)
    disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
    AdvectionEquationMod.evalResidual(mesh, sbp, eqn, opts)
    assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
    @fact eqn.res_vec --> roughly(zeros(mesh.numDof), atol=1e-12)
    =#
    fname = make_input_mms(opts, 3, dg=true)
    mesh, sbp, eqn, opts = solvePDE(ARGS[1])
    fill!(eqn.res, 0.0)
    disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
    AdvectionEquationMod.evalResidual(mesh, sbp, eqn, opts)
    assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
    @fact eqn.res_vec --> roughly(zeros(mesh.numDof), atol=1e-12)



    println("  -----testing degree 4 polynomial -----")
    #=
    fname = make_input_mms(opts, 4)
    ARGS[1] = fname
    mesh, sbp, eqn, opts = solvePDE(ARGS[1])
    fill!(eqn.res, 0.0)
    disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
  #  println("eqn.q_vec = ", eqn.q_vec)
    AdvectionEquationMod.evalResidual(mesh, sbp, eqn, opts)
    assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
    @fact eqn.res_vec --> roughly(zeros(mesh.numDof), atol=1e-12)
    =#

    fname = make_input_mms(opts, 4, dg=true)
    mesh, sbp, eqn, opts = solvePDE(ARGS[1])
    fill!(eqn.res, 0.0)
    disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
    AdvectionEquationMod.evalResidual(mesh, sbp, eqn, opts)
    assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
    @fact eqn.res_vec --> roughly(zeros(mesh.numDof), atol=1e-12)
=#
#=
    println("\n  -----testing 3D -----")
    ARGS[1] = "input_vals_3d.jl"
    mesh, sbp, eqn, opts = solvePDE(ARGS[1])

    println("  -----testing degree 1 polynomial-----")

    make_input_mms(opts, 1, dg=true)
    mesh, sbp, eqn, opts = solvePDE(ARGS[1])
    fill!(eqn.res, 0.0)
    disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
    AdvectionEquationMod.evalResidual(mesh, sbp, eqn, opts)
    assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
    @fact eqn.res_vec --> roughly(zeros(mesh.numDof), atol=1e-12)

    make_input_mms(opts, 1, dg=true, operator="SBPGamma")
    mesh, sbp, eqn, opts = solvePDE(ARGS[1])
    fill!(eqn.res, 0.0)
    disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
    AdvectionEquationMod.evalResidual(mesh, sbp, eqn, opts)
    assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
    @fact eqn.res_vec --> roughly(zeros(mesh.numDof), atol=1e-12)



    println("  -----testing degree 2 polynomial-----")
    make_input_mms(opts, 2, dg=true)
    mesh, sbp, eqn, opts = solvePDE(ARGS[1])
    fill!(eqn.res, 0.0)
    disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
    AdvectionEquationMod.evalResidual(mesh, sbp, eqn, opts)
    assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
    for i=1:mesh.numDof
      @fact eqn.res_vec[i] --> roughly(0.0, atol=1e-11)
    end

    println("testing gamma")
    make_input_mms(opts, 2, dg=true, operator="SBPGamma")
    mesh, sbp, eqn, opts = solvePDE(ARGS[1])
    fill!(eqn.res, 0.0)
    disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
    AdvectionEquationMod.evalResidual(mesh, sbp, eqn, opts)
    assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
    for i=1:length(eqn.res_vec)
      @fact eqn.res_vec[i] --> roughly(0.0, atol=1e-11)
    end

    println("  ----- testing degree 3 polynomial-----")
    println("testing gamma")
    make_input_mms(opts, 3, dg=true, operator="SBPGamma")
    mesh, sbp, eqn, opts = solvePDE(ARGS[1])
    fill!(eqn.res, 0.0)
    disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
    AdvectionEquationMod.evalResidual(mesh, sbp, eqn, opts)
    assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
    for i=1:length(eqn.res_vec)
      @fact eqn.res_vec[i] --> roughly(0.0, atol=1e-11)
    end

    println("  ----- testing degree 4 polynomial-----")
    println("testing gamma")
    make_input_mms(opts, 4, dg=true, operator="SBPGamma")
    mesh, sbp, eqn, opts = solvePDE(ARGS[1])
    fill!(eqn.res, 0.0)
    disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
    AdvectionEquationMod.evalResidual(mesh, sbp, eqn, opts)
    assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
    for i=1:length(eqn.res_vec)
      @fact eqn.res_vec[i] --> roughly(0.0, atol=1e-9)
    end
=#
  end  # end facts block

  return nothing
end  # end function

#test_mms()
add_func1!(AdvectionTests, test_mms, [TAG_MMS, TAG_SHORTTEST])


