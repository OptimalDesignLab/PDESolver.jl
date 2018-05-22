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
    array1DTo3D(mesh, sbp, eqn, opts, eqn.q_vec, eqn.q)

    # set BCs
    opts["BC1_name"] = string("p", d, "BC")
    AdvectionEquationMod.getBCFunctors(mesh, sbp, eqn, opts)

    # set source
    opts["SRCname"] = string("SRCp", d)
    AdvectionEquationMod.getSRCFunctors(mesh, sbp, eqn, opts)

    # compute residual
    fill!(eqn.res, 0.0)
    AdvectionEquationMod.evalResidual(mesh, sbp, eqn, opts)
    array3DTo1D(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)
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
    operator_names = ["SBPOmega", "SBPGamma", "SBPDiagonalE", "SBPDiagonalE2",
                      "SBPOmega3"]
    max_degree =     [         4,          4,              4,              2, 
                               4,]

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


  end  # end facts block

  return nothing
end  # end function

#test_mms()
add_func1!(AdvectionTests, test_mms, [TAG_MMS, TAG_SHORTTEST])


