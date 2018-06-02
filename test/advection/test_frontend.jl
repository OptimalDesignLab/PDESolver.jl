# define IC function, BC functor here
function TestIC(mesh, sbp, eqn, opts, q_vec)
  println("hello world")
end

type TestBCType <: BCType
end

function (obj::TestBCType)()
  println("hello edge of the world")
end

function test_frontend()

  facts("----- Testing PDESolver Module frontend -----") do
    physics_name = "Advection"
    @fact haskey(PDESolver.PhysicsModDict, physics_name) --> true
    mod, _createObjects, _checkOptions = retrieve_physics(physics_name)
    # FactCheck doesn't compare functions correctly, so use == instead
    @fact mod == AdvectionEquationMod --> true
    @fact _createObjects == AdvectionEquationMod.createObjects --> true
    @fact _checkOptions == AdvectionEquationMod.checkOptions --> true


    # test run_solver()
    fname = "input_vals_channel.jl"
    mesh, sbp, eqn, opts = solvePDE(fname)
    mesh, sbp, eqn2, opts = run_solver(fname)

    @fact eqn.q_vec --> roughly(eqn2.q_vec, atol=1e-13)

    ic_name = "TestIC"
    registerIC(AdvectionEquationMod, ic_name, TestIC)
    @fact haskey(AdvectionEquationMod.ICDict, ic_name) --> true
    @fact AdvectionEquationMod.ICDict[ic_name] == TestIC --> true

    bc_name = "TestBC"
    registerBC(AdvectionEquationMod, bc_name, TestBCType())
    @fact haskey(AdvectionEquationMod.BCDict, bc_name) --> true
    @fact AdvectionEquationMod.BCDict[bc_name] --> TestBCType()

  end  # end facts block

  return nothing
    
end

add_func1!(AdvectionTests, test_frontend, [TAG_FRONTEND, TAG_SHORTTEST])
