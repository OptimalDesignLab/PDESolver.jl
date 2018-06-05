# define IC function, BC functor here
function TestIC(mesh, sbp, eqn, opts, q_vec)
  println("hello world")
end

mutable struct TestBCType <: BCType
end

function (obj::TestBCType)()
  println("hello edge of the world")
end

function test_frontend()

  @testset "----- Testing PDESolver Module frontend -----" begin
    physics_name = "Advection"
    @test ( haskey(PDESolver.PhysicsModDict, physics_name) )== true
    mod, _createObjects, _checkOptions = retrieve_physics(physics_name)
    # FactCheck doesn't compare functions correctly, so use == instead
    @test ( mod == AdvectionEquationMod )== true
    @test ( _createObjects == AdvectionEquationMod.createObjects )== true
    @test ( _checkOptions == AdvectionEquationMod.checkOptions )== true


    # test run_solver()
    fname = "input_vals_channel.jl"
    mesh, sbp, eqn, opts = solvePDE(fname)
    mesh, sbp, eqn2, opts = run_solver(fname)

    @test isapprox( eqn.q_vec, eqn2.q_vec) atol=1e-13

    ic_name = "TestIC"
    registerIC(AdvectionEquationMod, ic_name, TestIC)
    @test ( haskey(AdvectionEquationMod.ICDict, ic_name) )== true
    @test ( AdvectionEquationMod.ICDict[ic_name] == TestIC )== true

    bc_name = "TestBC"
    registerBC(AdvectionEquationMod, bc_name, TestBCType())
    @test ( haskey(AdvectionEquationMod.BCDict, bc_name) )== true
    @test ( AdvectionEquationMod.BCDict[bc_name] )== TestBCType()

  end  # end facts block

  return nothing
    
end

add_func1!(AdvectionTests, test_frontend, [TAG_FRONTEND, TAG_SHORTTEST])
