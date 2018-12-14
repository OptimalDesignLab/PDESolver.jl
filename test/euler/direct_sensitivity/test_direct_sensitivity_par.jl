

#------------------------------------------------------------------------------
# must include this for calcCd use below
include("calcCd-for_CI.jl")

@doc """
Euler Equation -- test_direct_sensitivity

The function tests for the correctness of the direct sensitivity computations.
  
This is a serial test and uses the input files:
RK4_DS/input_naca0012_rk4_DS.jl
RK4_FDbase/input_naca0012_rk4_FDbase.jl
RK4_FDpos/input_naca0012_rk4_FDpos.jl

"""->

function test_direct_sensitivity()

  cd("direct_sensitivity")
  @testset "--- Testing Direct Sensitivities, Parallel 2 process ---" begin

    @test isdefined(:calcCd)      # checks to make sure calcCd is defined

    @testset "Checking RK4 Direct Sensitivity with Finite-Difference check" begin

      #------------------------------------------------------------------------------
      # DS
      cd("RK4_DS")
      fname = "input_naca0012_rk4_DS_par2.jl"
      mesh, sbp, eqn, opts = PDESolver.run_solver(fname)
      @assert mesh.isDG == true
      @assert opts["perturb_Ma"] == true

      fh = open("total_dCddM.dat")
      RK4_dCddM_DS = 99.0
      for line in eachline(fh)
        m = match(r" total dCd/dM:\s+(\S+)", line)
        if m != nothing
          RK4_dCddM_DS = float(m.captures[1])
        end
      end
      @test ! isapprox(RK4_dCddM_DS, 99.0, atol=1e-14)      # make sure Cd has been read and set properly
      cd("..")


      #=
      What we should see:
        RK4, 3 timesteps, 2 deg AoA, Ma = 0.25, FD pert = 1e-8:
        Finite-difference: -0.8263592449209156
        Direct-sensitivity: -0.8263592837861967

        FD, par2: -0.8263590728363468   (?)
      =#

      #------------------------------------------------------------------------------
      # FD
      cd("RK4_FDpos")
      fname = "input_naca0012_rk4_FDpos_par2.jl"
      mesh, sbp, eqn, opts = PDESolver.run_solver(fname)
      @assert mesh.isDG == true
      @assert opts["perturb_Ma"] == false
      cd("..")

      cd("RK4_FDbase")
      fname = "input_naca0012_rk4_FDbase_par2.jl"
      mesh, sbp, eqn, opts = PDESolver.run_solver(fname)
      @assert mesh.isDG == true
      @assert opts["perturb_Ma"] == false
      cd("..")

      FDbase_Ma, FDbase_Cd = calcCd("RK4_FDbase")
      FDpos_Ma, FDpos_Cd = calcCd("RK4_FDpos")
      RK4_dCddM_FD = (FDpos_Cd - FDbase_Cd)/(FDpos_Ma - FDbase_Ma)

      #------------------------------------------------------------------------------
      # The check
      @test isapprox(RK4_dCddM_FD, RK4_dCddM_DS, atol=1e-6)   # for some reason, parallel needs to be at 1e-6, not 1e-7

    end   # End @testset "Checking RK4 Direct Sensitivity with Finite-Difference check"

  end   # End @testset "--- Testing Direct Sensitivities, Parallel 2 process ---"

  cd("..")

end # End function test_direct_sensitivity

add_func1!(EulerTests, test_direct_sensitivity, [TAG_NLSOLVERS, TAG_DS, TAG_SHORTTEST])
