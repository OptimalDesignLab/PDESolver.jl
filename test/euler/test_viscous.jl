#Test functional Integrate and adjoint for euler equation.

@doc """
Euler Equation -- test_viscous

The function tests for the correctness of the viscous terms computation.
  
This is a serial test and uses the input file called

`input_viscous_polynomial_ser.jl`

"""->

function test_viscous()

  @testset "--- Testing Viscous Terms ---" begin

    ARGS[1] = "input_viscous_polynomial_ser.jl"
    # mesh, sbp, eqn, opts, pmesh = EulerEquationMod.createObjects(ARGS[1])
    mesh, sbp, eqn, opts = PDESolver.run_solver(ARGS[1])
    @assert mesh.isDG == true
    @assert opts["isViscous"] == true

    @testset "Checking Solution at Final Time Step" begin

      viscous_qvec_final_for_comparison = readdlm("viscous_files/solution_viscous_forcomparison_0.dat")

      for dof_ix = 1:length(eqn.q_vec)
        diff = eqn.q_vec[dof_ix] - viscous_qvec_final_for_comparison[dof_ix]
        @test isapprox( diff, 0.0) atol=1e-13
      end


    end # End  testset("Checking Solution at Final Time Step")

    @testset "Checking Diffusion Tensor" begin

      Tsol = opts["Tsol"]
      Tdim = mesh.dim
      GtL = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode, Tdim, Tdim, mesh.numNodesPerFace)
      GtR = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode, Tdim, Tdim, mesh.numNodesPerFace)

      # checking a couple interfaces
      interfaces = sview(mesh.interfaces, :)
      nfaces = length(mesh.interfaces)

      for interface = [3 4]
        fill!(GtL, 0.0)
        fill!(GtR, 0.0)

        q_faceL = Base.view(eqn.q_face, :, 1, :, interface)
        q_faceR = Base.view(eqn.q_face, :, 2, :, interface)

        EulerEquationMod.calcDiffusionTensor(eqn.params, q_faceL, GtL)
        EulerEquationMod.calcDiffusionTensor(eqn.params, q_faceR, GtR)

        GtL_vec = reshape(GtL, 128, )
        GtR_vec = reshape(GtR, 128, )

        GtL_file_to_check_against = string("viscous_files/viscous_test_face", interface, "_GtL_vec.dat")
        GtR_file_to_check_against = string("viscous_files/viscous_test_face", interface, "_GtR_vec.dat")

        GtL_vec_to_check_against = readdlm(GtL_file_to_check_against)
        GtR_vec_to_check_against = readdlm(GtR_file_to_check_against)

        for ii = 1:length(GtL_vec)

          diff = GtL_vec[ii] - GtL_vec_to_check_against[ii]
          @test isapprox( diff, 0.0) atol=1e-13

          diff = GtR_vec[ii] - GtR_vec_to_check_against[ii]
          @test isapprox( diff, 0.0) atol=1e-13

        end

      end   # end loop over interfaces

    end   # End  testset("Checking Diffusion Tensor")

  end   # End testset("--- Testing Viscous Terms ---")

end # End function test_viscous

add_func1!(EulerTests, test_viscous, [TAG_VISCOUS, TAG_LONGTEST])
