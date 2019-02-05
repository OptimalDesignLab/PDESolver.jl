function test_shockcapturing_parallel()

  arg_dict = Dict{String, Any}(
  "physics" => "Euler",
  "run_type" => 5,
  "jac_type" => 3,
  "jac_method" => 2,
  "dimensions" => 3,
  "operator_type" => "SBPOmega",
  "order" => 1,
  "IC_name" => "ICRho1E2U3",
  "numBC" => 1,
  "BC1" => [0,1,2,3,4,5],
  "BC1_name" => "Rho1E2U3BC",
  "delta_t" => 0.005,
  "t_max" => 500.000,
  "smb_name" => "SRCMESHES/tet8cube_parallel.smb",
  "res_abstol" => 1e-12,
  "res_reltol" => -1.0,
  "step_tol" => 1e-9,
  "itermax" => 20,
  "solve" => false,
  # make this DG
  "Flux_name" => "RoeFlux",
  "use_DG" => true, 
  "p_free" => 2.0, #2.0,
  "T_free" => 50.0, #50,
  "addShockCapturing" => true,
  "shock_capturing_name" => "SBPParabolic",
  )


  mesh, sbp, eqn, opts = solvePDE(arg_dict)

  @testset "Shock Capturing" begin
    test_shockmesh(mesh, sbp, eqn, opts)
  end

end


function getEntireMesh(mesh, sbp, eqn::EulerData{Tsol, Tres}, opts) where {Tsol, Tres}

  degree = sbp.degree

  # construct the shock mesh with all elements in it that are fully interior
  shockmesh = EulerEquationMod.ShockedElements{Tres}(mesh)

  for i=1:mesh.numEl
    push!(shockmesh, i, 1.0)
  end

  EulerEquationMod.completeShockElements(mesh, shockmesh)

  return shockmesh
end


function test_shockmesh(mesh, sbp, eqn, opts)

  shockmesh = getEntireMesh(mesh, sbp, eqn, opts)

  @test shockmesh.numShock == mesh.numEl
  @test shockmesh.numNeighbor == 0
  @test shockmesh.npeers == mesh.npeers

  for i=1:shockmesh.npeers
    @test shockmesh.peer_indices[i] == i
    @test shockmesh.numShared[i] == (mesh.shared_element_offsets[i+1] - mesh.shared_element_offsets[i])
    @test length(shockmesh.shared_els[i]) == shockmesh.numShared[i]
    @test shockmesh.numSharedInterfaces[i] == length(mesh.shared_interfaces[i])
  end


  return nothing
end


add_func1!(EulerTests, test_shockcapturing_parallel, [TAG_SHORTTEST, TAG_TMP]) 
