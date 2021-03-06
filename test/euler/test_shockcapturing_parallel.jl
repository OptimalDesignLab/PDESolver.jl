function test_shockcapturing_parallel()

  opts = Dict{String, Any}(
    "physics" => "Euler",
    "operator_type" => "SBPOmega",
    "dimensions" => 2,
    "run_type" => 5,
    "jac_method" => 2,
    "jac_type" => 3,
    "order" => 2,
    "IC_name" => "ICIsentropicVortex",
    "use_DG" => true,
    "volume_integral_type" => 2,
    "Volume_flux_name" => "IRFlux",
    "face_integral_type" => 2,
    "FaceElementIntegral_name" => "ESLFFaceIntegral",
    "Flux_name" => "IRFlux",
    "numBC" => 3,
    "BC1" => [0],
    "BC1_name" => "isentropicVortexBC",  # outlet
    "BC2" => [2],
    "BC2_name" => "isentropicVortexBC", # inlet
    "BC3" => [1, 3],
    "BC3_name" => "isentropicVortexBC",
    "aoa" => 0.0,
    "smb_name" => "SRCMESHES/vortex_3x3_p_.smb",
    "dmg_name" => ".null",
    "itermax" => 20,
    "res_abstol" => 1e-12,
    "res_reltol" => 1e-12,
    "do_postproc" => true,
    "exact_soln_func" => "ICIsentropicVortex",
    "solve" => false,
    )

  opts2 = copy(opts)
  opts2["operator_type"] = "SBPDiagonalE"
  opts2["face_integral_type"] = 1
  opts2["flux_name"] = "IRSLFFlux"
  opts2["use_lps"] = true



  mesh, sbp, eqn, opts = solvePDE(opts)
  mesh2, sbp2, eqn2, opts2 = solvePDE(opts2)

  @testset "Shock Capturing" begin
    test_shockmesh(mesh, sbp, eqn, opts)
    test_br2_parallelpart(mesh, sbp, eqn, opts)
    test_br2reduced_parallelpart(mesh2, sbp2, eqn2, opts2)
  end

end


function getEntireMesh(mesh, sbp, eqn::EulerData{Tsol, Tres}, opts) where {Tsol, Tres}

  degree = sbp.degree

  # construct the shock mesh with all elements in it that are fully interior
  shockmesh = EulerEquationMod.ShockedElements{Tres}(mesh)

  for i=1:mesh.numEl
    push!(shockmesh, i)
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


function test_br2_parallelpart(mesh, sbp, eqn::EulerData{Tsol, Tres},
                               _opts) where {Tsol, Tres}
# test consistency with the serial version

  opts = copy(_opts)
  sensor = EulerEquationMod.ShockSensorEverywhere{Tsol, Tres}(mesh, sbp, opts)
  capture = EulerEquationMod.SBPParabolicSC{Tsol, Tres}(mesh, sbp, eqn, opts, sensor)

  # solve the PDE to get a solution with non-zero jump between elements
  # that can be reproduced in parallel
  opts["solve"] = true
  solvePDE(mesh, sbp, eqn, opts)

  # the scheme is entropy stable when the dirichlet BC = 0
  for i=1:opts["numBC"]
    opts[string("BC", i, "_name")] = "zeroBC"
  end
  EulerEquationMod.getBCFunctors(mesh, sbp, eqn, opts)

  test_sc_parallelpart(mesh, sbp, eqn, opts, capture, "br2")

  return nothing
end

function test_br2reduced_parallelpart(mesh, sbp, eqn::EulerData{Tsol, Tres},
                                      _opts) where {Tsol, Tres}
# test consistency with the serial version

  opts = copy(_opts)
  sensor = EulerEquationMod.ShockSensorEverywhere{Tsol, Tres}(mesh, sbp, opts)
  capture = EulerEquationMod.SBPParabolicReducedSC{Tsol, Tres}(mesh, sbp, eqn, opts, sensor)

  # solve the PDE to get a solution with non-zero jump between elements
  # that can be reproduced in parallel
  opts["solve"] = true
  solvePDE(mesh, sbp, eqn, opts)

  test_sc_parallelpart(mesh, sbp, eqn, opts, capture, "br2reduced")

  return nothing
end




function test_sc_parallelpart(mesh, sbp, eqn::EulerData{Tsol, Tres},
                              opts, capture, prefix::String) where {Tsol, Tres}
  # need to communicate the final solution
  startSolutionExchange(mesh, sbp, eqn, opts, wait=true)

  w_vec = zeros(Tsol, mesh.numDof)
  copy!(w_vec, eqn.q_vec)
  EulerEquationMod.convertToIR(mesh, sbp, eqn, opts, w_vec)

  fill!(eqn.res, 0)
  EulerEquationMod.applyShockCapturing(mesh, sbp, eqn, opts, capture)

  val = dot(w_vec, eqn.res_vec)
  val = MPI.Allreduce(val, MPI.SUM, eqn.comm)

  println("val = ", val)
  @test val < 0

  # save this for parallel tests
  val2 = readdlm("$(prefix)_entropy_serial.dat")
  println("val2 = ", val2[1])
  println("diff = ", val - val2[1])

  @test abs(val - val2[1]) < 1e-8  # test the entropy dissipation is the same,
                                   # to the tolerance of the nonlinear solve

  saveSolutionToMesh(mesh, eqn.q_vec)
  writeVisFiles(mesh, "$(prefix)_parallel")

  return nothing
end


add_func1!(EulerTests, test_shockcapturing_parallel, [TAG_SHORTTEST, TAG_TMP]) 
