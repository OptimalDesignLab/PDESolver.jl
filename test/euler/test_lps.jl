# test local projection stabilization

function test_lps()

  opts = Dict{String, Any}(
    "physics" => "Euler",
    "operator_type" => "SBPDiagonalE",
    "dimensions" => 2,
    "run_type" => 5,
    "jac_method" => 2,
    "jac_type" => 2,
    "order" => 3,
    "IC_name" => "ICIsentropicVortex",
    "use_DG" => true,
    "volume_integral_type" => 2,
    "Volume_flux_name" => "IRFlux",
    "Flux_name" => "IRSLFFlux",
    "numBC" => 3,
    "BC1" => [0],
    "BC1_name" => "isentropicVortexBC",  # outlet
    "BC2" => [2],
    "BC2_name" => "isentropicVortexBC", # inlet
    "BC3" => [1, 3],
    "BC3_name" => "noPenetrationBC",  # was noPenetrationBC
    "aoa" => 0.0,
    "smb_name" => "SRCMESHES/vortex_3x3_.smb",
    "dmg_name" => ".null",
    "solve" => false,
    )

  mesh, sbp, eqn, opts = solvePDE(opts)

  @testset "LP Stabiization" begin
    # test the reference element
    degree = sbp.degree
    q = zeros(mesh.numNodesPerElement)
    coords = calcnodes(sbp, sbp.vtx)

    # the projection operator should not modify polynomials up to degree
    # p, and should decrease the energy of polynomials of higher degrees

    for d=0:(degree + 1)
      for i=1:mesh.numNodesPerElement
        x = coords[1, i]; y = coords[2, i]
        q[i] = x^(d) + 2*y^(d)
      end

      q2 = eqn.params.lps_data.P*q
      e1 = dot(q, sbp.w .* q)
      e2 = dot(q2, sbp.w .* q2)  # energy of higher order modes

      if d <= degree
        @test abs(e2) < 1e-15
      else
        @test e2 < 0.8*e1  # 0,8 is arbitrary factor
      end
    end   # end d

    # now test on non-reference elements
    for d=0:(degree + 1)
      for i=1:mesh.numEl
        for j=1:mesh.numNodesPerElement
          x = mesh.coords[1, j, i]; y = mesh.coords[2, j, i]
          q[j] = x^(d) + 2*y^(d)
        end

        q2 = eqn.params.lps_data.P*q
        e1 = dot(q, sbp.w .* q)
        e2 = dot(q2, sbp.w .* q2)  # energy of higher order modes

        if d <= degree
          @test abs(e2) < 1e-15
        else
          @test e2 < 0.8*e1  # 0,8 is arbitrary factor
        end
      end  # end i
    end  # end d

    # check entropy stability
    eqn.q .+= 0.01*rand(size(eqn.q))
    w = zeros(mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)
    entropy_vars = eqn.params.lps_data.entropy_vars
    for i=1:mesh.numEl
      for j=1:mesh.numNodesPerElement
        q_i = ro_sview(eqn.q, :, j, i)
        w_i = sview(w, :, j, i)
        EulerEquationMod.convertToEntropy(entropy_vars, eqn.params, q_i, w_i)
      end
    end

    fill!(eqn.res, 0)
    EulerEquationMod.applyLPStab(mesh, sbp, eqn, opts)
    val = sum(w .* eqn.res)
    println("val = ", val)
    @test val < 0
  end  # end testset

  return nothing
end


add_func1!(EulerTests, test_lps, [TAG_SHORTTEST, TAG_TMP])



