# test interpolation of fields from one SBP operator to another

function test_interp()

  @testset "----- Testing Field Interpolation -----" begin
    fname = "input_vals_channel_dg_large.jl"

    opts = read_input(fname)
    opts2 = copy(opts)
    opts2["order"] += 1

    degree = opts["order"]

    # solve problem on both grids

    mesh, sbp, eqn, opts = createObjects(opts)
    mesh2, sbp2, eqn2, opts2 = createObjects(opts2)
    #array1DTo3D(mesh2, sbp2, eqn2, opts2, eqn2.q_vec, eqn2.q)
    q_vec = copy(eqn.q_vec)
    q_vec2 = copy(eqn2.q_vec)

    # put degree p field into eqn
    for i=1:mesh.numEl
      for j=1:mesh.numNodesPerElement
        x = mesh.coords[1, j, i]
        y = mesh.coords[2, j, i]

        for k=1:mesh.numDofPerNode
          val = x^degree + 2*y^degree + k
          eqn.q[k, j, i] = val
          q_vec[mesh.dofs[k, j, i]] = val
        end
      end
    end

    interpField(sbp, eqn.q, sbp2, eqn2.q)
    interpField(mesh, sbp, q_vec, mesh2, sbp2, q_vec2)

    for i=1:mesh2.numEl
      for j=1:mesh2.numNodesPerElement
        x = mesh2.coords[1, j, i]
        y = mesh2.coords[2, j, i]

        for k=1:mesh2.numDofPerNode
          val = x^degree + 2*y^degree + k


          @test isapprox( abs(eqn2.q[k, j, i] - val), 0.0) atol=1e-13
          @test isapprox( abs(q_vec2[mesh2.dofs[k, j, i]] - val), 0.0) atol=1e-13
        end
      end
    end




  end  # end facts block

  return nothing
end


add_func1!(EulerTests, test_interp, [TAG_SHORTTEST, TAG_TMP])
