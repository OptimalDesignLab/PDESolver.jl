# test interpolation of fields from one SBP operator to another

function test_interp()

  facts("----- Testing Field Interpolation -----") do
    fname = "input_vals_channel_dg_large.jl"

    opts = read_input(fname)
    opts2 = copy(opts)
    opts2["order"] += 1

    degree = opts["order"]

    # solve problem on both grids

    mesh, sbp, eqn, opts = createObjects(opts)
    mesh2, sbp2, eqn2, opts2 = createObjects(opts2)
    disassembleSolution(mesh2, sbp2, eqn2, opts2, eqn2.q, eqn2.q_vec)

    # put degree p field into eqn
    for i=1:mesh.numEl
      for j=1:mesh.numNodesPerElement
        x = mesh.coords[1, j, i]
        y = mesh.coords[2, j, i]

        for k=1:mesh.numDofPerNode
          eqn.q[k, j, i] = x^degree + 2*y^degree + k
        end
      end
    end

    interpField(sbp, eqn.q, sbp2, eqn2.q)

    for i=1:mesh2.numEl
      for j=1:mesh2.numNodesPerElement
        x = mesh2.coords[1, j, i]
        y = mesh2.coords[2, j, i]

        for k=1:mesh2.numDofPerNode
          val = x^degree + 2*y^degree + k


          @fact abs(eqn2.q[k, j, i] - val) --> roughly(0.0, atol=1e-13)
        end
      end
    end

  end  # end facts block

  return nothing
end


add_func1!(EulerTests, test_interp, [TAG_SHORTTEST, TAG_TMP])
