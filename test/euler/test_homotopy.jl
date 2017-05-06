# tests for homotopy.jl

const test_homotopy_inputfile = "input_vals_channel.jl"
const test_homotopy_moddict = Dict{ASCIIString, Any}(
  "Flux_name" => "RoeFlux", 
  "use_DG" => true, 
  "IC_name" => "ICFreeStream",
  "BC1_name" => "FreeStreamBC",
  "new_fname" => "input_vals_channel_dg")


function test_homotopy(mesh, sbp, eqn, opts)

  # the initial condition is uniform flow, so the residual of the homotopy
  # should be zero

  res = zeros(eqn.res)
  fill!(eqn.res, 42)  # make sure eqn.res is not modified

  EulerEquationMod.calcHomotopyDiss(mesh, sbp, eqn, opts, res)

  for i=1:mesh.numEl
    @fact norm(res[:, :, i]) --> roughly(0.0, atol=1e-12)
  end

  for i=1:length(eqn.res)
    @fact eqn.res[i] --> 42
  end

  return nothing
end

add_func3!(EulerTests, test_homotopy, test_homotopy_inputfile, test_homotopy_moddict, [TAG_HOMOTOPY])
