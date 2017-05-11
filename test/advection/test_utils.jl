# test for the Utils module

"""
  Test calcBCNormal reverse mode.  This function works for any mesh dimension
"""
function test_calcBCNormal(mesh, eqn)
  dim = mesh.dim
  dxidx = rand(dim, dim)
  nrm = rand(dim)

  dxidx_c = zeros(Complex128, size(dxidx))
  copy!(dxidx_c, dxidx)
  nrm_out = zeros(Complex128, size(dxidx, 1))
  h = 1e-20
  pert = Complex128(0, h)

  jac = zeros(dim, dim*dim)
  for i=1:length(dxidx)
    dxidx_c[i] += pert
    calcBCNormal(eqn.params, dxidx_c, nrm, nrm_out)
    for j=1:length(nrm_out)
      jac[j, i] = imag(nrm_out[j])/h
    end
    dxidx_c[i] -= pert
  end

  nrm2_bar = zeros(dim)
  dxidx_bar = zeros(dim, dim)
  jac2 = zeros(jac)
  # compute reverse mode
  for i=1:length(nrm2_bar)
    nrm2_bar[i] = 1
    calcBCNormal_revm(eqn.params, dxidx, nrm, nrm2_bar, dxidx_bar)
    for j=1:length(dxidx_bar)
      jac2[i, j] = dxidx_bar[j]
    end

    fill!(dxidx_bar, 0.0)
    nrm2_bar[i] = 0
  end

  @fact norm(jac - jac2) --> roughly(0.0, atol=1e-13)
end


function test_utils2(mesh, sbp, eqn, opts)

  facts("----- Testing Utils 2D -----") do

    test_calcBCNormal(mesh, eqn)

    #TODO: test other things too
  end

end


function test_utils3(mesh, sbp, eqn, opts)

  facts("----- Testing Utils 3D -----") do

    test_calcBCNormal(mesh, eqn)
  end

end


add_func2!(AdvectionTests, test_utils3, test_3d_inputfile, [TAG_REVERSEMODE])
# it doesn't matter what mesh is loaded, so reuse the test_lowlevel one
add_func2!(AdvectionTests, test_utils2, test_lowlevel_inputfile, [TAG_REVERSEMODE])
