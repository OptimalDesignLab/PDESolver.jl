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

function test_area(mesh, sbp, eqn, opts)
  
  #----------------------------------------------------------------------------
  # test computeNormal
#  println("testing computeNormal")
  nout = mesh.dim*mesh.numNodesPerFace*mesh.numBoundaryFaces
  nin = mesh.dim*mesh.dim*mesh.numNodesPerFace*mesh.numBoundaryFaces
  jac = zeros(nout, nin)

  dxidx = rand(mesh.dim, mesh.dim, mesh.numNodesPerFace, mesh.numBoundaryFaces)

  nrm0 = Utils.computeNormal(mesh, eqn, mesh.bndryfaces, dxidx)
  pert = 1e-6
  for i=1:nin
    dxidx[i] += pert
    nrm = Utils.computeNormal(mesh, eqn, mesh.bndryfaces, dxidx)

    for j=1:nout
      jac[j, i] = (nrm[j] - nrm0[j])/pert
    end

    dxidx[i] -= pert
  end

  # reverse mode
  jac2 = zeros(jac)
  nrm_bar = zeros(nrm0)
  dxidx_bar = zeros(dxidx)
  for i=1:nout
    nrm_bar[i] = 1
    fill!(dxidx_bar, 0.0)
    Utils.computeNormal_rev(mesh, eqn, mesh.bndryfaces, nrm_bar, dxidx, dxidx_bar)

    for j=1:nin
      jac2[i, j] = dxidx_bar[j]
    end

    nrm_bar[i] = 0
  end

  for i=1:nin
    @fact norm(jac[:, i] - jac2[:, i])/size(jac, 1) --> roughly(0.0, atol=1e-5)
  end

i


  #----------------------------------------------------------------------------
  # test computeVolumeContribution
#  println("testing computeVolumecontribution")
  vol = calcVolumeContribution!(mesh, eqn, [1])
  @fact vol --> roughly(8, atol=1e-12)


  # test computeVolumeContrib
  nout = 1
  nin = mesh.dim*mesh.numNodesPerFace*mesh.numBoundaryFaces
  jac = zeros(nout, nin)


  v0 = calcVolumeContribution!(mesh, eqn, [1])
  pert = 1e-6
  for i=1:nin
    mesh.nrm_bndry[i] += pert
    v_i = calcVolumeContribution!(mesh, eqn, [1])

    jac[1, i] = (v_i - v0)/pert

    mesh.nrm_bndry[i] -= pert
  end

  # reverse mode
  fill!(mesh.nrm_bndry_bar, 0.0)
  Utils.calcVolumeContribution_rev!(mesh, eqn, [1], 1.0)

  @fact norm(jac - reshape(mesh.nrm_bndry_bar, 1, nin))/size(jac, 2) --> roughly(0.0, atol=1e-5)

  #----------------------------------------------------------------------------
  # test calcProjectedAreaContribution
#  println("testing calcProjectedAreaContribution")

  proj_area = calcProjectedAreaContribution!(mesh, eqn, [1], 1)
  @fact proj_area --> roughly(4, atol=1e-12)

  nout = 1
  nin = mesh.dim*mesh.numNodesPerFace*mesh.numBoundaryFaces
  jac = zeros(nout, nin)

  v0 = calcProjectedAreaContribution!(mesh, eqn, [1], 1)
  pert = 1e-6
  for i=1:nin
    mesh.nrm_bndry[i] += pert
    v_i = calcProjectedAreaContribution!(mesh, eqn, [1], 1)

    jac[1, i] = (v_i - v0)/pert

    mesh.nrm_bndry[i] -= pert
  end

  # reverse mode
  fill!(mesh.nrm_bndry_bar, 0.0)
  Utils.calcProjectedAreaContribution_rev!(mesh, eqn, [1], 1, 1.0)

  diffnorm = 0.0
  cnt = 0
  idx = 1
  for i=1:mesh.numBoundaryFaces
    for j=1:mesh.numNodesPerFace
      for k=1:mesh.dim
        if abs(mesh.nrm_bndry[k, j, i]) > 1e-8  # abs is not differentiable at 0
          diffnorm += abs(jac[1, idx] - mesh.nrm_bndry_bar[k, j, i])
          if abs(jac[1, idx] - mesh.nrm_bndry_bar[k, j, i]) > 1e-5
          end

          cnt += 1
        end
        idx += 1
      end
    end
  end

  diffnorm = sqrt(diffnorm)/cnt
  @fact diffnorm --> roughly(0.0, atol=1e-5)

  # check the other method
#  println("testing other method")
#  nrm0 = Utils.computeNormal(mesh, eqn, mesh.bndryfaces,)
  nrm0 = copy(mesh.nrm_bndry)
  coords = mesh.coords_bndry

  nrm = copy(nrm0)
  jac = zeros(1, length(nrm0))
  for i=1:length(nrm)
    nrm[i] += pert
    v_i = calcProjectedAreaContribution!(mesh.sbpface, nrm, 1)
    jac[1, i] = (v_i - v0)/pert
    nrm[i] -= pert
  end

  nrm_bar = zeros(nrm)
  calcProjectedAreaContribution_rev!(mesh.sbpface, nrm, nrm_bar, 1, 1.0)

  diffnorm = norm(vec(jac) - vec(nrm_bar))/length(jac)
  diffnorm = 0.0
  cnt = 0
  for i=1:length(jac)
    diff = jac[i] - nrm_bar[i]

    if abs(nrm[i]) > 1e-8
      cnt += 1
      diffnorm += diff*diff
    end
  end

  diffnorm = sqrt(diffnorm)/length(jac)
  @fact diffnorm --> roughly(0.0, atol=1e-5)

  return nothing
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
    test_area(mesh, sbp, eqn, opts)
  end

end


add_func2!(AdvectionTests, test_utils3, test_3d_inputfile, [TAG_REVERSEMODE, TAG_SHORTTEST])
# it doesn't matter what mesh is loaded, so reuse the test_lowlevel one
add_func2!(AdvectionTests, test_utils2, test_lowlevel_inputfile, [TAG_REVERSEMODE, TAG_SHORTTEST])
