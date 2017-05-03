# tests for curvilinear meshes

function test_curvilinear(mesh, sbp, eqn, opts)

  @assert mesh.sbpface.stencilsize == mesh.numNodesPerElement
  Tdim = mesh.dim
  sbpface = mesh.sbpface

  stencilsize = mesh.sbpface.stencilsize
#  ncontrib = mesh.numTypePerElement[end-1]

  # collect E for all the elements
  E_all = zeros(stencilsize, stencilsize, Tdim, mesh.numEl)
  contrib_count = zeros(Int, mesh.numEl)

  PL = zeros(stencilsize, stencilsize)
  PR = zeros(PL)

  Pnbr = zeros(mesh.numNodesPerFace, mesh.numNodesPerFace)

  println("checking 1^T E 1 = 0")
  # this is a periodic mesh, so we only need to look at interfaces
  for i=1:mesh.numInterfaces
    iface = mesh.interfaces[i]
    contrib_count[iface.elementL] += 1
    contrib_count[iface.elementR] += 1

    # get R_complete, B
    fill!(PL, 0.0)
    fill!(PR, 0.0)
    permMatrix!(sbpface.perm[:, iface.faceL], PL)
    permMatrix!(sbpface.perm[:, iface.faceR], PR)

    RL_complete = sbpface.interp.'*PL
    RR_complete = sbpface.interp.'*PR
    B = diagm(sbpface.wface)

    for dim = 1:Tdim  # loop over parametric directions

      # get N
      NL_dim = diagm(vec(mesh.nrm_face[dim, :, i]))

      # the node order is different for elementR
      permMatrix!(sbpface.nbrperm[:, iface.orient], Pnbr)
      NR_dim = -diagm(Pnbr*vec(mesh.nrm_face[dim, :, i]))

      # compute E in current direction and accumulate into E_all
      E_all[:, :, dim, iface.elementL] += RL_complete.'*NL_dim*B*RL_complete
      E_all[:, :, dim, iface.elementR] += RR_complete.'*NR_dim*B*RR_complete
    end
  end


  Ivec = ones(stencilsize)
  for i=1:mesh.numEl
    @fact contrib_count[i] --> 3
    for dim=1:Tdim
      E = E_all[:, :, dim, i]
      val = Ivec.'*E*Ivec

      @fact val[1] --> roughly(0.0, atol=1e-12)
      @fact norm(E - E.') --> roughly(0.0, atol=1e-12)
    end
  end

  # test S
  println("checking S and D")
  S = zeros(stencilsize, stencilsize, Tdim)

  for i=1:mesh.numEl
    calcSCurvilinear(sbp, mesh.dxidx[:, :, :, i], S)
    Qxi = sbp.Q[:, :, 1]
    Qeta = sbp.Q[:, :, 2]
    Hinv = inv(diagm(sbp.w./mesh.jac[:, i]))

    for dim=1:Tdim
      @fact norm(S[:, :, dim] + S[:, :, dim].') --> roughly(0.0, atol=1e-12)
      D = Hinv*(S[:, :, dim] + 0.5*E_all[:, :, dim, i])
      D1 = D*Ivec
      @fact norm(D1) --> roughly(0.0, atol=1e-12)
#=
      metrics_xi_x = diagm(vec(mesh.dxidx[1, dim, :, i]))
      metrics_eta_x = diagm(vec(mesh.dxidx[2, dim, :, i]))
      A_xi = metrics_xi_x*Qxi
      A_eta = metrics_eta_x*Qeta

      S2 = 0.5*A_xi - 0.5*A_xi.' + 0.5*A_eta - 0.5*A_eta.'
      println("S2 = \n", S2)
      D = Hinv*(S2 + 0.5*E_all[:, :, dim, i])
      D2 = D*Ivec
      @fact norm(D2) --> roughly(0.0, atol=1e-12)
=#
    end


  end

  println("checking theorem 6")
  D = zeros(stencilsize, stencilsize, Tdim)
  E = zeros(stencilsize, stencilsize, Tdim)
  for i=1:mesh.numEl
    Hinv = diagm(1./sbp.w)

    calcECurvilinear(sbp, mesh.dxidx[:, :, :, i], E)
    calcDCurvilinear(sbp, mesh.dxidx[:, :, :, i], D)
    for dim=1:Tdim
      lhs = D[:, :, dim]*Ivec
      rhs = Hinv*(E[:, :, dim] - E_all[:, :, dim, i])*Ivec
      @fact norm(lhs - rhs) --> roughly(0.0, atol=1e-12)
    end
  end
#=
  val = 0.0
  @time for i=1:mesh.numEl
    calcSCurvilinear(sbp, unsafe_view(mesh.dxidx, :, :, :, i), S)

    for dim=1:Tdim
      for i=1:stencilsize
        for j=1:stencilsize
          val += S[j, i, dim]
        end
      end
    end

  end
=#


  return nothing
end




infile = "input_vals_curve.jl"
add_func2!(EulerTests, test_curvilinear, infile, [TAG_CURVILINEAR, TAG_ESS])
