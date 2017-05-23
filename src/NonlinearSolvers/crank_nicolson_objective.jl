# import PDESolver.evalResidual
export calcdJdu_CS, calcObjectiveFn, obj_zero, calcVV, calcdRdA_CS

function calcdJdu_CS{Tmsh, Tsol}(mesh::AbstractMesh{Tmsh}, sbp::AbstractSBP, eqn::AbstractSolutionData{Tsol}, opts)

  # complex step it
  pert = complex(0, 1e-20)

  dJdu = zeros(Tsol, length(eqn.q_vec))

  for i = 1:length(eqn.q_vec)
    eqn.q_vec[i] += pert

    J_arr = calcObjectiveFn(mesh, sbp, eqn, opts)
    J = J_arr[1]
    # println("=== in dJdu_CS: typeof(J_arr): ", typeof(J_arr))
    # println("=== in dJdu_CS: typeof(J): ", typeof(J))
    # println("=== in dJdu_CS: typeof(integrand_deriv): ", typeof(integrand_deriv))
    # println("=== in dJdu_CS: typeof(pert): ", typeof(pert))
    dJdu[i] = imag(J)/norm(pert)
    eqn.q_vec[i] -= pert
  end

  return dJdu

end

# function calcObjectiveFn{Tmsh, Tsol}(mesh::AbstractMesh{Tmsh}, sbp::AbstractSBP, eqn::AbstractSolutionData{Tsol}, opts)
function calcObjectiveFn{Tmsh, Tsol}(mesh::AbstractMesh{Tmsh}, sbp::AbstractSBP, eqn::AbstractSolutionData{Tsol}, opts; isDeriv=false)

  eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
  if mesh.isDG
    boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, eqn.q, eqn.q_bndry)
  end

  # TODO: get functional edges in a non BS way
  functional_edges = 1
  nDof = 1

  local_functional_val = zeros(Tsol, nDof)
  # println("===dJdu=== size(local_functional_val): ", size(local_functional_val))
  # println("===dJdu=== size(eqn.q_vec): ", size(eqn.q_vec))
  # eqn.q_vec: (96,)
  # local_functional_val: (96,)


  for itr = 1:length(functional_edges)
    g_edge_number = functional_edges[itr]
    itr2 = 0
    for itr2 = 1:mesh.numBC
      if findfirst(mesh.bndry_geo_nums[itr2], g_edge_number) > 0
        break
      end
    end

    start_index = mesh.bndry_offsets[itr2]
    end_index = mesh.bndry_offsets[itr2+1]
    idx_range = start_index:(end_index-1)
    bndry_facenums = sview(mesh.bndryfaces, idx_range)

    nfaces = length(bndry_facenums)

    integrand = zeros(Tsol, 1, mesh.sbpface.numnodes, nfaces)
    # println("===dJdu=== size(integrand): ", size(integrand))
    # println("===dJdu=== nfaces: ", nfaces)
    # println("===dJdu=== mesh.sbpface.numnodes: ", mesh.sbpface.numnodes)

    for i = 1:nfaces
      bndry_i = bndry_facenums[i]
      global_facenum = idx_range[i]

      for j = 1:mesh.sbpface.numnodes
        #q = sview(eqn.q_bndry, :, j, global_facenum)
        q = eqn.q_bndry[:, j, global_facenum]
        # println("====== type of q: ", typeof(q))
        # println("====== size of q: ", size(q))
        # println("====== q: ", q)
        # convertToConservative(eqn.params, q, q2)

        # replaces calcBoundaryFunctionalIntegrand
        # integrand = zeros(Tsol, ndof, mesh.sbpface.numnodes, nfaces)    # dims?
        # integrand[1, j, i] = q.^2
        # TODO: figure out why [1]'s are required
        if isDeriv == false                 # calculates J = int(u^2)
          integrand[1, j, i] = q[1]*q[1]
        else                                # calculates dJdu = deriv(int(u^2)) = 2*u
          integrand[1, j, i] = 2*q[1]
        end

        # TODO: how to get the analytical derivative outside of integral


      end   # end of loop: j = 1:mesh.sbpfacenumnodes


      val_per_geom_edge = zeros(Tsol, 1)
      # println("===dJdu=== size(val_per_geom_edge): ", size(val_per_geom_edge))

      # use integratefunctional, not boundaryintegrate: why?
      integratefunctional!(mesh.sbpface, mesh.bndryfaces[idx_range], integrand, val_per_geom_edge)
#       boundaryintegrate!(mesh.sbpface, mesh.bndryfaces[idx_range], integrand, val_per_geom_edge)

#       println(" size of val_per_geom_edge: ", size(val_per_geom_edge))

      local_functional_val[:] += val_per_geom_edge[:]

      # TODO:
      # serial: print out local_functional_val, compare with analytical
      # parallel: mpi all reduce then do the same

    end   # end of loop: i = 1:nfaces

  end   # end of loop: itr = 1:length(functional_edges)

  return local_functional_val

end

"""
calculates v, which is dudA, for the advection adjoint test
"""
function calcVV(mesh, sbp, eqn, opts, t)

  # x = 2*pi
  x = 2.0       # make sure the mesh being run is u_adj_x0-2_y0-2_4x4_tri0.smb
  omega = 1.0

  v = sin(-x + omega*t)
  return v

end

"""
calculates dRdA for the advection adjoint test, using complex step
"""
function calcdRdA_CS(mesh, sbp, eqn, opts, t)

  println(" calculating dRdA using CS")

  # complex step perturbation
  pert = complex(0, 1e-20)

  eqn.params.sin_amplitude += pert
  # params.sin_amplitude += pert

  eqn_temp = eqn_deepcopy(eqn, mesh, sbp, opts)
  # eqn_temp.q = reshape(eqn_temp.q_vec, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)
  # eqn_temp.res = reshape(eqn_temp.res_vec, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)


  # J_arr = calcObjectiveFn(mesh, sbp, eqn, opts)
  # J = J_arr[1]
  # R = calcResidual
  evalResidual(mesh, sbp, eqn_temp, opts)
  assembleSolution(mesh, sbp, eqn_temp, opts, eqn_temp.res, eqn_temp.res_vec)

  dRdA = imag(eqn_temp.res_vec)/norm(pert)
  # println("dRdA: ", dRdA)

  eqn.params.sin_amplitude -= pert
  # params.sin_amplitude -= pert

  return dRdA
end

"""
calculates dRdA for the advection adjoint test, using finite difference
"""
function calcdRdA_FD(mesh, sbp, eqn, opts, t)

  println(" calculating dRdA using FD")

  # finite difference perturbation
  pert = 1e-6

  eqn.params.sin_amplitude += pert
  # params.sin_amplitude += pert

  eqn_temp = eqn_deepcopy(eqn, mesh, sbp, opts)
  # eqn_temp.q = reshape(eqn_temp.q_vec, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)
  # eqn_temp.res = reshape(eqn_temp.res_vec, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)

  # J_arr = calcObjectiveFn(mesh, sbp, eqn, opts)
  # J = J_arr[1]
  # R = calcResidual
  evalResidual(mesh, sbp, eqn_temp, opts)
  assembleSolution(mesh, sbp, eqn_temp, opts, eqn_temp.res, eqn_temp.res_vec)

  # this is required! without it, eqn.res_vec apparently has stale data.
  #   without this => difference between CS & FD is ~1e10
  evalResidual(mesh, sbp, eqn, opts)
  assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)

  # Ready to finish FD approx.
  #   f(x+h) is eqn_temp.res_vec
  #   f(x) is eqn.res_vec

  dRdA = (eqn_temp.res_vec - eqn.res_vec)/norm(pert)
  #println("dRdA: ", dRdA)

  eqn.params.sin_amplitude -= pert
  # params.sin_amplitude -= pert

  return dRdA
end

"""
  obj_zero

  Inputs: none
  Outputs: 0

  Zero-valued objective function
"""
function obj_zero()
  return 0.0
end
