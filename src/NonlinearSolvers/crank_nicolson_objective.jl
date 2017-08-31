export calcdJdu_CS, calcObjectiveFn, obj_zero, calcVV, calcdRdA_CS

"""
  calculates dJdu using the complex step method
"""
function calcdJdu_CS{Tmsh, Tsol}(mesh::AbstractMesh{Tmsh}, sbp::AbstractSBP, eqn::AbstractSolutionData{Tsol}, opts, h, t)

  # complex step it
  pert = complex(0, 1e-20)

  dJdu = zeros(Tsol, length(eqn.q_vec))

  for i = 1:length(eqn.q_vec)
    eqn.q_vec[i] += pert

    evalResidual(mesh, sbp, eqn, opts)
    assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)

    J_arr = calcObjectiveFn(mesh, sbp, eqn, opts, h, t)
    J = J_arr[1]
    # println("=== in dJdu_CS: typeof(J_arr): ", typeof(J_arr))
    # println("=== in dJdu_CS: typeof(J): ", typeof(J))
    # println("=== in dJdu_CS: typeof(integrand_deriv): ", typeof(integrand_deriv))
    # println("=== in dJdu_CS: typeof(pert): ", typeof(pert))
    dJdu[i] = imag(J)/norm(pert)
    eqn.q_vec[i] -= pert
  end

  evalResidual(mesh, sbp, eqn, opts)
  assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)

  return dJdu

end

"""
  calculates dJdu using the finite difference method
"""
function calcdJdu_FD{Tmsh, Tsol}(mesh::AbstractMesh{Tmsh}, sbp::AbstractSBP, eqn::AbstractSolutionData{Tsol}, opts, h, t)

  pert = 1e-6

  dJdu = zeros(Tsol, length(eqn.q_vec))

  # unperturbed
  J_unpert_arr = calcObjectiveFn(mesh, sbp, eqn, opts, h, t)
  J_unpert = J_unpert_arr[1]

  for i = 1:length(eqn.q_vec)

    eqn.q_vec[i] += pert

    evalResidual(mesh, sbp, eqn, opts)
    assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)

    J_arr = calcObjectiveFn(mesh, sbp, eqn, opts, h, t)
    J = J_arr[1]

    dJdu[i] = (J - J_unpert)/norm(pert)
    eqn.q_vec[i] -= pert

  end

  evalResidual(mesh, sbp, eqn, opts)
  assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)

  return dJdu

end

""" 
  calculates the time integral of the objective function, 
    storing it in a file each time as well as returning it
"""
function calcTimeIntObjFn(mesh, sbp, eqn, opts, h, t; isDeriv=false, isFirstOrLastStep=false)

  J_total = readdlm(filename)

  if isFirstOrLastStep
    fac = 1.0
  else
    fac = 2.0
  end

  J_i = calcObjectiveFn(mesh, sbp, eqn, opts, h, t; isDeriv)
  J_total += (h*0.5)*fac*J_i

  filename = "J_total.dat"
  writedlm(filename, J_total)

  return J_total

end

""" 
  calculates the objective function

  J = int(u^2)

  or optionally (if isDeriv is passed as true), calculates dJdu = int(2u)
"""
function calcObjectiveFn{Tmsh, Tsol}(mesh::AbstractMesh{Tmsh}, sbp::AbstractSBP, eqn::AbstractSolutionData{Tsol}, opts, h, t; isDeriv=false)

  eqn.disassembleSolution(mesh, sbp, eqn, opts, eqn.q, eqn.q_vec)
  if mesh.isDG
    boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, eqn.q, eqn.q_bndry)
  end

  # TODO: get functional edges in a non BS way
  # 1 corresponds to right side of square domain
  #   Numbering is CCW, with 0 on bottom
  functional_edges = 1
  nDof = 1

  if isDeriv == false
    J = zeros(Tsol, 1)
  else
    dJdu = zeros(Tsol, mesh.numNodesPerElement, mesh.numEl)
  end

  # mesh.bndry_geo_nums:
  #   [1]
  #   [0, 2, 3]

  for itr = 1:length(functional_edges)                  # 1:1
    g_edge_number = functional_edges[itr]               # 1
    itr2 = 0
    for itr2 = 1:mesh.numBC                             # 1:2
      if findfirst(mesh.bndry_geo_nums[itr2], g_edge_number) > 0
        break
      end
    end

    # itr2: 1
    # mesh.bndry_offsets:   [1; 3; 9]
    # mesh.bndryfaces:
    #   Boundary element, face = 7, 1
    #   Boundary element, face = 8, 1
    #   Boundary element, face = 1, 3
    #   Boundary element, face = 1, 1
    #   Boundary element, face = 2, 3
    #   Boundary element, face = 4, 1
    #   Boundary element, face = 6, 2
    #   Boundary element, face = 8, 2

    start_index = mesh.bndry_offsets[itr2]        # start_index = 1
    end_index = mesh.bndry_offsets[itr2+1]        # end_index = 3
    idx_range = start_index:(end_index-1)         # idx_range = 1:2
    bndry_facenums = sview(mesh.bndryfaces, idx_range)

    # bndry_facenums:
    #   Boundary element, face = 7, 1
    #   Boundary element, face = 8, 1

    nfaces = length(bndry_facenums)               # nfaces = 2

    # for scalar field version
    integrand = zeros(Tsol, mesh.sbpface.numnodes, nfaces)     # size: (2, 2)

    if isDeriv == false                 # calculates J = int(u^2)

      for i = 1:nfaces                # 1:2
        bndry_i = bndry_facenums[i]   
        global_facenum = idx_range[i]   # i of 1:2

        for j = 1:mesh.sbpface.numnodes       # 1:2
          q = eqn.q_bndry[:, j, global_facenum]

          integrand[j, i] = q[1]*q[1]

          # 20170622: need to scale integrand!
          # TODO: double and triple check jac_bndry index order here
          # TODO: get dimensionality in a non BS way
          dimensionality = 2
          scaling_factor = 1/(mesh.jac_bndry[j, i]^*(1/dimensionality))
          integrand[j, i] = integrand[j, i]*scaling_factor

        end   # end of loop: j = 1:mesh.sbpfacenumnodes

        # use integratefunctional, not boundaryintegrate: why?

      end   # end of loop: i = 1:nfaces
      J = integratefunctional!(mesh.sbpface, mesh.bndryfaces[idx_range], integrand)

    else                                # calculates dJdu = deriv(int(u^2)) = 2*u

      # for u_ix = 1:length(eqn.q_vec)
        for i = 1:nfaces                # 1:2
          bndry_i = bndry_facenums[i]   
          global_facenum = idx_range[i]   # i of 1:2

          for j = 1:mesh.sbpface.numnodes       # 1:2
            q = eqn.q_bndry[:, j, global_facenum]

            integrand[j, i] = 2*q[1]

            # 20170622: need to scale integrand!
            # TODO: double and triple check jac_bndry index order here
            # TODO: get dimensionality in a non BS way
            dimensionality = 2
            scaling_factor = 1/(mesh.jac_bndry[j, i]^*(1/dimensionality))
            integrand[j, i] = integrand[j, i]*scaling_factor

            # println(" ``````````````````````` j: $j  i: $i  q: ", q[1], "  integrand[j,i]: ", integrand[j,i])

          end   # end of loop: j = 1:mesh.sbpfacenumnodes

          # dJdu[u_ix] = mesh.sbpface.wface[j] * integrand[j,i]
          # use integratefunctional, not boundaryintegrate: why?
        end   # end of loop: i = 1:nfaces
      # end   # end of loop: u_ix = 1:length(eqn.q_vec)
        
        # dJdu = integratefunctional!(mesh.sbpface, mesh.bndryfaces[idx_range], integrand)
        # dJdu = integratefunctional!(mesh.sbpface, mesh.bndryfaces[idx_range], integrand)

        boundaryintegrate!(mesh.sbpface, mesh.bndryfaces[idx_range], integrand, dJdu)

    end   # end of isDeriv check around loops

  end   # end of loop: itr = 1:length(functional_edges)

  # 20170830 NOTE NOTE: this section is for setting the objective function to a constant (and dJdu to 0)
  #   the result of this is that the adjoint should be zero.
  #   Verified that the adjoint is zero when this is done. I'm leaving it commented out here for posterity,
  #   until the full adjoint works.
  #=
  if isDeriv == false
    println("type of J: ", typeof(J))
    println("size of J: ", size(J))
    J = 4.0
    return J
  else
    dJdu = zeros(dJdu)
    return dJdu
  end
  =#

  if isDeriv == false
    return J
  else
    return dJdu
  end


end

"""
calculates v, which is dudA, for the advection adjoint test
"""
function calcVV(mesh, sbp, eqn, opts, t)

  v = zeros(eqn.q_vec)

  for dof_ix = 1:mesh.numDof
    # st: subscript_tuple
    st = ind2sub((mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl), dof_ix)
    dofnum = st[1]    # dof number in this node
    nodenum = st[2]   # node number in this element
    elnum = st[3]     # element number

    coords_this_dof = getindex(mesh.coords, :, nodenum, elnum)

    x = coords_this_dof[1]

    omega = eqn.params.omega
    A = eqn.params.sin_amplitude

    # v is dudA. since u = A*sin(-x + omega*t), dudA = sin(-x + omega*t)
    v[dof_ix] = sin(-x + omega*t)

  end

  return v

end

"""
calculates dRdA for the advection adjoint test, using complex step
"""
function calcdRdA_CS(mesh, sbp, eqn, opts, t)

  println(" calculating dRdA using CS")

  # complex step perturbation
  pert = complex(0, 1e-20)

  # println(" <<< dRdA_CS: eqn.params.sin_amplitude: ", eqn.params.sin_amplitude)

  eqn_temp = eqn_deepcopy(eqn, mesh, sbp, opts)
  eqn_temp.params.sin_amplitude += pert
  # eqn_temp.params.omega += pert
  # println(" <<< dRdA_CS: eqn_temp.params.sin_amplitude: ", eqn_temp.params.sin_amplitude)

  evalResidual(mesh, sbp, eqn_temp, opts)
  assembleSolution(mesh, sbp, eqn_temp, opts, eqn_temp.res, eqn_temp.res_vec)

  dRdA = imag(eqn_temp.res_vec)/norm(pert)
  # println("dRdA: ", dRdA)

  eqn_temp.params.sin_amplitude -= pert
  # eqn_temp.params.omega -= pert

  return dRdA
end

"""
calculates dRdA for the advection adjoint test, using finite difference
"""
function calcdRdA_FD(mesh, sbp, eqn, opts, t)

  println(" calculating dRdA using FD")

  # finite difference perturbation
  pert = 1e-6

  # println(" <<< dRdA_FD: eqn.params.sin_amplitude: ", eqn.params.sin_amplitude)

  eqn_temp = eqn_deepcopy(eqn, mesh, sbp, opts)
  eqn_temp.params.sin_amplitude += pert
  # eqn_temp.params.omega += pert
  # println(" <<< dRdA_FD: eqn_temp.params.sin_amplitude: ", eqn_temp.params.sin_amplitude)

  println(" <<< dRdA_FD: using forward difference")

  evalResidual(mesh, sbp, eqn_temp, opts)
  assembleSolution(mesh, sbp, eqn_temp, opts, eqn_temp.res, eqn_temp.res_vec)

  eqn_temp.params.sin_amplitude -= pert
  # eqn_temp.params.omega -= pert

  # this is required! without it, eqn.res_vec apparently has stale data.
  #   without this => difference between CS & FD is ~1e10
  #   TODO: why is this required
  evalResidual(mesh, sbp, eqn, opts)
  assembleSolution(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)

  println(" <<< dRdA_FD: norm(eqn_temp.res_vec): ", norm(eqn_temp.res_vec))
  println(" <<< dRdA_FD: norm(eqn.res_vec): ", norm(eqn.res_vec))

  # Ready to finish FD approx.
  #   f(x+h) is eqn_temp.res_vec
  #   f(x) is eqn.res_vec

  dRdA = (eqn_temp.res_vec - eqn.res_vec)/norm(pert)
  #println("dRdA: ", dRdA)

  #=
  #----------------
  # central difference
  eqn_temp_fwd = eqn_deepcopy(eqn, mesh, sbp, opts)
  eqn_temp_fwd.params.sin_amplitude += pert

  eqn_temp_bck = eqn_deepcopy(eqn, mesh, sbp, opts)
  eqn_temp_bck.params.sin_amplitude -= pert

  println(" <<< dRdA_FD: using central difference")

  evalResidual(mesh, sbp, eqn_temp_fwd, opts)
  assembleSolution(mesh, sbp, eqn_temp_fwd, opts, eqn_temp_fwd.res, eqn_temp_fwd.res_vec)
  evalResidual(mesh, sbp, eqn_temp_bck, opts)
  assembleSolution(mesh, sbp, eqn_temp_bck, opts, eqn_temp_bck.res, eqn_temp_bck.res_vec)

  eqn_temp_fwd.params.sin_amplitude -= pert
  eqn_temp_bck.params.sin_amplitude += pert

  dRdA = (eqn_temp_fwd.res_vec - eqn_temp_bck.res_vec)/(2*norm(pert))
  =#

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
