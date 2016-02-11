

function applyGLS2{Tmsh, Tsol, Tres}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator, eqn::AdvectionData{Tsol, Tres, 2}, opts, src_func)

#  println("----- entered applyGLS2 -----")
  Dx = zeros(Float64, mesh.numNodesPerElement, mesh.numNodesPerElement)
  Dy = zeros(Dx)

  # working vectors that do the reduction
  red_vec1 = zeros(Tsol, mesh.numNodesPerElement)
  red_vec2 = zeros(Tsol, mesh.numNodesPerElement)
  red_vec3 = zeros(red_vec1)
  red_vec4 = zeros(red_vec2)
  u = zeros(Tsol, mesh.numNodesPerElement)
  dxidx = zeros(Tmsh, 2,2, mesh.numNodesPerElement)
  tau_vec = zeros(Tres, mesh.numNodesPerElement)
  # calculate Dxi and Deta
  Dxi = diagm(1./sbp.w)*sbp.Q[:, :, 1]
  Deta = diagm(1./sbp.w)*sbp.Q[:, :, 2]

  #DEBUGGING
  gls_res = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)
  gls_full = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)
  middle_term = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)
  res_before = copy(eqn.res)
  weighting_res = zeros(middle_term)
  weighting_vec = ones(mesh.numNodesPerElement)
  tau_sum = zero(Tres)
  tau_cnt = 0


  for i=1:mesh.numEl
    # calculate the quantities needed for this element
    dxidx_hat = view(mesh.dxidx[:, :, :, i])  
    jac = view(mesh.jac, :, i)
    res = view(eqn.res, :, :, i)
    # constant coefficient advection only!
    alpha_x = eqn.alpha_x[1, 1, i]
    alpha_y = eqn.alpha_y[1, 1, i]
    # calclate the true dxidx (not scaled by the mapping jacobian)
    for j=1:mesh.numNodesPerElement
      for k=1:2
        for p=1:2
          dxidx[p, k, j] = dxidx_hat[p, k, j]*jac[j]
        end
      end
    end

    for j=1:mesh.numNodesPerElement
      for k=1:mesh.numNodesPerElement
        Dx[k, j] = dxidx[1,1]*Dxi[k, j] + dxidx[2, 1]*Deta[k, j]
        Dy[k, j] = dxidx[1,2]*Dxi[k, j] + dxidx[2, 2]*Deta[k, j]
      end
    end

   for j=1:mesh.numNodesPerElement 
#     alpha_xj = eqn.alpha_x[1, j, i]
#     alpha_yj = eqn.alpha_y[1, j, i]
     dxidx_j = view(dxidx, :, :, j)
     p = 2
#     tau_vec[j] = getTau(alpha_xj, alpha_yj, dxidx_j, p)

     tau_vec[j] = getTau(alpha_x, alpha_y, jac[j], mesh.min_node_dist)
#     println("tau_vec[$j] = ", tau_vec[j])
     #DEBUGGING
     tau_sum += tau_vec[j]
     tau_cnt += 1
   end

   for j=1:mesh.numNodesPerElement
     u[j] = eqn.q[1, j, i]
   end


   # now compute the GLS term
   # trial space
   smallmatvec!(Dx, u, red_vec1)
   smallmatvec!(Dy, u, red_vec2)

   for j=1:mesh.numNodesPerElement
     red_vec1[j] = alpha_x*red_vec1[j] + alpha_y*red_vec2[j]
   end


   # add source term
   for j=1:mesh.numNodesPerElement
     coords_j = view(mesh.coords, :, j, i)
     red_vec1[j] -= src_func(coords_j, alpha_x, alpha_y, eqn.t)
   end

   #DEBUGGING
   for j=1:mesh.numNodesPerElement
     gls_res[1, j, i] = red_vec1[j]
   end

   
#   println("weighting space term = ", red_vec1)
   # middle terms
   # also copy pinto red_vec2 at same time
   for j=1:mesh.numNodesPerElement
     red_vec1[j] *= tau_vec[j]*sbp.w[j]/jac[j]
     red_vec2[j] = red_vec1[j]
     #DEBUGGING
     middle_term[1, j, i] = red_vec1[j]
   end

#   println("after middle term = ", red_vec1)

   smallmatTvec!(Dx, red_vec1, red_vec3)
   smallmatTvec!(Dy, red_vec2, red_vec4)

   for j=1:mesh.numNodesPerElement
     red_vec3[j] = alpha_x*red_vec3[j] + alpha_y*red_vec4[j]
   end


   weight_vals = alpha_x*smallmatTvec(Dx, weighting_vec)
   weight_vals2 =alpha_y*smallmatTvec(Dy, weighting_vec)


   # add source term to weighting term
   for j=1:mesh.numNodesPerElement
#     coords_j = view(mesh.coords, :, j, i)
#     q_j = u[j]
#     red_vec3[j] -= (src_func(coords_j, alpha_x, alpha_y, eqn.t)/q_j)*red_vec1[j]
#      red_vec3[j] = red_vec1[j]
      #DEBUGGING
      weighting_res[1, j, i] = weight_vals[j] + weight_vals2[j]
   end


   # update res
   for j=1:mesh.numNodesPerElement
     res[j] -= red_vec3[j]
     gls_full[1, j, i] -= red_vec3[j]
   end


  end  # end loop over elements

  gls_resvec = zeros(Tsol, mesh.numDof)
  full_resvec = zeros(Tsol, mesh.numDof)
  middle_vec = zeros(Tsol, mesh.numDof)
  old_resvec = zeros(Tsol, mesh.numDof)
  weight_resvec = zeros(Tsol, mesh.numDof)
  assembleSolution(mesh, sbp, eqn, opts, gls_res, gls_resvec)
  assembleSolution(mesh, sbp, eqn, opts, gls_full, full_resvec)
  assembleSolution(mesh, sbp, eqn, opts, middle_term, middle_vec)
  assembleSolution(mesh, sbp, eqn, opts, weighting_res, weight_resvec)
#  writedlm("gls_full.dat", real(gls_full))
#  writedlm("gls_fullvec.dat", full_resvec)
  assembleSolution(mesh, sbp, eqn, opts, res_before,  old_resvec)
  gls_norm = norm(gls_resvec, Inf)
  full_norm = norm(full_resvec, Inf)
  middle_norm = norm(middle_vec, Inf)
  old_norm = norm(old_resvec, Inf)
  weight_norm = norm(weight_resvec, Inf)
#  gls_norm = calcNorm(eqn, gls_resvec)
#  full_norm = calcNorm(eqn, full_resvec, strongres=true)
#  middle_norm = calcNorm(eqn, middle_vec, strongres=false)
#  old_norm = calcNorm(eqn, old_resvec, strongres=true)
#  weight_norm = calcNorm(eqn, weight_resvec, strongres=false)
  tau_avg = tau_sum/tau_cnt
  rmfile("gls_norm.dat")
  f = open("gls_norm.dat", "w")
  println(f, gls_norm, " ", old_norm, " ", full_norm, " ", real(tau_avg), " ", middle_norm, " ", weight_norm)
  close(f)




#  println("----- Finished applyGLS2 -----")
  return nothing

end  # end function




function getTau(alpha_x, alpha_y, dxidx::AbstractMatrix, p)
  fac = 2.5
  b1 = dxidx[1,1]*alpha_x + dxidx[2,1]*alpha_y
  b2 = dxidx[2,1]*alpha_x + dxidx[2,2]*alpha_y
  bp = absvalue(b1)^p + absvalue(b2)^p
  return fac*(bp^(-1/p))
end

function getTau{Tmsh}(alpha_x, alpha_y, jac::Tmsh, min_node_dist)

  fac = 1
  h = (1/sqrt(jac))/2  # /2 because reference element is -1 to 1
  alpha_nrm = sqrt(alpha_x*alpha_x + alpha_y*alpha_y)
  return fac*0.5*h/alpha_nrm
end
