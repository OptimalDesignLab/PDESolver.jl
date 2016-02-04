

function applyGLS2{Tmsh, Tsol, Tres}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator, eqn::AdvectionData{Tsol, Tres, 2}, opts)

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

  for i=1:1  # DEBUGGING: only do first element
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
     alpha_xj = eqn.alpha_x[1, j, i]
     alpha_yj = eqn.alpha_y[1, j, i]
     dxidx_j = view(dxidx, :, :, j)
     p = 2
     tau_vec[j] = getTau(alpha_xj, alpha_yj, dxidx_j, p)
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

#   println("weighting space term = ", red_vec1)
   # middle terms
   # also copy pinto red_vec2 at same time
   for j=1:mesh.numNodesPerElement
     red_vec1[j] *= tau_vec[j]*sbp.w[j]/jac[j]
     red_vec2[j] = red_vec1[j]
   end

#   println("after middle term = ", red_vec1)

   smallmatTvec!(Dx, red_vec1, red_vec3)
   smallmatTvec!(Dy, red_vec2, red_vec4)

   for j=1:mesh.numNodesPerElement
     res[j] -= alpha_x*red_vec3[j] + alpha_y*red_vec4[j]
   end


  end  # end loop over elements

#  println("----- Finished applyGLS2 -----")
  return nothing

end  # end function



function getTau(alpha_x, alpha_y, dxidx::AbstractMatrix, p)

  b1 = dxidx[1,1]*alpha_x + dxidx[2,1]*alpha_y
  b2 = dxidx[2,1]*alpha_x + dxidx[2,2]*alpha_y
  bp = absvalue(b1)^p + absvalue(b2)^p
  return bp^(1/p)
end
