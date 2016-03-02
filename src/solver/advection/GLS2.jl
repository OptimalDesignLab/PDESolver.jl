@doc """
### AdvectionEquationMod.applyGLS2

  This function updates the residual with the GLS term, using the definition
  of the parameter tau from Hughes Part III.

  Inputs:
    mesh
    sbp
    eqn
    opts
    src_func:  The functor that calculates the source term at a node

  Aliasing restrictions: none
"""->
function applyGLS2{Tmsh, Tsol, Tres}(mesh::AbstractMesh{Tmsh}, sbp::AbstractSBP, eqn::AdvectionData{Tsol, Tres, 2}, opts, src_func::SRCType)

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

   # middle terms
   # also copy pinto red_vec2 at same time
   for j=1:mesh.numNodesPerElement
     red_vec1[j] *= tau_vec[j]*sbp.w[j]/jac[j]
     red_vec2[j] = red_vec1[j]
   end

   smallmatTvec!(Dx, red_vec1, red_vec3)
   smallmatTvec!(Dy, red_vec2, red_vec4)

   for j=1:mesh.numNodesPerElement
     red_vec3[j] = alpha_x*red_vec3[j] + alpha_y*red_vec4[j]
   end


   # update res
   for j=1:mesh.numNodesPerElement
     res[j] -= red_vec3[j]
   end


  end  # end loop over elements


  #  println("----- Finished applyGLS2 -----")
  return nothing

end  # end function



@doc """
### AdvectionEquationMod.getTau

  This function calculates the stabilization parameter tau for GLS at a node, 
  using the definition given in Hughes Part III

  Inputs:
    alpha_x: advection coefficient in the x direction at the node
    alpha_y: advection coefficient in the y direction at the node
    dxidx: a 2x2 matrix containing the mapping jacobian at the node
    p : which p norm to use (typically either 1 or 2)

  Outputs:
    tau: the value of tau at the node

"""->
function getTau(alpha_x, alpha_y, dxidx::AbstractMatrix, p)
  fac = 1.0
  b1 = dxidx[1,1]*alpha_x + dxidx[1,2]*alpha_y
  b2 = dxidx[2,1]*alpha_x + dxidx[2,2]*alpha_y
  bp = absvalue(b1)^p + absvalue(b2)^p
  return fac*(bp^(-1/p))
end

@doc """
### AdvectionEquationMod.getTau

  This function calculate the stabilization parameter tau at a node using 
  the (heuristic) definition 0.5*h/sqrt(alpha_x^2 + alpha_y^2)

  Inputs:
    alpha_x: advection coefficient in the x direction at the node
    alpha_y: advection coefficient in the y direction at the node
    jac: determinant of the mapping jacobian at the node

  Outputs:
    tau: the value of tau at the node
"""->
function getTau{Tmsh}(alpha_x, alpha_y, jac::Tmsh, min_node_dist)

  fac = 1
  h = (1/sqrt(jac))/2  # /2 because reference element is -1 to 1
  alpha_nrm = sqrt(alpha_x*alpha_x + alpha_y*alpha_y)
  return fac*0.5*h/alpha_nrm
end
