# this function applies the form of GLS stabilization described in:
# A new finite element formulation for computational fluid dynamics:
# X The compressible Euler and Navier-Stokes Equations
# by Shakib, Hughes, and Johan

# this *should* work for both conservative and entropy variables, I think
function applyGLS2{Tmsh, Tsol, Tres}(mesh::AbstractMesh{Tmsh}, 
                   sbp::SBPOperator, eqn::EulerData{Tsol, Tres, Tdim}, opts)

  # reusable storage

  # flux jacobians: a n x n matrix for each coordinate direction for each node
  A_mats = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode, Tdim, 
                      mesh.numNodesPerElement)
  dxidx_hat = zeros(Tmsh, Tdim, Tdim, mesh.numNodesPerElement)

  # temporarily hold the transposed q variables for the element
  qtranspose = zeros(Tsol, mesh.numNodesPerElement, mesh.numDofPerNode)
  # temporarily hold the result of D*qtranspose
  qxtranspose = zeros(Tsol, mesh.numNodesPerElement, mesh.numDofPerNode, Tdim)
  # hold un-transposed qx
  qx = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerElement, Tdim)
  # hold tau, a n x n matrix for each node
  tau = zeros(Tres, mesh.numDofPerNode, mesh.numDofPerNode, 
                    mesh.numNodesPerElement)

  # temporary vectors 1 through 4, stored in a matrix so they can be
  # access programaticaly
  tmps = zeros(Tres, mesh.numDofPerNode, 2*Tdim)

  for el = 1:mesh.numEl
    # get all the quantities for this element






  end


end  # end function

function getGLSVars{Tmsh, Tsol, Tres, Tdim}(params::ParamType{Tdim},
                    q::AbstractArray{Tsol, 2}, aux_vars::AbstractArray{Tsol, 2},
                    dxidx_hat::AbstractArray{Tmsh, 3}, 
                    jac::AbstractArray{Tmsh, 1},
                    D::AbstractArray{Float64, 3}  
                    A_mats::AbstractArray{Tsol, 4}, # begin outputs
                    qtranspose::AbstractArray{Tsol, 2}, 
                    qxtranspose::AbstractArray{Tsol, 3}, 
                    qx::AbstractArray{Tsol, 3},
                    dxidx::AbstractArray{Tmsh, 3}, 
                    tau::AbstractArray{Tres, 3})


  numNodesPerElement = size(q, 2)
  numDofPerNode = size(q, 1)
  # move q into q_transpose
  for i=1:numDofPerNode
    for j=1:numNodesPerElement
      q_transpose[j, i] = q[i, j]
    end
  end

  # multiply by Dxi, Deta ...
  for d=1:Tdim
    smallmatmatT!(view(D, :, :, d), q_transpose, view(qxtranspose, :, :, d))
  end

  # un-transpose while putting into qx
  for d=1:Tdim
    for j=1:numNodesPerElement
      for i=1:numDofPerNode
        qx[i, j, d] = qxtranpose[j, i, d]
      end
    end
  end

  
  # get flux jacobian
  for k=1:numNodesPerElement
    q_k = view(q, :, k)
    A1_k = view(A_mats, :, :, 1, k)
    calcA1(params, q_k, A1_k)
    A2_k = view(A_mats, :, :, 2, k)
    calcA2(params, q_k, A2_k)

    if Tdim == 3  # three cheers for static analysis
      A3_k = view(A_mats, :, :, 3, k)
      calcA3(params, q_k, A3_k)
    end
  end

  # get dxidx - not scaled by 1/|J|
  for k = 1:numNodesPerElement
    for j=1:Tdim
      for i=1:Tdim
        dxidx[i, j, k] = dxidx_hat*jac[k]
      end
    end
  end

  # get tau
  for k=1:numNodesPerElement
    tau_k = view(tau, :, :, k)
    getTau
  end

  return nothing
end  # end function

function getTau{Tdim, var_type, Tsol, Tres, Tmsh}(
                params::ParamType{Tdim, var_type, Tsol, Tres, Tmsh}, 
                A_mat::AbstractArray{Tsol, 3}, dxidx::AbstractArray{Tmsh, 3}, 
                tau::AbstractArray{T, 2})

  numDofPerNode = size(A_mat, 1)
  tmp = params.A1
  tmp2 = params.A2

  for k=1:Tdim
    for j=1:Tdim
      Aj = view(A_mat, :, :, j)
      Ak = view(A_mat, :, :, k)
      smallmatmat(Aj, Ak, tmp)
      fac = zero(Tmsh)
      for i=1:Tdim
        fac += dxidx[i, j]*dxidx[i, k]
      end

      # accumulate in tmp2
      for p=1:numDofPerNode
        for q=1:numDofPerNode
          tmp2[p, q] += fac*tmp[p, q]
        end
      end

    end
  end
      # add a source term here


      # now take negative square root of the matrix
      # there are more efficient ways of doing this
      D, V = eig(tmp2)
      fill!(tmp, 0.0)
      for i=1:numDofPerNode
        D[i] = D[i]^(-0.5)
      end
      
      Vinv = inv(V)
      # reconstruct M = (V*D)*inv(V)
      for i=1:numDofPerNode
        entry = D[i]
        for j=1:numDofPerNode
          V[j, i] *= entry
        end
      end

      smallmatmat(V, Vinv, tau)

      
