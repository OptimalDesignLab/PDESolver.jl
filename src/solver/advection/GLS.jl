# GLS.jl for advection equation

@doc """
### AdvectionEquationMod.GLS

Add Galerkin Least-Squares stabilization to the weak residual. This
implementation is only for steady problems and conservative variables

**Inputs**

*  `mesh`: AbstractMesh type
*  `sbp` : Summation-by-parts operator
*  `eqn` : Equation object used elsewhere

**Outputs**

* None

"""->
function GLS{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator, 
                               eqn::AdvectionData{Tsol, Tres, Tdim})

  tau = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode, 
              mesh.numNodesPerElement, mesh.numEl)
  # calculate tau

  # Calculate the shape function derivatives
  Hinv = 1./sbp.w
  shapefuncderiv = zeros(Tsol, sbp.numnodes, sbp.numnodes, Tdim)
  for k = 1:Tdim
    for i = 1:sbp.numnodes
      for j = 1:sbp.numnodes
        shapefuncderiv[i,j,k] = Hinv[i]*sbp.Q[i,j,k]
      end
    end
  end 
  
  gls_res = zeros(eqn.res)

  for i = 1:mesh.numEl
    u = zeros(Tsol, mesh.numNodesPerElement)
    for j = 1:mesh.numNodesPerElement
      u[j] = eqn.q[1,j,i]  # get the element level q vector
    end
    AxiDxi = zeros(Tsol, mesh.numNodesPerElement, numNodesPerElement)
    alpha_x = view(eqn.alpha_x, 1, :, i)
    alpha_y = view(eqn.alpha_y, 1, :, i)
    dxidx = view(mesh.dxidx, :, :, :, i)
    calcAxiDxi(mesh, AxiDxi, dxidx, alpha_x, alpha_y, shapefuncderiv)

    intArr = zeros(AxiDxi) # intermediate array for storing tau*AxiDxi
    for j = 1:mesh.numNodesPerElement
      intArr = sbp.w[j]*tau[j,i]*AxiDxi
    end
    intvec = Axidxi.'*intArr*u
    gls_res[1,:,i] = intvec[:]
  end  # end for i = 1:mesh.numEl

  # Add to eqn.res
  for i = 1:mesh.numEl
    for j = 1:mesh.numNodesPerElement
      for k = 1:mesh.numDofPerNode
        eqn.res[k,j,i] -= gls_res[k,j,i]  # -ve sign because of right hand side.
      end                                 # Fluxes are computed explicitly and
    end                                   # not from weak redidual.
  end

  return nothing
end

function calcAxiDxi{Tmsh, Tsol}(mesh::AbstractMesh{Tmsh}, 
	                            AxiDxi::AbstractArray{Tsol,2},  
	                            dxidx::AbstractArray{Tmsh,3},
	                            alpha_x::AbstractArray{Tsol,1}, 
	                            alpha_y{Tsol,1}, shapefuncderiv),
                                shapefuncderiv::AbstractArray{Tsol,3})

  alpha_xi = eye(mesh.numNodesPerElement)
  alpha_eta = eye(mesh.numNodesPerElement)

  # Populate alpha_xi & alpha_eta
  for i = 1:mesh.numNodesPerElement
    alpha_xi[i,i] = dxidx[1, 1, i]*alpha_x[i] + dxidx[1, 2, i]*alpha_y[i]
    alpha_eta[i,i] = dxidx[2, 1, i]*alpha_x[i] + dxidx[2, 2, i]*alpha_y[i]
  end

  AxiDxi = alpha_xi*shapefuncderiv[:,:,1] + alpha_eta*shapefuncderiv[:,:,2]  

  return nothing
end