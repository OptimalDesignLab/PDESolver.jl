# GLS.jl for advection equation

@doc """
### AdvectionEquationMod.GLS

Add Galerkin Least-Squares stabilization to the weak residual. This
implementation is only for steady problems and conservative variables

**Inputs**

*  `mesh`: AbstractMesh type
*  `sbp` : Summation-by-parts operator
*  `eqn` : Advection equation object

**Outputs**

* None

"""->
function GLS{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh}, sbp::AbstractSBP, 
                                     eqn::AdvectionData{Tsol, Tres, Tdim})
  

  gls_res = zeros(eqn.res)
  tau = zeros(Tsol,mesh.numNodesPerElement, mesh.numEl)
  calcTau(mesh, sbp, eqn, tau)

  # Calculate the shape function derivatives
  shapefuncderiv = zeros(Tsol, sbp.numnodes, sbp.numnodes, Tdim)
  calcShapefuncDeriv(sbp, shapefuncderiv)
  
  for i = 1:mesh.numEl
    u = zeros(Tsol, mesh.numNodesPerElement)
    # get the element level q vector
    for j = 1:mesh.numNodesPerElement
      u[j] = eqn.q[1,j,i]  
    end
    alpha_x = sview(eqn.alpha_x, 1, :, i)
    alpha_y = sview(eqn.alpha_y, 1, :, i)
    dxidx = sview(mesh.dxidx, :, :, :, i)
    AxiDxi = calcAxiDxi(mesh, dxidx, alpha_x, alpha_y, shapefuncderiv)
    intArr = zeros(AxiDxi) # intermediate array for storing H*tau*AxiDxi
    for j = 1:mesh.numNodesPerElement
      intArr[j,:,i] = sbp.w[j]*tau[j]*AxiDxi[j,:,i]
    end
    
    intvec = AxiDxi.'*intArr*u
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


@doc """
### AdvectionEquationMod.calcAxidxi

Calculates Axi*Dxi + Aeta*Deta at the element level

**Inputs**

*  `mesh`   : Abstract mesh type
*  `Axidxi` : Product of flux jacobian with shape function derivative
              Axi*Dxi + Aeta*Deta
*  `dxidx`  : Mapping jacobian for all the nodes in an element
*  `alpha_x` & `alpha_y`: Velocities in the X & Y direction for all the nodes
                          in an element.
*  `shapefuncderiv`: Shape function derivative (Dxi and Deta above)

**Outputs**

*  None

"""->
function calcAxiDxi{Tmsh, Tsol}(mesh::AbstractMesh{Tmsh}, 
	                              dxidx::AbstractArray{Tmsh,3},
	                              alpha_x::AbstractArray{Tsol}, 
	                              alpha_y::AbstractArray{Tsol}, 
                                shapefuncderiv::AbstractArray{Tsol,3})

  alpha_xi = zeros(Tsol, mesh.numNodesPerElement, mesh.numNodesPerElement)
  alpha_eta = zeros(alpha_xi)

  # Populate alpha_xi & alpha_eta
  for i = 1:mesh.numNodesPerElement
    alpha_xi[i,i] = dxidx[1, 1, i]*alpha_x[i] + dxidx[1, 2, i]*alpha_y[i]
    alpha_eta[i,i] = dxidx[2, 1, i]*alpha_x[i] + dxidx[2, 2, i]*alpha_y[i]
  end
  AxiDxi = alpha_xi*shapefuncderiv[:,:,1] + alpha_eta*shapefuncderiv[:,:,2]
  
  return AxiDxi
end


@doc """
### AdvectionEquationMod.calcTau

Calculates the stabilization term tau for GLS. It operates at the global level.

**Inputs**

*  `mesh` : Abstract mesh type
*  `sbp`  : Summation-by-parts operator
*  `eqn`  : Advection equation object
*  `tau`  : Stabilization term. It is an element level array with size equal 
            to number of nodes in an element.

**Outputs**

* None
"""->
function calcTau{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh}, 
                 sbp::AbstractSBP, eqn::AdvectionData{Tsol, Tres, Tdim},
                 tau::AbstractArray{Tsol,2})
  
  # Using Glasby's implementation for advection equation

  # Get shape function derivatives at the nodes
  shapefuncderiv = zeros(Tsol, sbp.numnodes, sbp.numnodes, Tdim)
  calcShapefuncDeriv(sbp, shapefuncderiv)

  for i = 1:mesh.numEl
    alpha_xi = zeros(Tsol, mesh.numNodesPerElement)
    alpha_eta = zeros(alpha_xi)
  	for j = 1:mesh.numNodesPerElement
      alpha_x = eqn.alpha_x[1, j, i]
      alpha_y = eqn.alpha_y[1, j, i]
      alpha_xi[j] = mesh.dxidx[1, 1, j, i]*alpha_x + mesh.dxidx[1, 2, j, i]*alpha_y
      alpha_eta[j] = mesh.dxidx[2, 1, j, i]*alpha_x + mesh.dxidx[2, 2, j, i]*alpha_y
    end  # end for j = 1:mesh.numNodesPerElement
    
    for j = 1:mesh.numNodesPerElement
      for k = 1:mesh.numNodesPerElement
        tau[j,i] += abs(alpha_xi[j]*shapefuncderiv[j,k,1]) + 
                    abs(alpha_eta[j]*shapefuncderiv[j,k,2])
      end
      tau[j,i] = 1/tau[j,i]
    end
  
  end    # end for i = 1:mesh.numEl

  return nothing
end


@doc """
### AdvectionEquationMod.calcShapefuncDeriv

Calculates the shape function derivatives for a 2D problem

**Inputs**

*  `sbp` : Summation-by-parts operator
*  `shapefuncderiv` : Array containing the shape functions

**Outputs**

*  None

"""->
function calcShapefuncDeriv{Tsol}(sbp::AbstractSBP, 
                                  shapefuncderiv::AbstractArray{Tsol,3})

  Tdim = 2  # For a 2D problem
  Hinv = 1./sbp.w
  for k = 1:Tdim
    for i = 1:sbp.numnodes
      for j = 1:sbp.numnodes
        shapefuncderiv[i,j,k] = Hinv[i]*sbp.Q[i,j,k]
      end
    end
  end  

  return nothing
end