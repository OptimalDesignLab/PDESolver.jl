# This file has implementations of Galerkin Least-Squares (GLS) method and 
# Streamwise Upwind Petrov-Galerkin method in it.

@doc """
### EulerEquationMod.GLS

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
                               eqn::EulerData{Tsol, Tres, Tdim})
  
  #  println("Entered GLS")
  parametricFluxJacobian(mesh, sbp, eqn) # Calculate the euler flux jacobian  
  tau = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode, 
              mesh.numNodesPerElement, mesh.numEl) # Stabilization term
  calcTau(mesh, sbp, eqn, tau, eqn.params.order)
  # println("Tau = \n", tau)
  
  gls_res = zeros(eqn.res) # input array for volumeintegrate!
  
  # Get shape function derivatives
  Hinv = 1./sbp.w
  shapefuncderiv = zeros(Tsol, sbp.numnodes, sbp.numnodes, Tdim)
  for k = 1:Tdim
    for i = 1:sbp.numnodes
      for j = 1:sbp.numnodes
        shapefuncderiv[i,j,k] = Hinv[i]*sbp.Q[i,j,k]
      end
    end
  end

  # Populate gls_res
  for i = 1:mesh.numEl
    endof = mesh.numDofPerNode*mesh.numNodesPerElement # dofs in an element
    ndof = mesh.numDofPerNode
    u = zeros(Tsol,endof) # element level q_vec
    # Populate u
    for j = 1:mesh.numNodesPerElement
      for k = 1:mesh.numDofPerNode
        u[(j-1)*ndof+k] = eqn.q[k,j,i]
      end
    end
    # println("eqn.q = \n", round(eqn.q,2))
    # println("u = \n", round(u,2))
    Axidxi = zeros(Tsol, endof, endof)
    Axi = view(eqn.Axi,:,:,:,i) # Get flux jacobians for all nodes in an element
    Aeta = view(eqn.Aeta,:,:,:,i)
    calcAxidxi(Axidxi, shapefuncderiv, Axi, Aeta, ndof, mesh.numNodesPerElement)
    # println("Axidxi = \n", round(Axidxi,2))
    
    intArr = zeros(Axidxi) # intermediate array for storing tau*Axidxi
    for j = 1:mesh.numNodesPerElement
      m = (j-1)*ndof + 1
      intArr[m:m+ndof-1,:] = sbp.w[j]*mesh.jac[j,i]*tau[:,:,j,i]*Axidxi[m:m+ndof-1,:]
    end
    
    intvec = Axidxi.'*intArr*u
    for j = 1:mesh.numNodesPerElement
      for k = 1:mesh.numDofPerNode
        gls_res[k,j,i] = intvec[(j-1)*ndof+k]
      end
    end
  end   # end for i = 1:mesh.numEl

  # volumeintegrate!(sbp, gls_res, supg_res)
  
  # println(supg_res)
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
### EulerEquationMod.parametricFluxJacobian

Claculates the flux jacobian in the parametric space for all the nodes in the
mesh. It uses conservative variables

**Inputs**

*  `mesh` : Abstract mesh object
*  `sbp`  : Summation-by-parts operator
*  `eqn`  : Equation object

**Outputs**

*  None
"""->

function parametricFluxJacobian{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh}, 
                                sbp::SBPOperator, 
                                eqn::EulerData{Tsol, Tres, Tdim})

# Calculates the Flux Jacobian in the parametric space

  fill!(eqn.Axi, 0.0)   # Reset Axi & Aeta to zero before calculating flux-
  fill!(eqn.Aeta, 0.0)  # jacobaian

  for i = 1:mesh.numEl
    for j = 1:mesh.numNodesPerElement
      q = view(eqn.q,:, j, i)
      dxidx = view(mesh.dxidx,:,:,j,i)
      # Calculate the flux jacobian get Ax & Ay
      Ax= zeros(Tsol, 4, 4) # Flux jacobian in the x direction
      Ay = zeros(Ax)        # Flux jacobian in the y direction
      calcFluxJacobian(eqn, q, Ax, Ay)
      # get the fluz jacobian in xi & eta direction
      eqn.Axi[:,:,j,i] = Ax*dxidx[1,1] + Ay*dxidx[1,2]
      eqn.Aeta[:,:,j,i] = Ax*dxidx[2,1] + Ay*dxidx[2,2]
    end  # end for j = 1:mesh.numNodesPerElement
  end    # end for i = 1:mesh.numEL

  return nothing
end

@doc """
### EulerEquationMod.calcFluxJacobian

It computes the Euler flux Jacobian at the nodal level in te physical space 
(X & Y axes). It uses conservative variables.

**Inputs**

*  `eqn`  : Equation object
*  `q`    : Conservative variables at the node
*  `Ax`   : Flux Jacobian at a node in the X-direction
*  `Ay`   : Flux Jacobian at a node in the Y-direction

**Outputs**

*  None 
"""->
function calcFluxJacobian{Tsol, Tres, Tdim}(eqn::EulerData{Tsol, Tres, Tdim},
                          q::AbstractArray{Tsol,1}, Ax::AbstractArray{Tsol,2}, 
                          Ay::AbstractArray{Tsol,2})

  gamma_1 = eqn.params.gamma_1
  gamma = eqn.params.gamma
  R = eqn.params.R     # Gas constant
  cv = eqn.params.cv   # Specific heat at constant volume
  u = q[2]/q[1] # Get velocity in the x-direction 
  v = q[3]/q[1] # Get velocity in the x-direction
  intvar = (R/cv)*(q[4]/q[1] - 0.5*(u*u + v*v)) # intermediate variable
  
  # Populating Ax
  Ax[1,1] = 0
  Ax[1,2] = 1
  Ax[1,3] = 0
  Ax[1,4] = 0
  Ax[2,1] = -u*u + 0.5*R*(u*u + v*v)/cv 
  Ax[2,2] = 2*u - R*u/cv
  Ax[2,3] = -R*v/cv
  Ax[2,4] = R/cv
  Ax[3,1] = -u*v
  Ax[3,2] = v
  Ax[3,3] = u
  Ax[3,4] = 0
  Ax[4,1] = -q[2]*q[4]/(q[1]*q[1]) - u*intvar + 0.5*u*R*(u*u + v*v)/cv
  Ax[4,2] = q[4]/q[1] + intvar - R*u*u/cv
  Ax[4,3] = -R*u*v/cv
  Ax[4,4] = u + R*u/cv

  # Populating Ay
  Ay[1,1] = 0
  Ay[1,2] = 0
  Ay[1,3] = 1
  Ay[1,4] = 0
  Ay[2,1] = -v*u 
  Ay[2,2] = v
  Ay[2,3] = u
  Ay[2,4] = 0
  Ay[3,1] = -v*v + 0.5*R*(u*u + v*v)/cv
  Ay[3,2] = -R*u/cv
  Ay[3,3] = 2*v - R*v/cv
  Ay[3,4] = R/cv
  Ay[4,1] = -q[3]*q[4]/(q[1]*q[1]) - v*intvar + v*(R/cv)*0.5*(u*u + v*v)
  Ay[4,2] = -R*v*u/cv
  Ay[4,3] = q[4]/q[1] + intvar - R*v*v/cv
  Ay[4,4] = v + R*v/cv
                                    
  return nothing
end # end function calcFluxJacobian


@doc """
### EulerEquationMod.calcTau

Calculates the stabilization term tau for GLS

**Inputs**

*  `mesh` : Abstract mesh type
*  `sbp`  : Summation-by-parts operator
*  `eqn`  : Euler equation object
*  `tau`  : Stabilization term

**Outputs**

* None
"""->
function calcTau{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh}, 
                 sbp::SBPOperator, eqn::EulerData{Tsol, Tres, Tdim},
                 tau::AbstractArray{Tsol,4}, order::Integer)

  # Using Glasby's Method Equation 26 & 27

  # Get shape function derivatives
  Hinv = 1./sbp.w
  shapefuncderiv = zeros(Tsol, sbp.numnodes, sbp.numnodes, Tdim)
  for k = 1:Tdim
    for i = 1:sbp.numnodes
      for j = 1:sbp.numnodes
        shapefuncderiv[i,j,k] = Hinv[i]*sbp.Q[i,j,k]
      end
    end
  end

  for i = 1:mesh.numEl
    for j = 1:mesh.numNodesPerElement
      q = view(eqn.q,:,j,i)
      dxidx = view(mesh.dxidx,:,:,j,i) # get mapping jacobian at node j
      for k = 1:mesh.numNodesPerElement
        for l = 1:Tdim
          jac_vector = Tsol[dxidx[l,1];dxidx[l,2]] # get a jacobian vector for eigen value factorization
          T = zeros(Tsol,4,4)
          Tinv = zeros(T)
          Lambda = eye(T)
          calcEigenFactorization(eqn, q, jac_vector, T, Tinv, Lambda)
          tau[:,:,j,i] += mesh.jac[j,i]*T*abs(shapefuncderiv[j,k,l]*Lambda)*Tinv
        end # end for l = 1:Tdim
      end # for k = 1:mesh.numNodesPerElement
      
      tau[:,:,j,i] = inv(tau[:,:,j,i])
    end # end for j = 1:mesh.numNodesPerElement
  end   # end for i = 1:mesh.numEl
  
  #= Using Masayuki's thesis and Barth's chapter. Using latter's convention
  
  for i = 1:mesh.numEl
    # Get the normal vectors at the nodes
    for j = 1:mesh.numNodesPerElement
      q = view(eqn.q,:,j,i)
      dxidx = view(mesh.dxidx,:,:,j,i)
      jac = view(mesh.jac,j,i)
      Ax= zeros(Tsol, 4, 4) # Flux jacobian in the x direction
      Ay = zeros(Ax)        # Flux jacobian in the y direction
      calcFluxJacobian(eqn, q, Ax, Ay)
      n = zeros(Tsol,2,2)
      n[:,1] = jac.*[dxidx[1,1];dxidx[1,2]]
      n[:,2] = jac.*[dxidx[2,1];dxidx[2,2]]
      n[:,1] = n[:,1]/sqrt(n[1,1]*n[1,1] + n[2,1]*n[2,1])
      n[:,2] = n[:,2]/sqrt(n[1,2]*n[1,2] + n[2,2]*n[2,2])
      for k = 1:Tdim
        Bn = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode)
        Bn = (n[1,k]*Ax + n[2,k]*Ay)*sqrt(n[1,k]*n[1,k] + n[2,k]*n[2,k])
        E = eigfact(Bn)
        # println("E.values = \n", E.values)
        # println("size of E.values = ", size(E.values))
        absLambda = eye(Bn)
        for l = 1:mesh.numDofPerNode
          absLambda[l,l] = abs(E.values[l])*absLambda[l,l]
        end
        tau[:,:,j,i] += E.vectors*absLambda*inv(E.vectors)
      end
      tau[:,:,j,i] = inv(tau[:,:,j,i])
    end # end for j = 1:mesh.numNodesPerElement
  end   # end for i = 1:mesh.numEl

  =#

  #=  ORIGINAL NOT WORKING TAU
  
  for i = 1:mesh.numEl
    for j = 1:mesh.numNodesPerElement
      q = view(eqn.q,:,j,i)
      T = (q[4] - 0.5*(q[2]*q[2] + q[3]*q[3])/q[1])*(1/(q[1]*eqn.params.cv))
      c = sqrt(eqn.params.gamma*eqn.params.R*T)  # Speed of sound
      ux = q[2]/q[1]
      uy = q[3]/q[1]
      h = circumcircleDiameter(mesh.coords[:,:,i])
      # println("h = ", h)
      tau[j,i] = h/(c + sqrt(ux*ux + uy*uy))
    end
  end  =#

  return nothing
end


@doc """
### EulerEquationMod.calcAxidxi

Calculates Axi*Dxi + Aeta*Deta at the element level

**Inputs**

*  `Axidxi` : Product of flux jacobian with shape function derivative
              Axi*Dxi + Aeta*Deta
*  `shapefuncderiv`: Shape function derivative (Dxi and Deta above)
*  `Axi` : Flux jacobian in the xi direction
*  `Aeta` : Flux jacobian in the eta direction

**Outputs**

*  None
"""->
function calcAxidxi{Tsol}(Axidxi::AbstractArray{Tsol, 2}, 
                          shapefuncderiv::AbstractArray{Tsol,3},
                          Axi::AbstractArray{Tsol,3}, 
                          Aeta::AbstractArray{Tsol,3}, ndof, nnpe)
  
  # ndof = mesh.numDofPerNode
  # nnpe = mesh.numNodesPerElement

  for i = 1:nnpe
    for j = 1:nnpe
      m = (i-1)*ndof+1
      n = (j-1)*ndof+1
      nAxi = view(Axi,:,:,i) # Flux Jacobian at a node in xi direction
      nAeta = view(Aeta,:,:,i)
      Axidxi[m:(m+ndof-1),n:(n+ndof-1)] = nAxi*shapefuncderiv[i,j,1] + 
                                          nAeta*shapefuncderiv[i,j,2] 
    end # end for j = 1:nnpe
  end   # end for i = 1:nnpe


  return nothing
end

@doc """
### EulerEquationMod. calcEigenFactorization

Computes the eigen value facotrization of the flux jacobian in either xi or eta
direction. This is a nodal level function. It is only for 2D conservative 
variables.

Reference: ``Efficient Solution Methods for the Navier–Stokes Equations``, T.H,
           Pulliam, NASA Ames Research Center

**Inputs**

*  `eqn`  : Euler equation object
*  `q`    : Conservative variables at a node
*  `dxidx`: mapping jacobian vector. It can be either [dξ/dx dξ/dy] or 
            [dη/dx dη/dy] depending on what Axi or Aeta being factorized.
*  `T`    : Matrix of eigen vectors.
*  `Tinv` : Inverse of T
*  `Lambda` : Diagonal matrix of eigen values

**Outputs**

*  None

"""->
function calcEigenFactorization{Tsol, Tres, Tdim}(eqn::EulerData{Tsol, Tres, Tdim},
                          q::AbstractArray{Tsol,1}, dxidx::AbstractArray{Tsol,1},
                          T::AbstractArray{Tsol,2}, Tinv::AbstractArray{Tsol,2},
                          Lambda::AbstractArray{Tsol,2})

  kappax = dxidx[1]
  kappay = dxidx[2]
  kappa = sqrt(kappax*kappax + kappay*kappay)
  u = q[2]/q[1] # velocity in the x-direction
  v = q[3]/q[1] # velocity in the y-direction
  Temp = (q[4] - 0.5*(q[2]*q[2] + q[3]*q[3])/q[1])*(1/(q[1]*eqn.params.cv)) # Temperature
  c = sqrt(eqn.params.gamma*eqn.params.R*Temp)  # Speed of sound
  phi2 = 0.5*eqn.params.gamma_1*(u*u + v*v) # phi^2
  beta = 0.5/(c*c)
  thetacap = (kappax*u + kappay*v)/kappa
  U = kappax*u + kappay*v
  
  # Matrix of eigen vectors
  T[1,1] = 1.0
  T[2,1] = u
  T[3,1] = v
  T[4,1] = phi2/eqn.params.gamma_1
  T[1,2] = 0.0
  T[2,2] = kappay/kappa
  T[3,2] = -kappax/kappa
  T[4,2] = (kappay*u - kappax*v)/kappa
  T[1,3] = 1.0
  T[2,3] = u + kappax*c/kappa
  T[3,3] = v + kappay*c/kappa
  T[4,3] = (phi2 + c*c)/eqn.params.gamma_1 + c*thetacap
  T[1,4] = 1.0
  T[2,4] = u - kappax*c/kappa
  T[3,4] = v - kappay*c/kappa
  T[4,4] = (phi2 + c*c)/eqn.params.gamma_1 - c*thetacap

  # Inverse of matrix of eigen vectors
  Tinv[1,1] = 1 - phi2/(c*c)
  Tinv[2,1] = -(kappay*u - kappax*v)/kappa
  Tinv[3,1] = beta*(phi2 - c*thetacap)
  Tinv[4,1] = beta*(phi2 + c*thetacap)
  Tinv[1,2] = eqn.params.gamma_1*u/(c*c)
  Tinv[2,2] = kappay/kappa
  Tinv[3,2] = beta*(kappax*c/kappa - eqn.params.gamma_1*u)
  Tinv[4,2] = -beta*(kappax*c/kappa + eqn.params.gamma_1*u)
  Tinv[1,3] = eqn.params.gamma_1*v/(c*c)
  Tinv[2,3] = -kappax/kappa
  Tinv[3,3] = beta*(kappay*c/kappa - eqn.params.gamma_1*v)
  Tinv[4,3] = -beta*(kappay*c/kappa + eqn.params.gamma_1*v)
  Tinv[1,4] = -eqn.params.gamma_1/(c*c)
  Tinv[2,4] = 0.0
  Tinv[3,4] = beta*eqn.params.gamma_1
  Tinv[4,4] = beta*eqn.params.gamma_1
  
  # Get the eigen values
  Lambda[1,1] = U
  Lambda[2,2] = U
  Lambda[3,3] = U + c*kappa
  Lambda[4,4] = U - c*kappa

  return nothing
end


@doc """
### EulerEquationMod.calcElementArea

Calculates the element area of a 2D triangular element

**Inputs**

*  `coords` : Vertex coordinates of the element

**Outputs**

*  `area` : Area of the triangular element
"""->
function calcElementArea{Tmsh}(coords::AbstractArray{Tmsh, 2})
  # Calculates the element area using coordinates
  # 2D function for linear mapping

  A = coords[:,1]
  B = coords[:,2]
  C = coords[:,3]
  area = 0.5*(A[1]*(B[2] - C[2]) + B[1]*(C[2] - A[2]) + C[1]*(A[2] - B[2]))
  
  return area
end # end function calcElementArea


@doc """
### EulerEquationMod.circumcircleDiameter

Calculates the circumcircle diameter of a triangular element

**Inputs**

*  `coords` : Vertex coordinates of the element

**Outputs**

*  `dia` : Diameter of the circumcircle
"""->
function circumcircleDiameter{Tmsh}(coords::AbstractArray{Tmsh, 2})
  # Calculates the circumcircle diameter of the triangular element
  # 2D function with linear mapping
  # Reference: http://geometryatlas.com/entries/109
  
  x0 = coords[1,1]
  y0 = coords[2,1]
  x1 = coords[1,2]
  y1 = coords[2,2]
  x2 = coords[1,3]
  y2 = coords[2,3]

  a = sqrt((x0-x1)^2 + (y0-y1)^2)
  b = sqrt((x0-x2)^2 + (y0-y2)^2)
  c = sqrt((x1-x2)^2 + (y1-y2)^2)

  dia = a*b*c/sqrt((a+b+c)*(b+c-a)*(c+a-b)*(a+b-c))

  return dia
end  # end function circumcircleDiameter



@doc """
### EulerEquationMod.SUPG

Add Streamwise Upwind Petrov-Galerkin stabilization to the weak residual. This
implementation is only for steady problems

**Inputs**

*  `mesh`: AbstractMesh type
*  `sbp` : Summation-by-parts operator
*  `eqn` : Equation object used elsewhere

**Outputs**

* None

"""->

# SUPG implementation
function SUPG{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator, 
                              eqn::EulerData{Tsol, Tres, Tdim})

  FluxJacobian(mesh, sbp, eqn) # Calculate the euler flux jacobian  
  tau = zeros(Tsol, mesh.numNodesPerElement, mesh.numEl) # Stabilization term
  calcStabilizationTerm(mesh, sbp, eqn, tau)
  # println("tau = \n", tau)
  
  #=
  # Calculate strong residual
  strong_res = zeros(eqn.res)
  for i = 1:Tdim
    flux_parametric_i = view(eqn.flux_parametric,:,:,:,i)
    differentiate!(sbp, i, flux_parametric_i, strong_res)
  end
  =#
  supg_res = zeros(eqn.res)
  intvec = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl, 
                 Tdim) # intermediate vector for calculating integral 

  # SUPG works at the element level so we need to do do a loop over elements
  # get the strong residual from the weak residual. since it also includes the 
  # boundary conditions
  
  for i = 1:mesh.numEl  
    for j = 1:mesh.numNodesPerElement
      strong_res = zeros(Tsol, mesh.numDofPerNode)
      JHinverse = mesh.jac[j,i]/sbp.w[j]
      for k = 1:mesh.numDofPerNode
        strong_res[k] = JHinverse*eqn.res[k,j,i]
      end
      Axi = view(eqn.Axi,:,:,j,i)
      Aeta = view(eqn.Aeta,:,:,j,i)
      intvec[:,j,i,1] = (tau[j,i]*Axi).'*strong_res # [:,j,i]
      intvec[:,j,i,2] = (tau[j,i]*Aeta).'*strong_res # [:,j,i]
    end # end for j = 1:mesh.numNodesPerElement
  end   # end for i = 1:mesh.numEl
    
  # calculate the SUPG residual  
  for i = 1:Tdim
    intvec_i = view(intvec,:,:,:,i)
    weakdifferentiate!(sbp, i, intvec_i,supg_res, trans=true)
  end
  
  
  #=
  supg_res_vec = zeros(eqn.res_vec)
  for i=1:mesh.numEl  # loop over elements
    for j=1:mesh.numNodesPerElement
      for k=1:4  # loop over dofs on the node
        dofnum_k = mesh.dofs[k, j, i]
        supg_res_vec[dofnum_k] += supg_res[k,j,i]
      end
    end
  end
  =# 
  # Add tthe SUPG residual to the weak residual
  for i = 1:mesh.numEl
    for j = 1:mesh.numNodesPerElement
      for k = 1:mesh.numDofPerNode
        eqn.res[k,j,i] += supg_res[k,j,i] # because the negative sign is 
      end                                 # already incorporated in the weak 
    end                                   # residual.
  end
  
  # innerprod_supg = eqn.q_vec.'*supg_res_vec
  # println("innerprod_SUPG = ", innerprod_supg)
  
  #  println("eqn.res = \n", eqn.res)
  #  println("eqn.q = \n", eqn.q)
  return nothing
end # end function SUPG

@doc """
### EulerEquationMod.calcStabilizationTerm

Calculates the stabilization term tau for all the nodes in the mesh

**Inputs**

*  `mesh` : Abstract mesh type
*  `sbp`  : Summation-by-parts operator
*  `eqn`  : Equation object within EulerEquationMod
*  `tau`  : Stabilization term being calculated

**Outputs**

*  None
"""->
# Stabilization Term 3
function calcStabilizationTerm{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh}, 
                               sbp::SBPOperator, eqn::EulerData{Tsol, Tres, Tdim},
                               tau::AbstractArray{Tsol,2})
  
  # Reference: http://enu.kz/repository/2010/AIAA-2010-1183.pdf, eqn 15

  # Get shape function derivatives
  Hinv = 1./sbp.w
  # println(Hinv)
  
  shapefuncderiv = zeros(sbp.numnodes, sbp.numnodes, Tdim)
  for k = 1:Tdim
    for i = 1:sbp.numnodes
      for j = 1:sbp.numnodes
        shapefuncderiv[j,i,k] = Hinv[i]*sbp.Q[j,i,k]
      end
    end
  end
  # println(eqn.Axi)
  
  
  for i = 1:mesh.numEl
    for j = 1:mesh.numNodesPerElement
      q = view(eqn.q,:,j,i)
      T = (q[4] - 0.5*(q[2]*q[2] + q[3]*q[3])/q[1])*(1/(q[1]*eqn.params.cv))
      c = sqrt(eqn.params.gamma*eqn.params.R*T)  # Speed of sound
      ux = q[2]/q[1]
      uy = q[3]/q[1]
      h_supg = 0.0
      for k = 1:sbp.numnodes
        h_supg += abs(ux*shapefuncderiv[j,k,1] + uy*shapefuncderiv[j,k,2])
      end
      h_supg = 2/h_supg
      tau[j,i] = 0.5*h_supg/(c + sqrt(ux*ux + uy*uy))
    end
  end


  return nothing
end # end calcStabilizationTerm

#=
# Stabilization Term 1
function calcStabilizationTerm{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh}, 
                               sbp::SBPOperator, eqn::EulerData{Tsol, Tres, Tdim},
                               tau::AbstractArray{Tsol,2})
  
  # q in the parametric space. Since everything happens in this space
  q_param = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)
  
  for i = 1:mesh.numEl
    for j = 1:mesh.numNodesPerElement
      q = view(eqn.q,:,j,i)
      dxidx = view(mesh.dxidx,:,:,j,i)
      rhoe = q[4] -0.5*(q[2]*q[2] + q[3]*q[3])/(q[1]*q[1]) 
      uxi = (q[2]*dxidx[1,1] + q[3]*dxidx[1,2])/q[1]
      ueta = (q[2]*dxidx[2,1] + q[3]*dxidx[2,2])/q[1]
      q_param[1,j,i] = q[1]
      q_param[2,j,i] = q[1]*uxi
      q_param[3,j,i] = q[1]*ueta
      q_param[4,j,i] = rhoe
    end # end for j = 1:mesh.numNodesPerElement
  end   # end for i = 1:mesh.numEl

  beta = zeros(Tsol,2, mesh.numNodesPerElement, mesh.numEl)
  
  for k = 1:Tdim
    res = zeros(q_param)
    differentiate!(sbp, k, q_param, res)
    for i = 1:mesh.numEl
      for j = 1:mesh.numNodesPerElement
        for l = 1:mesh.numDofPerNode
          beta[k,j,i] += q_param[l,j,i]*res[l,j,i]
        end # end for l = 1:mesh.numDofPerNode
      end   # end for j = 1:mesh.numNodesPerElement
    end     # end for i = 1:mesh.numEl
  end       # end for k = 1:Tdim

  for i = 1:mesh.numEl
    elem_area = calcElementArea(mesh.coords[:,:,i])
    h = sqrt(2*elem_area)
    for j = 1:mesh.numNodesPerElement
      beta[:,j,i] = beta[:,j,i]/norm(beta[:,j,i],2)
      q = view(eqn.q,:,j,i)
      T = (q[4] - 0.5*(q[2]*q[2] + q[3]*q[3])/q[1])*(1/(q[1]*eqn.params.cv))
      c = sqrt(eqn.params.gamma*eqn.params.R*T)  # Speed of sound
      uxi = zeros(Tsol,2) # Array of velocities in the xi & eta direction
      uxi[1] = q_param[2,j,i]/q_param[1,j,i]
      uxi[2] = q_param[3,j,i]/q_param[1,j,i]
      # Advective stabilization term
      tau_a = 0.5*h/(c + norm(uxi.'*beta[:,j,i], 1))
      # tau[j,i] = 0.005*max(0.0, tau_a)
      tau[j,i] = max(0.0, tau_a)
    end # end for j = 1:mesh.numNodesPerElement
  end   # for i = 1:mesh.numEl

  return nothing
end # end calcStabilizationTerm
=#

#=
# Stabilization term 2
function calcStabilizationTerm{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh}, 
                               sbp::SBPOperator, eqn::EulerData{Tsol, Tres, Tdim},
                               tau::AbstractArray{Tsol,2})
  
  # Reference for stabilization: http://enu.kz/repository/2010/AIAA-2010-1183.pdf

  for i = 1:mesh.numEl
    for j = 1:mesh.numNodesPerElement
      C = 2 # for 2D quad with ref. coord [-1,1] 
            # for 2D tri with ref coord [0,1], C = 1
      q = view(eqn.q,:,j,i)
      T = (q[4] - 0.5*(q[2]*q[2] + q[3]*q[3])/q[1])*(1/(q[1]*eqn.params.cv))
      c = sqrt(eqn.params.gamma*eqn.params.R*T)  # Speed of sound
      ux = q[2]/q[1]
      uy = q[3]/q[1]
      dxidx = view(mesh.dxidx,:,:, j, i)
      if abs(ux) < 1e-13 || abs(uy) < 1e-13
        tau[j,i] = 0.0
      else
        g_ij = dxidx[1,1]*dxidx[1,2] + dxidx[2,1]*dxidx[2,2]
        if abs(g_ij) < 1e-13
          tau[j,i] = 0.0
        else
          elem_area = calcElementArea(mesh.coords[:,:,i])
          # h_supg = sqrt(2*elem_area)
          h_supg = C*sqrt((ux*ux + uy*uy)/(ux*g_ij*uy))
          tau[j,i] = 0.5*h_supg/(c + sqrt(ux*ux + uy*uy))
        end # end if
      end   # end if
    end # end for j = 1:mesh.numNodesPerElement
  end   # end for i = 1:mesh.numEl

  return nothing
end # end calcStabilizationTerm
=#
#=
# Stabilization term 4
function calcStabilizationTerm{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh}, 
                               sbp::SBPOperator, eqn::EulerData{Tsol, Tres, Tdim},
                               tau::AbstractArray{Tsol,2})

  # Reference: Three-Dimensional Stabilized Finite Elements for Compressible 
  #            Navier–Stokes, T. Taylor Erwin, AIAA Journal Vol 51, No. 6,
  #            June 2013

  # Get shape function derivatives
  Hinv = 1./sbp.w
  # println(Hinv)
  
  shapefuncderiv = zeros(sbp.numnodes, sbp.numnodes, Tdim)
  for k = 1:Tdim
    for i = 1:sbp.numnodes
      for j = 1:sbp.numnodes
        shapefuncderiv[j,i,k] = Hinv[i]*sbp.Q[j,i,k]
      end
    end
  end

  for i = 1:mesh.numEl
    for j = 1:mesh.numNodesPerElement
      Axi = view(eqn.Axi,:,:,j,i)
      Aeta = view(eqn.Aeta,:,:,j,i)
      invtau = zeros(mesh.numDofPerNode, mesh.numDofPerNode)
      for k = 1:sbp.numnodes
        invtau += shapefuncderiv[j,k,1]*Axi + shapefuncderiv[j,k,2]*Aeta
        #println("invtau = ", invtau)
      end
      T = inv(invtau)
      println("T = \n", T)
    end
  end
  

  return nothing
end # end calcStabilizationTerm
=#

#------------------------------------------------------------------------------
# Debugging code
# calculate the boundary integral using the actual euler flux

function getPhysBCFluxes(mesh, sbp, eqn, opts, bndryfluxPhysical)
  # Calculate the physical BC flux

  #println("mesh.bndry_funcs = ", mesh.bndry_funcs)
  
  functor_i = BCDict["isentropicVortexBC_physical"]
  for i=1:mesh.numBC
  #  println("computing flux for boundary condition ", i)
    # functor_i = mesh.bndry_funcs[i]
    start_index = mesh.bndry_offsets[i]
    end_index = mesh.bndry_offsets[i+1]
    bndry_facenums_i = view(mesh.bndryfaces, start_index:(end_index - 1))
    bndryflux_i = view(bndryfluxPhysical, :, :, start_index:(end_index - 1))
    #functor_i(q, flux_parametric, aux_vars, x, dxidx, nrm, bndryflux_i, eqn.params)

    calcBoundaryFlux(mesh, sbp, eqn, functor_i, bndry_facenums_i, bndryflux_i)
  end


  return nothing
end


function residualComparison(mesh, sbp, eqn, opts)

 # Get strong residual by differentiating
  differentiation_strong_res = zeros(eqn.res)
  
  for i = 1:Tdim
    flux_parametric_i = view(eqn.flux_parametric,:,:,:,i)
    differentiate!(sbp, i, flux_parametric_i, differentiation_strong_res)
  end
  

  # Get strong residual from the weak form
  strong_res_from_weak = zeros(eqn.res)
  for i = 1:mesh.numEl  
    for j = 1:mesh.numNodesPerElement
      JHinverse = mesh.jac[j,i]/sbp.w[j]
      for k = 1:mesh.numDofPerNode
        strong_res_from_weak[k,j,i] = JHinverse*eqn.res[k,j,i]
      end # end for k = 1:mesh.numDofPerNode
    end # end for i = 1:mesh.numEl 
  end   # end for j = 1:mesh.numNodesPerElement

  println("differentiation_strong_res = \n", differentiation_strong_res)
  println("strong_res_from_weak = \n", strong_res_from_weak)
  # println("\n eqn.res = \n", eqn.res)
  Error = zeros(eqn.res)
  for i = 1:mesh.numEl
    Error[:,:,i] = differentiation_strong_res[:,:,i] + strong_res_from_weak[:,:,i]
  end
  println("\nError = \n", Error)

  return nothing
end # end function residualComparison


#----------------------------------------------------------------------------
