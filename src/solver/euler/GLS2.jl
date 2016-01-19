# this function applies the form of GLS stabilization described in:
# A new finite element formulation for computational fluid dynamics:
# X The compressible Euler and Navier-Stokes Equations
# by Shakib, Hughes, and Johan

# this *should* work for both conservative and entropy variables, I think
function applyGLS2{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh}, 
                   sbp::SBPOperator, eqn::EulerData{Tsol, Tres, Tdim}, opts)

  println("----- Entered applyGLS2 -----")
  # extract some constants
  numDofPerNode = mesh.numDofPerNode
  numNodesPerElement = mesh.numNodesPerElement
  w = sbp.w
  # reusable storage

  # flux jacobians: a n x n matrix for each coordinate direction for each node
  A_mats = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode, Tdim, 
                      mesh.numNodesPerElement)
#  dxidx_hat = zeros(Tmsh, Tdim, Tdim, mesh.numNodesPerElement)

  # temporarily hold the transposed q variables for the element
  qtranspose = zeros(Tsol, mesh.numNodesPerElement, mesh.numDofPerNode)
  # temporarily hold the result of D*qtranspose
  qxitranspose = zeros(Tsol, mesh.numNodesPerElement, mesh.numDofPerNode, Tdim)
  # hold un-transposed qxi (temporarily)
  qxi = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerElement, Tdim)
  # hold non-scaled mapping jacobian
  dxidx = zeros(Tmsh, Tdim, Tdim, numNodesPerElement)
  # hold tau, a n x n matrix for each node
  taus = zeros(Tres, mesh.numDofPerNode, mesh.numDofPerNode, 
                    mesh.numNodesPerElement)

  # hold differentiation matrices D
  D = zeros(Float64, numNodesPerElement, numNodesPerElement, Tdim)
  Dtranspose = zeros(D)  # transpose D because we access it in that order
  # temporary vectors 1 through 4, stored in a matrix so they can be
  # access programaticaly
  tmps = zeros(Tres, mesh.numDofPerNode, 2*Tdim)
  qx = zeros(Tres, numDofPerNode, Tdim)  # accumulation vectors

  # calculate D
  for d=1:Tdim
    smallmatmat!(diagm(1./sbp.w), view(sbp.Q, :, :, d), view(D, :, :, d))
    Dtranspose[:, :, d] = D[:, :, d].'
  end

  
   for el = 1:1  # DEBUGGING: do only first element
#  for el = 1:mesh.numEl
    # get all the quantities for this element
    dxidx_hat = view(mesh.dxidx, :, :, :, el)
    aux_vars = view(eqn.aux_vars, :, :, el)
    getGLSVars(eqn.params, view(eqn.q, :, :, el), aux_vars, dxidx_hat, view(mesh.jac, :, el),
               D, A_mats, qtranspose, qxitranspose, qxi, dxidx, taus)
#=
    if el == 1
      println("q = ", view(eqn.q, :, :, el))
      println("dxidx_hat = ", dxidx_hat)
      println("D = ", D)
      println("A_mats = ", A_mats)
      println("qxi = ", qxi)
      println("dxidx = ", dxidx)
    end
=#
    for i=1:numNodesPerElement
      res_i = view(eqn.res, :, i, el)
      for j=1:numNodesPerElement

        # zero  out some things
        fill!(qx, 0.0)
        fill!(tmps, 0.0)
        # trial space part
        # rotate qxi, qeta to qx, qy
        for d1=1:Tdim  # the x-y coordinate direction
          for d2=1:Tdim  # the xi-eta coordinate direction
            for n=1:numDofPerNode
              qx[n, d1] += dxidx_hat[d2, d1, j]*qxi[n, j, d2]
            end
          end
        end

#        println("\n  qx = \n", qx)

        # multiply qx, qy  flux jacobians Ax, Ay
        for d1=1:Tdim
          A_d1 = view(A_mats, :, :, d1, j)
          q_d1 = view(qx, :, d1)
          tmp_d1= view(tmps, :, d1)
          smallmatvec!(A_d1, q_d1, tmp_d1)
        end

#        println("\n  trial space terms = \n", tmps[:, 1:2])

        # now add them together
        tmp1 = view(tmps, :, 1)
        for d1=2:Tdim
          for n=1:numDofPerNode
            tmp1[n] += tmps[n, d1]
          end
        end

#        println("\n trial space term = \n", tmp1)
        # now multiply by tau
        tau_j = view(taus, :, :, j)
#        println("  \ntau = \n", tau_j) 
        tmp2 = view(tmps, :, 2)  # store result of multiplication here
        smallmatvec!(tau_j, tmp1, tmp2)  # this overwrites tmp2

#        println("\n  after multiplication by tau, reduction vector = \n", tmp2)
        # copy tmp2 into tmp1 so the next phase of reductions can proceed in
        # parallel
        # multiply by integration weight (from sbp.w) at the same time
        for n=1:numDofPerNode
          tmp2[n] *= w[j]
          tmp1[n] = tmp2[n]
        end

#        println("\n  after multiplication by w, reduction vector = \n", tmp2)

        # now do weighting space 
         # now multiply by flux jacobians transposed
         for d1=1:Tdim
           A_d1 = view(A_mats, :, :, d1, j)
           x_d1 = view(tmps, :, d1)
           b_d1 = view(tmps, :, d1+Tdim)  # now we use the second half of tmps
           smallmatTvec!(A_d1, x_d1, b_d1)
         end

#         println("\n  after multiplication by flux jacobians, reduction vector = \n", tmps[:, 3:4])

         # rotate the differentiation matrices D into x-y
        for d1=1:Tdim  # the x-y coordinate direction
#          println("d1 = ", d1)
          fac = zero(Tmsh)  # variable to accumulate the factor in
          for d2=1:Tdim  # the xi-eta coordinate direction
#              println("d2 = ", d2)
#              println("dxidx = ", dxidx[d2, d1, i])
#              println("d value = ", Dtranspose[i, j, d2])
              fac += dxidx[d2, d1, i]*Dtranspose[i, j, d2]
          end

#          println("fac = ", fac)

          # now update upper half of tmps with the factor
          for n=1:numDofPerNode
            tmps[n, d1 + Tdim] *= fac
          end
        end


#       println("\n  after multiplication by Dx, Dy, reduction vector = \n", tmps[:, 3:4])

       # now update res
       for d1=(Tdim+1):(2*Tdim)  # loop over the second half of tmps
         @simd for n=1:numDofPerNode
           res_i[n] += tmps[n, d1]
         end
       end




      end  # end loop over i
    end  # end loop over j

  end  # end loop over elements

  println("----- Finished applyGLS2 -----")
  return nothing

end  # end function

@doc """
### EulerEquationMod.getGLSVars

  This function calculates all the quantities needed to calculate the GLS
  stabilization term for an element.

  Although I only worked out the math for 2D, I am pretty sure this will
  work in 3D as well

  Inputs:
    params: ParamType, var_type can be anything
    q: array of solution variables for the element, 
       numDofPerNode x numNodesPerElement
    aux_vars: array of auxiliary variables for the element, 
              numAuxVars numNodesPerElement
    dxidx_hat: an array of the mapping jacobian scaled by 1/|J|, ie.
               (dxi/dx)/|J| for the entire element, 
               (Tdim x Tdim) x numNodesPerElement
    D : an array of the SBP differentiation matrices in the parametric 
        coordinate directions, (numNodesPerElement x numNodesPerElement) x Tdim

  Inputs/Outputs:
    A_mats: array the be populated with flux jacobian for the entire element,
            (numDofPerNode x NumDofPerNode) x Tdim x numNodesPerElement
    qtranspose: array to hold transpose(q) (temporarily)
    qxitranspose: array to hold [Dxi*transpose(q), Deta*transpose(q) ...]
                 numNodesPerElement x numDofPerNode x Tdim  (temporarily)
    qxi: transpose of qxtranspose (reversing the first two dimensions)
    dxidx: dxidx_hat * |J| (so the scaling by 1/|J| is undone)
    tau:  array of tau matrices for each node, 
          numDofPerNode x numDofPerNode x numNodesPerElement

            
"""->
function getGLSVars{Tmsh, Tsol, Tres, Tdim}(params::ParamType{Tdim},
                    q::AbstractArray{Tsol, 2}, aux_vars::AbstractArray{Tsol, 2},
                    dxidx_hat::AbstractArray{Tmsh, 3}, 
                    jac::AbstractArray{Tmsh, 1},
                    D::AbstractArray{Float64, 3},  
                    A_mats::AbstractArray{Tsol, 4}, # begin outputs
                    qtranspose::AbstractArray{Tsol, 2}, 
                    qxitranspose::AbstractArray{Tsol, 3}, 
                    qxi::AbstractArray{Tsol, 3},
                    dxidx::AbstractArray{Tmsh, 3}, 
                    tau::AbstractArray{Tres, 3})

#  println("----- Entered getGLSVars -----")

  numNodesPerElement = size(q, 2)
  numDofPerNode = size(q, 1)
  # move q into q_transpose
#  for i=1:numDofPerNode
#    for j=1:numNodesPerElement
#      qtranspose[j, i] = q[i, j]
#    end
#  end

  # multiply by Dxi, Deta ...
  for d=1:Tdim
    smallmatmatT!(view(D, :, :, d), q, view(qxitranspose, :, :, d))
  end

  # un-transpose while putting into qx
  for d=1:Tdim
    for j=1:numNodesPerElement
      for i=1:numDofPerNode
        qxi[i, j, d] = qxitranspose[j, i, d]
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
        dxidx[i, j, k] = dxidx_hat[i,j,k]*jac[k]
      end
    end
  end

  # get tau for each node
  for k=1:numNodesPerElement
    tau_k = view(tau, :, :, k)
    A_mat_k = view(A_mats, :, :, :, k)
    dxidx_k = view(dxidx_hat, :, :, k)
    getTau(params, A_mat_k, dxidx_k, tau_k)
  end

  return nothing
end  # end function

@doc """
### EulerEquationMod.getTau

  This function computes the tau matrix for a single node, for a steady problem
  as described in Hughes part X.

  Inputs
    params:  a ParamType, var_type can be anything
    A_mat: a matrix containing the flux jacobians at the node
           (numDofPerNode x numDofPerNode) x Tdim
    dxidx: a matrix the mapping jacobian for the node (dxi_i/dx_j), Tdim x Tdim
           This should *not* be scaled by 1/|J|

  Inputs/Outputs:
    tau: the numDofPerNode x numDofPerNode matrix to be populated


"""->
function getTau{Tdim, var_type, Tsol, Tres, Tmsh}(
                params::ParamType{Tdim, var_type, Tsol, Tres, Tmsh}, 
                A_mat::AbstractArray{Tsol, 3}, dxidx::AbstractArray{Tmsh, 2}, 
                tau::AbstractArray{Tres, 2})

#  println("----- Entered getTau -----")

  numDofPerNode = size(A_mat, 1)
  AjAk = params.A1
  flux_term = params.A2
  fill!(AjAk, 0.0)
  fill!(flux_term, 0.0)
 
  for k=1:Tdim
    for j=1:Tdim
      Aj = view(A_mat, :, :, j)
      Ak = view(A_mat, :, :, k)
      smallmatmat!(Aj, Ak, AjAk)  # Aj*Ak

      # calculate factor of dxidx*dxidx + dxidx*dxidx ...
      jacobian_fac = zero(Tmsh)
      for i=1:Tdim
        jacobian_fac += dxidx[i, j]*dxidx[i, k]
      end

      # accumulate dxidx*dxidx*Aj*Ak for all j, k in tmp2
      for p=1:numDofPerNode
        for q=1:numDofPerNode
          flux_term[p, q] += jacobian_fac*AjAk[p, q]
        end
      end

    end  # end loop over j
  end  # end loop over k

  # add a source term here


  # now take negative square root of the matrix
  # there are more efficient ways of doing this
  D, V = eig(flux_term)
  # check that D contains only positive numbers
  # this should be inside a @debug1
  for i=1:numDofPerNode
    if real(D[i]) < 0.0
      println("warning, D[$i] = ", D[i])
    end
  end

  for i=1:numDofPerNode
    D[i] = D[i]^(-0.5)
  end
  
  # reconstruct M = (V*D)*inv(V)
  Vinv = inv(V)
  for i=1:numDofPerNode
    entry = D[i]
    for j=1:numDofPerNode
      V[j, i] *= entry
    end
  end

  smallmatmat!(V, Vinv, tau)

#  println("----- Finished getTau -----")
  return nothing
end  # end function
