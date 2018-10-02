# this function applies the form of GLS stabilization described in:
# A new finite element formulation for computational fluid dynamics:
# X The compressible Euler and Navier-Stokes Equations
# by Shakib, Hughes, and Johan

# this *should* work for both conservative and entropy variables, I think
function applyGLS2(mesh::AbstractMesh{Tmsh}, 
sbp::AbstractSBP, eqn::EulerData{Tsol, Tres, Tdim}, opts) where {Tmsh, Tsol, Tres, Tdim}

#  println("----- Entered applyGLS2 -----")
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

  # DEBUGGING

  gls_res = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerElement, 
                  mesh.numEl)

  gls_full = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerElement, 
                  mesh.numEl)
  tau_eval_sum = 0.0
  tau_eval_cnt = 0
  res_before = copy(eqn.res)

  # calculate D
  for d=1:Tdim
    smallmatmat!(diagm(1./sbp.w), sview(sbp.Q, :, :, d), sview(D, :, :, d))
    Dtranspose[:, :, d] = D[:, :, d].'
  end

#   for el =1:1  #DEBUGGING 
  for el = 1:mesh.numEl
    # get all the quantities for this element
    dxidx_hat = sview(mesh.dxidx, :, :, :, el)
    jac = sview(mesh.jac, :, el)
    aux_vars = sview(eqn.aux_vars, :, :, el)
    getGLSVars(eqn.params, eqn.params_conservative, sview(eqn.q, :, :, el), aux_vars, dxidx_hat, sview(mesh.jac, :, el), D, A_mats, qtranspose, qxitranspose, qxi, dxidx, taus)
#=
    if el == 1
      println("q = ", sview(eqn.q, :, :, el))
      println("dxidx_hat = ", dxidx_hat)
      println("D = ", D)
      println("A_mats = ", A_mats)
      println("qxi = ", qxi)
      println("dxidx = ", dxidx)
    end
=#
    for i=1:numNodesPerElement
      res_i = sview(eqn.res, :, i, el)
      #DEBUGGING
      gls_res_i = sview(gls_res, :, i, el)
      gls_full_i = sview(gls_full, :, i, el)
      for j=1:numNodesPerElement

        # zero  out some things
        fill!(qx, 0.0)
        fill!(tmps, 0.0)
        # trial space part
        # rotate qxi, qeta to qx, qy
        for d1=1:Tdim  # the x-y coordinate direction
          for d2=1:Tdim  # the xi-eta coordinate direction
            for n=1:numDofPerNode
              qx[n, d1] += dxidx[d2, d1, j]*qxi[n, j, d2]
            end
          end
        end

#        println("\n  qx = \n", qx)

        # multiply qx, qy  flux jacobians Ax, Ay
        for d1=1:Tdim
          A_d1 = sview(A_mats, :, :, d1, j)
          q_d1 = sview(qx, :, d1)
          tmp_d1= sview(tmps, :, d1)
          smallmatvec!(A_d1, q_d1, tmp_d1)
        end

#        println("\n  trial space terms = \n", tmps[:, 1:2])

        # now add them together
        tmp1 = sview(tmps, :, 1)
        for d1=2:Tdim
          for n=1:numDofPerNode
            tmp1[n] += tmps[n, d1]
          end
        end
        
        #DEBUGGING
       
        for n=1:numDofPerNode
          gls_res_i[n] += tmp1[n]
        end
        
        
#        println("\n trial space term = \n", tmp1)
        # now multiply by tau
        tau_j = sview(taus, :, :, j)
#        println("  \ntau = \n", tau_j) 
        tmp2 = sview(tmps, :, 2)  # store result of multiplication here
        smallmatvec!(tau_j, tmp1, tmp2)  # this overwrites tmp2

#        println("\n  after multiplication by tau, reduction vector = \n", tmp2)
        # copy tmp2 into tmp1 so the next phase of reductions can proceed in
        # parallel
        # multiply by integration weight (from sbp.w) at the same time
        # the indices should be i not j ???
        for n=1:numDofPerNode
          tmp2[n] *= w[i]/jac[i]
          tmp1[n] = tmp2[n]
        end

#        println("\n  after multiplication by w, reduction vector = \n", tmp2)

        # now do weighting space 
         # now multiply by flux jacobians transposed
         for d1=1:Tdim
           A_d1 = sview(A_mats, :, :, d1, j)
           x_d1 = sview(tmps, :, d1)
           b_d1 = sview(tmps, :, d1+Tdim)  # now we use the second half of tmps
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
           gls_full_i[n] -= tmps[n, d1]  # DEBUGGING
           res_i[n] -= tmps[n, d1]
         end
       end




      end  # end loop over j


      # DEBUGGING
      Dvals, V = eig(sview(taus, :, :, i))
      tau_eval_sum += real(Dvals[1])
      tau_eval_cnt += 1

    end  # end loop over i

  end  # end loop over elements

  #DEBUGGING
 #= 
  gls_resvec = zeros(Tsol, mesh.numDof)
  full_resvec = zeros(Tsol, mesh.numDof)
  old_resvec = zeros(Tsol, mesh.numDof)
  array3DTo1D(mesh, sbp, eqn, opts, gls_res, gls_resvec)
  array3DTo1D(mesh, sbp, eqn, opts, gls_full, full_resvec)
  writedlm("gls_full.dat", real(gls_full))
  writedlm("gls_fullvec.dat", full_resvec)
  array3DTo1D(mesh, sbp, eqn, opts, res_before, old_resvec)
  gls_norm = calcNorm(eqn, gls_resvec)
  full_norm = calcNorm(eqn, full_resvec)
  old_norm = calcNorm(eqn, old_resvec)
  tau_eval_avg = tau_eval_sum/tau_eval_cnt
  rmfile("gls_norm.dat")
#  println("printing gls_norm.dat")
  f = open("gls_norm.dat", "w")
  println(f, gls_norm, " ", old_norm, " ", full_norm, " ", tau_eval_avg)
  close(f)
#  println("----- Finished applyGLS2 -----")
 =#
  return nothing

end  # end function

@doc """
### EulerEquationMod.getGLSVars

  This function calculates all the quantities needed to calculate the GLS
  stabilization term for an element.

  Although I only worked out the math for 2D, I am pretty sure this will
  work in 3D as well

  Inputs:
    params: ParamType, var_type can be anythinga
    params_c: ParamType, var_type must be :conservative
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
function getGLSVars(params::ParamType{Tdim},
params_c::ParamType{Tdim, :conservative},
q::AbstractArray{Tsol, 2}, aux_vars::AbstractArray{Tsol, 2},
dxidx_hat::AbstractArray{Tmsh, 3}, 
jac::AbstractArray{Tmsh, 1},
D::AbstractArray{Float64, 3},  
A_mats::AbstractArray{Tsol, 4}, # begin outputs
qtranspose::AbstractArray{Tsol, 2}, 
qxitranspose::AbstractArray{Tsol, 3}, 
qxi::AbstractArray{Tsol, 3},
dxidx::AbstractArray{Tmsh, 3}, 
tau::AbstractArray{Tres, 3}) where {Tmsh, Tsol, Tres, Tdim}


  # who cares about performance?
  tau_type = params.tau_type



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
    smallmatmatT!(sview(D, :, :, d), q, sview(qxitranspose, :, :, d))
  end

  # un-transpose while putting into qx
  for d=1:Tdim
    for j=1:numNodesPerElement
      for i=1:numDofPerNode
        qxi[i, j, d] = qxitranspose[j, i, d]
      end
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

 
  # get conservative variable flux jacobian 
  if tau_type == 1
#    println("getting conservative variable flux jacobians")
    q_c = params_c.q_vals
    fill!(q_c, 0.0)
    for k=1:numNodesPerElement
      q_k = sview(q, :, k)
      convertToConservative(params, q_k, q_c)
      A1_k = sview(A_mats, :, :, 1, k)
      calcA1(params_c, q_c, A1_k)
      A2_k = sview(A_mats, :, :, 2, k)
      calcA2(params_c, q_c, A2_k)

      if Tdim == 3  # three cheers for static analysis
        A3_k = sview(A_mats, :, :, 3, k)
        calcA3(params_c, q_c, A3_k)
      end
    end
  else
#    println("getting entropy variable flux jacobians")
    # get entropy variable  flux jacobian 
    for k=1:numNodesPerElement
      q_k = sview(q, :, k)
#      println("q_k code = ", q_k)
      A1_k = sview(A_mats, :, :, 1, k)
      calcA1(params, q_k, A1_k)
#      println("A1 code = \n", A1_k)
      A2_k = sview(A_mats, :, :, 2, k)
      calcA2(params, q_k, A2_k)
#      println("A2 code = \n", A2_k)

      if Tdim == 3  # three cheers for static analysis
        A3_k = sview(A_mats, :, :, 3, k)
        calcA3(params, q_k, A3_k)
      end
    end
  end




  # get tau for each node - using conservative variable flux jacobian
  for k=1:numNodesPerElement
    q_k = sview(q, :, k)
    tau_k = sview(tau, :, :, k)
    A_mat_k = sview(A_mats, :, :, :, k)
    dxidx_k = sview(dxidx, :, :, k)

    # dance branch prediction, dance
    if tau_type == 1
     getTau(params, q_k, A_mat_k, dxidx_k, tau_k)
    elseif tau_type == 2
      getTau(params, jac[k], tau_k)
    elseif tau_type == 3
      p = 1
#      println("\n arguments passed to getTau from applyGLS2:")
#      println("q_k = ", q_k)
#      println("A_mat_k = \n", A_mat_k)
#      println("dxidx_k = \n", dxidx_k)
#      println("p_val = ", p)
#      println("tau_k = \n", tau_k)
      getTau(params, q_k, A_mat_k, dxidx_k, p, tau_k)
    else
      println(STDERR, "Warning: unsupported Tau requested for GLS")
    end

  end

  if tau_type == 1
    # get entropy variable  flux jacobian 
    for k=1:numNodesPerElement
      q_k = sview(q, :, k)
      A1_k = sview(A_mats, :, :, 1, k)
      calcA1(params, q_k, A1_k)
      A2_k = sview(A_mats, :, :, 2, k)
      calcA2(params, q_k, A2_k)

      if Tdim == 3  # three cheers for static analysis
        A3_k = sview(A_mats, :, :, 3, k)
        calcA3(params, q_k, A3_k)
      end
    end
  end

  return nothing
end  # end function

#=
function getGLSVars(params::ParamType{Tdim},
                    params_c::ParamType{Tdim, :conservative},
                    q::AbstractArray{Tsol, 2}, aux_vars::AbstractArray{Tsol, 2},
                    dxidx_hat::AbstractArray{Tmsh, 3}, 
                    jac::AbstractArray{Tmsh, 1},
                    D::AbstractArray{Float64, 3},  
                    A_mats::AbstractArray{Tsol, 4}, # begin outputs
                    qtranspose::AbstractArray{Tsol, 2}, 
                    qxitranspose::AbstractArray{Tsol, 3}, 
                    qxi::AbstractArray{Tsol, 3},
                    dxidx::AbstractArray{Tmsh, 3}, 
                    tau::AbstractArray{Tres, 3}) where {Tmsh, Tsol, Tres, Tdim}

=#

function applyGLS3(mesh::AbstractMesh{Tmsh}, sbp, eqn::EulerData{Tsol, Tres, Tdim}, opts) where {Tsol, Tres, Tdim, Tmsh}

  numDofPerNode = mesh.numDofPerNode
  numNodesPerElement = mesh.numNodesPerElement
  w = sbp.w
  # reusable storage

  # flux jacobians: a n x n matrix for each coordinate direction for each node
  A_mats = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode, Tdim, 
                      mesh.numNodesPerElement)
#  dxidx_hat = zeros(Tmsh, Tdim, Tdim, mesh.numNodesPerElement)

  # hold non-scaled mapping jacobian
  dxidx = zeros(Tmsh, Tdim, Tdim, numNodesPerElement)
  # hold tau, a n x n matrix for each node
  taus = zeros(Tres, mesh.numDofPerNode, mesh.numDofPerNode, 
                    mesh.numNodesPerElement)

  # hold differentiation matrices Dxi (parametric directions) and Dx 
  # (physical directions)
  Dxi = zeros(Float64, numNodesPerElement, numNodesPerElement, Tdim)
  Dx = zeros(Dxi)

  # temporary vectors 1 through 4, stored in a matrix so they can be
  # access programaticaly
  tmps = zeros(Tres, mesh.numDofPerNode, 2*Tdim)
  trial_term = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement)
  middle_term = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement)
  complete_term = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement)

  # calculate Dxi
  # calculate D
  for d=1:Tdim
    smallmatmat!(diagm(1./sbp.w), sview(sbp.Q, :, :, d), sview(Dxi, :, :, d))
  end

  # DEBUGGING
  gls_res = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerElement, 
                  mesh.numEl)
  middle_terms = zeros(gls_res)
  gls_full = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerElement, 
                  mesh.numEl)
  tau_eval_sum = 0.0
  tau_eval_cnt = 0
  res_before = copy(eqn.res)




#   for el = 1:1  # DEBUGGING
  for el = 1:mesh.numEl

    q_el = sview(eqn.q, :, :, el)
    aux_vars_el = sview(eqn.aux_vars, :, :, el)
    dxidx_hat_el = sview(mesh.dxidx, :, :, :, el)
    jac_el = sview(mesh.jac, :, el)
   
    fill!(Dx, 0.0)
    fill!(taus, 0.0)
    # get all the variables needed for this element
    getGLSVars3(eqn.params, eqn.params_conservative, q_el, aux_vars_el, 
                dxidx_hat_el, jac_el, Dxi, Dx, A_mats, dxidx, taus)

    fill!(trial_term, 0.0)
    fill!(middle_term, 0.0)
    fill!(complete_term, 0.0)

    for i=1:mesh.numNodesPerElement  # row index  (free index)
      gls_res_debug_i = sview(gls_res, :, i, el)

      # variable to accumulate the trial term at node i
      trial_term_i = sview(trial_term, :, i)
      for j=1:mesh.numNodesPerElement  # column index  (summed index)
        fill!(tmps, 0.0)
        q_j = sview(eqn.q, :, j, el)

        ### Trial space term ###

        # multiply q_j by entries of the Dx and store in the lower half of
        # tmps
        for d1 = 1:Tdim
          tmp_d1 = sview(tmps, :, d1)
          d_val = Dx[i, j, d1]
          for n=1:numDofPerNode
            tmp_d1[n] = d_val*q_j[n]
          end
        end

        # now multiply by flux jacobian evaluated at i, storing in upper
        # half of tmps
        for d1 = 1:Tdim
          A_i = sview(A_mats, :, :, d1, i)
          tmp_1 = sview(tmps, :, d1)
          tmp_2 = sview(tmps, :, d1 + Tdim)
          smallmatvec!(A_i, tmp_1, tmp_2)
        end

        # accumulate into trial term vector
        # Perf Note: do parallel reads of tmps vectors
        # is there an AVX instruction for scaling by a constant?
        for n=1:numDofPerNode
          for d1 = (Tdim+1):2*Tdim
            trial_term_i[n] += tmps[n, d1]
            #DEBUGGING
            gls_res_debug_i[n] += tmps[n, d1]
          end
        end

      end  # end sum over j in the trial term
    end  # end loop over i for the trial term

#    println("trial_term = \n", trial_term)

    ### now do middle terms ###

    for i=1:mesh.numNodesPerElement
      # multiply by integration weight, jacobian factor
      trial_term_i = sview(trial_term, :, i)
      middle_term_i = sview(middle_term, :, i)
      middle_terms_debug_i = sview(middle_terms, :, i, el)

      fac = w[i]/jac_el[i]
      for n=1:numDofPerNode
        trial_term_i[n] *= fac
      end

      # multiply by tau
      tau_i = sview(taus, :, :, i)
      smallmatvec!(tau_i, trial_term_i, middle_term_i)

      # DEBUGGING
      for n=1:numDofPerNode
        middle_terms_debug_i[n] += middle_term[n]
      end
    
    end  # end loop over i for middle terms

#    println("middle_term = ", middle_term)

    ### now do weighting space term ###

    for i=1:mesh.numNodesPerElement  # free index
      complete_term_i = sview(complete_term, :, i)
      res_i = sview(eqn.res, :, i, el)
      gls_full_debug_i = sview(gls_full, :, i, el)

      for j=1:mesh.numNodesPerElement  # summed index
        fill!(tmps, 0.0)
        middle_term_j = sview(middle_term, :, j)


        # multiply by transposed flux jacobian, store in lower half of tmps
        for d1 = 1:Tdim
          A_d1 = sview(A_mats, :, :, d1, j)
          tmp_d1 = sview(tmps, :, d1)
          smallmatTvec!(A_d1, middle_term_j, tmp_d1)
        end

        # multiply by entries of the Dx
        for d1 = 1:Tdim
          d_d1 = Dx[j, i, d1]
          tmp_d1 = sview(tmps, :, d1)
          for n=1:numDofPerNode
            tmp_d1[n] *= d_d1
          end
        end

        # sum into complete_term
        for n=1:numDofPerNode
          for d1=1:Tdim
            tmp_d1 = sview(tmps, :, d1)
            complete_term_i[n] += tmp_d1[n]
          end
        end

      end  # end loop over j

      # update res
      for n=1:numDofPerNode
        res_i[n] -= complete_term_i[n]
        #DEBUGGING
        gls_full_debug_i[n] -= complete_term_i[n]
      end


    end  # end loop over i

#    println("complete_term = \n", complete_term)
#    println("res = \n", eqn.res[:, :, el])

  end  # end loop over elements

#=
  println("element 1 gls_res = \n", gls_res[:, :, 1])
  println("element 1 middle_terms = \n", middle_terms[:, :, 1])
  println("element 1 gls_full = \n", gls_full[:, :, 1])
=#
  #=
  gls_res_vec = zeros(Tres, mesh.numDof)
  middle_vec = zeros(gls_res_vec)
  gls_fullvec = zeros(gls_res_vec)
  array3DTo1D(mesh, sbp, eqn, opts, gls_res_vec
  =#

  return nothing

end


function getGLSVars3(params::ParamType{Tdim}, 
params_c::ParamType{Tdim, :conservative}, 
q::AbstractArray{Tsol, 2}, 
aux_vars::AbstractArray{Tsol, 2}, 
dxidx_hat::AbstractArray{Tmsh, 3}, 
jac::AbstractArray{Tmsh, 1}, D::AbstractArray{Float64, 3},
# begin output parameters
Dx::AbstractArray{Tmsh, 3}, A_mats::AbstractArray{Tsol, 4},
dxidx::AbstractArray{Tmsh,3}, tau::AbstractArray{Tres, 3}) where {Tmsh, Tsol, Tres, Tdim}


  numDofPerNode = size(q, 1)
  numNodesPerElement = size(q, 2)

   # get dxidx - not scaled by 1/|J|
  for k = 1:numNodesPerElement
    for j=1:Tdim
      for i=1:Tdim
        dxidx[i, j, k] = dxidx_hat[i,j,k]*jac[k]
      end
    end
  end

  # calculate Dx, Dy, Dz
  for d_phys = 1:Tdim  # loop over x, y, z directions
    D_phys = sview(Dx, :, :, d_phys)
    for d_param = 1:Tdim  # loop for xi, eta, ...
      D_param = sview(D, :, :, d_param)
      for j=1:numNodesPerElement
        for i=1:numNodesPerElement
          fac = dxidx[d_param, d_phys, i]  # get jacobian term
          D_phys[i, j] += fac*D_param[i, j]
        end
      end
    end
  end

  tau_type = params.tau_type
  if tau_type == 1
#    println("getting conservative variable flux jacobians")
    q_c = params_c.q_vals
    fill!(q_c, 0.0)
    for k=1:numNodesPerElement
      q_k = sview(q, :, k)
      convertToConservative(params, q_k, q_c)
      A1_k = sview(A_mats, :, :, 1, k)
      calcA1(params_c, q_c, A1_k)
      A2_k = sview(A_mats, :, :, 2, k)
      calcA2(params_c, q_c, A2_k)

      if Tdim == 3  # three cheers for static analysis
        A3_k = sview(A_mats, :, :, 3, k)
        calcA3(params_c, q_c, A3_k)
      end
    end
  else
#    println("getting entropy variable flux jacobians")
    # get entropy variable  flux jacobian 
    for k=1:numNodesPerElement
      q_k = sview(q, :, k)
#      println("q_k code = ", q_k)
      A1_k = sview(A_mats, :, :, 1, k)
      calcA1(params, q_k, A1_k)
#      println("A1 code = \n", A1_k)
      A2_k = sview(A_mats, :, :, 2, k)
      calcA2(params, q_k, A2_k)
#      println("A2 code = \n", A2_k)

      if Tdim == 3  # three cheers for static analysis
        A3_k = sview(A_mats, :, :, 3, k)
        calcA3(params, q_k, A3_k)
      end
    end
  end

  # get tau for each node - using conservative variable flux jacobian
  for k=1:numNodesPerElement
    q_k = sview(q, :, k)
    tau_k = sview(tau, :, :, k)
    A_mat_k = sview(A_mats, :, :, :, k)
    dxidx_k = sview(dxidx, :, :, k)

    # dance branch prediction, dance
    if tau_type == 1
     getTau(params, q_k, A_mat_k, dxidx_k, tau_k)
    elseif tau_type == 2
      getTau(params, jac[k], tau_k)
    elseif tau_type == 3
      p = 1
#      println("\n arguments passed to getTau from applyGLS2:")
#      println("q_k = ", q_k)
#      println("A_mat_k = \n", A_mat_k)
#      println("dxidx_k = \n", dxidx_k)
#      println("p_val = ", p)
#      println("tau_k = \n", tau_k)
      getTau(params, q_k, A_mat_k, dxidx_k, p, tau_k)
    else
      println(STDERR, "Warning: unsupported Tau requested for GLS")
    end

  end

  if tau_type == 1
    # get entropy variable  flux jacobian 
    for k=1:numNodesPerElement
      q_k = sview(q, :, k)
      A1_k = sview(A_mats, :, :, 1, k)
      calcA1(params, q_k, A1_k)
      A2_k = sview(A_mats, :, :, 2, k)
      calcA2(params, q_k, A2_k)

      if Tdim == 3  # three cheers for static analysis
        A3_k = sview(A_mats, :, :, 3, k)
        calcA3(params, q_k, A3_k)
      end
    end
  end



  return nothing

end  # end function


function test_GLS(mesh::AbstractMesh{Tmsh}, sbp, eqn::AbstractSolutionData{Tsol, Tres}, opts) where {Tsol, Tres, Tmsh}

  eqn.params.tau_type = 2
  Dxi = diagm(1./sbp.w)*sbp.Q[:, :, 1]
  Deta = diagm(1./sbp.w)*sbp.Q[:, :, 2]

  # create indices
  idx_range = Array{UnitRange{Int64}}(mesh.numNodesPerElement)
  for i=1:mesh.numNodesPerElement
    start_idx = (i-1)*mesh.numDofPerNode + 1
    end_idx = i*mesh.numDofPerNode
    idx_range[i] = copy(start_idx:end_idx)
  end

#    println("idx_range = ", idx_range)

    size_block = mesh.numNodesPerElement*mesh.numDofPerNode
#    println("size_block = ", size_block)


  # testing: only do one element
  for el =1:mesh.numEl
#    println("testing element ", el)
    res_el = sview(eqn.res, :, :, el)

    # constant mapping elements only
    dxidx = zeros(2,2)
    dxidx_hat_el = sview(mesh.dxidx, :, :, 1, el)
    jac_el = sview(mesh.jac, :, el)
    q_el = reshape(copy(eqn.q[:, :, el]), size_block)

    for i=1:2
      for j=1:2
        dxidx[i,j] = dxidx_hat_el[i, j, 1]*jac_el[1]
      end
    end

    # calculate Dx, Dy
    Dx = dxidx[1,1]*Dxi + dxidx[2, 1]*Deta
    Dy = dxidx[1,2]*Dxi + dxidx[2, 2]*Deta

#    println("Dx = \n", Dx)
#    println("Dy = \n", Dy)

    # create block Dx, Dy
    Dx_tilde = zeros(Tmsh, size_block, size_block)
    Dy_tilde = zeros(Dx_tilde)

   
    for i=1:mesh.numNodesPerElement
      idx_i = idx_range[i]
#      println("idx_i = ", idx_i)
      for j=1:mesh.numNodesPerElement
        idx_j = idx_range[j]
#        println("  idx_j = ", idx_j)
        Dx_tilde[idx_i, idx_j] = Dx[i, j]*eye(mesh.numDofPerNode)
        Dy_tilde[idx_i, idx_j] = Dy[i, j]*eye(mesh.numDofPerNode)
      end
    end

#    println("Dx_tilde = \n", Dx_tilde)
#    println("Dy_tilde = \n", Dy_tilde)

    # create A1 tilde and A2 tilde
    A1_tilde = zeros(Tsol, size_block, size_block)
    A2_tilde = zeros(Tsol, size_block, size_block)

    for i=1:mesh.numNodesPerElement
      idx_i = idx_range[i]
      q_i = q_el[idx_i]

      A1 = sview(A1_tilde, idx_i, idx_i)
      EulerEquationMod.calcA1(eqn.params, q_i, A1)

      A2 = sview(A2_tilde, idx_i, idx_i)
      EulerEquationMod.calcA2(eqn.params, q_i, A2)
    end

#    println("A1 = \n", A1_tilde)
#    println("A2 = \n", A2_tilde)

    # create middle terms, including tau
    middle_tilde = zeros(Tres, size_block, size_block)

    for i=1:mesh.numNodesPerElement

      idx_i = idx_range[i]
      tau = zeros(Tres, mesh.numDofPerNode, mesh.numDofPerNode)
      EulerEquationMod.getTau(eqn.params, jac_el[i], tau)  # tau number 2

      middle_tilde[idx_i, idx_i] = (sbp.w[i]/jac_el[i])*tau
    end


    # create the operator
    fancy_L  = A1_tilde*Dx_tilde + A2_tilde*Dy_tilde
    gls_operator = fancy_L.'*middle_tilde*fancy_L
    
    @assert isSymmetric(gls_operator, 1e-12)
    #println("max asymmetry = ", maximum(abs(gls_operator - gls_operator.')))

    gls_test = -gls_operator*q_el

    for i=1:size_block
      res_el[i] += gls_test[i]
    end
  end  # end loop over elements

  return nothing

end  # end function





@doc """
### EulerEquationMod.getTau

  This function computes the tau matrix for a single node, for a steady problem
  as described in Hughes part X.

  Inputs
    params:  a ParamType, var_type can be anything
    A_mat: a matrix containing the flux jacobians of the conservative variables 
           at the node
           (numDofPerNode x numDofPerNode) x Tdim
    dxidx: a matrix the mapping jacobian for the node (dxi_i/dx_j), Tdim x Tdim
           This should *not* be scaled by 1/|J|

  Inputs/Outputs:
    tau: the numDofPerNode x numDofPerNode matrix to be populated


"""->
function getTau(params::ParamType{Tdim, var_type, Tsol, Tres, Tmsh}, 
                q::AbstractVector{Tsol}, A_mat::AbstractArray{Tsol, 3}, 
                dxidx::AbstractArray{Tmsh, 2}, tau::AbstractArray{Tres, 2}) where {Tdim, var_type, Tsol, Tres, Tmsh}

#  println("----- Entered getTau original-----")

  numDofPerNode = size(A_mat, 1)

  data = params.get_tau_data
  @unpack data AjAk flux_term A0inv tmp_mat

  fill!(AjAk, 0.0)  # unneeded?
  fill!(flux_term, 0.0)
 
  for k=1:Tdim
    for j=1:Tdim
      Aj = sview(A_mat, :, :, j)
      Ak = sview(A_mat, :, :, k)
      smallmatmat!(Aj, Ak, AjAk)  # Aj*Ak

      # calculate factor of dxidx*dxidx + dxidx*dxidx ...
      jacobian_fac = zero(Tmsh)
      for i=1:Tdim
        jacobian_fac += dxidx[i, j]*dxidx[i, k]
      end

      # accumulate dxidx*dxidx*Aj*Ak for all j, k in tmp2
      for p=1:numDofPerNode
        for r=1:numDofPerNode
          flux_term[p, r] += jacobian_fac*AjAk[p, r]
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

  smallmatmat!(V, Vinv, tmp_mat)

  calcA0Inv(params, q, A0inv)
  smallmatmat!(A0inv, tmp_mat, tau)
#  println("----- Finished getTau -----")
  return nothing
end  # end function

function getTau(params::ParamType, jac::Number, tau::AbstractArray{Tres, 2}) where Tres

  fac = 1.0
  for i=1:size(tau, 1)
    tau[i,i] = fac*1/(jac^(1/2))
  end


end

# not sure if this works with conservative variables
# A_mats must be the entropy variable flux jacobians
function getTau(params::ParamType{Tdim, :entropy}, 
q::AbstractVector{Tsol}, A_mat::AbstractArray{Tsol, 3}, 
dxidx::AbstractArray{Tmsh, 2}, p::Integer, tau::AbstractArray{Tres, 2}) where {Tsol, Tres, Tmsh, Tdim}

  data = params.get_tau_data
  @unpack data B_d B_p A_mat_hat A0 tmp_mat2
  fill!(B_p, 0.0)
  calcA0(params, q, A0)
  L = chol(A0)' # is the hermitian transpose right here?

  Linv = full(inv(L))  # convert to full matrices because I think it is more
                       # efficient for small matrices

  # calculate the A_i hats = inv(L)*flux_jacobian_i*inv(L).'
  for d1=1:Tdim
    A_hat_d1 = sview(A_mat_hat, :, :, d1)
    A_d1 = sview(A_mat, :, :, d1)
    smallmatmat!(Linv, A_d1, tmp_mat2)
    smallmatmatT!(tmp_mat2, Linv, A_hat_d1)
  end

#  println("A_mat_hat = \n", A_mat_hat)

  for d1=1:Tdim  # loop over B_d
    fill!(B_d, 0.0)

    for d2=1:Tdim  # summed index dxi_d1/dx_d2 A_d2
      A_hat_d2 = sview(A_mat_hat, :, :, d2)
      dxidx_d2 = dxidx[d1, d2]
      # accumulate into B_d
      for i=1:size(B_d, 1)
        for j=1:size(B_d, 2)
          B_d[i, j] += dxidx_d2*A_hat_d2[i, j]
        end
      end

    end  # end loop d2

 #   println("before make_symmetric, symmetry norm = ", vecnorm(B_d - B_d.'))
    make_symmetric!(B_d)  # make sure it is symmetric so we get real eigenvalues
                          # when using real variables

#    println("B$d1 = \n", B_d)
    D, V = eig(B_d)
#    println("D$d1 = \n", D)
#    println("V$d1 = \n", V)
#    println("p = ", p)
 #   println("before update, B_p = \n", B_p)
    # take absolute value, raise to  power p while accumulating into B_p
    for k=1:length(D)  # for each eigenvalue, do  outer product
      v_k = sview(V, :, k)
      val_k = absvalue(D[k])^p
      for i=1:size(B_p, 1)
        for j=1:size(B_p, 2)
          B_p[i, j] += val_k*v_k[i]*v_k[j]
        end
      end
    end

#    println("after adding b$d1, B_p = \n", B_p)

  end  # end loop d1

#  println("B_p = \n", B_p)
  D2, V2 = eig(B_p)

  # now calculate tau: invert and take the pth root of D_p at the same time
  fac = 1.0
  fill!(tau, 0.0)
  for k=1:length(D2)
    v_k = sview(V2, :, k)
    val_k = D2[k]^(-1/p)
    for i=1:size(B_p, 1)
      for j=1:size(B_p, 2)
        tau[i, j] += fac*val_k*v_k[i]*v_k[j]
      end
    end
  end

#  println("tau = \n", tau)

#  println("----- finished getTau -----")
  return nothing
end
 




