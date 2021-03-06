# differentiated version of functions in euler_funcs.jl

"""
  Applies the element level inverse mass matrix to an element-level jacobian
  contribution

  **Inputs**

   * jac_el: vector, length `numNodesPerElement`, containing the determinant
              of the mapping jacobian at each node of the element
   * weights: the SBP integration weights for each node of the element,
               length `numNodesPerElement

  **Inputs/Outputs**

   * res_jac: array containing element level Jacobian, to be multiplied by
              the inverse mass matrix, size `numDofPerNode` x `numDofPerNode`
              x `numNodesPerElement` x `numNodesPerElement`
"""
function applyMinvElement(jac_el::AbstractVector, weights::AbstractVector,
                          res_jac::AbstractArray{T, 4}) where {T}

  # multiply by Minv if needed
  for q=1:size(res_jac, 4)
    for p=1:size(res_jac, 3)
      val = jac_el[p]/weights[p]  # entry in Minv
      @simd for m=1:size(res_jac, 2)
        @simd for n=1:size(res_jac, 1)
          res_jac[n, m, p, q] *= val
        end
      end
    end
  end

  return nothing
end


"""
  Computes the derivataive of the volume term of the Roe scheme with respect
  to `q`, ie

  d/dq (Q^T * f)

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * assembler: used to assemble the contribution of each element into the
                Jacobian
"""
function calcVolumeIntegrals_nopre_diff(
                                   mesh::AbstractMesh{Tmsh},
                                   sbp::AbstractOperator,
                                   eqn::EulerData{Tsol, Tres, Tdim},
                                   opts,
                                   assembler::AssembleElementData) where {Tmsh, Tsol, Tres, Tdim}



  data = eqn.params.calc_volume_integrals_data
  @unpack data flux_jac res_jac nrm

  for i=1:mesh.numEl
    fill!(flux_jac, 0)
    fill!(res_jac, 0)
    for j=1:mesh.numNodesPerElement
      q_j = ro_sview(eqn.q, :, j, i)
      aux_vars_j = ro_sview(eqn.aux_vars, :, j, i)

      # compute dF/dq
      for k=1:Tdim
        fluxjac_k = sview(flux_jac, :, :, j, k)

        # get the direction vector
        for p=1:Tdim
          nrm[p] = mesh.dxidx[k, p, j, i]
        end
        calcEulerFlux_diff(eqn.params, q_j, aux_vars_j, nrm, fluxjac_k)
      end  # end loop k
    end  # end loop j

    # compute dR/dq
    for k=1:Tdim
      weakDifferentiateElement_jac!(sbp, k, sview(flux_jac, :, :, :, k), res_jac, SummationByParts.Add(), true)
    end

    
    if eqn.params.use_Minv == 1
      jac_el = sview(mesh.jac, :, i)
      applyMinvElement(jac_el, sbp.w, res_jac)
    end

    # assemble element level jacobian into the residual
    assembleElement(assembler, mesh, i, res_jac)
#    fill!(res_jac, 0.0)
    # flux_jac gets overwritten, so no need to zero it 

  end  # end loop i

  return nothing
end  # end function


"""
  Reverse mode wrt metrics of [`calcVolumeIntegrals_nopre`](@ref)
"""
function calcVolumeIntegrals_nopre_revm(
           mesh::AbstractMesh{Tmsh},
           sbp::AbstractOperator,
           eqn::EulerData{Tsol, Tres, Tdim},
           opts) where {Tmsh, Tsol, Tres, Tdim}


  # flux in the parametric directions for a given element
  flux_el_bar = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, Tdim)
  nrm = zeros(Tmsh, mesh.dim)
  nrm_bar = zeros(Tmsh, mesh.dim)
#=
  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      q_j = ro_sview(eqn.q, :, j, i)
      aux_vars_j = ro_sview(eqn.aux_vars, :, j, i)

      for k=1:Tdim
        flux_k = sview(flux_el, :, j, k)

        # get the direction vector
        for p=1:Tdim
          nrm[p] = mesh.dxidx[k, p, j, i]
        end
        # consider calculating all directions at once
        # not sure if that will help because the data dependencies are
        # really simple
        calcEulerFlux(eqn.params, q_j, aux_vars_j, nrm, flux_k)
      end  # end loop k
    end  # end loop j

    res_i = sview(eqn.res, :, :, i)
    for k=1:Tdim
      weakDifferentiateElement!(sbp, k, sview(flux_el, :, :, k), res_i, SummationByParts.Add(), true)
    end
=#
  #-----------------------------
  # reverse sweep
  for i=1:mesh.numEl
    res_bar_i = sview(eqn.res_bar, :, :, i)
    fill!(flux_el_bar, 0)
    for k=1:Tdim
      flux_bar_k = sview(flux_el_bar, :, :, k)
      weakDifferentiateElement_rev!(sbp, k, flux_bar_k, res_bar_i,
                                    SummationByParts.Add(), true)
    end

    for j=1:mesh.numNodesPerElement
      q_j = ro_sview(eqn.q, :, j, i)
      aux_vars_j = ro_sview(eqn.aux_vars, :, j, i)

      for k=1:Tdim
        fill!(nrm_bar, 0)
        flux_bar_k = sview(flux_el_bar, :, j, k)

        for p=1:Tdim
          nrm[p] = mesh.dxidx[k, p, j, i]
        end

        calcEulerFlux_revm(eqn.params, q_j, aux_vars_j, nrm, nrm_bar, flux_bar_k)
        for p=1:Tdim
          mesh.dxidx_bar[k, p, j, i] += nrm_bar[p]
        end
      end  # end loop k
    end  # end loop j
  end  # end loop i

  return nothing
end  # end function

"""
  Reverse mode wrt q of [`calcVolumeIntegrals_nopre`](@ref)
"""
function calcVolumeIntegrals_nopre_revq(
           mesh::AbstractMesh{Tmsh},
           sbp::AbstractOperator,
           eqn::EulerData{Tsol, Tres, Tdim},
           opts) where {Tmsh, Tsol, Tres, Tdim}


  # flux in the parametric directions for a given element
  flux_el_bar = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, Tdim)
  nrm = zeros(Tmsh, mesh.dim)
#=
  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      q_j = ro_sview(eqn.q, :, j, i)
      aux_vars_j = ro_sview(eqn.aux_vars, :, j, i)

      for k=1:Tdim
        flux_k = sview(flux_el, :, j, k)

        # get the direction vector
        for p=1:Tdim
          nrm[p] = mesh.dxidx[k, p, j, i]
        end
        # consider calculating all directions at once
        # not sure if that will help because the data dependencies are
        # really simple
        calcEulerFlux(eqn.params, q_j, aux_vars_j, nrm, flux_k)
      end  # end loop k
    end  # end loop j

    res_i = sview(eqn.res, :, :, i)
    for k=1:Tdim
      weakDifferentiateElement!(sbp, k, sview(flux_el, :, :, k), res_i, SummationByParts.Add(), true)
    end
=#
  #-----------------------------
  # reverse sweep
  for i=1:mesh.numEl
    res_bar_i = sview(eqn.res_bar, :, :, i)
    fill!(flux_el_bar, 0)
    for k=1:Tdim
      flux_bar_k = sview(flux_el_bar, :, :, k)
      weakDifferentiateElement_rev!(sbp, k, flux_bar_k, res_bar_i,
                                    SummationByParts.Add(), true)
    end

    for j=1:mesh.numNodesPerElement
      q_j = ro_sview(eqn.q, :, j, i)
      q_bar_j = sview(eqn.q_bar, :, j, i)
      aux_vars_j = ro_sview(eqn.aux_vars, :, j, i)

      for k=1:Tdim
        flux_bar_k = sview(flux_el_bar, :, j, k)

        for p=1:Tdim
          nrm[p] = mesh.dxidx[k, p, j, i]
        end

        calcEulerFlux_revq(eqn.params, q_j, q_bar_j, aux_vars_j, nrm, flux_bar_k)
      end  # end loop k
    end  # end loop j
  end  # end loop i

  return nothing
end  # end function



"""
  Computes the derivative of the strong form volume terms with
  respect to `q`, ie.

  d/dq (-Q * f)

  but only the mesh.numDofPerNode x mesh.numDofPerNode diagonal block

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * assembler: used to assemble the contribution of each element into the
                Jacobian
  
"""
function calcVolumeIntegralsStrong_nopre_diff(
                                   mesh::AbstractMesh{Tmsh},
                                   sbp::AbstractOperator,
                                   eqn::EulerData{Tsol, Tres, Tdim},
                                   opts,
                                   assembler::AssembleElementData) where {Tmsh, Tsol, Tres, Tdim}


  @assert eqn.params.use_Minv != 1  # use_Minv not supported (TODO: why?)

  data = eqn.params.calc_volume_integrals_data
  @unpack data flux_jac res_jac nrm

  for i=1:mesh.numEl
    fill!(flux_jac, 0)
    fill!(res_jac, 0)
    for j=1:mesh.numNodesPerElement
      q_j = ro_sview(eqn.q, :, j, i)
      aux_vars_j = ro_sview(eqn.aux_vars, :, j, i)

      # compute dF/dq
      for k=1:Tdim
        fluxjac_k = sview(flux_jac, :, :, j, k)

        # get the direction vector
        for p=1:Tdim
          nrm[p] = mesh.dxidx[k, p, j, i]
        end
        calcEulerFlux_diff(eqn.params, q_j, aux_vars_j, nrm, fluxjac_k)
      end  # end loop k
    end  # end loop j

    # compute dR/dq
    for k=1:Tdim
      weakDifferentiateElement_jac!(sbp, k, sview(flux_jac, :, :, :, k), res_jac, SummationByParts.Subtract(), false)
    end

    
    if eqn.params.use_Minv == 1
      jac_el = sview(mesh.jac, :, i)
      applyMinvElement(jac_el, sbp.w, res_jac)
    end

    # assemble element level jacobian into the residual
    assembleElement(assembler, mesh, i, res_jac)

  end  # end loop i

  return nothing
end  # end function

"""
  Differentiated version of [`calcVolumeIntegralsSplitForm`](@ref)

  **Inputs**
   * mesh
   * sbp
   * eqn
   * opts
   * functor: the numerical flux function F, of type FluxType_diff
"""
function calcVolumeIntegralsSplitForm_diff(
                mesh::AbstractMesh{Tmsh},
                sbp::AbstractOperator,
                eqn::EulerData{Tsol, Tres, Tdim}, opts,
                functor::FluxType_diff,
                assembler::AssembleElementData) where {Tmsh, Tsol, Tres, Tdim}

  if opts["use_staggered_grid"]
    error("Jacobian of staggered grid volume integrals not supported")
  else
    # do the curvilinear version even if msh is linear
    calcVolumeIntegralsSplitFormCurvilinear_diff(mesh, sbp, eqn, opts, functor,
                                                 assembler)
  end

  return nothing
end



"""
  Computes the Jacobian contribution from [`calcVolumeIntegralsSplitFormCurvilinear`](@ref).
  
  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * functor: [`FluxType_diff`](@ref). 
   * assembler: used to assemble the contribution of each element into the
                Jacobian
 
"""
function calcVolumeIntegralsSplitFormCurvilinear_diff(
                mesh::AbstractMesh{Tmsh}, sbp::AbstractOperator,
                eqn::EulerData{Tsol, Tres, Tdim}, opts,
                functor::FluxType_diff,
                assembler::AssembleElementData) where {Tmsh, Tsol, Tres, Tdim}

  dxidx = mesh.dxidx
  res = eqn.res
  aux_vars = eqn.aux_vars
  params = eqn.params

  data = params.calc_volume_integrals_data
  @unpack data nrmD F_d Sx res_jac

  flux_jacL = zeros(Tres, mesh.numDofPerNode, mesh.numDofPerNode, Tdim)
  flux_jacR = zeros(Tres, mesh.numDofPerNode, mesh.numDofPerNode, Tdim)

  # S is calculated in x-y-z, so the normal vectors should be the unit normals
  fill!(nrmD, 0.0)
  for d=1:Tdim
    nrmD[d, d] = 1
  end

  for i=1:mesh.numEl
    # get S for this element
    dxidx_i = ro_sview(dxidx, :, :, :, i)
    calcSCurvilinear(sbp, dxidx_i, Sx)

    fill!(res_jac, 0.0)
    for j=1:mesh.numNodesPerElement
      q_j = ro_sview(eqn.q, :, j, i)
      aux_vars_j = ro_sview(aux_vars, :, j, i)
      for k=1:(j-1)  # loop over lower triangle of S
        q_k = ro_sview(eqn.q, :, k, i)

        # calculate the numerical flux functions in all Tdim
        # directions at once
        fill!(flux_jacL, 0.0)
        fill!(flux_jacR, 0.0)
        functor(params, q_j, q_k, aux_vars_j, nrmD, flux_jacL, flux_jacR)

        @simd for d=1:Tdim
          @simd for p=1:mesh.numDofPerNode
            @simd for q=1:mesh.numDofPerNode
              res_jac[q, p, j, j] -= 2*Sx[j, k, d]*flux_jacL[q, p, d]
              res_jac[q, p, j, k] -= 2*Sx[j, k, d]*flux_jacR[q, p, d]

              res_jac[q, p, k, j] += 2*Sx[j, k, d]*flux_jacL[q, p, d]
              res_jac[q, p, k, k] += 2*Sx[j, k, d]*flux_jacR[q, p, d]
            end
          end
        end

      end  # end k loop
    end  # end j loop

    if eqn.params.use_Minv == 1
      jac_el = sview(mesh.jac, :, i)
      applyMinvElement(jac_el, sbp.w, res_jac)
    end

    # assemble element level jacobian into the residual
    assembleElement(assembler, mesh, i, res_jac)
  end  # end i loop

  return nothing
end


function calcVolumeIntegralsSplitForm_revm(
                mesh::AbstractMesh{Tmsh},
                sbp::AbstractOperator,
                eqn::EulerData{Tsol, Tres, Tdim}, opts,
                functor::FluxType,
                functor_revm::FluxType_revm) where {Tmsh, Tsol, Tres, Tdim}

  if opts["use_staggered_grid"]
    error("staggered grid not supported in reverse mode wrt metrics")
  else
    if mesh.coord_order == 1
      calcVolumeIntegralsSplitFormLinear_revm(mesh, sbp, eqn, opts, functor,
                                                   functor_revm)
    else
      calcVolumeIntegralsSplitFormCurvilinear_revm(mesh, sbp, eqn, opts, functor,
                                                   functor_revm)
    end
  end

  return nothing
end


function calcVolumeIntegralsSplitForm_revq(
                mesh::AbstractMesh{Tmsh},
                sbp::AbstractOperator,
                eqn::EulerData{Tsol, Tres, Tdim}, opts,
                functor::FluxType,
                functor_revq::FluxType_revq) where {Tmsh, Tsol, Tres, Tdim}

  if opts["use_staggered_grid"]
    error("staggered grid not supported in reverse mode wrt metrics")
  else
    if mesh.coord_order == 1
      calcVolumeIntegralsSplitFormLinear_revq(mesh, sbp, eqn, opts, functor,
                                                   functor_revq)
    else
      calcVolumeIntegralsSplitFormCurvilinear_revq(mesh, sbp, eqn, opts, functor,
                                                   functor_revq)
    end
  end

  return nothing
end




"""
  Reverse mode of [`calcVolumeIntegralsSplitFormLinear`](@ref).  The _bar fields
  of the mesh are updated.

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * functor: the two point flux function (`FluxType`)
   * functor_revm: the reverse mode wrt metrics flux function (`FluxType_revm`)
"""
function calcVolumeIntegralsSplitFormLinear_revm(
                                        mesh::AbstractMesh{Tmsh},
                                        sbp::AbstractOperator,
                                        eqn::EulerData{Tsol, Tres, Tdim}, opts,
                                        functor::FluxType,
                                        functor_revm::FluxType_revm) where {Tmsh, Tsol, Tres, Tdim}

#  println("----- entered calcVolumeIntegralsSplitForm -----")
  dxidx = mesh.dxidx
  res = eqn.res
  res_bar = eqn.res_bar
  q = eqn.q
  aux_vars = eqn.aux_vars
  params = eqn.params
  data = params.calc_volume_integrals_data
  @unpack data nrmD F_d S F_d_bar nrmD_bar

  F_d_bar = zeros(F_d)
  nrmD_bar = zeros(Tmsh, mesh.dim, mesh.dim)


  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      q_j = ro_sview(q, :, j, i)
      aux_vars_j = ro_sview(aux_vars, :, j, i)
      for k=1:(j-1)  # loop over lower triangle of S
        q_k = ro_sview(q, :, k, i)
        # calcaulate the normal vector in each parametric directions
        for d=1:Tdim
          # get the normal vector
          for p=1:Tdim
            nrmD[p, d] = dxidx[d, p, j, i]
          end
        end

        # calculate the numerical flux functions in all Tdim
        # directions at once
#        functor(params, q_j, q_k, aux_vars_j, nrmD, F_d)

        fill!(F_d_bar, 0)
        @simd for d=1:Tdim
          # update residual
          @simd for p=1:(Tdim+2)
            #res[p, j, i] -= 2*S[j, k, d]*F_d[p, d]
            #res[p, k, i] += 2*S[j, k, d]*F_d[p, d]

            #--------------------
            # reverse sweep
            F_d_bar[p, d] -= 2*S[j, k, d]*res_bar[p, j, i]
            F_d_bar[p, d] += 2*S[j, k, d]*res_bar[p, k, i]
          end  # end p loop
        end  # end d loop

        fill!(nrmD_bar, 0)
        functor_revm(params, q_j, q_k, aux_vars_j, nrmD, nrmD_bar, F_d_bar)

        for d=1:Tdim
          for p=1:Tdim
            mesh.dxidx_bar[d, p, j, i] += nrmD_bar[p, d]
          end
        end

      end  # end k loop
    end  # end j loop
  end  # end i loop

  return nothing
end


"""
  Reverse mode of [`calcVolumeIntegralsSplitFormLinear`](@ref) wrt q.
  The _bar fields
  of the mesh are updated.

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * functor: the two point flux function (`FluxType`)
   * functor_revq: the reverse mode wrt metrics flux function (`FluxType_revq`)
"""
function calcVolumeIntegralsSplitFormLinear_revq(
                                        mesh::AbstractMesh{Tmsh},
                                        sbp::AbstractOperator,
                                        eqn::EulerData{Tsol, Tres, Tdim}, opts,
                                        functor::FluxType,
                                        functor_revq::FluxType_revq) where {Tmsh, Tsol, Tres, Tdim}

#  println("----- entered calcVolumeIntegralsSplitForm -----")
  dxidx = mesh.dxidx
  res = eqn.res
  res_bar = eqn.res_bar
  q = eqn.q
  aux_vars = eqn.aux_vars
  params = eqn.params
  data = params.calc_volume_integrals_data
  @unpack data nrmD F_d S F_d_bar

  F_d_bar = zeros(F_d)
  nrmD_bar = zeros(Tmsh, mesh.dim, mesh.dim)


  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      q_j = ro_sview(q, :, j, i)
      q_bar_j = sview(eqn.q_bar, :, j, i)
      aux_vars_j = ro_sview(aux_vars, :, j, i)
      for k=1:(j-1)  # loop over lower triangle of S
        q_k = ro_sview(q, :, k, i)
        q_bar_k = sview(eqn.q_bar, :, k, i)
        # calcaulate the normal vector in each parametric directions
        for d=1:Tdim
          # get the normal vector
          for p=1:Tdim
            nrmD[p, d] = dxidx[d, p, j, i]
          end
        end

        # calculate the numerical flux functions in all Tdim
        # directions at once
#        functor(params, q_j, q_k, aux_vars_j, nrmD, F_d)

        fill!(F_d_bar, 0)
        @simd for d=1:Tdim
          # update residual
          @simd for p=1:(Tdim+2)
            #res[p, j, i] -= 2*S[j, k, d]*F_d[p, d]
            #res[p, k, i] += 2*S[j, k, d]*F_d[p, d]

            #--------------------
            # reverse sweep
            F_d_bar[p, d] -= 2*S[j, k, d]*res_bar[p, j, i]
            F_d_bar[p, d] += 2*S[j, k, d]*res_bar[p, k, i]
          end  # end p loop
        end  # end d loop

        functor_revq(params, q_j, q_bar_j, q_k, q_bar_k, aux_vars_j, nrmD, 
                     F_d_bar)
      end  # end k loop
    end  # end j loop
  end  # end i loop

  return nothing
end





"""
  Reverse mode wrt metrics of [`calcVolumeIntegralsSplitFormCurvilinear`](@ref)
"""
function calcVolumeIntegralsSplitFormCurvilinear_revm(
                      mesh::AbstractMesh{Tmsh},
                      sbp::AbstractOperator,
                      eqn::EulerData{Tsol, Tres, Tdim}, opts,
                      functor::FluxType,
                      functor_revm::FluxType_revm) where {Tmsh, Tsol, Tres, Tdim}


#  println("\nentered calcVolumeIntegralsSplitFormCurvilinear_revm")
  dxidx = mesh.dxidx
  res = eqn.res
  res_bar = eqn.res_bar
  q = eqn.q
  aux_vars = eqn.aux_vars
  params = eqn.params

  data = params.calc_volume_integrals_data
  @unpack data nrmD F_d Sx Sx_bar F_d_bar

  # S is calculated in x-y-z, so the normal vectors should be the unit normals
  fill!(nrmD, 0.0)
  for d=1:Tdim
    nrmD[d, d] = 1
  end

  for i=1:mesh.numEl
    # get S for this element
    dxidx_i = ro_sview(dxidx, :, :, :, i)
    dxidx_i_bar = sview(mesh.dxidx_bar, :, :, :, i)
    calcSCurvilinear(sbp, dxidx_i, Sx)
    fill!(Sx_bar, 0)

    for j=1:mesh.numNodesPerElement
      q_j = ro_sview(q, :, j, i)
      aux_vars_j = ro_sview(aux_vars, :, j, i)
      for k=1:(j-1)  # loop over lower triangle of S
        q_k = ro_sview(q, :, k, i)

        # calculate the numerical flux functions in all Tdim
        # directions at once
        functor(params, q_j, q_k, aux_vars_j, nrmD, F_d)

        @simd for d=1:Tdim
          # update residual
          @simd for p=1:(Tdim+2)
#            res[p, j, i] -= 2*Sx[j, k, d]*F_d[p, d]
#            res[p, k, i] += 2*Sx[j, k, d]*F_d[p, d]

            # because the flux is computed in the unit normal directions, it
            # has no connection to anything in mesh, so flux_bar is not needed
            Sx_bar[j, k, d] -= 2*F_d[p, d]*res_bar[p, j, i]
            Sx_bar[j, k, d] += 2*F_d[p, d]*res_bar[p, k, i]
          end
        end  # end d loop

      end  # end k loop
    end  # end j loop

    calcSCurvilinear_rev(sbp, dxidx_i, dxidx_i_bar, Sx, Sx_bar)

  end  # end i loop

  return nothing
end


"""
  Reverse mode wrt q of [`calcVolumeIntegralsSplitFormCurvilinear`](@ref)
"""
function calcVolumeIntegralsSplitFormCurvilinear_revq(
                      mesh::AbstractMesh{Tmsh},
                      sbp::AbstractOperator,
                      eqn::EulerData{Tsol, Tres, Tdim}, opts,
                      functor::FluxType,
                      functor_revq::FluxType_revq) where {Tmsh, Tsol, Tres, Tdim}


#  println("\nentered calcVolumeIntegralsSplitFormCurvilinear_revq")
  dxidx = mesh.dxidx
  res = eqn.res
  res_bar = eqn.res_bar
  q = eqn.q
  aux_vars = eqn.aux_vars
  params = eqn.params

  data = params.calc_volume_integrals_data
  @unpack data nrmD F_d Sx F_d_bar

  # S is calculated in x-y-z, so the normal vectors should be the unit normals
  fill!(nrmD, 0.0)
  for d=1:Tdim
    nrmD[d, d] = 1
  end

  for i=1:mesh.numEl
    # get S for this element
    dxidx_i = ro_sview(dxidx, :, :, :, i)
    dxidx_i_bar = sview(mesh.dxidx_bar, :, :, :, i)
    calcSCurvilinear(sbp, dxidx_i, Sx)

    for j=1:mesh.numNodesPerElement
      q_j = ro_sview(q, :, j, i)
      q_bar_j = sview(eqn.q_bar, :, j, i)
      aux_vars_j = ro_sview(aux_vars, :, j, i)
      for k=1:(j-1)  # loop over lower triangle of S
        q_k = ro_sview(q, :, k, i)
        q_bar_k = sview(eqn.q_bar, :, k, i)

        # calculate the numerical flux functions in all Tdim
        # directions at once
        #functor(params, q_j, q_k, aux_vars_j, nrmD, F_d)

        fill!(F_d_bar, 0)
        @simd for d=1:Tdim
          # update residual
          @simd for p=1:(Tdim+2)
#            res[p, j, i] -= 2*Sx[j, k, d]*F_d[p, d]
#            res[p, k, i] += 2*Sx[j, k, d]*F_d[p, d]

            F_d_bar[p, d] -= 2*Sx[j, k, d]*res_bar[p, j, i]
            F_d_bar[p, d] += 2*Sx[j, k, d]*res_bar[p, k, i]
          end
        end  # end d loop

        functor_revq(params, q_j, q_bar_j, q_k, q_bar_k, aux_vars_j, nrmD,
                     F_d_bar)

      end  # end k loop
    end  # end j loop

  end  # end i loop

  return nothing
end





"""
  Computes the jacobian of [`calcEulerFlux`](@ref) with respect to `q`.
  Methods are available for 2D and 3D

  The caller must zero out the output array (if required)

  **Inputs**

   * params: ParamType, conservative variables only
   * q: vector of conservative variables at node
   * aux_vars: auxiliary variables at the node
   * dir: direction vector (possibly scaled) to compute the flux jacobian in

  **Inputs/Outputs**

   * Fjac: flux jacobian, numDofPerNode x numDofPerNode, summed into

"""
function calcEulerFlux_diff(params::ParamType{2, :conservative},
                      q::AbstractArray{Tsol,1},
                      aux_vars::AbstractArray{Tres, 1},
                      dir::AbstractArray{Tmsh},  Fjac::AbstractArray{Tsol,2}) where {Tmsh, Tsol, Tres}
# calculates the Euler flux in a particular direction at a point
# eqn is the equation type
# q is the vector (of length 4), of the conservative variables at the point
# aux_vars is the vector of auxiliary variables at the point
# dir is a vector of length 2 that specifies the direction
# F is populated with the flux Jacobian
# 2D  only


  p_dot = params.eulerfluxdata.p_dot
  press = calcPressure_diff(params, q, p_dot)
#  press = getPressure(aux_vars)
#  press = @getPressure(aux_vars)
  fac = 1/q[1]
  U = (q[2]*dir[1] + q[3]*dir[2])*fac
  U_dot1 = -(q[2]*dir[1] + q[3]*dir[2])*fac*fac
  U_dot2 = dir[1]*fac
  U_dot3 = dir[2]*fac

  # F[1] = q[1]*U
  # F[2] = q[2]*U + dir[1]*press
  # F[3] = q[3]*U + dir[2]*press
  # F[4] = (q[4] + press)*U
  Fjac[1, 1] += U + q[1]*U_dot1
  Fjac[2, 1] +=     q[2]*U_dot1 + dir[1]*p_dot[1]
  Fjac[3, 1] +=     q[3]*U_dot1 + dir[2]*p_dot[1]
  Fjac[4, 1] +=     q[4]*U_dot1 + press*U_dot1 + U*p_dot[1]

  Fjac[1, 2] +=     q[1]*U_dot2
  Fjac[2, 2] += U + q[2]*U_dot2 + dir[1]*p_dot[2]
  Fjac[3, 2] +=     q[3]*U_dot2 + dir[2]*p_dot[2]
  Fjac[4, 2] +=     q[4]*U_dot2 + press*U_dot2 + U*p_dot[2]

  Fjac[1, 3] +=     q[1]*U_dot3
  Fjac[2, 3] +=     q[2]*U_dot3 + dir[1]*p_dot[3]
  Fjac[3, 3] += U + q[3]*U_dot3 + dir[2]*p_dot[3]
  Fjac[4, 3] +=     q[4]*U_dot3 + press*U_dot3 + U*p_dot[3]

  Fjac[1, 4] += 0
  Fjac[2, 4] += dir[1]*p_dot[4]
  Fjac[3, 4] += dir[2]*p_dot[4]
  Fjac[4, 4] += U + U*p_dot[4]

  return nothing

end


function calcEulerFlux_diff(params::ParamType{3},
                      q::AbstractArray{Tsol,1},
                      aux_vars::AbstractArray{Tres, 1},
                      dir::AbstractArray{Tmsh},  Fjac::AbstractArray{Tsol,2}) where {Tmsh, Tsol, Tres}
# calculates the Euler flux in a particular direction at a point
# eqn is the equation type
# q is the vector (of length 4), of the conservative variables at the point
# aux_vars is the vector of auxiliary variables at the point
# dir is a vector of length 2 that specifies the direction
# F is populated with the flux Jacobian
# 2D  only


  p_dot = params.eulerfluxdata.p_dot
  press = calcPressure_diff(params, q, p_dot)
#  press = getPressure(aux_vars)
#  press = @getPressure(aux_vars)
  fac = 1/q[1]
  U = (q[2]*dir[1] + q[3]*dir[2] + q[4]*dir[3])*fac
  U_dot1 = -(q[2]*dir[1] + q[3]*dir[2] + q[4]*dir[3])*fac*fac
  U_dot2 = dir[1]*fac
  U_dot3 = dir[2]*fac
  U_dot4 = dir[3]*fac

  # F[1] = q[1]*U
  # F[2] = q[2]*U + dir[1]*press
  # F[3] = q[3]*U + dir[2]*press
  # F[4] = q[4]*U + dir[3]*press
  # F[4] = (q[5] + press)*U
  Fjac[1, 1] += U + q[1]*U_dot1
  Fjac[2, 1] +=     q[2]*U_dot1 + dir[1]*p_dot[1]
  Fjac[3, 1] +=     q[3]*U_dot1 + dir[2]*p_dot[1]
  Fjac[4, 1] +=     q[4]*U_dot1 + dir[3]*p_dot[1]
  Fjac[5, 1] +=     q[5]*U_dot1 + press*U_dot1 + U*p_dot[1]

  Fjac[1, 2] +=     q[1]*U_dot2
  Fjac[2, 2] += U + q[2]*U_dot2 + dir[1]*p_dot[2]
  Fjac[3, 2] +=     q[3]*U_dot2 + dir[2]*p_dot[2]
  Fjac[4, 2] +=     q[4]*U_dot2 + dir[3]*p_dot[2]
  Fjac[5, 2] +=     q[5]*U_dot2 + press*U_dot2 + U*p_dot[2]

  Fjac[1, 3] +=     q[1]*U_dot3
  Fjac[2, 3] +=     q[2]*U_dot3 + dir[1]*p_dot[3]
  Fjac[3, 3] += U + q[3]*U_dot3 + dir[2]*p_dot[3]
  Fjac[4, 3] +=     q[4]*U_dot3 + dir[3]*p_dot[3]
  Fjac[5, 3] +=     q[5]*U_dot3 + press*U_dot3 + U*p_dot[3]

  Fjac[1, 4] +=     q[1]*U_dot4
  Fjac[2, 4] +=     q[2]*U_dot4 + dir[1]*p_dot[4]
  Fjac[3, 4] +=     q[3]*U_dot4 + dir[2]*p_dot[4]
  Fjac[4, 4] += U + q[4]*U_dot4 + dir[3]*p_dot[4]
  Fjac[5, 4] +=     q[5]*U_dot4 + press*U_dot4 + U*p_dot[4]

  Fjac[1, 5] += 0
  Fjac[2, 5] += dir[1]*p_dot[5]
  Fjac[3, 5] += dir[2]*p_dot[5]
  Fjac[4, 5] += dir[3]*p_dot[5]
  Fjac[5, 5] += U + U*p_dot[5]

  return nothing

end

"""
  Computes the gradient of pressure with respect to `q` at a node.
  Methods are available in 2D and 3D

  **Inputs**

   * params: ParamType, conservative variables only
   * q: vector of conservative variables at the node

  **Inputs/Outputs**

   * pdot: vector of length numDofPerNode, overwritten with derivative of `p` 
           wrt `q` (overwritten)
"""
function calcPressure_diff(params::ParamType{2, :conservative},
                      q::AbstractArray{Tsol,1}, p_dot::AbstractVector{Tsol} ) where Tsol
  # calculate pressure for a node
  # q is a vector of length 4 of the conservative variables

  t1 = 1/(q[1]*q[1])
  t2 = q[2]*q[2]
  t3 = q[3]*q[3]

  p_dot[1] = (params.gamma_1)*( 0.5*(t2*t1 + t3*t1))
  p_dot[2] = -(params.gamma_1)*(q[2]/q[1])
  p_dot[3] = -(params.gamma_1)*(q[3]/q[1])
  p_dot[4] = params.gamma_1

  return  (params.gamma_1)*(q[4] - 0.5*(t2 + t3)/q[1])
end



function calcPressure_diff(params::ParamType{3, :conservative},
                      q::AbstractArray{Tsol,1}, p_dot::AbstractVector{Tsol} ) where Tsol
  # calculate pressure for a node

  t1 = 1/(q[1]*q[1])
  t2 = q[2]*q[2]
  t3 = q[3]*q[3]
  t4 = q[4]*q[4]

  p_dot[1] =  params.gamma_1*( 0.5*(t2 + t3 + t4)*t1)
  p_dot[2] = -params.gamma_1*(q[2]/q[1])
  p_dot[3] = -params.gamma_1*(q[3]/q[1])
  p_dot[4] = -params.gamma_1*(q[4]/q[1])
  p_dot[5] =  params.gamma_1

  return (params.gamma_1)*(q[5] - 0.5*(t2 + t3 + t4)/q[1])
end

function calcPressure_hess(params::ParamType{2, :conservative},
                      q::AbstractArray{Tsol,1}, p_dot::AbstractMatrix{Tsol} ) where Tsol
  # calculate pressure for a node
  # q is a vector of length 4 of the conservative variables

  t1 = 1/(q[1]*q[1]); t1_dot1 = -2*t1/q[1]
  t2 = q[2]*q[2]; t2_dot2 = 2*q[2]
  t3 = q[3]*q[3]; t3_dot3 = 2*q[3]

  #p_dot[1] = (params.gamma_1)*( 0.5*(t2*t1 + t3*t1))
  #p_dot[2] = -(params.gamma_1)*(q[2]/q[1])
  #p_dot[3] = -(params.gamma_1)*(q[3]/q[1])
  #p_dot[4] = params.gamma_1

  p_dot[1, 1] = params.gamma_1*0.5*t1_dot1*(t2 + t3)
  p_dot[2, 1] = params.gamma_1*0.5*(t2_dot2*t1)
  p_dot[3, 1] = params.gamma_1*0.5*(t3_dot3*t1)
  p_dot[4, 1] = 0

  p_dot[1, 2] = params.gamma_1*q[2]/(q[1]*q[1])
  p_dot[2, 2] = -params.gamma_1/q[1]
  p_dot[3, 2] = 0
  p_dot[4, 2] = 0

  p_dot[1, 3] = p_dot[3, 1]
  p_dot[2, 3] = 0
  p_dot[3, 3] = -params.gamma_1/q[1]
  p_dot[4, 3] = 0

  p_dot[1, 4] = 0
  p_dot[2, 4] = 0
  p_dot[3, 4] = 0
  p_dot[4, 4] = 0

  return  (params.gamma_1)*(q[4] - 0.5*(t2 + t3)/q[1])
end


function calcPressure_hess(params::ParamType{3, :conservative},
                      q::AbstractArray{Tsol,1}, p_dot::AbstractMatrix{Tsol} ) where Tsol
  # calculate pressure for a node
  # q is a vector of length 4 of the conservative variables

  t1 = 1/(q[1]*q[1]); t1_dot1 = -2*t1/q[1]
  t2 = q[2]*q[2]; t2_dot2 = 2*q[2]
  t3 = q[3]*q[3]; t3_dot3 = 2*q[3]
  t4 = q[4]*q[4]; t4_dot4 = 2*q[4]

  #p_dot[1] =  params.gamma_1*( 0.5*(t2 + t3 + t4)*t1)
  #p_dot[2] = -params.gamma_1*(q[2]/q[1])
  #p_dot[3] = -params.gamma_1*(q[3]/q[1])
  #p_dot[4] = -params.gamma_1*(q[4]/q[1])
  #p_dot[5] =  params.gamma_1


  p_dot[1, 1] = params.gamma_1*0.5*t1_dot1*(t2 + t3 + t4)
  p_dot[2, 1] = params.gamma_1*0.5*(t2_dot2*t1)
  p_dot[3, 1] = params.gamma_1*0.5*(t3_dot3*t1)
  p_dot[4, 1] = params.gamma_1*0.5*(t4_dot4*t1)
  p_dot[5, 1] = 0

  p_dot[1, 2] = params.gamma_1*q[2]/(q[1]*q[1])
  p_dot[2, 2] = -params.gamma_1/q[1]
  p_dot[3, 2] = 0
  p_dot[4, 2] = 0
  p_dot[5, 2] = 0

  p_dot[1, 3] = p_dot[3, 1]
  p_dot[2, 3] = 0
  p_dot[3, 3] = -params.gamma_1/q[1]
  p_dot[4, 3] = 0
  p_dot[5, 3] = 0

  p_dot[1, 4] = p_dot[4, 1]
  p_dot[2, 4] = 0
  p_dot[3, 4] = 0
  p_dot[4, 4] = -params.gamma_1/q[1]
  p_dot[5, 4] = 0

  p_dot[1, 5] = 0
  p_dot[2, 5] = 0
  p_dot[3, 5] = 0
  p_dot[5, 5] = 0

  return  (params.gamma_1)*(q[4] - 0.5*(t2 + t3)/q[1])
end


"""
  Reverse mode of the pressure calculation

  **Inputs**

   * params: ParamType object
   * q: vector of solution variables at the node (length numDofPerNode)
   * p_bar: seed value for pressure

  **Inputs/Outputs**

   * q_bar: vector to be updated (not overwritten) with the result
"""
function calcPressure_revq(params::ParamType{2, :conservative},
                           q::AbstractArray{Tsol, 1}, q_bar::AbstractArray{Tsol, 1},
                           p_bar::Number) where {Tsol}

  q_bar[1] +=  params.gamma_1*0.5*(q[2]*q[2] + q[3]*q[3])/(q[1]*q[1])*p_bar
  q_bar[2] += -params.gamma_1*q[2]*p_bar/q[1]
  q_bar[3] += -params.gamma_1*q[3]*p_bar/q[1]
  q_bar[4] +=  params.gamma_1*p_bar

  return nothing
end


function calcPressure_revq(params::ParamType{3, :conservative},
                           q::AbstractArray{Tsol, 1}, q_bar::AbstractArray{Tsol, 1},
                           p_bar::Number) where {Tsol}

  q_bar[1] +=  params.gamma_1*0.5*(q[2]*q[2] + q[3]*q[3] + q[4]*q[4])/(q[1]*q[1])*p_bar
  q_bar[2] += -params.gamma_1*q[2]*p_bar/q[1]
  q_bar[3] += -params.gamma_1*q[3]*p_bar/q[1]
  q_bar[4] += -params.gamma_1*q[4]*p_bar/q[1]
  q_bar[5] +=  params.gamma_1*p_bar

  return nothing
end

function calcPressure_diff_revq(params::ParamType{2, :conservative},
                      q::AbstractArray{Tsol,1}, q_bar::AbstractVector{Tres},
                      p_dot_bar::AbstractVector, press_bar::Number ) where {Tsol, Tres}
  # calculate pressure for a node
  # q is a vector of length 4 of the conservative variables

  t1 = 1/(q[1]*q[1])
  t2 = q[2]*q[2]
  t3 = q[3]*q[3]
 

#=
  press = (params.gamma_1)*(q[4] - 0.5*(t2 + t3)/q[1])
  p_dot[1] = (params.gamma_1)*( 0.5*(t2*t1 + t3*t1))
  p_dot[2] = -(params.gamma_1)*(q[2]/q[1])
  p_dot[3] = -(params.gamma_1)*(q[3]/q[1])
  p_dot[4] = params.gamma_1
=#
  # reverse sweep

  # pressure
  q_bar[1] += (params.gamma_1)*( 0.5*(t2*t1 + t3*t1))*press_bar
  q_bar[2] += -(params.gamma_1)*(q[2]/q[1])*press_bar
  q_bar[3] += -(params.gamma_1)*(q[3]/q[1])*press_bar
  q_bar[4] += params.gamma_1*press_bar

  # p_dot
  t1_bar = zero(Tres)
  t2_bar = zero(Tres)
  t3_bar = zero(Tres)

  t1_bar += params.gamma_1*(0.5*(t2 + t3))*p_dot_bar[1]
  t2_bar += params.gamma_1*0.5*t1*p_dot_bar[1]
  t3_bar += params.gamma_1*0.5*t1*p_dot_bar[1]

  q_bar[1] += params.gamma_1*q[2]*t1*p_dot_bar[2]
  q_bar[2] -= params.gamma_1*p_dot_bar[2]/q[1]

  q_bar[1] += params.gamma_1*q[3]*t1*p_dot_bar[3]
  q_bar[3] -= params.gamma_1*p_dot_bar[3]/q[1]

  # nothing to do for q_bar[4]

  q_bar[3] += 2*q[3]*t3_bar
  q_bar[2] += 2*q[2]*t2_bar
  q_bar[1] -= 2*t1_bar/(q[1]^3)

  return nothing
end


function calcPressure_diff_revq(params::ParamType{3, :conservative},
                      q::AbstractArray{Tsol,1}, q_bar::AbstractVector{Tres},
                      p_dot_bar::AbstractVector, press_bar::Number ) where {Tsol, Tres}

  # calculate pressure for a node

  t1 = 1/(q[1]*q[1])
  t2 = q[2]*q[2]
  t3 = q[3]*q[3]
  t4 = q[4]*q[4]
#=
  (params.gamma_1)*(q[5] - 0.5*(t2 + t3 + t4)/q[1])
  p_dot[1] =  params.gamma_1*( 0.5*(t2 + t3 + t4)*t1)
  p_dot[2] = -params.gamma_1*(q[2]/q[1])
  p_dot[3] = -params.gamma_1*(q[3]/q[1])
  p_dot[4] = -params.gamma_1*(q[4]/q[1])
  p_dot[5] =  params.gamma_1
=#

  # pressure
  q_bar[1] += (params.gamma_1)*( 0.5*(t2*t1 + t3*t1 + t4*t1))*press_bar
  q_bar[2] += -(params.gamma_1)*(q[2]/q[1])*press_bar
  q_bar[3] += -(params.gamma_1)*(q[3]/q[1])*press_bar
  q_bar[4] += -(params.gamma_1)*(q[4]/q[1])*press_bar
  q_bar[5] += params.gamma_1*press_bar

  # p_dot
  t1_bar = zero(Tres)
  t2_bar = zero(Tres)
  t3_bar = zero(Tres)
  t4_bar = zero(Tres)

  t1_bar += params.gamma_1*(0.5*(t2 + t3 + t4))*p_dot_bar[1]
  t2_bar += params.gamma_1*0.5*t1*p_dot_bar[1]
  t3_bar += params.gamma_1*0.5*t1*p_dot_bar[1]
  t4_bar += params.gamma_1*0.5*t1*p_dot_bar[1]

  q_bar[1] += params.gamma_1*q[2]*t1*p_dot_bar[2]
  q_bar[2] -= params.gamma_1*p_dot_bar[2]/q[1]

  q_bar[1] += params.gamma_1*q[3]*t1*p_dot_bar[3]
  q_bar[3] -= params.gamma_1*p_dot_bar[3]/q[1]

  q_bar[1] += params.gamma_1*q[4]*t1*p_dot_bar[4]
  q_bar[4] -= params.gamma_1*p_dot_bar[4]/q[1]

  # nothing to do for q_bar[5]

  q_bar[4] += 2*q[4]*t4_bar
  q_bar[3] += 2*q[3]*t3_bar
  q_bar[2] += 2*q[2]*t2_bar
  q_bar[1] -= 2*t1_bar/(q[1]^3)

  return nothing
end




"""
  Differentiated version of [`getLambdaMax`](@ref)

  **Inputs**

   * params: ParamType
   * qL: vector of conservative variables at a node
   * dir: direction vector (can be scaled)

  **Inputs/Outputs**

   * lambda_dot: derivative of lambda max wrt qL (overwritten)

  **Outputs**

   * lambda_max: maximum eigenvalue
"""
function getLambdaMax_diff(params::ParamType{2},
                      qL::AbstractVector{Tsol},
                      dir::AbstractVector{Tmsh},
                      lambda_dot::AbstractVector{Tres}) where {Tsol, Tres, Tmsh}

  gamma = params.gamma
  Un = zero(Tres)
  dA = zero(Tmsh)
  rhoLinv = 1/qL[1]
  rhoLinv_dotL1 = -rhoLinv*rhoLinv

  p_dot = params.get_lambda_max_data.p_dot
  pL = calcPressure_diff(params, qL, p_dot)
  aL = sqrt(gamma*pL*rhoLinv)  # speed of sound
  t1 = gamma*rhoLinv/(2*aL)
  t2 = gamma*pL/(2*aL)
  aL_dotL1 = t1*p_dot[1] + t2*rhoLinv_dotL1
  aL_dotL2 = t1*p_dot[2]
  aL_dotL3 = t1*p_dot[3]
  aL_dotL4 = t1*p_dot[4]


  Un_dotL1 = dir[1]*qL[2]*rhoLinv_dotL1
  Un_dotL2 = dir[1]*rhoLinv
  Un += dir[1]*qL[2]*rhoLinv

  Un_dotL1 += dir[2]*qL[3]*rhoLinv_dotL1
  Un_dotL3 = dir[2]*rhoLinv
  Un += dir[2]*qL[3]*rhoLinv

  for i=1:2
    dA += dir[i]*dir[i]
  end

  dA = sqrt(dA)

  lambda_max = absvalue3(Un) + dA*aL
  lambda_dot[1] = dA*aL_dotL1
  lambda_dot[2] = dA*aL_dotL2
  lambda_dot[3] = dA*aL_dotL3
  lambda_dot[4] = dA*aL_dotL4

  fac = absvalue3_deriv(Un)
  lambda_dot[1] += fac*Un_dotL1
  lambda_dot[2] += fac*Un_dotL2
  lambda_dot[3] += fac*Un_dotL3

  return lambda_max
end



function getLambdaMax_diff(params::ParamType{3},
                      qL::AbstractVector{Tsol},
                      dir::AbstractVector{Tmsh},
                      lambda_dot::AbstractVector{Tres}) where {Tsol, Tres, Tmsh}

  gamma = params.gamma
  Un = zero(Tres)
  dA = zero(Tmsh)
  rhoLinv = 1/qL[1]
  rhoLinv_dotL1 = -rhoLinv*rhoLinv

  p_dot = params.get_lambda_max_data.p_dot
  pL = calcPressure_diff(params, qL, p_dot)
  aL = sqrt(gamma*pL*rhoLinv)  # speed of sound
  t1 = gamma*rhoLinv/(2*aL)
  t2 = gamma*pL/(2*aL)
  aL_dotL1 = t1*p_dot[1] + t2*rhoLinv_dotL1
  aL_dotL2 = t1*p_dot[2]
  aL_dotL3 = t1*p_dot[3]
  aL_dotL4 = t1*p_dot[4]
  aL_dotL5 = t1*p_dot[5]


  Un_dotL1 = dir[1]*qL[2]*rhoLinv_dotL1
  Un_dotL2 = dir[1]*rhoLinv
  Un += dir[1]*qL[2]*rhoLinv

  Un_dotL1 += dir[2]*qL[3]*rhoLinv_dotL1
  Un_dotL3 = dir[2]*rhoLinv
  Un += dir[2]*qL[3]*rhoLinv

  Un_dotL1 += dir[3]*qL[4]*rhoLinv_dotL1
  Un_dotL4 = dir[3]*rhoLinv
  Un += dir[3]*qL[4]*rhoLinv


  for i=1:3
    dA += dir[i]*dir[i]
  end

  dA = sqrt(dA)

  lambda_max = absvalue(Un) + dA*aL
  lambda_dot[1] = dA*aL_dotL1
  lambda_dot[2] = dA*aL_dotL2
  lambda_dot[3] = dA*aL_dotL3
  lambda_dot[4] = dA*aL_dotL4
  lambda_dot[5] = dA*aL_dotL5

  fac = absvalue3_deriv(Un)
  lambda_dot[1] += fac*Un_dotL1
  lambda_dot[2] += fac*Un_dotL2
  lambda_dot[3] += fac*Un_dotL3
  lambda_dot[4] += fac*Un_dotL4

  return lambda_max
end

"""
  Method that computes the maximum eigenvalue in the volume (not in any
  particular direction

  **Inputs**

   * params: ParamType
   * qL: vector of conservative variables at a node
   * dir: direction vector (can be scaled)

  **Inputs/Outputs**

   * lambda_dot: derivative of lambda max wrt qL (overwritten)

  **Outputs**

   * lambda_max: maximum eigenvalue

"""
function getLambdaMax_diff(params::ParamType{2},
                      qL::AbstractVector{Tsol},
                      lambda_dot::AbstractVector{Tres}) where {Tsol, Tres}

  gamma = params.gamma
  Un = zero(Tres)
  rhoLinv = 1/qL[1]
  rhoLinv_dotL1 = -rhoLinv*rhoLinv

  p_dot = params.get_lambda_max_data.p_dot
  pL = calcPressure_diff(params, qL, p_dot)
  aL = sqrt(gamma*pL*rhoLinv)  # speed of sound
  t1 = gamma*rhoLinv/(2*aL)
  t2 = gamma*pL/(2*aL)
  aL_dotL1 = t1*p_dot[1] + t2*rhoLinv_dotL1
  aL_dotL2 = t1*p_dot[2]
  aL_dotL3 = t1*p_dot[3]
  aL_dotL4 = t1*p_dot[4]


  u_i = qL[2]*rhoLinv
  Un_dotL1 = 2*u_i*qL[2]*rhoLinv_dotL1
  Un_dotL2 = 2*u_i*rhoLinv
  Un += u_i*u_i

  u_i = qL[3]*rhoLinv
  Un_dotL1 += 2*u_i*qL[3]*rhoLinv_dotL1
  Un_dotL3  = 2*u_i*rhoLinv
  Un += u_i*u_i

  Un1 = sqrt(Un)
  Un1_dotL1 = (0.5/Un1)*Un_dotL1
  Un1_dotL2 = (0.5/Un1)*Un_dotL2
  Un1_dotL3 = (0.5/Un1)*Un_dotL3

  lambda_max = Un1 + aL
  lambda_dot[1] = Un1_dotL1 + aL_dotL1
  lambda_dot[2] = Un1_dotL2 + aL_dotL2
  lambda_dot[3] = Un1_dotL3 + aL_dotL3
  lambda_dot[4] =             aL_dotL4

  return lambda_max
end


function getLambdaMax_diff(params::ParamType{3},
                      qL::AbstractVector{Tsol},
                      lambda_dot::AbstractVector{Tres}) where {Tsol, Tres}

  gamma = params.gamma
  Un = zero(Tres)
  rhoLinv = 1/qL[1]
  rhoLinv_dotL1 = -rhoLinv*rhoLinv

  p_dot = params.get_lambda_max_data.p_dot
  pL = calcPressure_diff(params, qL, p_dot)
  aL = sqrt(gamma*pL*rhoLinv)  # speed of sound
  t1 = gamma*rhoLinv/(2*aL)
  t2 = gamma*pL/(2*aL)
  aL_dotL1 = t1*p_dot[1] + t2*rhoLinv_dotL1
  aL_dotL2 = t1*p_dot[2]
  aL_dotL3 = t1*p_dot[3]
  aL_dotL4 = t1*p_dot[4]
  aL_dotL5 = t1*p_dot[5]


  u_i = qL[2]*rhoLinv
  Un_dotL1 = 2*u_i*qL[2]*rhoLinv_dotL1
  Un_dotL2 = 2*u_i*rhoLinv
  Un += u_i*u_i

  u_i = qL[3]*rhoLinv
  Un_dotL1 += 2*u_i*qL[3]*rhoLinv_dotL1
  Un_dotL3  = 2*u_i*rhoLinv
  Un += u_i*u_i

  u_i = qL[4]*rhoLinv
  Un_dotL1 += 2*u_i*qL[4]*rhoLinv_dotL1
  Un_dotL4  = 2*u_i*rhoLinv
  Un += u_i*u_i


  Un1 = sqrt(Un)
  Un1_dotL1 = (0.5/Un1)*Un_dotL1
  Un1_dotL2 = (0.5/Un1)*Un_dotL2
  Un1_dotL3 = (0.5/Un1)*Un_dotL3
  Un1_dotL4 = (0.5/Un1)*Un_dotL4

  lambda_max = Un1 + aL
  lambda_dot[1] = Un1_dotL1 + aL_dotL1
  lambda_dot[2] = Un1_dotL2 + aL_dotL2
  lambda_dot[3] = Un1_dotL3 + aL_dotL3
  lambda_dot[4] = Un1_dotL4 + aL_dotL4
  lambda_dot[5] =             aL_dotL5

  return lambda_max
end



"""
  Reverse mode of [`getLambdaMax`](@ref) with respect to the metrics

  **Inputs**

   * params
   * qL
   * dir
   * lambda_bar: seed value for reverse mode

  **Inputs/Outputs**

   * nrm_bar: vector to be updated (not overwritten) with the result
"""
function getLambdaMax_revm(params::ParamType{Tdim}, 
                      qL::AbstractVector{Tsol}, 
                      dir::AbstractVector{Tmsh},
                      dir_bar::AbstractVector{Tmsh}, lambda_bar::Number) where {Tsol, Tmsh, Tdim}

  Tres = promote_type(Tsol, Tmsh)
  gamma = params.gamma
  Un = zero(Tres)
  dA = zero(Tmsh)
  rhoLinv = 1/qL[1]

  
  pL = calcPressure(params, qL)
  aL = sqrt(gamma*pL*rhoLinv)  # speed of sound

  for i=1:Tdim
    Un += dir[i]*qL[i+1]*rhoLinv
    dA += dir[i]*dir[i]
  end

  dA2 = sqrt(dA)

  lambda_max = absvalue(Un) + dA2*aL

  # reverse sweep
  fac = Un > 0 ? 1 : -1
  Un_bar = fac*lambda_bar
  dA2_bar = aL*lambda_bar
#  aL_bar = dA2*lambda_bar

  dA_bar = dA2_bar/(2*dA2)

  for i=1:Tdim
    # Un
    dir_bar[i] += qL[i+1]*rhoLinv*Un_bar
    # dA
    dir_bar[i] += 2*dir[i]*dA_bar
  end


  return lambda_max
end


"""
  Reverse mode wrt q of [`getLambdaMax`](@ref)

  **Inputs**

   * params
   * qL
   * dir
   * lambda_bar: reverse mode seed value

  **Inputs/Outputs**

   * qL_bar: adjoint part of qL, will be updated (not overwritten)
"""
function getLambdaMax_revq(params::ParamType{Tdim}, 
                           qL::AbstractVector{Tsol}, 
                           qL_bar::AbstractVector{Tsol},
                           dir::AbstractVector{Tmsh},
                           lambda_bar::Number) where {Tsol, Tmsh, Tdim}

  Tres = promote_type(Tsol, Tmsh)
  gamma = params.gamma
  Un = zero(Tres)
  dA = zero(Tmsh)
  rhoLinv = 1/qL[1]

  pL = calcPressure(params, qL)
  aL = sqrt(gamma*pL*rhoLinv)  # speed of sound

  for i=1:Tdim
    Un += dir[i]*qL[i+1]*rhoLinv
    dA += dir[i]*dir[i]
  end

  dA2 = sqrt(dA)

  lambda_max = absvalue(Un) + dA2*aL
  rhoLinv_bar = zero(Tsol)
 

  # reverse sweep
  fac = Un > 0 ? 1 : -1
  Un_bar = fac*lambda_bar
#  dA2_bar = aL*lambda_bar
  aL_bar = dA2*lambda_bar

  
  # dA2 = sqrt(dA)
#  dA_bar = dA2_bar/(2*dA2)

  rhoLinv_bar = zero(Tsol)
  for i=1:Tdim
    qL_bar[i+1] += dir[i]*rhoLinv*Un_bar
    rhoLinv_bar += dir[i]*qL[i+1]*Un_bar
    # dir is not being differentated here
  end

  # aL
  pL_bar = gamma*rhoLinv*aL_bar/(2*aL)
  rhoLinv_bar += gamma*pL*aL_bar/(2*aL)

  calcPressure_revq(params, qL, qL_bar, pL_bar)

  qL_bar[1] += -rhoLinv*rhoLinv*rhoLinv_bar
  
  return lambda_max
end


function getLambdaMax_revq(params::ParamType{Tdim}, 
                           qL::AbstractVector{Tsol}, 
                           qL_bar::AbstractVector{Tsol},
                           lambda_bar::Number) where {Tsol, Tdim}

  gamma = params.gamma
  Un = zero(Tsol)
  rhoLinv = 1/qL[1]

  pL = calcPressure(params, qL)
  aL = sqrt(gamma*pL*rhoLinv)  # speed of sound

  for i=1:Tdim
    u_i = qL[i+1]*rhoLinv
    Un += u_i*u_i
  end

  Un2 = sqrt(Un)
  lambda_max = Un2 + aL
 

  # reverse sweep
  fac = Un > 0 ? 1 : -1
  Un_bar = lambda_bar/(2*Un2)
  aL_bar = lambda_bar

  
  rhoLinv_bar = zero(Tsol)
  for i=1:Tdim
    u_i = qL[i+1]*rhoLinv

    u_i_bar = 2*u_i*Un_bar
    qL_bar[i+1] += u_i_bar*rhoLinv
    rhoLinv_bar += qL[i+1]*u_i_bar
  end

  # aL
  pL_bar = gamma*rhoLinv*aL_bar/(2*aL)
  rhoLinv_bar += gamma*pL*aL_bar/(2*aL)

  calcPressure_revq(params, qL, qL_bar, pL_bar)

  qL_bar[1] += -rhoLinv*rhoLinv*rhoLinv_bar
  
  return lambda_max
end

"""
  Differentiated version of [`calc2RWaveSpeeds`](@ref)

  **Inputs**

   * params
   * qL
   * qR
   * nrm: normal vector

  **Inputs/Outputs**

   * sL_dot: derivative of sL wrt qL and qR, `numDofPerNode` x 2
   * sR_dot: derivative of sR wrt qL and qR, `numDofPerNode` x 2

  **Outputs**

   * sL: slowest wave speed
   * sR: fastest wave speed
"""
function calc2RWaveSpeeds_diff(params::ParamType{Tdim}, qL::AbstractVector{Tsol},
                          qR::AbstractVector, nrm::AbstractVector{Tmsh},
                          sL_dot::AbstractMatrix, sR_dot::AbstractMatrix
                         ) where {Tdim, Tsol, Tmsh}

  @unpack params.tworwavespeeddata pL_dot pR_dot aL_dot aR_dot u_nrmL_dot u_nrmR_dot p_tr_dot qfL_dot qfR_dot

  Tres = promote_type(Tsol, Tmsh)
  numDofPerNode = length(qL)

  # compute pressure and speed of sound
  pL = calcPressure_diff(params, qL, pL_dot)
  pR = calcPressure_diff(params, qR, pR_dot)

  facL = params.gamma/qL[1]; facR = params.gamma/qR[1]
  aL = sqrt(facL*pL); aR = sqrt(facR*pR)

  facL /= 2*aL; facR /= 2*aR
  for i=1:numDofPerNode
    aL_dot[i] = facL*pL_dot[i]
    aR_dot[i] = facR*pR_dot[i]
  end
  aL_dot[1] -= facL*pL/qL[1]; aR_dot[1] -= facR*pR/qR[1]
  z = params.gamma_1/(2*params.gamma)

  # compute velocity in face normal direction
  u_nrmL = zero(Tres); u_nrmR = zero(Tres)
  fac = calcLength(params, nrm)
  for i=1:Tdim
    u_nrmL += qL[i+1]*nrm[i]; u_nrmL_dot[i+1] = nrm[i]/(fac*qL[1])
    u_nrmR += qR[i+1]*nrm[i]; u_nrmR_dot[i+1] = nrm[i]/(fac*qR[1])
  end
  u_nrmL_dot[1] = -u_nrmL/(fac*qL[1]*qL[1])
  u_nrmR_dot[1] = -u_nrmR/(fac*qR[1]*qR[1])

  u_nrmL /= fac*qL[1]; u_nrmR /= fac*qR[1]

  #qfL = 1.3
  #qfR = 1.02

  # compute p_tr
  pLz = pL^z; pRz = pR^z
  num = aL + aR - 0.5*params.gamma_1*(u_nrmR - u_nrmL)
  den = aL/pLz + aR/pRz

  t1 = num/den
  p_tr = t1^(1/z)

  # derivatives of num, den, t1, and p_tr
  for i=1:numDofPerNode
    num_dotL_i = aL_dot[i] + 0.5*params.gamma_1*u_nrmL_dot[i]
    num_dotR_i = aR_dot[i] - 0.5*params.gamma_1*u_nrmR_dot[i]

    den_dotL_i = aL_dot[i]/pLz + -z*aL*pL_dot[i]/(pLz*pL)
    den_dotR_i = aR_dot[i]/pRz + -z*aR*pR_dot[i]/(pRz*pR)

    t1_dotL_i = (num_dotL_i*den - num*den_dotL_i)/(den*den)
    t1_dotR_i = (num_dotR_i*den - num*den_dotR_i)/(den*den)

    fac_i = p_tr/(z*t1)
    p_tr_dot[1, i] = fac_i*t1_dotL_i
    p_tr_dot[2, i] = fac_i*t1_dotR_i
  end

  # compute q values
  if p_tr <= pL
    qfL = Tres(1.0)
    for i=1:numDofPerNode
      qfL_dot[1, i] = 0 ; qfL_dot[2, i] = 0
    end
  else
    qfL = sqrt(1 + (params.gamma + 1)*(p_tr/pL - 1)/(2*params.gamma))
    for i=1:numDofPerNode
      t1_dotL_i = (p_tr_dot[1, i]*pL - p_tr*pL_dot[i])/(pL*pL)
      t1_dotR_i = p_tr_dot[2, i]/pL

      qfL_dot[1, i] = (params.gamma + 1)*t1_dotL_i/(4*qfL*params.gamma)
      qfL_dot[2, i] = (params.gamma + 1)*t1_dotR_i/(4*qfL*params.gamma)
    end
  end

  if p_tr <= pR
    qfR = Tres(1.0)
    for i=1:numDofPerNode
      qfR_dot[1, i] = 0; qfR_dot[2, i] = 0
    end
  else
    qfR = sqrt(1 + (params.gamma + 1)*(p_tr/pR - 1)/(2*params.gamma))
    for i=1:numDofPerNode
      t1_dotL_i = p_tr_dot[1, i]/pR
      t1_dotR_i = (p_tr_dot[2, i]*pR - p_tr*pR_dot[i])/(pR*pR)

      qfR_dot[1, i] = (params.gamma + 1)*t1_dotL_i/(4*qfR*params.gamma)
      qfR_dot[2, i] = (params.gamma + 1)*t1_dotR_i/(4*qfR*params.gamma)
    end
  end

  # compute sL and sR
  sL = u_nrmL - aL*qfL
  sR = u_nrmR + aR*qfR

  for i=1:numDofPerNode
    sL_dot[i, 1] = u_nrmL_dot[i] - aL_dot[i]*qfL - aL*qfL_dot[1, i]
    sL_dot[i, 2] =                               - aL*qfL_dot[2, i]

    sR_dot[i, 1] =                                 aR*qfR_dot[1, i]
    sR_dot[i, 2] = u_nrmR_dot[i] + aR_dot[i]*qfR + aR*qfR_dot[2, i]
  end

#=
  # compute sL and sR
  sL = u_nrmL - aL*qfL
  sR = u_nrmR + aR*qfR
  for i=1:numDofPerNode
    sL_dot[i, 1] = u_nrmL_dot[i] - aL_dot[i]*qfL
    sR_dot[i, 2] = u_nrmR_dot[i] + aR_dot[i]*qfR 
  end
=#
  return sL, sR
end


function calc2RWaveSpeeds_revq(params::ParamType{Tdim},
                              qL::AbstractVector{Tsol}, qL_bar::AbstractVector,
                              qR::AbstractVector, qR_bar::AbstractVector,
                              nrm::AbstractVector{Tmsh},
                              sL_bar::Number, sR_bar::Number,
                              ) where {Tdim, Tsol, Tmsh}

  Tres = promote_type(Tsol, Tmsh)

  pL = calcPressure(params, qL); pR = calcPressure(params, qR)
  aL = sqrt(params.gamma*pL/qL[1]); aR = sqrt(params.gamma*pR/qR[1])
  z = params.gamma_1/(2*params.gamma)

  # compute velocity in face normal direction
  u_nrmL = zero(Tres); u_nrmR = zero(Tres)
  fac = calcLength(params, nrm)
  for i=1:Tdim
    u_nrmL += qL[i+1]*nrm[i]
    u_nrmR += qR[i+1]*nrm[i]
  end
  u_nrmL_orig = u_nrmL; u_nrmR_orig = u_nrmR
  u_nrmL /= fac*qL[1]; u_nrmR /= fac*qR[1]

  #qfL = 1.3
  #qfR = 1.02


  num = aL + aR - 0.5*params.gamma_1*(u_nrmR - u_nrmL)
  pLz = pL^z; pRz = pR^z
  den = aL/pLz + aR/pRz

  frac = num/den
  p_tr = frac^(1/z)

  # compute q values
  if p_tr <= pL
    qfL = Tres(1.0)
  else
    qfL = sqrt(1 + (params.gamma + 1)*(p_tr/pL - 1)/(2*params.gamma))
  end

  if p_tr <= pR
    qfR = Tres(1.0)
  else
    qfR = sqrt(1 + (params.gamma + 1)*(p_tr/pR - 1)/(2*params.gamma))
  end

  sL = u_nrmL - aL*qfL
  sR = u_nrmR + aR*qfR

  #------------------------
  # reverse sweep

  u_nrmL_bar =  sL_bar
  aL_bar     = -qfL*sL_bar
  pL_bar     = zero(Tres)
  qfL_bar    = -aL*sL_bar

  u_nrmR_bar = sR_bar
  aR_bar     = qfR*sR_bar
  qfR_bar    = aR*sR_bar
  pR_bar     = zero(Tres)

  p_tr_bar = zero(Tres)
  if p_tr <= pR
    pR_bar   = zero(Tres)
  else
    p_tr_bar += Tres( qfR_bar*(params.gamma + 1)/(4*params.gamma*qfR*pR) )
    pR_bar    = Tres( qfR_bar*( -(params.gamma + 1)*p_tr)/(4*qfR*params.gamma*pR*pR))
  end

  if p_tr <= pL
    pL_bar = zero(Tres)
  else
    p_tr_bar += Tres( qfL_bar*(params.gamma + 1)/(4*params.gamma*qfL*pL) )
    pL_bar    = Tres( qfL_bar*( -(params.gamma + 1)*p_tr)/(4*qfL*params.gamma*pL*pL))
  end

  #TODO: precompute things raised to fractional powers
  # TODO: if p_tr_bar = 0, is any of this necessary?
  num_bar = p_tr_bar*(1/(z*den))*p_tr/frac
  den_bar = -p_tr_bar*p_tr/(z*den)

  # compute num and den
  aL_bar += den_bar/pLz; aR_bar += den_bar/pRz
  pL_bar += -den_bar*z*aL/(pLz*pL); pR_bar += -den_bar*z*aR/(pRz*pR)

  aL_bar += num_bar; aR_bar += num_bar
  u_nrmL_bar += 0.5*params.gamma_1*num_bar
  u_nrmR_bar -= 0.5*params.gamma_1*num_bar

  # compute velocity in face normal direction
  qL_bar[1] += -u_nrmL_bar*u_nrmL_orig/(fac*qL[1]*qL[1])
  qR_bar[1] += -u_nrmR_bar*u_nrmR_orig/(fac*qR[1]*qR[1])
  u_nrmL_bar = u_nrmL_bar/(fac*qL[1]); u_nrmR_bar = u_nrmR_bar/(fac*qR[1])
  for i=1:Tdim
    qL_bar[i+1] += u_nrmL_bar*nrm[i]
    qR_bar[i+1] += u_nrmR_bar*nrm[i]
  end


  pL_bar += aL_bar*params.gamma/(2*aL*qL[1])
  pR_bar += aR_bar*params.gamma/(2*aR*qR[1])
  qL_bar[1] += aL_bar*(-params.gamma*pL/(qL[1]*qL[1]))/(2*aL)
  qR_bar[1] += aR_bar*(-params.gamma*pR/(qR[1]*qR[1]))/(2*aR)

  calcPressure_revq(params, qL, qL_bar, pL_bar)
  calcPressure_revq(params, qR, qR_bar, pR_bar)
end


function calc2RWaveSpeeds_revm(params::ParamType{Tdim},
                            qL::AbstractVector{Tsol}, qR::AbstractVector,
                            nrm::AbstractVector{Tmsh}, nrm_bar::AbstractVector,
                            sL_bar::Number, sR_bar::Number,
                            ) where {Tdim, Tsol, Tmsh}

  Tres = promote_type(Tsol, Tmsh)

  pL = calcPressure(params, qL); pR = calcPressure(params, qR)
  aL = sqrt(params.gamma*pL/qL[1]); aR = sqrt(params.gamma*pR/qR[1])
  z = params.gamma_1/(2*params.gamma)

  # compute velocity in face normal direction
  u_nrmL = zero(Tres); u_nrmR = zero(Tres)
  fac = calcLength(params, nrm)
  for i=1:Tdim
    u_nrmL += qL[i+1]*nrm[i]
    u_nrmR += qR[i+1]*nrm[i]
  end
  u_nrmL_orig = u_nrmL; u_nrmR_orig = u_nrmR
  u_nrmL /= fac*qL[1]; u_nrmR /= fac*qR[1]
#  qfL = 1.3
#  qfR = 1.02


  num = aL + aR - 0.5*params.gamma_1*(u_nrmR - u_nrmL)
  pLz = pL^z; pRz = pR^z
  den = aL/pLz + aR/pRz

  frac = num/den
  p_tr = frac^(1/z)

  # compute q values
  if p_tr <= pL
    qfL = Tres(1.0)
  else
    qfL = sqrt(1 + (params.gamma + 1)*(p_tr/pL - 1)/(2*params.gamma))
  end

  if p_tr <= pR
    qfR = Tres(1.0)
  else
    qfR = sqrt(1 + (params.gamma + 1)*(p_tr/pR - 1)/(2*params.gamma))
  end

  sL = u_nrmL - aL*qfL
  sR = u_nrmR + aR*qfR

  #------------------------
  # reverse sweep

  u_nrmL_bar =  sL_bar
  qfL_bar    = -aL*sL_bar

  u_nrmR_bar = sR_bar
  qfR_bar    = aR*sR_bar

  p_tr_bar = zero(Tres)
  if p_tr > pR
    p_tr_bar += Tres( qfR_bar*(params.gamma + 1)/(4*params.gamma*qfR*pR) )
  end

  if p_tr > pL
    p_tr_bar += Tres( qfL_bar*(params.gamma + 1)/(4*params.gamma*qfL*pL) )

  end
  num_bar = p_tr_bar*(1/(z*den))*p_tr/frac

  u_nrmL_bar += 0.5*params.gamma_1*num_bar
  u_nrmR_bar -= 0.5*params.gamma_1*num_bar

  # compute velocity in face normal direction
  fac_bar  = -u_nrmL_bar*u_nrmL_orig/(qL[1]*fac*fac)
  fac_bar += -u_nrmR_bar*u_nrmR_orig/(qR[1]*fac*fac)
  u_nrmL_bar = u_nrmL_bar/(fac*qL[1]); u_nrmR_bar = u_nrmR_bar/(fac*qR[1])
  for i=1:Tdim
    nrm_bar[i] += qL[i+1]*u_nrmL_bar
    nrm_bar[i] += qR[i+1]*u_nrmR_bar
  end
  calcLength_rev(params, nrm, nrm_bar, fac_bar)
end


