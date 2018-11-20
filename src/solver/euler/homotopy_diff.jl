# differentiated version of homotopy.jl
import PDESolver.evalHomotopyJacobian

function evalHomotopyJacobian(mesh::AbstractMesh, sbp::AbstractSBP,
                              eqn::EulerData, opts::Dict, 
                              assembler::AssembleElementData, lambda::Number)

  calcHomotopyDiss_jac(mesh, sbp, eqn, opts, assembler, lambda)
end

function calcHomotopyDiss_jac(mesh::AbstractDGMesh{Tmsh}, sbp, 
                          eqn::EulerData{Tsol, Tres}, opts, assembler, lambda) where {Tsol, Tres, Tmsh}

  # some checks for when parallelism is enabled
  @assert opts["parallel_data"] == PARALLEL_DATA_ELEMENT
  for i=1:mesh.npeers
    @assert eqn.shared_data[i].recv_waited
  end

  params = eqn.params
  #----------------------------------------------------------------------------
  # volume dissipation

  # compute the D operator in each direction
  D = zeros(mesh.numNodesPerElement, mesh.numNodesPerElement, mesh.dim)
  for d=1:mesh.dim
    D[:, :, d] = inv(diagm(sbp.w))*sbp.Q[:, :, d]
  end

  res_jac = eqn.params.calc_face_integrals_data.res_jacLL
  t2_dot = eqn.params.calc_face_integrals_data.res_jacLR  # work array
  t1 = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerElement)
  lambda_dot = zeros(Tres, mesh.numDofPerNode)

  nrm = zeros(Tmsh, mesh.dim)
  for el=1:mesh.numEl
    q_el = sview(eqn.q, :, :, el)
    fill!(res_jac, 0.0)
    for d1=1:mesh.dim
      fill!(t1, 0.0)
      fill!(t2_dot, 0.0)

      differentiateElement!(sbp, d1, q_el, t1)

      # t2 = t1*lambda_p
      # contribution t2_dot = lambda_p*t1_dot
      for p=1:mesh.numNodesPerElement
        q_p = sview(eqn.q, :, p, el)

        # get vector in xi direction defined by dim
        for k=1:mesh.dim
          nrm[k] = mesh.dxidx[d1, k, p, el]
        end

        lambda_max = getLambdaMax(eqn.params, q_p, nrm)

        for q=1:mesh.numNodesPerElement
          for j=1:mesh.numDofPerNode
#            for i=1:mesh.numDofPerNode
              t2_dot[j, j, p, q] += lambda_max*D[p, q, d1]
#            end
          end
        end
      end   # end loop p

      # contribution t2_dot += t1*lambda_p_dot
      for p=1:mesh.numNodesPerElement
        q_p = sview(eqn.q, :, p, el)

        # get vector in xi direction defined by dim
        for k=1:mesh.dim
          nrm[k] = mesh.dxidx[d1, k, p, el]
        end

        getLambdaMax_diff(eqn.params, q_p, nrm, lambda_dot)

        for j=1:mesh.numDofPerNode
          for i=1:mesh.numDofPerNode
            t2_dot[i, j, p, p] += t1[i, p]*lambda_dot[j]
          end
        end
      end  # end loop p


      # apply Q^T
#      for d2=1:mesh.dim
        for p=1:mesh.numNodesPerElement
          for q=1:mesh.numNodesPerElement
            for c=1:mesh.numNodesPerElement
              for j=1:mesh.numDofPerNode
                for i=1:mesh.numDofPerNode
                  res_jac[i, j, p, q] += sbp.Q[c, p, d1]*t2_dot[i, j, c, q]
                end
              end
            end
          end
        end
#      end  # end loop d2

    end  # end loop d1
 
    # negate res_jac for consistency with physics module
    for i=1:length(res_jac)
      res_jac[i] = -lambda*res_jac[i]
    end

     
    assembleElement(assembler, mesh, el, res_jac)
  end  # end loop el


  fill!(res_jac, 0.0)
  fill!(t2_dot, 0.0)


  #----------------------------------------------------------------------------
  # interface terms

  @unpack params.calc_face_integrals_data q_faceL q_faceR flux_dotL flux_dotR res_jacLL res_jacLR res_jacRL res_jacRR

  #=
  q_faceL = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
  q_faceR = zeros(q_faceL)
#  nrm2 = eqn.params.nrm2
  flux_dotL = zeros(Tres, mesh.numDofPerNode, mesh.numDofPerNode, mesh.numNodesPerFace)
  flux_dotR = zeros(flux_dotL)

  res_jacLL = params.res_jacLL
  res_jacLR = params.res_jacLR
  res_jacRL = params.res_jacRL
  res_jacRR = params.res_jacRR
  =#
  lambda_dotL = zeros(Tres, mesh.numDofPerNode)
  lambda_dotR = zeros(Tres, mesh.numDofPerNode)

  for i=1:mesh.numInterfaces
    iface_i = mesh.interfaces[i]
    qL = sview(eqn.q, :, :, iface_i.elementL)
    qR = sview(eqn.q, :, :, iface_i.elementR)
    fill!(res_jacLL, 0.0)
    fill!(res_jacLR, 0.0)
    fill!(res_jacRL, 0.0)
    fill!(res_jacRR, 0.0)

    interiorFaceInterpolate!(mesh.sbpface, iface_i, qL, qR, q_faceL, q_faceR)

    # calculate the flux jacobian at each face node
    for j=1:mesh.numNodesPerFace
      qL_j = sview(q_faceL, :, j)
      qR_j = sview(q_faceR, :, j)

      # get the face normal
      nrm2 = sview(mesh.nrm_face, :, j, i)

      lambda_max = getLambdaMaxSimple_diff(eqn.params, qL_j, qR_j, nrm2,
                                           lambda_dotL, lambda_dotR)

      for k=1:mesh.numDofPerNode
        # flux[k, j] = 0.5*lambda_max*(qL_j[k] - qR_j[k])
        for m=1:mesh.numDofPerNode
          flux_dotL[m, k, j] = 0.5*lambda_dotL[k]*(qL_j[m] - qR_j[m])
          flux_dotR[m, k, j] = 0.5*lambda_dotR[k]*(qL_j[m] - qR_j[m])
        end
        flux_dotL[k, k, j] += 0.5*lambda_max
        flux_dotR[k, k, j] -= 0.5*lambda_max
      end
    end  # end loop j

    # multiply by lambda here and it will get carried through
    # interiorFaceIntegrate_jac
    scale!(flux_dotL, lambda)
    scale!(flux_dotR, lambda)

    # compute dR/dq
    interiorFaceIntegrate_jac!(mesh.sbpface, iface_i, flux_dotL, flux_dotR,
                             res_jacLL, res_jacLR, res_jacRL, res_jacRR,
                             SummationByParts.Subtract())
    # assemble into the Jacobian
    assembleInterface(assembler, mesh.sbpface, mesh, iface_i, res_jacLL, res_jacLR,
                                                res_jacRL, res_jacRR)

  end  # end loop i

  fill!(res_jacLL, 0.0)
  fill!(res_jacLR, 0.0)
  fill!(res_jacRL, 0.0)
  fill!(res_jacRR, 0.0)


  #----------------------------------------------------------------------------
  # skipping boundary integrals
  # use nrm2, flux_jfacL from interface terms above
  if opts["homotopy_addBoundaryIntegrals"]
    qg = zeros(Complex128, mesh.numDofPerNode)
    q_faceLc = zeros(Complex128, mesh.numDofPerNode, mesh.numNodesPerFace)
    h = 1e-20
    pert = Complex128(0, h)
    for i=1:mesh.numBoundaryFaces
      bndry_i = mesh.bndryfaces[i]
      qL = sview(eqn.q, :, :, bndry_i.element)
#      resL = sview(res, :, :, bndry_i.element)
      fill!(q_faceLc, 0.0)
      fill!(res_jac, 0.0)

      boundaryFaceInterpolate!(mesh.sbpface, bndry_i.face, qL, q_faceLc)

      # compute flux jacobian at each node
      for j=1:mesh.numNodesPerFace
        q_j = sview(q_faceLc, :, j)
        for m=1:mesh.numDofPerNode
          q_j[m] += pert
    #      dxidx_j = sview(mesh.dxidx_bndry, :, :, j, i)

          # calculate boundary state
          coords = sview(mesh.coords_bndry, :, j, i)
          calcFreeStream(eqn.params_complex, coords, qg)

          # calculate face normal
          nrm2 = sview(mesh.nrm_bndry, :, j, i)

          # calculate lambda_max
          lambda_max = getLambdaMaxSimple(eqn.params_complex, q_j, qg, nrm2)

          # calculate dissipation
          for k=1:mesh.numDofPerNode
            flux_dotL[k, m, j] = lambda*imag(0.5*lambda_max*(q_j[k] - qg[k]))/h
          end

          q_j[m] -= pert
        end  # end loop m
      end  # end loop j

      
      boundaryFaceIntegrate_jac!(mesh.sbpface, bndry_i.face, flux_dotL, res_jac,
                               SummationByParts.Subtract())

      assembleBoundary(assembler, mesh.sbpface, mesh, bndry_i, res_jac)
    end  # end loop i

    fill!(res_jac, 0.0)
  end


  #---------------------------------------------------------------------------- 
  # shared face integrals
  # use q_faceL, q_faceR, lambda_dotL, lambda_dotR, flux_dotL, flux_dotR
  # from above

  workarr = zeros(q_faceR)
  for peer=1:mesh.npeers
    # get data for this peer
    interfaces_peer = mesh.shared_interfaces[peer]

    qR_peer = eqn.shared_data[peer].q_recv
#    dxidx_peer = mesh.dxidx_sharedface[peer]
    nrm_peer = mesh.nrm_sharedface[peer]
    start_elnum = mesh.shared_element_offsets[peer]

    for i=1:length(interfaces_peer)
      iface_i = interfaces_peer[i]
      qL_i = sview(eqn.q, :, :, iface_i.elementL)
      qR_i = sview(qR_peer, :, :, iface_i.elementR - start_elnum + 1)
      fill!(res_jacLL, 0.0)
      fill!(res_jacLR, 0.0)

      # interpolate to face
      interiorFaceInterpolate!(mesh.sbpface, iface_i, qL_i, qR_i, q_faceL, q_faceR)
      # compute flux at every face node
      for j=1:mesh.numNodesPerFace
        qL_j = sview(q_faceL, :, j)
        qR_j = sview(q_faceR, :, j)
        nrm2 = sview(nrm_peer, :, j, i)

        # get max wave speed
        lambda_max = getLambdaMaxSimple_diff(eqn.params, qL_j, qR_j, nrm2,
                                             lambda_dotL, lambda_dotR)

        # calculate flux
        for k=1:mesh.numDofPerNode
          # flux[k, j] = 0.5*lambda_max*(qL_j[k] - qR_j[k])
          for m=1:mesh.numDofPerNode
            flux_dotL[m, k, j] = 0.5*lambda_dotL[k]*(qL_j[m] - qR_j[m])
            flux_dotR[m, k, j] = 0.5*lambda_dotR[k]*(qL_j[m] - qR_j[m])
          end
          flux_dotL[k, k, j] += 0.5*lambda_max
          flux_dotR[k, k, j] -= 0.5*lambda_max
        end
      end  # end loop j

      # multiply by lambda here and it will get carried through
      # interiorFaceIntegrate_jac
      scale!(flux_dotL, lambda)
      scale!(flux_dotR, lambda)

      # compute dR/dq
      interiorFaceIntegrate_jac!(mesh.sbpface, iface_i, flux_dotL, flux_dotR,
                                res_jacLL, res_jacLR, res_jacRL, res_jacRR,
                                SummationByParts.Subtract())


     assembleSharedFace(assembler, mesh.sbpface, mesh, iface_i, res_jacLL, res_jacLR)
    end  # end loop i
  end  # end loop peer
  
  fill!(res_jacLL, 0.0)
  fill!(res_jacLR, 0.0)
  fill!(res_jacRL, 0.0)
  fill!(res_jacRR, 0.0)




  return nothing
end



"""
  Differentiated version of [`getLambdaMaxSimple`](@ref)

  **Inputs**

   * params
   * qL
   * qR
   * dir

  **Inputs/Outputs**

   * lambda_dotL: derivative of lambda wrt. qL
   * lambda_dotR: derivative of lambda wrt. qR
"""
function getLambdaMaxSimple_diff(params::ParamType{Tdim}, 
                      qL::AbstractVector{Tsol}, qR::AbstractVector{Tsol}, 
                      dir::AbstractVector{Tmsh},
                      lambda_dotL::AbstractVector{Tres},
                      lambda_dotR::AbstractVector{Tres}) where {Tsol, Tres, Tmsh, Tdim}

  q_avg = params.get_lambda_max_simple_data.q_avg

  for i=1:length(q_avg)
    q_avg[i] = 0.5*(qL[i] + qR[i])
  end

  lambda_max = getLambdaMax_diff(params, q_avg, dir, lambda_dotL)

  for i=1:length(lambda_dotL)
    lambda_dotL[i] *= 0.5
    lambda_dotR[i] = lambda_dotL[i]
  end

  return lambda_max
end


