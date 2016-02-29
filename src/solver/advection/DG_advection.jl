# DG.jl

function evalInteriorFlux{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator, 
                         eqn::AdvectionData{Tsol, Tres, Tdim}, opts)

  interfaces = view(mesh.interfaces,:)
  fluxface = zeros(Tsol, 1, sbp.numfacenodes, size(interfaces,1))

  for (findex, face) in enumerate(interfaces)
  	for j = 1:sbp.numfacenodes

      # Work on the Left side of interface
	    iL = sbp.facenodes[i, face.faceL]::Int
      q = eqn.q[1, iL, face.elementL]
      alpha_x = eqn.alpha_x[1, iL, face.elementL]
      alpha_y = eqn.alpha_y[ 1, iL, face.elementL]
      nrm = view(sbp.facenormal,:,face.faceL)
      dxidx = view(mesh.dxidx,:,:,iL, face.elementL)
      FL = calcNorm(q, alpha_x, alpha_y, nrm, dxidx)
      
      # Work on the Right side of interface
      iR = getnbrnodeindex(sbp, face, iL)
      q = eqn.q[1, iR, face.elementR]
      alpha_x = eqn.alpha_x[1, iR, face.elementR]
      alpha_y = eqn.alpha_y[1, iR, face.elementR]
      nrm = view(sbp.facenormal, :, face.faceR)
      dxidx = view(mesh.dxidx,:,:,iR,face.elementR)
      FR = calcNorm(q, alpha_x, alpha_y, nrm, dxidx)

      # Calculate the average flux across the interface
      fluxface[1,j,findex] = -0.5*(FL + FR) # negative sign because the way residual is calculated
                                            # volumeintegrate! has a +ve sign

    end # end for j = 1:sbp.numfacenodes
  end # end for (findex, face) in enumerate(interfaces)

  interiorfaceintegrate!(sbp, interfaces, fluxface, eqn.res)

  return nothing
end

function calcNormFlux(q, alpha_x, alpha_y, nrm, dxidx)

  alpha_xi = dxidx[1,1]*alpha_x + dxidx[1,2]*alpha_y
  alpha_eta = dxidx[2,1]*alpha_x + dxidx[2,2]*alpha_y
  alpha_n  = alpha_xi*nrm[1] + alpha_eta*nrm[2]

  flux = alpha_n*q

  return flux
end