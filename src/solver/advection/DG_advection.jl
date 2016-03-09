# DG.jl

function evalInteriorFlux{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh}, sbp::AbstractSBP, 
                         eqn::AdvectionData{Tsol, Tres, Tdim}, opts)

  fluxface = zeros(Tsol, 1, sbp.numfacenodes, size(mesh.interfaces,1))
  nbrnodeindex = Array(sbp.numfacenodes:-1:1)

  for (findex, face) in enumerate(mesh.interfaces)
  	for j = 1:sbp.numfacenodes
      
      # Work on the Left side of interface
	    iL = sbp.facenodes[j, face.faceL]::Int
      qL = eqn.q[1, iL, face.elementL]
      alpha_x = eqn.alpha_x[1, iL, face.elementL]
      alpha_y = eqn.alpha_y[ 1, iL, face.elementL]
      nrm = view(sbp.facenormal,:,face.faceL)
      dxidx = view(mesh.dxidx,:,:,iL, face.elementL)
      
      # Work on the Right side of interface
      iR = sbp.facenodes[nbrnodeindex[j], face.faceR]::Int
      qR = eqn.q[1, iR, face.elementR]

      # Calculate the upwind flux across the interface. negative sign because 
      # the way residual is calculated volumeintegrate! has a +ve sign
      fluxface[1,j,findex] = RoeSolver(qR, qL, alpha_x, alpha_y, nrm, dxidx)
      # println("fluxface[1,j,findex] = ", fluxface[1,j,findex]) 
    
    end # end for j = 1:sbp.numfacenodes
  end # end for (findex, face) in enumerate(interfaces)
  
  # println("fluxface = " fluxface)
  interiorfaceintegrate!(sbp, mesh.interfaces, fluxface, eqn.res)

  return nothing
end


#=
function calcNormFlux(q, alpha_x, alpha_y, nrm, dxidx)

  alpha_xi = dxidx[1,1]*alpha_x + dxidx[1,2]*alpha_y
  alpha_eta = dxidx[2,1]*alpha_x + dxidx[2,2]*alpha_y
  alpha_n  = alpha_xi*nrm[1] + alpha_eta*nrm[2]

  flux = alpha_n*q

  return flux
end =#