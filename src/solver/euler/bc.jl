export getBCFunctors


# Euler flux calculator used by isentropicVortexBC ONLY!!

#export rho1Energy2BC, isentropicVortexBC

#=
# this function no longer works
function rho1Energy2BC(q, x, dxidx, nrm)

  E1dq = zeros(Float64, 4)
  E2dq = zeros(Float64, 4)

  #=
  println("q: ",q)
  println("x: ",x)
  println("dxidx: ",dxidx)
  println("nrm: ",nrm)
  =#

  # getting qg
  qg = zeros(Float64, 4)
#  calcRho1Energy2(x, eqn, qg)
   calcRho1Energy2U3(x, eqn, qg)
  
#   println("qg: ",qg)

  # Declaring constants 
  d1_0 = 1.0
  d0_0 = 0.0
  d0_5 = 0.5
  tau = 1.0
  sgn = 1.0
  gamma = 1.4
  gami = gamma - 1
  sat_Vn = 0.025
  sat_Vl = 0.025

  # Begin main executuion
  nx = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
  ny = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]

#  println("nrm: ",nrm)
#  println("nx: ",nx)
#  println("ny: ",ny)

  dA = sqrt(nx*nx + ny*ny)
  
  fac = d1_0/q[1]
#   println(typeof(fac))
#   println(typeof(q[4]))
  uL = q[2]*fac; vL = q[3]*fac;
  phi = d0_5*(uL*uL + vL*vL)

  HL = gamma*q[4]*fac - gami*phi
  
  fac = d1_0/qg[1]
  uR = qg[2]*fac; vR = qg[3]*fac;
  phi = d0_5*(uR*uR + vR*vR)
  HR = gamma*qg[4]*fac - gami*phi

  sqL = sqrt(q[1]); sqR = sqrt(qg[1])
  fac = d1_0/(sqL + sqR)
  u = (sqL*uL + sqR*uR)*fac
  v = (sqL*vL + sqR*vR)*fac
  
  H = (sqL*HL + sqR*HR)*fac
  phi = d0_5*(u*u + v*v)
  
  a = sqrt(gami*(H - phi))
  Un = u*nx + v*ny

  lambda1 = Un + dA*a
  lambda2 = Un - dA*a
  lambda3 = Un
  rhoA = absvalue(Un) + dA*a
  lambda1 = d0_5*(tau*max(absvalue(lambda1),sat_Vn *rhoA) + sgn*lambda1)
  lambda2 = d0_5*(tau*max(absvalue(lambda2),sat_Vn *rhoA) + sgn*lambda2)
  lambda3 = d0_5*(tau*max(absvalue(lambda3),sat_Vl *rhoA) + sgn*lambda3)

  dq1 = q[1] - qg[1] 
  dq2 = q[2] - qg[2]
  dq3 = q[3] - qg[3]
  dq4 = q[4] - qg[4]

  #-- diagonal matrix multiply
  sat = zeros(Float64, 4)
  sat[1] = lambda3*dq1
  sat[2] = lambda3*dq2
  sat[3] = lambda3*dq3
  sat[4] = lambda3*dq4

#   println("sat 1: ",sat)

  #-- get E1*dq
  E1dq[1] = phi*dq1 - u*dq2 - v*dq3 + dq3
  E1dq[2] = E1dq[1]*u
  E1dq[3] = E1dq[1]*v
  E1dq[4] = E1dq[1]*H

  #-- get E2*dq
  E2dq[1] = d0_0
  E2dq[2] = -Un*dq1 + nx*dq2 + ny*dq3
  E2dq[3] = E2dq[2]*ny
  E2dq[4] = E2dq[2]*Un
  E2dq[2] = E2dq[2]*nx

  #-- add to sat
  tmp1 = d0_5*(lambda1 + lambda2) - lambda3
  tmp2 = gami/(a*a)
  tmp3 = d1_0/(dA*dA)
  sat[:] = sat[:] + tmp1*(tmp2*E1dq[:] + tmp3*E2dq[:])
#   println("sat 2: ",sat)
  
  #-- get E3*dq
  E1dq[1] = -Un*dq1 + nx*dq2 + ny*dq3
  E1dq[2] = E1dq[1]*u
  E1dq[3] = E1dq[1]*v
  E1dq[4] = E1dq[1]*H

  #-- get E4*dq
  E2dq[1] = d0_0
  E2dq[2] = phi*dq1 - u*dq2 - v*dq3 + dq4
  E2dq[3] = E2dq[2]*ny
  E2dq[4] = E2dq[2]*Un
  E2dq[2] = E2dq[2]*nx

  #-- add to sat
  tmp1 = d0_5*(lambda1 - lambda2)/(dA*a)
  sat[:] = sat[:] + tmp1*(E1dq[:] + gami*E2dq[:])
#   println("sat 3: ",sat)
  
  returnval = sat + getEulerFlux(q, nx, ny)
#   println("returnval: ", returnval)
  return returnval

end # ends the function eulerRoeSAT

# Euler Roe Solver for boundary integrate
#function isentropicVortexBC{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator{Tsbp}, eqn::EulerEquation{Tsol})
function isentropicVortexBC{T}(q::AbstractArray{T,1}, x::AbstractArray{T,1}, dxidx::AbstractArray{T,2}, nrm::AbstractArray{T,1}, flux::AbstractArray{T, 1}, mesh::AbstractMesh, eqn::EulerEquation)

  E1dq = zeros(Float64, 4)
  E2dq = zeros(Float64, 4)

  # getting qg
  qg = zeros(Float64, 4)
  calcIsentropicVortex(x, eqn, qg)
#  calcVortex(x, eqn, qg)

  # Declaring constants 
  d1_0 = 1.0
  d0_0 = 0.0
  d0_5 = 0.5
  tau = 1.0
  sgn = -1.0
  gamma = 1.4
  gami = gamma - 1
  sat_Vn = 0.025
  sat_Vl = 0.025

  # Begin main executuion
  nx = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
  ny = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]

  dA = sqrt(nx*nx + ny*ny)
  
  fac = d1_0/q[1]
#   println(typeof(fac))
#   println(typeof(q[4]))
  uL = q[2]*fac; vL = q[3]*fac;
  phi = d0_5*(uL*uL + vL*vL)

  HL = gamma*q[4]*fac - gami*phi
  
  fac = d1_0/qg[1]
  uR = qg[2]*fac; vR = qg[3]*fac;
  phi = d0_5*(uR*uR + vR*vR)
  HR = gamma*qg[4]*fac - gami*phi

  sqL = sqrt(q[1]); sqR = sqrt(qg[1])
  fac = d1_0/(sqL + sqR)
  u = (sqL*uL + sqR*uR)*fac
  v = (sqL*vL + sqR*vR)*fac
  
  H = (sqL*HL + sqR*HR)*fac
  phi = d0_5*(u*u + v*v)
  
  a = sqrt(gami*(H - phi))
  Un = u*nx + v*ny

  lambda1 = Un + dA*a
  lambda2 = Un - dA*a
  lambda3 = Un
  rhoA = absvalue(Un) + dA*a
  lambda1 = d0_5*(tau*max(absvalue(lambda1),sat_Vn *rhoA) + sgn*lambda1)
  lambda2 = d0_5*(tau*max(absvalue(lambda2),sat_Vn *rhoA) + sgn*lambda2)
  lambda3 = d0_5*(tau*max(absvalue(lambda3),sat_Vl *rhoA) + sgn*lambda3)

  dq1 = q[1] - qg[1] 
  dq2 = q[2] - qg[2]
  dq3 = q[3] - qg[3]
  dq4 = q[4] - qg[4]

  #-- diagonal matrix multiply
  sat = zeros(Float64, 4)
  sat[1] = lambda3*dq1
  sat[2] = lambda3*dq2
  sat[3] = lambda3*dq3
  sat[4] = lambda3*dq4

  #-- get E1*dq
  E1dq[1] = phi*dq1 - u*dq2 - v*dq3 + dq4
  E1dq[2] = E1dq[1]*u
  E1dq[3] = E1dq[1]*v
  E1dq[4] = E1dq[1]*H

  #-- get E2*dq
  E2dq[1] = d0_0
  E2dq[2] = -Un*dq1 + nx*dq2 + ny*dq3
  E2dq[3] = E2dq[2]*ny
  E2dq[4] = E2dq[2]*Un
  E2dq[2] = E2dq[2]*nx

  #-- add to sat
  tmp1 = d0_5*(lambda1 + lambda2) - lambda3
  tmp2 = gami/(a*a)
  tmp3 = d1_0/(dA*dA)
  sat[:] = sat[:] + tmp1*(tmp2*E1dq[:] + tmp3*E2dq[:])
  
  #-- get E3*dq
  E1dq[1] = -Un*dq1 + nx*dq2 + ny*dq3
  E1dq[2] = E1dq[1]*u
  E1dq[3] = E1dq[1]*v
  E1dq[4] = E1dq[1]*H

  #-- get E4*dq
  E2dq[1] = d0_0
  E2dq[2] = phi*dq1 - u*dq2 - v*dq3 + dq4
  E2dq[3] = E2dq[2]*ny
  E2dq[4] = E2dq[2]*Un
  E2dq[2] = E2dq[2]*nx

  #-- add to sat
  tmp1 = d0_5*(lambda1 - lambda2)/(dA*a)
  sat[:] = sat[:] + tmp1*(E1dq[:] + gami*E2dq[:])

  euler_flux = zeros(4)
  calcEulerFlux(eqn, q, [nx, ny], euler_flux)

#  flux[:] = sat + getEulerFlux(q, nx, ny, eqn)
  flux[:] = -(sat + euler_flux)
 
#  return sat + getEulerFlux(q, nx, ny)
   return nothing

end # ends the function eulerRoeSAT


function isentropicVortexBC{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator{Tsbp}, eqn::EulerEquation{Tsol})


#function isentropicVortexBC{Tmsh, Tsbp}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator{Tsbp}, q::AbstractArray, coords::AbstractArray, dxidx::AbstractArray, flux::AbstractArray)


  # q = conserv. vars (vec)
  # x = coords (vec)
  # dxidx = 2x2 matrix
  # nrm = vew(sbp.facenormal, :, bndry.face
  # flux = returned value
  E1dq = zeros(Tsol, 4)
  E2dq = zeros(Tsol, 4)
  sat = zeros(Tsol, 4)

  # getting qg
  qg = zeros(Tsol, 4)
  euler_flux = zeros(Tsol, 4)
#  calcIsentropicVortex(x, eqn, qg)
#  calcVortex(x, eqn, qg)

  # Declaring constants 
  d1_0 = 1.0
  d0_0 = 0.0
  d0_5 = 0.5
  tau = 1.0
  sgn = -1.0
  gamma = 1.4
  gami = gamma - 1
  sat_Vn = 0.025
  sat_Vl = 0.025

  for i=1:mesh.numBoundaryEdges
    bndry_i = mesh.bndryfaces[i]
    for j = 1:sbp.numfacenodes
      k = sbp.facenodes[j, bndry_i.face]

      q = view(eqn.q, :, k, bndry_i.element)
      x = view(mesh.coords, :, k, bndry_i.element)
      dxidx = view(mesh.dxidx, :, :, k, bndry_i.element)
      nrm = view(sbp.facenormal, :, bndry_i.face)
      #println("eqn.bndryflux = ", eqn.bndryflux)
      flux = view(eqn.bndryflux, :, i, j)

      # Begin main executuion
      # get qg
      calcIsentropicVortex(x, eqn, qg)

      nx = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
      ny = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]

      dA = sqrt(nx*nx + ny*ny)
      
      fac = d1_0/q[1]
    #   println(typeof(fac))
    #   println(typeof(q[4]))
      uL = q[2]*fac; vL = q[3]*fac;
      phi = d0_5*(uL*uL + vL*vL)

      HL = gamma*q[4]*fac - gami*phi
      
      fac = d1_0/qg[1]
      uR = qg[2]*fac; vR = qg[3]*fac;
      phi = d0_5*(uR*uR + vR*vR)
      HR = gamma*qg[4]*fac - gami*phi

      sqL = sqrt(q[1]); sqR = sqrt(qg[1])
      fac = d1_0/(sqL + sqR)
      u = (sqL*uL + sqR*uR)*fac
      v = (sqL*vL + sqR*vR)*fac
      
      H = (sqL*HL + sqR*HR)*fac
      phi = d0_5*(u*u + v*v)
      
      a = sqrt(gami*(H - phi))
      Un = u*nx + v*ny

      lambda1 = Un + dA*a
      lambda2 = Un - dA*a
      lambda3 = Un
      rhoA = absvalue(Un) + dA*a
      lambda1 = d0_5*(tau*max(absvalue(lambda1),sat_Vn *rhoA) + sgn*lambda1)
      lambda2 = d0_5*(tau*max(absvalue(lambda2),sat_Vn *rhoA) + sgn*lambda2)
      lambda3 = d0_5*(tau*max(absvalue(lambda3),sat_Vl *rhoA) + sgn*lambda3)

      dq1 = q[1] - qg[1] 
      dq2 = q[2] - qg[2]
      dq3 = q[3] - qg[3]
      dq4 = q[4] - qg[4]

      #-- diagonal matrix multiply
      sat[1] = lambda3*dq1
      sat[2] = lambda3*dq2
      sat[3] = lambda3*dq3
      sat[4] = lambda3*dq4

      #-- get E1*dq
      E1dq[1] = phi*dq1 - u*dq2 - v*dq3 + dq4
      E1dq[2] = E1dq[1]*u
      E1dq[3] = E1dq[1]*v
      E1dq[4] = E1dq[1]*H

      #-- get E2*dq
      E2dq[1] = d0_0
      E2dq[2] = -Un*dq1 + nx*dq2 + ny*dq3
      E2dq[3] = E2dq[2]*ny
      E2dq[4] = E2dq[2]*Un
      E2dq[2] = E2dq[2]*nx

      #-- add to sat
      tmp1 = d0_5*(lambda1 + lambda2) - lambda3
      tmp2 = gami/(a*a)
      tmp3 = d1_0/(dA*dA)
      sat[:] = sat[:] + tmp1*(tmp2*E1dq[:] + tmp3*E2dq[:])
      
      #-- get E3*dq
      E1dq[1] = -Un*dq1 + nx*dq2 + ny*dq3
      E1dq[2] = E1dq[1]*u
      E1dq[3] = E1dq[1]*v
      E1dq[4] = E1dq[1]*H

      #-- get E4*dq
      E2dq[1] = d0_0
      E2dq[2] = phi*dq1 - u*dq2 - v*dq3 + dq4
      E2dq[3] = E2dq[2]*ny
      E2dq[4] = E2dq[2]*Un
      E2dq[2] = E2dq[2]*nx

      #-- add to sat
      tmp1 = d0_5*(lambda1 - lambda2)/(dA*a)
      sat[:] = sat[:] + tmp1*(E1dq[:] + gami*E2dq[:])

      calcEulerFlux(eqn, q, [nx, ny], euler_flux)

    #  flux[:] = sat + calcEulerFlux(q, nx, ny, eqn)a
      # make flux negative until boundaryintegrate! supports -=
      flux = -(sat + euler_flux)
    end
  end
     
    #  return sat + getEulerFlux(q, nx, ny)
       return nothing

end # ends the function eulerRoeSAT

=#



# mid level function
function calcBoundaryFlux{Tmsh, Tsbp, Tsol, Tres}( mesh::AbstractMesh{Tmsh}, sbp::SBPOperator{Tsbp}, eqn::EulerEquation{Tsol}, functor::BCType, bndry_facenums::AbstractArray{Int,1}, flux::AbstractArray{Tres, 3})
# calculate the boundary flux for the boundary condition evaluated by the functor

  nfaces = length(bndry_facenums)

  for i=1:n  # loop over faces with this BC
    bndry_i = bndry_facenums[i]
    for j = 1:sbp.numfacenodes
      k = sbp.facenodes[j, bndry_i.face]

      # get components
      q = view(eqn.q, :, k, bndry_i.element)
      x = view(mesh.coords, :, k, bndry_i.element)
      dxidx = view(mesh.dxidx, :, :, k, bndry_i.element)
      nrm = view(sbp.facenormal, :, bndry_i.face)
      #println("eqn.bndryflux = ", eqn.bndryflux)
      flux = view(flux, :, j, i)

      functor(q, x, dxidx, nrm, flux, eqn)
    end
  end


  return nothing
end





# mid level function
# no longer needed
#=
function getIsentropicVortexBoundaryFlux{Tmsh, Tsbp, Tsol}(mesh::AbstractMesh{Tmsh}, sbp::SBPOperator{Tsbp}, eqn::EulerEquation{Tsol})


  for i=1:mesh.numBoundaryEdges
    bndry_i = mesh.bndryfaces[i]
    for j = 1:sbp.numfacenodes
      k = sbp.facenodes[j, bndry_i.face]

      # get components
      q = view(eqn.q, :, k, bndry_i.element)
      x = view(mesh.coords, :, k, bndry_i.element)
      dxidx = view(mesh.dxidx, :, :, k, bndry_i.element)
      nrm = view(sbp.facenormal, :, bndry_i.face)
      #println("eqn.bndryflux = ", eqn.bndryflux)
      flux = view(eqn.bndryflux, :, j, i)

      isentropicVortexBC(q, x, dxidx, nrm, flux, eqn)
    end
  end

  return nothing
end
=#



type isentropicVortexBC <: BCType
end

# low level function
function call{Tmsh, Tsol, Tres}(obj::isentropicVortexBC, q::AbstractArray{Tsol,1}, x::AbstractArray{Tmsh,1}, dxidx::AbstractArray{Tmsh,2}, nrm::AbstractArray{Tmsh,1}, flux::AbstractArray{Tres, 1}, eqn::EulerEquation{Tsol, 2})

  E1dq = zeros(Tres, 4)
  E2dq = zeros(Tres, 4)

#  println("entered isentropicOvrtexBC (low level)")
#  println("Tsol = ", Tsol)

  # getting qg
  qg = zeros(Tsol, 4)
  calcIsentropicVortex(x, eqn, qg)
#  calcVortex(x, eqn, qg)

  # Declaring constants 
  d1_0 = 1.0
  d0_0 = 0.0
  d0_5 = 0.5
  tau = 1.0
  sgn = -1.0
  gamma = 1.4
  gami = gamma - 1
  sat_Vn = convert(Tsol, 0.025)
  sat_Vl = convert(Tsol, 0.025)

  # Begin main executuion
  nx = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
  ny = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]

  dA = sqrt(nx*nx + ny*ny)
  
  fac = d1_0/q[1]
#   println(typeof(fac))
#   println(typeof(q[4]))
  uL = q[2]*fac; vL = q[3]*fac;
  phi = d0_5*(uL*uL + vL*vL)

  HL = gamma*q[4]*fac - gami*phi
  
  fac = d1_0/qg[1]
  uR = qg[2]*fac; vR = qg[3]*fac;
  phi = d0_5*(uR*uR + vR*vR)
  HR = gamma*qg[4]*fac - gami*phi

  sqL = sqrt(q[1]); sqR = sqrt(qg[1])
  fac = d1_0/(sqL + sqR)
  u = (sqL*uL + sqR*uR)*fac
  v = (sqL*vL + sqR*vR)*fac
  
  H = (sqL*HL + sqR*HR)*fac
  phi = d0_5*(u*u + v*v)
  
  a = sqrt(gami*(H - phi))
  Un = u*nx + v*ny

  lambda1 = Un + dA*a
  lambda2 = Un - dA*a
  lambda3 = Un
  rhoA = absvalue(Un) + dA*a

#  println("sat_Vn = ", sat_Vn)
#  println("lambda1 = ", lambda1)
#  println("absvalue(lambda1) = ", absvalue(lambda1))
  lambda1 = d0_5*(tau*max(absvalue(lambda1),sat_Vn *rhoA) + sgn*lambda1)
  lambda2 = d0_5*(tau*max(absvalue(lambda2),sat_Vn *rhoA) + sgn*lambda2)
  lambda3 = d0_5*(tau*max(absvalue(lambda3),sat_Vl *rhoA) + sgn*lambda3)

  dq1 = q[1] - qg[1] 
  dq2 = q[2] - qg[2]
  dq3 = q[3] - qg[3]
  dq4 = q[4] - qg[4]

  #-- diagonal matrix multiply
  sat = zeros(Tres, 4)
  sat[1] = lambda3*dq1
  sat[2] = lambda3*dq2
  sat[3] = lambda3*dq3
  sat[4] = lambda3*dq4

  #-- get E1*dq
  E1dq[1] = phi*dq1 - u*dq2 - v*dq3 + dq4
  E1dq[2] = E1dq[1]*u
  E1dq[3] = E1dq[1]*v
  E1dq[4] = E1dq[1]*H

  #-- get E2*dq
  E2dq[1] = d0_0
  E2dq[2] = -Un*dq1 + nx*dq2 + ny*dq3
  E2dq[3] = E2dq[2]*ny
  E2dq[4] = E2dq[2]*Un
  E2dq[2] = E2dq[2]*nx

  #-- add to sat
  tmp1 = d0_5*(lambda1 + lambda2) - lambda3
  tmp2 = gami/(a*a)
  tmp3 = d1_0/(dA*dA)
  sat[:] = sat[:] + tmp1*(tmp2*E1dq[:] + tmp3*E2dq[:])
  
  #-- get E3*dq
  E1dq[1] = -Un*dq1 + nx*dq2 + ny*dq3
  E1dq[2] = E1dq[1]*u
  E1dq[3] = E1dq[1]*v
  E1dq[4] = E1dq[1]*H

  #-- get E4*dq
  E2dq[1] = d0_0
  E2dq[2] = phi*dq1 - u*dq2 - v*dq3 + dq4
  E2dq[3] = E2dq[2]*ny
  E2dq[4] = E2dq[2]*Un
  E2dq[2] = E2dq[2]*nx

  #-- add to sat
  tmp1 = d0_5*(lambda1 - lambda2)/(dA*a)
  sat[:] = sat[:] + tmp1*(E1dq[:] + gami*E2dq[:])

  euler_flux = zeros(Tsol, 4)
  calcEulerFlux(eqn, q, [nx, ny], euler_flux)

#  flux[:] = sat + getEulerFlux(q, nx, ny, eqn)
#  flux[:] = -(sat + euler_flux)
  for i=1:4  # ArrayViews does not support flux[:] = ...
    flux[i] = -(sat[i] + euler_flux[i])
  end
 
#  return sat + getEulerFlux(q, nx, ny)
   return nothing

end # ends the function eulerRoeSAT



# every time a new boundary condition is created,
# add it to the dictionary
const isentropicVortexBC_ = isentropicVortexBC()
global const BCDict = Dict{ASCIIString, BCType} (
"isentropicVortexBC" => isentropicVortexBC_
)


function getBCFunctors(mesh::PumiMesh, sbp::SBPOperator, eqn::EulerEquation, opts)
# populate the array mesh.bndry_funcs with the functors for the boundary condition types

println("Entered getBCFunctors")
println("BCDict = ", BCDict)

for i=1:mesh.numBC
  key_i = string("BC", i, "_name")
  val = opts[key_i]
  println("BCDict[val] = ", BCDict[val])
  mesh.bndry_funcs[i] = BCDict[val]
end

return nothing

end

