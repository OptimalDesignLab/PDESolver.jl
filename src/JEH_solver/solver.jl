module Solver
# simple Euler solver

include("mesh2d.jl")

using SummationByParts
using .Mesh2D

export EulerSolver

const γ = 1.4

type EulerSolver{T}
  Nx::Int
  Ny::Int
  mesh::Mesh{T}
  sbp::TriSBP{T}
  massinv::Array{T}
  x::Array{T}
  dξdx::Array{T}
  α::Array{T}
  jac::Array{T}
  edgestab::Bool
  exact::Function

  function EulerSolver(degree::Int, Nx::Int, Ny::Int; bubble::Int=-1, 
                           edgestab::Bool=true)
    sbp = TriSBP{T}(degree=degree, bubble=bubble)
    mesh = Mesh{T}(Nx, Nx, sbp, 1.0, 1.0)
    dξdx = zeros(T, (2, 2, sbp.numnodes, mesh.numelem))
    α = zeros(dξdx)
    jac = zeros(T, (sbp.numnodes, mesh.numelem))
    x = zeros(T, (2, sbp.numnodes, mesh.numelem))
    for k = 1:mesh.numelem
      for i = 1:sbp.numnodes
        x[:,i,k] = mesh.x[:,mesh.elemnodes[i,k]]
      end
    end
    mappingjacobian!(sbp, x, dξdx, jac)
    if edgestab
      for k = 1:mesh.numelem
        for i = 1:sbp.numnodes
          for di1 = 1:2
            for di2 = 1:2
              α[di1,di2,i,k] = (dξdx[di1,1,i,k].*dξdx[di2,1,i,k] + 
                                dξdx[di1,2,i,k].*dξdx[di2,2,i,k])*jac[i,k]
            end
          end
        end
      end
    end
    massinv = Solver.calcinvmassmatrix!(sbp, mesh, jac, Nx, Ny)
    new(Nx, Ny, mesh, sbp, massinv, x, dξdx, α, jac, edgestab)
  end
end

@doc """
Builds the inverse of the mass matrix and returns it as a 1D array

**Inputs**

* `sbp`: An triangular SBP type
* `mesh`: a simple 2D mesh 
* `jac`: the determinant of the mapping Jacobian, stored in (node,elem) format
* `Nx`,`Ny`: number of edges in the x and y directions

**Returns**
* `massinv` : the inverse mass matrix, stored as a 1D array

"""->
function calcinvmassmatrix!{T}(sbp::TriSBP{T}, mesh::Mesh{T}, jac::Array{T},
                               Nx::Int, Ny::Int)
  massinv = zeros(T, (mesh.numnodes))
  masselem = zeros(T, (sbp.numnodes, mesh.numelem))
  for k = 1:mesh.numelem
    for i = 1:sbp.numnodes
      masselem[i,k] = sbp.w[i]/jac[i,k]
    end
  end

  # scatter element mass matrices to global mass matrix
  fill!(massinv, 0.0)
  for k = 1:mesh.numelem
    for i = 1:sbp.numnodes
      massinv[mesh.elemnodes[i,k]] += masselem[i,k]
    end
  end

  # invert the matrix
  for k = 1:mesh.numnodes
    massinv[k] = 1.0/massinv[k]
  end
  return massinv
end

@doc """
Returns the Euler flux in the direction dξ
"""->
function EulerFlux(q, dξ)
  F = zeros(q)
  press = (γ-1)*(q[4] - 0.5*(q[2]^2 + q[3]^2)/q[1])
  U = (q[2]*dξ[1] + q[3]*dξ[2])/q[1]
  F[1] = q[1]*U
  F[2] = q[2]*U + dξ[1]*press
  F[3] = q[3]*U + dξ[2]*press
  F[4] = (q[4] + press)*U
  return F
end

function calcIsentropicVortex!{T}(coords::Array{T}, sol::Array{T})
  x = coords[1]
  y = coords[2]

  gamma = γ
  gami = gamma-1.0

  ri = 1.0
  Mai = 0.95 
  rhoi = 2.0
  prsi = 1.0/gamma    

  rinv = ri/sqrt(x*x + y*y)
  rho = rhoi*(1.0 + 0.5*gami*Mai*Mai*(1.0 - rinv*rinv))^(1.0/gami)
  Ma = sqrt((2.0/gami)*( ( (rhoi/rho)^(gami) ) * 
                         (1.0 + 0.5*gami*Mai*Mai) - 1.0 ) )
  theta = atan2(y,x)
  press = prsi* ( (1.0 + 0.5*gami*Mai*Mai) / 
                 (1.0 + 0.5*gami*Ma*Ma) )^(gamma/gami)
  a = sqrt(gamma*press/rho)
        
  sol[1] = rho
  sol[2] = -rho*a*Ma*sin(theta)
  sol[3] = rho*a*Ma*cos(theta)
  sol[4] = press/gami + 0.5*rho*a*a*Ma*Ma
end

function isentropicVortexBC(q, x, dxidx, nrm)
  E1dq = zeros(Float64, 4)
  E2dq = zeros(Float64, 4)

  # getting qg
  qg = zeros(Float64, 4)
  calcIsentropicVortex!(x, qg)

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
  rhoA = abs(Un) + dA*a
  lambda1 = d0_5*(max(abs(lambda1),sat_Vn *rhoA) + sgn*lambda1)
  lambda2 = d0_5*(max(abs(lambda2),sat_Vn *rhoA) + sgn*lambda2)
  lambda3 = d0_5*(max(abs(lambda3),sat_Vl *rhoA) + sgn*lambda3)

  dq1 = q[1] - qg[1] 
  dq2 = q[2] - qg[2]
  dq3 = q[3] - qg[3]
  dq4 = q[4] - qg[4]

  # diagonal matrix multiply
  sat = zeros(Float64, 4)
  sat[1] = lambda3*dq1
  sat[2] = lambda3*dq2
  sat[3] = lambda3*dq3
  sat[4] = lambda3*dq4

  # get E1*dq
  E1dq[1] = phi*dq1 - u*dq2 - v*dq3 + dq4
  E1dq[2] = E1dq[1]*u
  E1dq[3] = E1dq[1]*v
  E1dq[4] = E1dq[1]*H

  # get E2*dq
  E2dq[1] = d0_0
  E2dq[2] = -Un*dq1 + nx*dq2 + ny*dq3
  E2dq[3] = E2dq[2]*ny
  E2dq[4] = E2dq[2]*Un
  E2dq[2] = E2dq[2]*nx

  # add to sat
  tmp1 = d0_5*(lambda1 + lambda2) - lambda3
  tmp2 = gami/(a*a)
  tmp3 = d1_0/(dA*dA)
  sat[:] = sat[:] + tmp1*(tmp2*E1dq[:] + tmp3*E2dq[:])
  
  # get E3*dq
  E1dq[1] = -Un*dq1 + nx*dq2 + ny*dq3
  E1dq[2] = E1dq[1]*u
  E1dq[3] = E1dq[1]*v
  E1dq[4] = E1dq[1]*H

  # get E4*dq
  E2dq[1] = d0_0
  E2dq[2] = phi*dq1 - u*dq2 - v*dq3 + dq4
  E2dq[3] = E2dq[2]*ny
  E2dq[4] = E2dq[2]*Un
  E2dq[2] = E2dq[2]*nx

  # add to sat
  tmp1 = d0_5*(lambda1 - lambda2)/(dA*a)
  sat[:] = sat[:] + tmp1*(E1dq[:] + gami*E2dq[:])
  
  return sat + EulerFlux(q, [nx; ny])
  #return EulerFlux(q, [nx; ny])
end

@doc """
Computes the Euler equation residual.

**Inputs**

* `solver`: the EulerSolver type that defines the problem
* `u`: the current solution, stored as a 2D array

**Outputs**
* `res` : the residual, stored as a 2D array (same size as u)

"""->
function calcweakresidual!{T}(solver::EulerSolver{T}, u::Array{T}, res::Array{T})
  @assert( size(u,2) == solver.mesh.numnodes )
  @assert( length(u) == length(res) )

  function stabscale(u, dξdx, nrm)
    h = 1.0/solver.mesh.Nx # <--- this is not general
#     return 0.01*h^(solver.sbp.degree+1)
    return 0.025*h^(solver.sbp.degree+1)
  end

  # gather u onto elements
  reselem = zeros(T, (4, solver.sbp.numnodes, solver.mesh.numelem))
  uelem = zeros(reselem)
  for k = 1:solver.mesh.numelem
    for i = 1:solver.sbp.numnodes
      uelem[:,i,k] = u[:,solver.mesh.elemnodes[i,k]]
    end
  end

  # calc fluxes
  Fξ = zeros(T, (4, solver.sbp.numnodes, solver.mesh.numelem))
  Fη = zeros(Fξ)
  for k = 1:solver.mesh.numelem
    for i = 1:solver.sbp.numnodes
      unode = u[:,solver.mesh.elemnodes[i,k]]
      dξ = solver.dξdx[:,:,i,k]
      Fξ[:,i,k] = EulerFlux(unode, dξ[1,:])
      Fη[:,i,k] = EulerFlux(unode, dξ[2,:])
    end
  end

  # apply the stiffness matrix operators to Fx and Fy
  weakdifferentiate!(solver.sbp, 1, Fξ, reselem, trans=true)
  weakdifferentiate!(solver.sbp, 2, Fη, reselem, trans=true)
  reselem *= -1.0

  # apply the boundary conditions weakly
  boundaryintegrate!(solver.sbp, solver.mesh.bndryfaces, uelem, solver.x,
                     solver.dξdx, isentropicVortexBC, reselem)

  # apply edge stabilization
  edgestabilize!(solver.sbp, solver.mesh.ifaces, uelem, solver.x, solver.dξdx,
                 solver.jac, solver.α, stabscale, reselem)

  # scatter element data back to global residual
  fill!(res, 0.0)
  for k = 1:solver.mesh.numelem
    for i = 1:solver.sbp.numnodes
      res[:,solver.mesh.elemnodes[i,k]] += reselem[:,i,k]
    end
  end

  # apply the inverse of the global mass matrix
  for k = 1:solver.mesh.numnodes
    res[:,k] *= solver.massinv[k]
  end
end

@doc """
Apply one step of RK4 to a vector based on the problem defined in solver

**Inputs**

* `solver`: the AdvectionSolver type that defines the problem and mesh
* `dt`: the time step size

**InOuts**

* `u`: solution vector being updated

**Returns**

* `resnorm`: the residual norm of the first stage

"""->
function timemarchRK!{T}(solver::EulerSolver{T}, u::Array{T}, dt::T)
  @assert(dt > 0.0)

  # 1st stage
  k1 = zeros(u)
  Solver.calcweakresidual!(solver, u, k1)
  k1 *= -1.0
  u1 = deepcopy(u)
  u1 += (0.5*dt).*k1

  # 2nd stage
  k2 = zeros(u)
  Solver.calcweakresidual!(solver, u1, k2)
  k2 *= -1.0
  u2 = deepcopy(u)
  u2 += (0.5*dt).*k2

  # 3rd stage
  k3 = zeros(u)
  Solver.calcweakresidual!(solver, u2, k3)
  k3 *= -1.0
  u3 = deepcopy(u)
  u3 += dt.*k3

  # 4th stage
  k4 = zeros(u)
  Solver.calcweakresidual!(solver, u3, k4)
  k4 *= -1.0

  # update solution
  for i = 1:solver.mesh.numnodes
    u[:,i] += (dt/6.0)*(k1[:,i] + 2.0*k2[:,i] + 2.0*k3[:,i] + k4[:,i])
  end
  return norm(k1[1,:])
end

@doc """
Prints the solution and nodes to a file

**Inputs**

* `solver`: the AdvectionSolver type that defines the problem and mesh
* `u`: solution vector being printed to file
* `filename` (optional): the name of the file to write to

"""->
function printsolution{T}(solver::EulerSolver{T}, u::Array,
                          filename::AbstractString="euler.dat")
  uelem = zeros(T, (4,solver.sbp.numnodes, solver.mesh.numelem))
  for k = 1:solver.mesh.numelem
    for i = 1:solver.sbp.numnodes
      uelem[:,i,k] = u[:,solver.mesh.elemnodes[i,k]]
    end
  end
  f = open(filename, "w")
  for x in solver.x[1,:,:]
    print(f, x, " ")
  end
  println(f)
  for y in solver.x[2,:,:]
    print(f, y, " ")
  end
  println(f)
  for ρ in uelem[1,:,:]
    print(f, ρ, " ")
  end
  println(f)
  for ρu in uelem[2,:,:]
    print(f, ρu, " ")
  end
  println(f)
  close(f)
end

@doc """
Computes the L2 and max error

**Inputs**

* `solver`: the EulerSolver type that defines the problem and mesh
* `u`: solution vector whose error is begin computed

**Returns**

* `L2error`: the L2 error computed using the SBP norm
* `maxerror`: the max error over the nodes

"""->
function calcerror{T}(solver::EulerSolver{T}, u::Array)
  L2error = 0.0
  maxerror = 0.0
  for k = 1:solver.mesh.numelem
    for i = 1:solver.sbp.numnodes
      qexact = zeros(Float64, 4)
      calcIsentropicVortex!(solver.x[:,i,k], qexact)
      nodeerror = norm(u[:,solver.mesh.elemnodes[i,k]] - qexact)
      nodeerror > maxerror ? maxerror = nodeerror : nothing
      nodeerror *= nodeerror
      L2error += nodeerror*solver.sbp.w[i]/solver.jac[i,k]
    end
  end
  return sqrt(L2error), maxerror
end

end
