# this file contains the functions to evaluate the right hand side of the weak form in the pdf
# SL stands for semi-linear form, contents inside the inv(M)(...) parenthesis on pg. 7 of Jared's derivation


using SummationByParts
using PdePumiInterface


function dataPrep(mesh::AbstractMesh, sbp::SBPOperator, eqn::EulerEquation, SL::AbstractVector, SL0::AbstractVector)
# gather up all the data needed to do vectorized operatinos on the mesh
# linear elements only

# need: u (previous timestep solution), x (coordinates), dxidx, jac, res, array of interfaces

  u = Array(Float64, mesh.numDofPerNode, sbp.numnodes, mesh.numEl)  # hold previous timestep solution
  x = Array(Float64, 2, sbp.numnodes, mesh.numEl)  # hold x y coordinates of nodes
  dxidx = Array(Float64, 2, 2, sbp.numnodes, mesh.numEl)  
  jac = Array(Float64, sbp.numnodes, mesh.numEl)
  res = zeros(Float64, mesh.numDofPerNode, sbp.numnodes, mesh.numEl)  # hold result of computation

  # reused loop variables
  dofnums = zeros(mesh.numDofPerNode, sbp.numnodes)

  for i=1:mesh.numEl  # loop over elements
    dofnums = getGlobalNodeNumbers(mesh, i)
    u[:, :, i] = SL0[dofnums]

    x[:,:,i] = getElementVertCoords(mesh, [i])[1:2, :, :]

  end

  # get dxidx, jac using x
  mappingjacobian!(sbp, x, dxidx, jac)  

  # only need internal boundaries (not external)
  num_ext_edges = size(getBoundaryEdgeNums(mesh))[1]  # bad memory efficiency
  num_int_edges = getNumEdges(mesh) - num_ext_edges

  new_bndry = Boundary(2, 3)
#   println("new_bndry = ", new_bndry)

  new_interface = Interface(1, 2, 3, 4)
#   println("new_interface = ", new_interface)

#   println("num_int_edges = ", num_int_edges)

  interfaces = Array(typeof(new_interface), num_int_edges)

  pos = 1 # current position in interfaces
  for i=1:getNumEdges(mesh)
#     println("i = ", i)
#     println("pos = ", pos)
    # get number of elements using the edge
    adjacent_nums, num_adjacent = getAdjacentEntityNums(mesh, i, 1, 2)
#     println("num_adjacent = ", num_adjacent)
#     println("adjacent_nums = ", adjacent_nums)
    if num_adjacent > 1  # internal edge
#       println("this is an internal edge")
      element1 = adjacent_nums[1]
      element2 = adjacent_nums[2]

      coords_1 = x[:, :, element1]
      coords_2 = x[:, :, element2]

      # calculate centroid
      centroid1 = sum(coords_1, 2)
      centroid2 = sum(coords_2, 2)

      if abs(centroid1[1] - centroid2[2]) < 1e-10  # if big enough difference
        if centroid1[1] < centroid2[1]
    elementL = element1
    elementR = element2
        else
    elementL = element2
    elementR = element1
        end
      else  # use y coordinate to decide
        if centroid1[2] < centroid2[2]
    elementL = element1
    elementR = element2
        else
    elementL = element2
    elementR = element1
        end
      end

      edgeL = getEdgeLocalNum(mesh, i, elementL)
      edgeR = getEdgeLocalNum(mesh, i, elementR)

      interfaces[pos] = Interface(elementL, elementR, edgeL, edgeR)
#       println("updating pos")
      pos += 1

  #    print("\n")

    end  # end if internal edge

#     print("\n")
  end  # end loop over edges

  return u, x, dxidx, jac, res, interfaces

end # end function dataPrep

function evalVolumeIntegrals(mesh::AbstractMesh, sbp::SBPOperator, eqn::EulerEquation, SL::AbstractVector, SL0::AbstractVector)
# evaluate all the integrals over the elements (but not the boundary)
# does not do boundary integrals
# mesh : a mesh type, used to access information about the mesh
# sbp : an SBP operator, used to get stiffness matricies and stuff
# eqn : a type that holds some constants needed to evaluate the equation
#	: also used for multiple dispatch
# SL : solution vector to be populated (mesh.numDof entries)
# SL0 : solution vector at previous timesteps (mesh.numDof entries)

  numEl = getNumEl(mesh)
  
  p = 1           # linear elements
  
  # this is commented out because it gets zeroed out earlier, and we don't want to throw out other work
#   SL = zeros(Float64,mesh.numDof)

  println("==== in evalVolumeIntegrals ====")
  
  for element = 1:numEl
  
    dofnums = getGlobalNodeNumbers(mesh, element)

  #   sbp = TriSBP{Float64}(degree=p)

    # From PUMI, PdePumiInterface.jl, line ~142
    x_3D = getElementVertCoords(mesh,[element])

    # need x_3D -> x_2D
    # required by mappingjacobian, dimensions:
    #   1: coord
    #   2: node
    #   3: elem
    x = sub(x_3D, 1:2, 1:3, :)

    # required by mappingjacobian, dimensions:
    #   1: ref coord
    #   2: phys coord
    #   3: node
    #   4: elem
    dxidx = zeros(Float64, 2, 2, 3, 1)

    # required by mappingjacobian, dimensions:
    #   1: node
    #   2: elem
    jac = zeros(Float64, 3, 1)
    
    mappingjacobian!(sbp, x, dxidx, jac)
  
    # cheating a little, constant across all nodes since linear, so use node #1 for all
    # dxidx dimensions:
    #   1: ref coord
    #   2: phys coord
    #   3: node
    #   4: elem
    dxi_dx = dxidx[1,1,1,1]
    dxi_dy = dxidx[1,2,1,1]
    deta_dx = dxidx[2,1,1,1]
    deta_dy = dxidx[2,2,1,1]
  
    F1 = zeros(Float64, 12)
    F2 = zeros(Float64, 12)

    getF1(mesh, sbp, eqn, SL0, element, F1)
    getF2(mesh, sbp, eqn, SL0, element, F2)
  
    # du/dt = inv(M)*(volumeintegrate!(f1/jac) + volumeintegrate!(f2/jac)
  
  
    # ------------- setting source terms
    src_rho = 0.0
    src_rhou = 0.0
    src_rhov = 0.0
    src_E = 0.0

    # f's need 2 indices for 2D, and 3 entries corresponding to the vertices of a tri element
    # little f's: source terms; f1, for ex, source @ all 3 nodes for 1st conserved variable
    f1 = zeros(Float64,3,1)
    f2 = zeros(Float64,3,1)
    f3 = zeros(Float64,3,1)
    f4 = zeros(Float64,3,1)

    # ------------- populating source terms
    for node = 1:3
      f1[node,1] = src_rho
      f2[node,1] = src_rhou
      f3[node,1] = src_rhov
      f4[node,1] = src_E
    end

    result1 = zeros(size(f1))
    result2 = zeros(size(f2))
    result3 = zeros(size(f3))
    result4 = zeros(size(f4))

    # volumeintegrate! only does += to the result, so no need for result1, result2, etc
    # it is assumed that `u` is a rank-2 array, with the first dimension for the local-node index, and the second dimension for the element index.
    volumeintegrate!(sbp, f1/jac[1,1], result1)
    volumeintegrate!(sbp, f2/jac[1,1], result2)
    volumeintegrate!(sbp, f3/jac[1,1], result3)
    volumeintegrate!(sbp, f4/jac[1,1], result4)
  
    source_result = zeros(Float64,12)
    i = 1
    for node = 1:3
      source_result[i] = result1[node,1]
      source_result[i+1] = result2[node,1]
      source_result[i+2] = result3[node,1]
      source_result[i+3] = result4[node,1]
      i = i+4
    end
     
    flux_result = eqn.bigQT_xi*(F1*dxi_dx + F2*dxi_dy) + eqn.bigQT_eta*(F1*deta_dx + F2*deta_dy)

    #=
    println("F1hat: \n",F1)
    println("F2hat: \n",F2)
#     println("bigQT_xi: \n",eqn.bigQT_xi)
#     println("bigQT_eta: \n",eqn.bigQT_eta)
    println("dofnums: \n",dofnums)
    println("dxi_dx: \n",dxi_dx)
    println("dxi_dy: \n",dxi_dy)
    println("deta_dx: \n",deta_dx)
    println("deta_dy: \n",deta_dy)
    =#
  
    F1 = zeros(Float64, 12)
    F2 = zeros(Float64, 12)
  
    SL_el = source_result + flux_result
    println("source_result: ",source_result)
    println("flux_result: ",flux_result)
  
    # function assembleU(vec::AbstractVector, element::Integer, u::AbstractVector)
    # assembles a vector vec (of size 12, coresponding to solution values for an element), into the global solution vector u
    # element specifies which element number the number in vec belong to

    println("- element #: ",element,"   SL_el: ",SL_el)

    assembleSL(SL_el, element, SL)
  
  end
  
  println("==== end of evalVolumeIntegrals ====")
  return nothing 

end
#------------- end of evalVolumeIntegrals

# Euler flux calculator used by isentropicVortexBC ONLY!!
function getEulerFlux(u_vals::AbstractArray, nx, ny)

  eqn = EulerEquation(sbp)
  f1 = zeros(Float64, 4)
  f2 = zeros(f1)

  pressure = calcPressure(u_vals, eqn)
  f1[1] = u_vals[2]
  f1[2] = (u_vals[2]^2)/u_vals[1] + pressure
  f1[3] = (u_vals[2]*u_vals[3])/u_vals[1]
  f1[4] = (u_vals[4] + pressure)*u_vals[2]/u_vals[1]

  f2[1] = u_vals[3]
  f2[2] = (u_vals[2]*u_vals[3])/u_vals[1]
  f2[3] = (u_vals[3]^2)/u_vals[1] + pressure
  f2[4] = (u_vals[4] + pressure)*u_vals[3]/u_vals[1]

  eulerflux = f1*nx + f2*ny
  return eulerflux

end


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
  calcRho1Energy2(x, eqn, qg)
  
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
  rhoA = abs(Un) + dA*a
  lambda1 = d0_5*(tau*max(abs(lambda1),sat_Vn *rhoA) + sgn*lambda1)
  lambda2 = d0_5*(tau*max(abs(lambda2),sat_Vn *rhoA) + sgn*lambda2)
  lambda3 = d0_5*(tau*max(abs(lambda3),sat_Vl *rhoA) + sgn*lambda3)

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
function isentropicVortexBC(q, x, dxidx, nrm)

  E1dq = zeros(Float64, 4)
  E2dq = zeros(Float64, 4)

  # getting qg
  qg = zeros(Float64, 4)
  calcIsentropicVortex(x, eqn, qg)

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
  rhoA = abs(Un) + dA*a
  lambda1 = d0_5*(tau*max(abs(lambda1),sat_Vn *rhoA) + sgn*lambda1)
  lambda2 = d0_5*(tau*max(abs(lambda2),sat_Vn *rhoA) + sgn*lambda2)
  lambda3 = d0_5*(tau*max(abs(lambda3),sat_Vl *rhoA) + sgn*lambda3)

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
  
  return sat + getEulerFlux(q, nx, ny)

end # ends the function eulerRoeSAT


function evalBoundaryIntegrals(mesh::AbstractMesh, sbp::SBPOperator, eqn::EulerEquation, SL::AbstractVector, SL0::AbstractVector)
# evaluate all the integrals over the boundary
# does not do boundary integrals
# mesh : a mesh type, used to access information about the mesh
# operator : an SBP operator, used to get stiffness matricies and stuff
# eqn : a type that holds some constants needed to evaluate the equation
#	: also used for multiple dispatch
# u : solution vector to be populated (mesh.numDof entries), partially populated by evalVolumeIntegrals
# u0 : solution vector at previous timesteps (mesh.numDof entries)

  println("====== start of evalBoundaryIntegrals ======")

# Nodal Coordinates
x = zeros(Float64,(2,sbp.numnodes,getNumEl(mesh))); # nodal Coordinates of the marix
for i = 1:getNumEl(mesh)
  vtxcoord = getElementVertCoords(mesh, [i]);
  vtxcoord = squeeze(vtxcoord,3);
  vtxcoord = vtxcoord[1:2,:]
  vtxcoord = vtxcoord'
  # println(vtxcoord)
  x[:,:,i] = calcnodes(sbp, vtxcoord);
end

dxidx = zeros(Float64, (2,2,sbp.numnodes,getNumEl(mesh))); # Jacobian Matrix
dofJacobian = zeros(Float64, (2,2,4*sbp.numnodes,getNumEl(mesh))) # Jacobian Matrix needed for boundary integrate
jac = zeros(Float64, (sbp.numnodes,getNumEl(mesh))); # Determinant of the Jacobian Matrix

mappingjacobian!(sbp, x, dxidx, jac) # Get the Jocabian for transformation between actual and iso-parametric space

bndryfaces = Array(Boundary, mesh.numBoundaryEdges)
bndryfaces = getBoundaryArray(mesh)

numBoundaryElements = getNumBoundaryElements(mesh)
numEl = getNumEl(mesh)

# Calculate the intital condition
# u0 = ones(Float64, 4, sbp.numnodes, numEl)      # FOR TESTING

# res dimensions
#   1: num dofs
#   2: num nodes
#   3: num elem
# res = zeros(SL0) # stores the result of boundary integrate


u, x, dxidx, jac, result, interfaces = dataPrep(mesh, sbp, eqn, SL, SL0)

#=
println("type of SL0: ", typeof(SL0))
println("size of dxidx: ", size(dxidx,4))
println("size of SL0: ", size(SL0,3))
println("size of result: ", size(result,3))
println("size of x: ", size(x,3))
=#

boundaryintegrate!(sbp, bndryfaces, u, x, dxidx, isentropicVortexBC, result)
# boundaryintegrate!(sbp, bndryfaces, u, x, dxidx, rho1Energy2BC, result)

result = (-1)*result

# assembling into global SL vector
for element = 1:numEl
  for node = 1:sbp.numnodes
    vec = result[:,node,element]
    assembleSLNode(vec, element, node, SL)
  end
  println("- element #: ",element,"   result[:,:,element]:",result[:,:,element])
end

  println("==== end of evalBoundaryIntegrals ====")


  return nothing

end
#------------- end of evalBoundaryIntegrals

# This function adds edge stabilization to a residual using Prof. Hicken's edgestabilize! in SBP
function addEdgeStabilize(mesh::AbstractMesh, sbp::SBPOperator, eqn::EulerEquation, SL::AbstractVector, SL0::AbstractVector)

  println("==== start of addEdgeStabilize ====")
  # alpha calculated like in edgestabilize! documentation
  # stabscale (U+a)*gamma*h^2 where U=u*n, where u is the velocity 
  #   (remember to scale by rho) and n is the unit normal vector, from nrm->dxidx, then scaled by length
  # ifaces needs to be calculated
  # x needs to be passed
  # need to clarify u vs res. maybe change the u variable name to semilinear 

#   println("====== Entering edge stabilize ======")
  

    # dxidx dimensions:
    #   1: ref coord
    #   2: phys coord
    #   3: node
    #   4: elem

  u, x, dxidx, jac, result, interfaces = dataPrep(mesh, sbp, eqn, SL, SL0)
  numEl = getNumEl(mesh)

  alpha = zeros(Float64,2,2,sbp.numnodes,numEl)

  # calculating alpha, required by edgestabilize!
  for k = 1:numEl
    for i = 1:sbp.numnodes                                                                              
      for di1 = 1:2                                                                                     
        for di2 = 1:2                                                                                   
          alpha[di1,di2,i,k] = (dxidx[di1,1,i,k].*dxidx[di2,1,i,k] + dxidx[di1,2,i,k].*dxidx[di2,2,i,k])*jac[i,k]
        end                                                                                             
      end                                                                                               
    end                                                                                                 
  end     

  function stabscale(u, dxidx, nrm)

#     println("==== entering stabscale ====")

    # grabbing conserved variables
    rho = u[1]
    vel_x = u[2]/rho
    vel_y = u[3]/rho
    Energy = u[4]

    # from JC's code below, eqn should still be in scope
    pressure = calcPressure(u, eqn)

    # solved eqn for e: E = rho*e + (1/2)*rho*u^2
    vel_squared = vel_x^2 + vel_y^2
    energy = Energy/rho - (1/2)*vel_squared

    # gamma stored in EulerEquation type
    gamma = eqn.gamma

    println("pressure: ",pressure)
    println("gamma: ",gamma)
    println("rho: ",rho)
    # ideal gas law
    speed_sound = sqrt((gamma*pressure)/rho)

    # choice for edge stabilization constant: 
    #   refer to email from JH, 20150504:
    #   Anthony: there is little guidance in the literature regarding 
    #     gamma for the Euler equations.  I suggest starting with 
    #     gamma = 0.01.  If that fails (with a cfl of 1.0), then decrease 
    #     it by an order of magnitude at at time until RK is stable.  
    #     Once you find a value that works, try increasing it slowly.
#     edge_stab_gamma = 0.01  # default
    edge_stab_gamma = 0.0 
#     edge_stab_gamma = 0.00001

    # edge lengths component wise
    h_x = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
    h_y = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]

    # edge length
    h = sqrt(h_x^2 + h_y^2)

    # scaled velocity scalar
    U = vel_x*(nrm[1]/h) + vel_y*(nrm[2]/h)

    return (U + speed_sound)*edge_stab_gamma*h^2

  end

  # u argument here is SL in a different format
  edgestabilize!(sbp, interfaces, u, x, dxidx, jac, alpha, stabscale, result)

  # assembling into global SL vector
  for element = 1:numEl
    for node = 1:sbp.numnodes
      vec = result[:,node,element]
      assembleSLNode(vec, element, node, SL)
    end
    println("- element #: ",element,"   result[:,:,element]:",result[:,:,element])
  end
  println("==== end of addEdgeStabilize ====")


#   println("- element #: ",element,"   result[:,:,element]:",result[:,:,element])
#   println("==== end of evalBoundaryIntegrals ====")

  return nothing

end





function applyMassMatrixInverse(mesh::AbstractMesh, sbp::SBPOperator, eqn::EulerEquation, SL::AbstractVector, SL0::AbstractVector)
# apply the inverse of the mass matrix to the entire solution vector
# this is a good, memory efficient implimentation

  numEl = getNumEl(mesh)
  nnodes = sbp.numnodes
  dofpernode = getNumDofPerNode(mesh)
  for i=1:numEl
    dofnums_i = getGlobalNodeNumbers(mesh, i)  # get dof nums for this element
    for j=1:nnodes
      for k=1:dofpernode
        dofnum_k = dofnums_i[k,j]
        SL[dofnum_k] /= sbp.w[j]
      end
    end
  end
    
  return nothing 

end

# some helper functions

function getF1(mesh::AbstractMesh, sbp::SBPOperator, eqn::EulerEquation, SL0::AbstractVector, element::Integer, f1::AbstractVector)
# gets the vector F1 (see weak form derivation) for a particular element
# for linear triangles, the size of F1 is 3 nodes * 4 dof per node = 12 entries
# element : number of element to fetch F1 for
# f1 : vector (of length 12) to populate with F1.  This vector is overwritten

#   println("entered getF1")
#   println("element number = ", element)

  dofnums = getGlobalNodeNumbers(mesh, element)
#   println("dofnums = \n", dofnums)
  SL_vals = zeros(4)  # hold SL0 values for a single node

  for i=1:3  # loop over nodes
#     println("at node ", i)
    SL_vals = SL0[dofnums[:,i]]  # get the SL0 values
#     println("SL_vals = \n", SL_vals)

    # calculate pressure
  #  internal_energy = SL_vals[4]/SL_vals[1] - 0.5*(SL_vals[2]^2 + SL_vals[3]^2)/(SL_vals[1]^2)
  #  pressure = SL_vals[1]*eqn.R*internal_energy/eqn.cv
  #  println("internal_energy = ", internal_energy, " , pressure = ", pressure)
    pressure = calcPressure(SL_vals, eqn)

    # calculate F1 for this node
    start_index = 4*(i-1) + 1
    f1[start_index] = SL_vals[2]  # f1_1 (density equation)
    f1[start_index + 1] = (SL_vals[2]^2)/SL_vals[1] + pressure  # f1_2 (momentum-x equation)
    f1[start_index + 2] = (SL_vals[2]*SL_vals[3])/SL_vals[1]  # f1_3 (momentum-y equation)
    f1[start_index + 3] = (SL_vals[4] + pressure)*SL_vals[2]/SL_vals[1] # f1_4 (energy equation)

#     print("\n")
  end

#   println("F1 = \n", f1)

  return nothing

end

function getF2(mesh::AbstractMesh, sbp::SBPOperator, eqn::EulerEquation, SL0::AbstractVector, element::Integer, f2::AbstractVector)
# gets the vector F2 (see weak form derivation) for a particular element
# for linear triangles, the size of F2 is 3 nodes * 4 dof per node = 12 entries
# element : number of element to fetch F2 for
# f2 : vector (of length 12) to populate with F2.  This vector is overwritten

#   println("entered getF2")
#   println("element number = ", element)

  dofnums = getGlobalNodeNumbers(mesh, element)
#   println("dofnums = \n", dofnums)
  SL_vals = zeros(4)  # hold SL0 values for a single node

  for i=1:3  # loop over nodes
#     println("at node ", i)
    SL_vals = SL0[dofnums[:,i]]  # get the SL0 values
#     println("SL_vals = \n", SL_vals)

    # calculate pressure
  #  internal_energy = SL_vals[4]/SL_vals[1] - 0.5*(SL_vals[2]^2 + SL_vals[3]^2)/(SL_vals[1]^2)
  #  pressure = SL_vals[1]*eqn.R*internal_energy/eqn.cv
  #  println("internal_energy = ", internal_energy, " , pressure = ", pressure)
    pressure = calcPressure(SL_vals, eqn)

    # calculate F1 for this node
    start_index = 4*(i-1) + 1
    f2[start_index] = SL_vals[3]  # f2_1 (density equation)
    f2[start_index + 1] = (SL_vals[2]*SL_vals[3])/SL_vals[1] # f2_2 (momentum-x equation)
    f2[start_index + 2] = (SL_vals[3]^2)/SL_vals[1] + pressure  # f2_3 (momentum-y equation)
    f2[start_index + 3] = (SL_vals[4] + pressure)*SL_vals[3]/SL_vals[1] # f2_4 (energy equation)

#     print("\n")
  end

#   println("F2 = \n", f2)

  return nothing

end

function calcPressure(SL_vals::AbstractVector, eqn::EulerEquation)
  # calculate pressure for a node
  # SL is a vector of length 4

  internal_energy = SL_vals[4]/SL_vals[1] - 0.5*(SL_vals[2]^2 + SL_vals[3]^2)/(SL_vals[1]^2)
  pressure = SL_vals[1]*eqn.R*internal_energy/eqn.cv
  println("internal_energy = ", internal_energy, " , pressure = ", pressure)

  return pressure

end



function assembleSL(vec::AbstractVector, element::Integer, SL::AbstractVector)
# assembles a vector vec (of size 12, coresponding to solution values for an element), into the global solution vector SL
# element specifies which element number the number in vec belong to
#   println("entered assembleU")
#   println("element = ", element)
#   println("vec = \n", vec)
  dofnums = getGlobalNodeNumbers(mesh, element)

  for i=1:3  # loop over nodes
    dofnums_i = dofnums[:,i]
    start_index = 4*(i-1) + 1
    
    SL[dofnums_i] += vec[start_index:(start_index+3)]
  end

#   println("SL = \n", SL)

  return nothing

end

function assembleSL(vec::AbstractVector, element::Integer, component::Integer, SL::AbstractVector)
# assembles a vector vec (of size 3, corresponding to the solution for one degree of freedom at each node) into the global solution vector SL
# element specifies which element number the numbers in vec belong to
#  component specifies which dof of each node (1,2,3, or 4)

#   println("entered assembleU")
#   println("element = ", element)
#   println("component = ", component)
#   println("vec = \n", vec)
  dofnums = getGlobalNodeNumbers(mesh, element)
#   println("dofnums = ", dofnums)

  dofnums_comp = dofnums[component,:]
  SL[dofnums_comp.'] += vec

#   println("SL = \n", SL)

  return nothing

end

function assembleSLNode(vec::AbstractVector, element::Integer, node::Integer, SL::Vector)
# vec is vector of length 4 (conservative variables for a single node
# element is the element number that the node belongs to
# node is the node on the element (1, 2 or 3 for linear triangles)
# SL is global solution vector

  dofnums = getGlobalNodeNumbers(mesh, element)
  dofnums_n = dofnums[:, node]

  SL[dofnums_n] = vec
  
  return nothing

end  # end function


