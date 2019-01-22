# tests for shock capturing

import EulerEquationMod: AbstractShockSensor, AbstractShockCapturing

function test_shocksensorPP()

  @testset "Shock sensor PP" begin
    opts = read_input_file("input_vals_jac2d.jl")
    opts["order"] = 2
    delete!(opts, "calc_jac_explicit")
    opts["force_solution_complex"] = true
    mesh, sbp, eqn, opts = solvePDE(opts)

    Tsol = eltype(eqn.q); Tres = eltype(eqn.res)
    q = eqn.q[:, :, 1]
    jac = ones(Float64, mesh.numNodesPerElement)
    res = zeros(eltype(eqn.res), mesh.numDofPerNode, mesh.numNodesPerElement)

    sensor = eqn.params.sensor_pp
    capture = eqn.params.projection_shock_capturing
    # initial condition is constant, check the sensor reports no shock
    Se, ee = EulerEquationMod.getShockSensor(eqn.params, sbp, sensor, q, jac)

    @test abs(Se) < 1e-12
    @test ee == 0

    fill!(res, 0)
    EulerEquationMod.applyShockCapturing(eqn.params, sbp, sensor, capture, q, jac, res)
    @test maximum(abs.(res)) < 1e-13

    # test when a shock is present
    q[1, 3] += 5
    Se, ee = EulerEquationMod.getShockSensor(eqn.params, sbp, sensor, q, jac)

    @test abs(Se) > 1e-12
    @test ee > 0.99

    fill!(res, 0)
    w = copy(q)
    for i=1:mesh.numNodesPerElement
      w_i = sview(w, :, i)
      q_i = sview(q, :, i)
      EulerEquationMod.convertToIR(eqn.params, q_i, w_i)
    end
    EulerEquationMod.applyShockCapturing(eqn.params, sbp, sensor, capture, q, jac, res)

    @test sum(res .* w) < 0  # the term is negative definite

    # case 3: ee = 1
    test_shocksensor_diff(eqn.params, sbp, sensor, q, jac)
    test_shockcapturing_diff(eqn.params, sbp, sensor, capture, q, jac)

    # case 2: ee on sin wave
    q[1, 3] = 1.0105
    for i=1:mesh.numNodesPerElement
      for j=2:mesh.numDofPerNode
        q[j, i] += 0.1*(i + j)
      end
    end

    test_shocksensor_diff(eqn.params, sbp, sensor, q, jac)
    test_shockcapturing_diff(eqn.params, sbp, sensor, capture, q, jac)

  end  # end testset


  return nothing
end


"""
  Tests derivative of the shock sensor at a given state
"""
function test_shocksensor_diff(params, sbp, sensor::AbstractShockSensor, _q, jac)

  numDofPerNode, numNodesPerElement = size(_q)
  q = zeros(Complex128, numDofPerNode, numNodesPerElement)
  copy!(q, _q)
  Se_jac = zeros(q); Se_jac2 = zeros(q)
  ee_jac = zeros(q); ee_jac2 = zeros(q)

  h = 1e-20
  pert = Complex128(0, h)
  for i=1:numNodesPerElement
    for j=1:numDofPerNode
      q[j, i] += pert
      Se, ee = EulerEquationMod.getShockSensor(params, sbp, sensor, q, jac)
      Se_jac[j, i] = imag(Se)/h
      ee_jac[j, i] = imag(ee)/h
      q[j, i] -= pert
    end
  end

  EulerEquationMod.getShockSensor_diff(params, sbp, sensor, q, jac,
                                       Se_jac2, ee_jac2)

  @test maximum(abs.(Se_jac - Se_jac2)) < 1e-12
  @test maximum(abs.(ee_jac - ee_jac2)) < 1e-12

  # test vector mode
  q_dot = rand_realpart(size(q))
  q .+= pert*q_dot
  Se, ee = EulerEquationMod.getShockSensor(params, sbp, sensor, q, jac)
  Se_dot = imag(Se)/h
  ee_dot = imag(ee)/h
  q .-= pert*q_dot

  # run again to make sure intermediate arrays are zeroed out
  EulerEquationMod.getShockSensor_diff(params, sbp, sensor, q, jac,
                                       Se_jac2, ee_jac2)

  Se_dot2 = sum(Se_jac2 .* q_dot)
  ee_dot2 = sum(ee_jac2 .* q_dot)

  @test abs(Se_dot - Se_dot2) < 1e-12
  @test abs(ee_dot - ee_dot2) < 1e-12

  
  return nothing
end

"""
  Tests Jacobian of shock capturing scheme
"""
function test_shockcapturing_diff(params, sbp, sensor::AbstractShockSensor,
                                  capture::AbstractShockCapturing,
                                   _q, jac)

  numDofPerNode, numNodesPerElement = size(_q)
  q = zeros(Complex128, size(_q))
  copy!(q, _q)
  q_dot = rand_realpart(size(q))
  #q_dot = zeros(numDofPerNode, numNodesPerElement)
  #q_dot[:, 5] = 2
  res = zeros(Complex128, size(q))
  
  h = 1e-20
  pert = Complex128(0, h)

  # complex step
  q .+= pert*q_dot
  EulerEquationMod.applyShockCapturing(params, sbp, sensor, capture, q, jac, res)
  res_dot = imag(res)./h
  q .-= pert*q_dot


  # AD
  res_jac = zeros(numDofPerNode, numDofPerNode, numNodesPerElement, numNodesPerElement)
  EulerEquationMod.applyShockCapturing_diff(params, sbp, sensor, capture, q,
                                            jac, res_jac)

  res_dot2 = zeros(Complex128, numDofPerNode, numNodesPerElement)
  for i=1:numNodesPerElement
    for j=1:numDofPerNode
      res_dot2 .+= res_jac[:, j, :, i] * q_dot[j, i]
    end
  end

  @test maximum(abs.(res_dot - res_dot2)) < 1e-12

  return nothing
end

add_func1!(EulerTests, test_shocksensorPP, [TAG_SHORTTEST])


#------------------------------------------------------------------------------
# test LDG shock capturing


function test_ldg()
  opts = Dict{String, Any}(
    "physics" => "Euler",
    "operator_type" => "SBPOmega",
    "dimensions" => 2,
    "run_type" => 5,
    "jac_method" => 2,
    "jac_type" => 2,
    "order" => 2,
    "IC_name" => "ICIsentropicVortex",
    "use_DG" => true,
    "volume_integral_type" => 2,
    "Volume_flux_name" => "IRFlux",
    "face_integral_type" => 2,
    "FaceElementIntegral_name" => "ESLFFaceIntegral",
    "Flux_name" => "IRFlux",
    "numBC" => 3,
    "BC1" => [0],
    "BC1_name" => "isentropicVortexBC",  # outlet
    "BC2" => [2],
    "BC2_name" => "isentropicVortexBC", # inlet
    "BC3" => [1, 3],
    "BC3_name" => "noPenetrationBC",  # was noPenetrationBC
    "aoa" => 0.0,
    "smb_name" => "SRCMESHES/vortex_3x3_.smb",
    "dmg_name" => ".null",
    "itermax" => 20,
    "res_abstol" => 1e-9,
    "res_reltol" => 1e-9,
    "do_postproc" => true,
    "exact_soln_func" => "ICIsentropicVortex",
    "force_solution_complex" => true,
    "force_mesh_complex" => true,
    "solve" => false,
    )


  @testset "Local DG shock capturing" begin

    mesh, sbp, eqn, opts = solvePDE(opts)

    testQx(mesh, sbp, eqn, opts)
    test_shockmesh(mesh, sbp, eqn, opts)
    test_thetaface(mesh, sbp, eqn, opts)

  end

  return nothing
end


add_func1!(EulerTests, test_ldg, [TAG_SHORTTEST, TAG_TMP])


"""
  Set q to be a polynomial of specified degree
"""
function setPoly(mesh, q::Abstract3DArray, degree::Int)

  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      x = mesh.coords[1, j, i]
      y = mesh.coords[2, j, i]
      val = x^degree + 2*(y^degree)
      if mesh.dim == 3
        z = mesh.coords[3, j, i]
        val += 3*(z^degree)
      end

      for k=1:mesh.numDofPerNode
        q[k, j, i] = val + k
      end
    end
  end

  return nothing
end

function setPolyDeriv(mesh, qx, degree::Int)
# qx should be mesh.numDofPerNode x mesh.dim x mesh.numNodesPerElement x
# mesh.numEl

  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      x = mesh.coords[1, j, i]
      y = mesh.coords[2, j, i]
      #val = x^degree + 2*(y^degree)
      valx = degree*x^(degree-1)
      valy = 2*degree*y^(degree-1)
      if mesh.dim == 3
        z = mesh.coords[3, j, i]
        #val += 3*(z^degree)
        valz = 3*degree*z^(degree-1)
      end

      for k=1:mesh.numDofPerNode
        qx[k, 1, j, i] = valx
        qx[k, 2, j, i] = valy
        if mesh.dim == 3
          qx[k, 3, j, i] = valz
        end
      end
    end
  end

  return nothing
end



"""
  Returns array: `iface_idx`, `numFacesPerElement`
  x `numEl`, containing the indices of the faces that compose
  this element.
"""
function getInterfaceList(mesh)

  iface_idx = zeros(Int, mesh.dim + 1, mesh.numEl)  # interface index

  for i=1:mesh.numInterfaces
    iface_i = mesh.interfaces[i]
    for k=1:(mesh.dim + 1)
      if iface_idx[k, iface_i.elementL] == 0
        iface_idx[k, iface_i.elementL] = i
        break
      end
    end

    for k=1:(mesh.dim + 1)
      if iface_idx[k, iface_i.elementR] == 0
        iface_idx[k, iface_i.elementR] = i
        break
      end
    end
  end

  return iface_idx
end

"""
  Computes [Ex, Ey, Ez] * q

  **Inputs**

   * mesh
   * i: element number
   * iface_idx: the indices of the faces of element `i` `in mesh.interfaces`
   * q_i: mesh.numDofPerNode x mesh.numNodesPerElement
   * qface: mesh.numDofPerNode x mesh.numNodesPerFace work array
   * work2: mesh.numDofPerNode x mesh.numNodesPerFace x mesh.dim work array
   * E_term: mesh.numDofPerNode x mesh.numNodesPerElement x meshh.dim output
             array
"""
function applyE(mesh, i::Integer, iface_idx::AbstractVector, 
                q_i::AbstractMatrix, qface::AbstractMatrix,
                work2::Abstract3DArray, E_term::Abstract3DArray)

  fill!(E_term, 0)
  for f=1:size(iface_idx, 1)
    idx_f = iface_idx[f]
    iface_f = mesh.interfaces[idx_f]
    if iface_f.elementL == i
      face = iface_f.faceL
    else
      @assert iface_f.elementR == i
      face = iface_f.faceR
    end
   
    boundaryFaceInterpolate!(mesh.sbpface, face, q_i, qface)
    for j=1:mesh.numNodesPerFace
      for d=1:mesh.dim
        # figure out the right normal vector
        if  iface_f.elementL == i
          nj = mesh.nrm_face[d, j, idx_f]
        else
          nj = -mesh.nrm_face[d, mesh.sbpface.nbrperm[j, iface_f.orient], idx_f]
        end

        for k=1:mesh.numDofPerNode
          work2[k, j, d] = nj*qface[k, j]
        end

      end   # end d
    end  # end j

    for d=1:mesh.dim
      work2_d = sview(work2, :, :, d)
      E_d = sview(E_term, :, :, d)
      boundaryFaceIntegrate!(mesh.sbpface, face, work2_d, E_d)
    end
  end  # end f

  return nothing
end


function testQx(mesh, sbp, eqn::EulerData{Tsol, Tres}, opts) where {Tsol, Tres}

  degree = sbp.degree
#  degree = 1
  setPoly(mesh, eqn.q, degree)
  qderiv = zeros(eltype(eqn.q), mesh.numDofPerNode, mesh.dim, mesh.numNodesPerElement, mesh.numEl)
  setPolyDeriv(mesh, qderiv, degree)

  iface_idx = getInterfaceList(mesh)

  # test that Dx = -M * Qx^T + M*Ex, where Dx is exact for polynomials
  qx_term = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
  work = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
  op = SummationByParts.Subtract()

  qface = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
  work2 = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace, mesh.dim)
  E_term = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
  nel = 0  # number of elements tests
  for i=1:mesh.numEl
    if iface_idx[end, i] == 0  # non-interior element
      continue
    end
    nel += 1

    q_i = ro_sview(eqn.q, :, :, i)
    dxidx_i = ro_sview(mesh.dxidx, :, :, :, i)
    
    # do Qx
    fill!(qx_term, 0)
#    fill!(work, 0)
    EulerEquationMod.applyQxTransposed(sbp, q_i, dxidx_i, work, qx_term, op)

    # Do Ex: interpolate to face, apply normal vector (apply nbrperm and fac),
    #        reverse interpolate
    #TODO: this only works for fully interior elements

    applyE(mesh, i, sview(iface_idx, :, i), q_i, qface,  work2, E_term)


    # check against analytical derivative
    for j=1:mesh.numNodesPerElement
      for d=1:mesh.dim
        fac = mesh.jac[j, i]/sbp.w[j]
        for k=1:mesh.numDofPerNode
          val = fac*(qx_term[k, j, d] + E_term[k, j, d])
          @test abs(val - qderiv[k, d, j, i]) < 1e-12
        end
      end
    end


  end  # end i

  # make sure some elements were done
  @test nel > 0

  testQx2(mesh, sbp, eqn, opts)

  return nothing
end


function testQx2(mesh, sbp, eqn::EulerData{Tsol, Tres}, opts)  where {Tsol, Tres}
  # now that the first method of applyQxTransposed is verified, use it to
  # test the second

  degree = sbp.degree
  wx = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
  wxi = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
  res1_tmp = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
  res1 = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement)
  res2 = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement)

  for i=1:mesh.numEl

    # setup polynomial
    for j=1:mesh.numNodesPerElement
      x = mesh.coords[1, j, i]
      y = mesh.coords[2, j, i]
      for d=1:mesh.dim
        for k=1:mesh.numDofPerNode
          # make the polynomials and their derivatives different in each direction
          facx = d + k
          facy = 2*d + k  
          facz = 3*d + k

          wx[k, j, d] = facx*x^degree + facy*y^degree
          if mesh.dim == 3
            z = mesh.coords[3, j, k]
            wx[k, j, d] += facz*z^degree
          end
        end
      end  # end 
    end   # end j

    dxidx_i = ro_sview(mesh.dxidx, :, :, :, i)

    # call first method
    # This computes [Qx, Qy, Qz] * w_d, sum only Q_x * w_d into res1
    fill!(res1, 0)
    for d=1:mesh.dim
      fill!(res1_tmp, 0)
      EulerEquationMod.applyQxTransposed(sbp, sview(wx, :, :, d), dxidx_i,
                                         wxi, res1_tmp)
      for j=1:mesh.numNodesPerElement
        for k=1:mesh.numDofPerNode
          res1[k, j] += res1_tmp[k, j, d]
        end
      end
    end  # end d

    # second method
    fill!(res2, 0)
    EulerEquationMod.applyQxTransposed(sbp, wx, dxidx_i, wxi, res2)

    @test maximum(abs.(res2 - res1)) < 1e-13
  end  # end i

  return nothing
end


function test_shockmesh(mesh, sbp, eqn::EulerData{Tsol, Tres}, opts) where {Tsol, Tres}
   
  # construct the shock mesh with all elements in it that are fully interior
  iface_idx = getInterfaceList(mesh)
  shockmesh = EulerEquationMod.ShockedElements{Tres}(mesh)

  for i=1:mesh.numEl
    if iface_idx[end, i] != 0
      push!(shockmesh, i, 0.0)
    end
  end

  EulerEquationMod.completeShockElements(mesh, shockmesh)

  # all elements that are fully interior should be listed as shocked, all
  # other elements should be on the boundary

  elnums_shock = sview(shockmesh.elnums_all, 1:shockmesh.numShock)
  elnums_bndry = sview(shockmesh.elnums_all, (shockmesh.numShock+1):(shockmesh.numEl))
  for i=1:mesh.numEl
    if iface_idx[end, i] != 0
      @test i in elnums_shock
    else
      @test i in elnums_bndry
    end
  end

  # get all elements on boundary
  boundary_els = Array{Int}(mesh.numBoundaryFaces)
  for i=1:mesh.numBoundaryFaces
    boundary_els[i] = mesh.bndryfaces[i].element
  end
  boundary_els = unique(boundary_els)
  sort!(boundary_els)

  boundary_els_shock = sort!(shockmesh.elnums_all[(shockmesh.numShock+1):end])

  @test length(boundary_els) == length(boundary_els_shock)
  @test maximum(boundary_els - boundary_els_shock) == 0

  @test shockmesh.numEl == mesh.numEl
  
  # make sure no elements are double-counted
  @test length(unique(shockmesh.elnums_all[1:shockmesh.numEl])) == mesh.numEl

  # check interfaces
  for i=1:shockmesh.numInterfaces
    iface_red = shockmesh.ifaces[i]
    idx_orig = iface_red.idx_orig
    elnum_fullL = shockmesh.elnums_all[iface_red.iface.elementL]
    elnum_fullR = shockmesh.elnums_all[iface_red.iface.elementR]
    @test elnum_fullL == Int(mesh.interfaces[idx_orig].elementL)
    @test elnum_fullR == Int(mesh.interfaces[idx_orig].elementR)
  end
end


function test_thetaface(mesh, sbp, eqn::EulerData{Tsol, Tres}, opts) where {Tsol, Tres}

  degree = sbp.degree

  # construct the shock mesh with all elements in it that are fully interior
  iface_idx = getInterfaceList(mesh)
  shockmesh = EulerEquationMod.ShockedElements{Tres}(mesh)

  for i=1:mesh.numEl
    if iface_idx[end, i] != 0
      push!(shockmesh, i, 0.0)
    end
  end

  EulerEquationMod.completeShockElements(mesh, shockmesh)

  capture = EulerEquationMod.LDGShockCapturing{Tsol, Tres}()
  EulerEquationMod.allocateArrays(capture, mesh, shockmesh)
  flux = EulerEquationMod.LDG_ESFlux()

  # set w_el to be polynomial
  for i=1:shockmesh.numEl
    i_full = shockmesh.elnums_all[i]
    for j=1:mesh.numNodesPerElement
      x = mesh.coords[1, j, i_full]
      y = mesh.coords[2, j, i_full]
      for k=1:mesh.numDofPerNode
        capture.w_el[k, j, i] = k*x^degree + (k+5)*y^degree
        if mesh.dim == 3
          z = mesh.coords[3, j, i_full]
          capture.w_el[k, j, i] += (k+10)*z^degree
        end
      end
    end
  end

  # use the LDG code
  fill!(capture.q_j, 0)
  EulerEquationMod.computeThetaFaceContribution(mesh, sbp, eqn, opts,
                                    capture, shockmesh, flux)

  # compute against Ex.  For polynomials, the face contribution
  # reduces to the E_i w (sum i) operator
  qface = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
  work = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace, mesh.dim)
  E_term = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)

  for i=1:shockmesh.numShock
    i_full = shockmesh.elnums_all[i]
    w_i = sview(capture.w_el, :, :, i)
    idx_i = sview(iface_idx, :, i_full)
    applyE(mesh, i_full,  idx_i, w_i, qface, work, E_term)

    for d=1:mesh.dim
      @test maximum(abs.(capture.q_j[:, :, d, i] - E_term[:, :, d])) < 1e-13
    end
  end

  return nothing
end

    


