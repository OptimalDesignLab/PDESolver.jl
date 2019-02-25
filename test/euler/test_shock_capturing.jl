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

    sensor = EulerEquationMod.ShockSensorPP{Tsol, Tres}(mesh, sbp, opts)
    capture = EulerEquationMod.ProjectionShockCapturing{Tsol, Tres}(mesh, sbp, eqn, opts)
    # initial condition is constant, check the sensor reports no shock
    Se, ee = EulerEquationMod.getShockSensor(eqn.params, sbp, sensor, q, jac)
    test_shockcapturing_diff(eqn.params, sbp, sensor, capture, q, jac)

    @test abs(Se) < 1e-12
    @test ee == 0

    fill!(res, 0)
    EulerEquationMod.calcShockCapturing(eqn.params, sbp, sensor, capture, q, jac, res)
    @test maximum(abs.(res)) < 1e-13

    # test when a shock is present
    q[1, 3] += 5
    Se, ee = EulerEquationMod.getShockSensor(eqn.params, sbp, sensor, q, jac)

    @test abs(Se) > 1e-12
    @test ee > 0.01

    fill!(res, 0)
    w = copy(q)
    for i=1:mesh.numNodesPerElement
      w_i = sview(w, :, i)
      q_i = sview(q, :, i)
      EulerEquationMod.convertToIR(eqn.params, q_i, w_i)
    end
    EulerEquationMod.calcShockCapturing(eqn.params, sbp, sensor, capture, q, jac, res)

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
  EulerEquationMod.calcShockCapturing(params, sbp, sensor, capture, q, jac, res)
  res_dot = imag(res)./h
  q .-= pert*q_dot


  # AD
  res_jac = zeros(numDofPerNode, numDofPerNode, numNodesPerElement, numNodesPerElement)
  EulerEquationMod.calcShockCapturing_diff(params, sbp, sensor, capture, q,
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

add_func1!(EulerTests, test_shocksensorPP, [TAG_SHORTTEST, TAG_TMP])


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
    "res_abstol" => 1e-12,
    "res_reltol" => 1e-12,
    "do_postproc" => true,
    "exact_soln_func" => "ICIsentropicVortex",
    "force_solution_complex" => true,
    "force_mesh_complex" => true,
    "solve" => false,
    "addVolumeIntegrals" => false,
    "addFaceIntegrals" => false,
    "addBoundaryIntegrals" => false,
    )


  @testset "Local DG shock capturing" begin

    mesh, sbp, eqn, opts = solvePDE(opts)

    testQx(mesh, sbp, eqn, opts)

    test_shockmesh(mesh, sbp, eqn, opts)
    
    test_thetaface(mesh, sbp, eqn, opts)
    test_qj(mesh, sbp, eqn, opts)
    test_qface(mesh, sbp, eqn, opts)
    test_q(mesh, sbp, eqn, opts)
    ic_func = EulerEquationMod.ICDict[opts["IC_name"]]
    ic_func(mesh, sbp, eqn, opts, eqn.q_vec)
    test_ldg_ESS(mesh, sbp, eqn, opts)

    test_br2_gradw(mesh, sbp, eqn, opts)
    test_br2_volume(mesh, sbp, eqn, opts)
    test_br2_face(mesh, sbp, eqn, opts)
    test_br2_Dgk(mesh, sbp, eqn, opts)

    for i=1:10
      ic_func = EulerEquationMod.ICDict[opts["IC_name"]]
      ic_func(mesh, sbp, eqn, opts, eqn.q_vec)
      test_br2_ESS(mesh, sbp, eqn, opts; fullmesh=true)
      
      ic_func = EulerEquationMod.ICDict[opts["IC_name"]]
      ic_func(mesh, sbp, eqn, opts, eqn.q_vec)
      test_br2_ESS(mesh, sbp, eqn, opts; fullmesh=false)
    end

    ic_func = EulerEquationMod.ICDict[opts["IC_name"]]
    ic_func(mesh, sbp, eqn, opts, eqn.q_vec)
    test_br2_serialpart(mesh, sbp, eqn, opts)

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
  Sets capture.w_el to be polynomial.  w_el is numDofPerNode x
  numNodesPerElement x shockmesh.numEl
"""
function setWPoly(mesh, shockmesh, w_el, degree::Int)

  # set w_el to be polynomial
  for i=1:shockmesh.numEl
    i_full = shockmesh.elnums_all[i]
    for j=1:mesh.numNodesPerElement
      x = mesh.coords[1, j, i_full]
      y = mesh.coords[2, j, i_full]
      for k=1:mesh.numDofPerNode
        w_el[k, j, i] = k*x^degree + (k+1)*y^degree
        if mesh.dim == 3
          z = mesh.coords[3, j, i_full]
          w_el[k, j, i] += (k+2)*z^degree
        end
      end
    end
  end

  return nothing
end

"""
  Get the xyz derivatives of setWPoly.  `w_elx` should be numDofPerNode x
  numNodesPerElement x dim x shockmesh.numEl
"""
function setWPolyDeriv(mesh, shockmesh, w_elx::AbstractArray{T, 4}, degree::Int) where {T}

  # set w_el to be polynomial
  for i=1:shockmesh.numEl
    i_full = shockmesh.elnums_all[i]
    for j=1:mesh.numNodesPerElement
      x = mesh.coords[1, j, i_full]
      y = mesh.coords[2, j, i_full]
      for k=1:mesh.numDofPerNode
        #capture.w_el[k, j, i] = k*x^degree + (k+1)*y^degree
        w_elx[k, j, 1, i] = degree*k*x^(degree-1)
        w_elx[k, j, 2, i] = degree*(k+1)*y^(degree-1)
        if mesh.dim == 3
          z = mesh.coords[3, j, i_full]
          #capture.w_el[k, j, i] += (k+2)*z^degree
          w_elx[k, j, 3, i] = degree*(k+2)*z^(degree-1)
        end
      end
    end
  end

  return nothing
end

"""
  Computes second derivative.  `w_elx` is numDofPerNode x numNodesPerElement
  x dim x dim x shockmesh.numEl, each dim x dim block contains d/dx_i dx_j.
"""
function setWPolyDeriv2(mesh, shockmesh, w_elx::AbstractArray{T, 5}, degree::Int) where {T}

  # set w_el to be polynomial
  fill(w_elx, 0)
  for i=1:shockmesh.numEl
    i_full = shockmesh.elnums_all[i]
    for j=1:mesh.numNodesPerElement
      x = mesh.coords[1, j, i_full]
      y = mesh.coords[2, j, i_full]
      for k=1:mesh.numDofPerNode
        #capture.w_el[k, j, i] = k*x^degree + (k+1)*y^degree
        w_elx[k, j, 1, 1, i] = degree*(degree-1)*k*x^(degree-2)
        w_elx[k, j, 2, 2, i] = degree*(degree-1)*(k+1)*y^(degree-2)
        if mesh.dim == 3
          z = mesh.coords[3, j, i_full]
          #capture.w_el[k, j, i] += (k+2)*z^degree
          w_elx[k, j, 3, 3, i] = degree*(degree-1)*(k+2)*z^(degree-2)
        end
      end
    end
  end

  return nothing
end



"""
  q_j is numDofPerNode x numNodesPerElement x dim x shockmesh.numEl
"""
function setQjPoly(mesh, shockmesh, q_j, degree::Int)

  # set w_el to be polynomial
  for i=1:shockmesh.numEl
    i_full = shockmesh.elnums_all[i]
    for j=1:mesh.numNodesPerElement
      x = mesh.coords[1, j, i_full]
      y = mesh.coords[2, j, i_full]
      for d=1:mesh.dim
        for k=1:mesh.numDofPerNode
          q_j[k, j, d, i] = (k+d)*x^degree + (k+1 + 2*d)*y^degree
          if mesh.dim == 3
            z = mesh.coords[3, j, i_full]
            q_j[k, j, d, i] += (k+2 + 3*d)*z^degree
          end
        end
      end
    end
  end

  return nothing
end


"""
  q_jx is numDofPerNode x numNodesPerElement x dim x dim x shockmesh.numEl
  the dim x dim block contains dq_i/dx_j
"""
function setQjPolyDeriv(mesh, shockmesh, q_jx, degree::Int)

  # set w_el to be polynomial
  for i=1:shockmesh.numEl
    i_full = shockmesh.elnums_all[i]
    for j=1:mesh.numNodesPerElement
      x = mesh.coords[1, j, i_full]
      y = mesh.coords[2, j, i_full]
      for d=1:mesh.dim
        for k=1:mesh.numDofPerNode
          #q_j[k, j, d, i] = (k+d)*x^degree + (k+1 + 2*d)*y^degree
          q_jx[k, j, d, 1, i] = degree*(k+d)*x^(degree-1)
          q_jx[k, j, d, 2, i] = degree*(k+1 + 2*d)*y^(degree-1)


          if mesh.dim == 3
            z = mesh.coords[3, j, i_full]
            #q_j[k, j, d, i] += (k+2 + 3*d)*z^degree
            q_jx[k, j, d, 3, i] += degree*(k+2 + 3*d)*z^degree
          end
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
  Gets a shockmesh consisting of only elements on the interior of the domain
  (so all faces of the elements are contained in mesh.interfaces)

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts

  **Outputs**

   * iface_idx: output of getInterfaceList
   * shockmesh: the shock mesh object
"""
function getInteriorMesh(mesh, sbp, eqn::EulerData{Tsol, Tres}, opts) where {Tsol, Tres}

  degree = sbp.degree

  # construct the shock mesh with all elements in it that are fully interior
  iface_idx = getInterfaceList(mesh)
  shockmesh = EulerEquationMod.ShockedElements{Tres}(mesh)

  for i=1:mesh.numEl
    if iface_idx[end, i] != 0
      push!(shockmesh, i, 1.0)
    end
  end

  EulerEquationMod.completeShockElements(mesh, shockmesh)

  return iface_idx, shockmesh
end

"""
  Like `getInteriorMesh`, but puts all elements in the shockmesh
"""
function getEntireMesh(mesh, sbp, eqn::EulerData{Tsol, Tres}, opts) where {Tsol, Tres}

  degree = sbp.degree

  # construct the shock mesh with all elements in it that are fully interior
  iface_idx = getInterfaceList(mesh)
  shockmesh = EulerEquationMod.ShockedElements{Tres}(mesh)

  for i=1:mesh.numEl
    push!(shockmesh, i, 1.0)
  end

  EulerEquationMod.completeShockElements(mesh, shockmesh)

  return iface_idx, shockmesh
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
  qxT_term = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
  work = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
  op = SummationByParts.Subtract()

  # test Dx
  dx_term = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)

  # test Qx
  qx_term = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)

  # test Dx^T
  dxT_term = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)

  # test calculating the operator matrices themselves
  Dx = zeros(Tres, mesh.numNodesPerElement, mesh.numNodesPerElement, mesh.dim)
  Qx = zeros(Tres, mesh.numNodesPerElement, mesh.numNodesPerElement, mesh.dim)
  DxT = zeros(Tres, mesh.numNodesPerElement, mesh.numNodesPerElement, mesh.dim)
  QxT = zeros(Tres, mesh.numNodesPerElement, mesh.numNodesPerElement, mesh.dim)


  dx_term2 = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
  qx_term2 = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)


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
    jac_i = ro_sview(mesh.jac, :, i)
    
    # do Qx^T
    fill!(qxT_term, 0)
#    fill!(work, 0)
    EulerEquationMod.applyQxTransposed(sbp, q_i, dxidx_i, work, qxT_term, op)

    # Do Ex: interpolate to face, apply normal vector (apply nbrperm and fac),
    #        reverse interpolate
    applyE(mesh, i, sview(iface_idx, :, i), q_i, qface,  work2, E_term)

    # test apply Dx
    fill!(dx_term, 0)
    EulerEquationMod.applyDx(sbp, q_i, dxidx_i, jac_i, work, dx_term)

    # test apply Qx
    fill!(qx_term, 0)
    EulerEquationMod.applyQx(sbp, q_i, dxidx_i, work, qx_term)

    # test apply Dx^T
    fill!(dxT_term, 0)
    EulerEquationMod.applyDxTransposed(sbp, q_i, dxidx_i, jac_i, work, dxT_term)

    # test explicitly computed operator matrices
    EulerEquationMod.calcDx(sbp, dxidx_i, jac_i, Dx)
    EulerEquationMod.calcQx(sbp, dxidx_i, Qx)
    EulerEquationMod.calcDxTransposed(sbp, dxidx_i, jac_i, DxT)
    EulerEquationMod.calcQxTransposed(sbp, dxidx_i, QxT)


    for d1=1:mesh.dim
      for k=1:mesh.numDofPerNode
        dx_term2[k, :, d1] = Dx[:, :, d1]*q_i[k, :]
        qx_term2[k, :, d1] = Qx[:, :, d1]*q_i[k, :]
      end
      @test maximum(abs.(DxT[:, :, d1] - Dx[:, :, d1].')) < 1e-13
      @test maximum(abs.(QxT[:, :, d1] - Qx[:, :, d1].')) < 1e-13
    end

    # check against analytical derivative
    for j=1:mesh.numNodesPerElement
      for d=1:mesh.dim
        fac = mesh.jac[j, i]/sbp.w[j]
        for k=1:mesh.numDofPerNode
          val = fac*(qxT_term[k, j, d] + E_term[k, j, d])
          val2 = dx_term[k, j, d]
          val3 = fac*qx_term[k, j, d]
          val4 = fac*qx_term2[k, j, d]
          @test abs(val - qderiv[k, d, j, i]) < 1e-12
          @test abs(val2 - qderiv[k, d, j, i]) < 1e-12
          @test abs(val3 - qderiv[k, d, j, i]) < 1e-12
          @test abs(dx_term2[k, j, d] - qderiv[k, d, j, i]) < 1e-12
          @test abs(val4 - qderiv[k, d, j, i]) < 1e-12
          
          @test abs(dx_term2[k, j, d] - dx_term[k, j, d]) < 1e-12
          @test abs(qx_term2[k, j, d] - qx_term[k, j, d]) < 1e-12

        end
      end
    end

    for d=1:mesh.dim
      for k=1:mesh.numDofPerNode
        # Dx^T is a bit weird because it isn't easily related to a polynomial.
        # Instead do: v^T Dx^T u = (v^T Dx^T) u = dv/dx^T u
        #                        = v^T (Dx^T u), where v is a polynomial
        val4 = dot(dxT_term[k, :, d], q_i[k, :])
        val5 = dot(qderiv[k, d, :, i], q_i[k, :])
        @test abs(val4 - val5) < 5e-12
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

  res1_qxT_tmp = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
  res1_qx_tmp = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
  res1_dx_tmp = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
  res1_dxT_tmp = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)

  res1_qxT = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement)
  res1_dx = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement)
  res1_qx = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement)
  res1_dxT = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement)
  
  res2_qxT = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement)
  res2_dx = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement)
  res2_qx = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement)
  res2_dxT = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement)

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
    jac_i = ro_sview(mesh.jac, :, i)

    # call first method
    # This computes [Qx, Qy, Qz] * w_d, sum only Q_x * w_d into res1
    fill!(res1_qxT, 0); fill!(res1_qx, 0); fill!(res1_dx, 0); fill!(res1_dxT, 0)
    for d=1:mesh.dim
      fill!(res1_qxT_tmp, 0); fill!(res1_qx_tmp, 0), fill!(res1_dx_tmp, 0)
      fill!(res1_dxT_tmp, 0)
      wx_d = sview(wx, :, :, d)
      EulerEquationMod.applyQxTransposed(sbp, wx_d, dxidx_i, wxi, res1_qxT_tmp)
      EulerEquationMod.applyQx(sbp, wx_d, dxidx_i, wxi, res1_qx_tmp)
      EulerEquationMod.applyDx(sbp, wx_d, dxidx_i, jac_i, wxi, res1_dx_tmp)
      EulerEquationMod.applyDxTransposed(sbp, wx_d, dxidx_i, jac_i, wxi, res1_dxT_tmp)
      for j=1:mesh.numNodesPerElement
        for k=1:mesh.numDofPerNode
          res1_qxT[k, j] += res1_qxT_tmp[k, j, d]
          res1_qx[k, j] += res1_qx_tmp[k, j, d]
          res1_dx[k, j] += res1_dx_tmp[k, j, d]
          res1_dxT[k, j] += res1_dxT_tmp[k, j, d]
        end
      end
    end  # end d

    # second method
    fill!(res2_qxT, 0); fill!(res2_qx, 0); fill!(res2_dx, 0); fill!(res2_dxT, 0)
    EulerEquationMod.applyQxTransposed(sbp, wx, dxidx_i, wxi, res2_qxT)
    EulerEquationMod.applyQx(sbp, wx, dxidx_i, wxi, res2_qx)
    EulerEquationMod.applyDx(sbp, wx, dxidx_i, jac_i, wxi, res2_dx)
    EulerEquationMod.applyDxTransposed(sbp, wx, dxidx_i, jac_i, wxi, res2_dxT)

    @test maximum(abs.(res2_qxT - res1_qxT)) < 1e-13
    @test maximum(abs.(res2_qx - res1_qx)) < 1e-13
    @test maximum(abs.(res2_dx - res1_dx)) < 1e-13
    @test maximum(abs.(res2_dxT - res1_dxT)) < 1e-12
  end  # end i

  return nothing
end


function test_shockmesh(mesh, sbp, eqn::EulerData{Tsol, Tres}, opts) where {Tsol, Tres}
   
  # construct the shock mesh with all elements in it that are fully interior
  iface_idx, shockmesh = getInteriorMesh(mesh, sbp, eqn, opts)

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

  boundary_els_shock = sort!(shockmesh.elnums_all[(shockmesh.numShock+1):shockmesh.numEl])

  @test length(boundary_els) == length(boundary_els_shock)
  @test maximum(boundary_els - boundary_els_shock) == 0

  @test shockmesh.numEl == mesh.numEl
  
  # make sure no elements are double-counted
  @test length(unique(shockmesh.elnums_all[1:shockmesh.numEl])) == mesh.numEl

  @test shockmesh.numBoundaryFaces == 0  # neighbor elements should not add
                                         # boundaries
  # check interfaces
  for i=1:shockmesh.numInterfaces
    iface_red = shockmesh.ifaces[i]
    idx_orig = iface_red.idx_orig
    elnum_fullL = shockmesh.elnums_all[iface_red.iface.elementL]
    elnum_fullR = shockmesh.elnums_all[iface_red.iface.elementR]
    @test elnum_fullL == Int(mesh.interfaces[idx_orig].elementL)
    @test elnum_fullR == Int(mesh.interfaces[idx_orig].elementR)
  end


  # add only boundary elements, to check that Boundaries are correctly added
  shockmesh2 = EulerEquationMod.ShockedElements{Tres}(mesh)

  for i in boundary_els
    push!(shockmesh2, i, 1.0)
  end
  EulerEquationMod.completeShockElements(mesh, shockmesh2)

  @test shockmesh2.numShock == length(boundary_els)
  println("length(shockmesh.bndryfaces) = ", length(shockmesh2.bndryfaces))
  println("shockmesh.bndryfaces = ", shockmesh.bndryfaces)
  @test shockmesh2.numBoundaryFaces == length(mesh.bndryfaces)
  
  # check original indices
  idx_orig = zeros(Int, 0)
  for i=1:shockmesh2.numBoundaryFaces
    idx_orig_i = shockmesh2.bndryfaces[i].idx_orig
    push!(idx_orig, idx_orig_i)
  end
  sort!(idx_orig)
  @test idx_orig == collect(1:mesh.numBoundaryFaces)


end


function test_thetaface(mesh, sbp, eqn::EulerData{Tsol, Tres}, opts) where {Tsol, Tres}

  degree = sbp.degree
  iface_idx, shockmesh = getInteriorMesh(mesh, sbp, eqn, opts)

  capture = EulerEquationMod.LDGShockCapturing{Tsol, Tres}()
  EulerEquationMod.allocateArrays(capture, mesh, shockmesh)
  flux = EulerEquationMod.LDG_ESFlux()

  # set w_el to be polynomial
  setWPoly(mesh, shockmesh, capture.w_el, degree)

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


function test_qj(mesh, sbp, eqn::EulerData{Tsol, Tres}, opts) where {Tsol, Tres}

  iface_idx, shockmesh = getInteriorMesh(mesh, sbp, eqn, opts)
  degree = sbp.degree

  capture = EulerEquationMod.LDGShockCapturing{Tsol, Tres}()
  EulerEquationMod.allocateArrays(capture, mesh, shockmesh)
  flux = EulerEquationMod.LDG_ESFlux()
  diffusion = EulerEquationMod.ShockDiffusion(shockmesh.ee)

  # q_j = D_j * w when interpolation is exact
  setWPoly(mesh, shockmesh, capture.w_el, degree)
  wx_el = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim, shockmesh.numEl)
  setWPolyDeriv(mesh, shockmesh, wx_el, degree)

  # use LDG code
  EulerEquationMod.computeThetaVolumeContribution(mesh, sbp, eqn, opts, capture, shockmesh)
  EulerEquationMod.computeThetaFaceContribution(mesh, sbp, eqn, opts, capture, shockmesh, flux)
  EulerEquationMod.computeQFromTheta(mesh, sbp, eqn, opts, capture, shockmesh, diffusion)

  # compare against analytical value
  for i=1:shockmesh.numShock
    for d=1:mesh.dim
      @test maximum(abs.(capture.q_j[:, :, d, i] - wx_el[:, :, d, i])) < 1e-12
    end
  end

  # all the neighbor elements should have q_j = 0 because epsilon = 0 there
  for i=(shockmesh.numShock+1):shockmesh.numEl
    for d=1:mesh.dim
      @test maximum(abs.(capture.q_j[:, :, d, i])) == 0
    end
  end

  return nothing
end


function test_qface(mesh, sbp, eqn::EulerData{Tsol, Tres}, opts) where {Tsol, Tres}

  degree = sbp.degree
  iface_idx, shockmesh = getInteriorMesh(mesh, sbp, eqn, opts)

  capture = EulerEquationMod.LDGShockCapturing{Tsol, Tres}()
  EulerEquationMod.allocateArrays(capture, mesh, shockmesh)
  flux = EulerEquationMod.LDG_ESFlux()


  setQjPoly(mesh, shockmesh, capture.q_j, degree)
  fill!(capture.w_el, 0)  # this only shows up in a jump term

  # use LDG code
  fill(eqn.res, 0)
  EulerEquationMod.computeQFaceTerm(mesh, sbp, eqn, opts, capture, shockmesh, flux)

  # for exact interpolation, the LDG term reduces to E_x * q_x + Ey*q_y
  qface = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
  work = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace, mesh.dim)
  E_term = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
  E_term2 = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement)

  for i=1:shockmesh.numShock
    i_full = shockmesh.elnums_all[i]
    idx_i = sview(iface_idx, :, i_full)

    # compute using E operator
    fill!(E_term2, 0)
    for d=1:mesh.dim
      q_j = sview(capture.q_j, :, :, d, i)
      # its a bit wasteful to compute E_y*q_x and E_x*q_y, but its only a test
      applyE(mesh, i_full,  idx_i, q_j, qface, work, E_term)
      for j=1:mesh.numNodesPerElement
        for k=1:mesh.numDofPerNode
          E_term2[k, j] += E_term[k, j, d]
        end
      end
    end

    # compare
    @test maximum(abs.(eqn.res[:, :, i_full] - E_term2)) < 1e-12
  end
end


function test_q(mesh, sbp, eqn::EulerData{Tsol, Tres}, opts) where {Tsol, Tres}

  # for exact interpolation, the q terms are -Q^T + E, which is equal to Q,
  # Thus M^-1 * Q = D, which we can compute analytically.
  # Because of the sum, we get, D_x * q_x + D_y * q_y

  degree = sbp.degree
  iface_idx, shockmesh = getInteriorMesh(mesh, sbp, eqn, opts)

  capture = EulerEquationMod.LDGShockCapturing{Tsol, Tres}()
  EulerEquationMod.allocateArrays(capture, mesh, shockmesh)
  flux = EulerEquationMod.LDG_ESFlux()


  q_jx = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim, mesh.dim, shockmesh.numEl)
  setQjPoly(mesh, shockmesh, capture.q_j, degree)
  setQjPolyDeriv(mesh, shockmesh, q_jx, degree)
  fill!(capture.w_el, 0)  # this only shows up in a jump term

  # use LDG code
  fill!(eqn.res, 0)
  EulerEquationMod.computeQFaceTerm(mesh, sbp, eqn, opts, capture, shockmesh, flux)
  EulerEquationMod.computeQVolumeTerm(mesh, sbp, eqn, opts, capture, shockmesh)

  # test against analytical derivative
  res2 = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement)
  for i=1:shockmesh.numShock
    i_full = shockmesh.elnums_all[i]

    fill!(res2, 0)
    for d=1:mesh.dim
      for j=1:mesh.numNodesPerElement
        for k=1:mesh.numDofPerNode
          res2[k, j] += q_jx[k, j, d, d, i]
        end
      end
    end

    # apply inverse mass matrix to LDG terms
    for j=1:mesh.numNodesPerElement
      fac = mesh.jac[j, i_full]/sbp.w[j]
      for k=1:mesh.numDofPerNode
        eqn.res[k, j, i_full] *= fac
      end
    end

    @test maximum(abs.(eqn.res[:, :, i_full] - res2)) < 1e-11
  end

  return nothing
end


function test_ldg_ESS(mesh, sbp, eqn::EulerData{Tsol, Tres}, opts) where {Tsol, Tres}

  # construct the shock mesh with all elements in it that are fully interior
  iface_idx, shockmesh = getInteriorMesh(mesh, sbp, eqn, opts)

  capture = EulerEquationMod.LDGShockCapturing{Tsol, Tres}()
  EulerEquationMod.allocateArrays(capture, mesh, shockmesh)
  flux = EulerEquationMod.LDG_ESFlux()

  # add random component to q
  # test entropy stability
  q_pert = 0.1*rand(size(eqn.q_vec))
  eqn.q_vec .+= q_pert
  array1DTo3D(mesh, sbp, eqn, opts, eqn.q_vec, eqn.q)
  fill!(eqn.res, 0)

  EulerEquationMod.calcShockCapturing(mesh, sbp, eqn, opts, capture, shockmesh)

  array3DTo1D(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)

  w_vec = zeros(Tsol, mesh.numDof)
  copy!(w_vec, eqn.q_vec)
  EulerEquationMod.convertToIR(mesh, sbp, eqn, opts, w_vec)

  val = dot(w_vec, eqn.res_vec)
  println("val = ", val)
  @test val < 0

  return nothing
end


#------------------------------------------------------------------------------
# test BR2

function test_br2_gradw(mesh, sbp, eqn::EulerData{Tsol, Tres}, opts) where {Tsol, Tres}

  iface_idx, shockmesh = getInteriorMesh(mesh, sbp, eqn, opts)

  # set q to be convertToConservative(polynomial)
  # Then the inverse process will produce known values, and have known
  # derivatives

  degree = sbp.degree
  setPoly(mesh, eqn.q, degree)
  q_deriv = zeros(Tres, mesh.numDofPerNode, mesh.dim, mesh.numNodesPerElement, mesh.numEl)
  setPolyDeriv(mesh, q_deriv, degree)

  w_vals = copy(eqn.q)

  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      q_j = sview(eqn.q, :, j, i)
      w_j = sview(w_vals, :, j, i)
      EulerEquationMod.convertToIR_(eqn.params, q_j, w_j)
    end
  end

  capture = EulerEquationMod.SBPParabolicSC{Tsol, Tres}(mesh, sbp, eqn, opts)
  EulerEquationMod.allocateArrays(capture, mesh, shockmesh)

  EulerEquationMod.computeGradW(mesh, sbp, eqn, opts, capture, shockmesh,
                                capture.entropy_vars, capture.diffusion)

  A0inv = zeros(Tsol, mesh.numDofPerNode, mesh.numDofPerNode)
  dw_dx = zeros(Tsol, mesh.numDofPerNode)
  for i=1:shockmesh.numEl
    i_full = shockmesh.elnums_all[i]
    @test maximum(abs.(capture.w_el[:, :, i] - w_vals[:, :, i_full])) < 1e-13
  end

  return nothing
end


function test_br2_volume(mesh, sbp, eqn::EulerData{Tsol, Tres}, opts) where {Tsol, Tres}

  degree = sbp.degree
  iface_idx, shockmesh = getInteriorMesh(mesh, sbp, eqn, opts)

  capture = EulerEquationMod.SBPParabolicSC{Tsol, Tres}(mesh, sbp, eqn, opts)
  EulerEquationMod.allocateArrays(capture, mesh, shockmesh)


  # set w_el to be polynomial
  setWPoly(mesh, shockmesh, capture.w_el, degree)
  setWPolyDeriv(mesh, shockmesh, capture.grad_w, degree)
  w_deriv2 = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim,
                   mesh.dim, shockmesh.numEl)
  setWPolyDeriv2(mesh, shockmesh, w_deriv2, degree)

  fill!(eqn.res, 0)
  EulerEquationMod.computeVolumeTerm(mesh, sbp, eqn, opts, capture, shockmesh)

  # test: w^T Dx^T H Dx w = wx^T H wx
  vals = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement)
  for i=1:shockmesh.numShock
    i_full = shockmesh.elnums_all[i]
    fill!(vals, 0)
    for j=1:mesh.numNodesPerElement
      for k=1:mesh.numDofPerNode
        eqn.res[k, j, i_full] *= mesh.jac[j, i_full]/sbp.w[j]
        for d1=1:mesh.dim
          vals[k, j] += w_deriv2[k, j, d1, d1, i]
        end
      end
    end

    @test maximum(abs.(vals - eqn.res[:, :, i_full])) < 1e-12
  end  # end i

  return nothing
end


function test_br2_face(mesh, sbp, eqn::EulerData{Tsol, Tres}, opts) where {Tsol, Tres}

  degree = sbp.degree
  iface_idx, shockmesh = getInteriorMesh(mesh, sbp, eqn, opts)

  capture = EulerEquationMod.SBPParabolicSC{Tsol, Tres}(mesh, sbp, eqn, opts)
  EulerEquationMod.allocateArrays(capture, mesh, shockmesh)

  # test alpha sums to 1
  alpha_sum = zeros(shockmesh.numEl)
  for i=1:shockmesh.numInterfaces
    iface_i = shockmesh.ifaces[i].iface
    alpha_sum[iface_i.elementL] += capture.alpha[1, i]
    alpha_sum[iface_i.elementR] += capture.alpha[2, i]
  end

  for i=1:shockmesh.numShock
    @test abs(alpha_sum[i] - 1) < 1e-15
  end

  if !shockmesh.isNeumann
    println("\nTesting alpha values")
    # in Dirichlet case, all alphas should be 1/numFacesPerElement (at least
    # for shocked elements)
    val = 1/(mesh.dim + 1)
    for i=1:shockmesh.numInterfaces
      iface_i = shockmesh.ifaces[i].iface
      if iface_i.elementL <= shockmesh.numShock
        @test abs(capture.alpha[1, i] - val) < 1e-13
      end
      if iface_i.elementR <= shockmesh.numShock
        @test abs(capture.alpha[2, i] - val) < 1e-13
      end
    end

    for i=1:shockmesh.numBoundaryFaces
      @test abs(capture.alpha_b[i] - val) < 1e-13
    end
  end  # end if




  # set w_el to be polynomial
  setWPoly(mesh, shockmesh, capture.w_el, degree)
  setWPolyDeriv(mesh, shockmesh, capture.grad_w, degree)

  fill!(eqn.res, 0)
  EulerEquationMod.computeFaceTerm(mesh, sbp, eqn, opts, capture, shockmesh,
                                   capture.diffusion, capture.penalty)

  # because interpolation and differentiation are exact for polynomials, this
  # term should be zero
  @test maximum(abs.(eqn.res)) < 1e-12

  # test applyPenalty
  delta_w = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
  theta = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)
  res_wL = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)
  res_wR = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)
  res_thetaL = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)
  res_thetaR = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace)

  tmp_el = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement)
  for i=1:shockmesh.numInterfaces
    iface_red = shockmesh.ifaces[i].iface
    idx_orig = shockmesh.ifaces[i].idx_orig
    elnumL = shockmesh.elnums_all[iface_red.elementL]
    elnumR = shockmesh.elnums_all[iface_red.elementL]

    wL = ro_sview(capture.w_el, :, :, iface_red.elementL)
    wR = ro_sview(capture.w_el, :, :, iface_red.elementR)
    nrm_face = ro_sview(mesh.nrm_face, :, :, idx_orig)
    alphas = ro_sview(capture.alpha, :, i)
    jacL = ro_sview(mesh.jac, :, elnumL)
    jacR = ro_sview(mesh.jac, :, elnumR)

    # test T3 by checking consistency with boundaryIntegrate
    rand_realpart!(delta_w); fill!(res_wL, 0); fill!(res_wR, 0)
    fill!(theta, 0); fill!(res_thetaL, 0); fill!(res_thetaR, 0)
    EulerEquationMod.applyPenalty(capture.penalty, sbp, mesh.sbpface,
                                  capture.diffusion, iface_red, delta_w,
                                  theta, wL, wR, nrm_face, alphas,
                                  jacL, jacR, res_wL, res_wR,
                                  res_thetaL, res_thetaR)
    fill!(tmp_el, 0)
    boundaryFaceIntegrate!(mesh.sbpface, iface_red.faceL, delta_w, tmp_el)
    for k=1:mesh.numDofPerNode
      @test maximum(abs.(2*sum(res_thetaL[k, :]) + sum(tmp_el[k, :]))) < 1e-12
      @test maximum(abs.(2*sum(res_thetaR[k, :]) - sum(tmp_el[k, :]))) < 1e-12
    end

    # test T2 by checking consistency with boundaryIntegrate
    fill!(delta_w, 0); fill!(res_wL, 0); fill!(res_wR, 0)
    rand_realpart!(theta); fill!(res_thetaL, 0); fill!(res_thetaR, 0)

    EulerEquationMod.applyPenalty(capture.penalty, sbp, mesh.sbpface,
                                  capture.diffusion, iface_red, delta_w,
                                  theta, wL, wR, nrm_face, alphas, 
                                  jacL, jacR, res_wL, res_wR,
                                  res_thetaL, res_thetaR)
    fill!(tmp_el, 0)
    boundaryFaceIntegrate!(mesh.sbpface, iface_red.faceL, theta, tmp_el)
    for k=1:mesh.numDofPerNode
      @test maximum(abs.(2*sum(res_wL[k, :]) - sum(tmp_el[k, :]))) < 1e-12
      @test maximum(abs.(2*sum(res_wR[k, :]) - sum(tmp_el[k, :]))) < 1e-12
    end

    # test entropy stability of T1
    rand_realpart!(delta_w); fill!(res_wL, 0); fill!(res_wR, 0)
    fill!(theta, 0); fill!(res_thetaL, 0); fill!(res_thetaR, 0)

    EulerEquationMod.applyPenalty(capture.penalty, sbp, mesh.sbpface,
                                  capture.diffusion, iface_red, delta_w,
                                  theta, wL, wR, nrm_face, alphas,
                                  jacL, jacR, res_wL, res_wR,
                                  res_thetaL, res_thetaR)
    for k=1:mesh.numDofPerNode
      @test dot(delta_w[k, :], res_wL[k, :]) > 0
      @test -dot(delta_w[k, :], res_wR[k, :]) > 0
    end

  end

  return nothing
end


function test_br2_Dgk(mesh, sbp, eqn::EulerData{Tsol, Tres}, opts) where {Tsol, Tres}

  degree = sbp.degree
  iface_idx, shockmesh = getInteriorMesh(mesh, sbp, eqn, opts)

  capture = EulerEquationMod.SBPParabolicSC{Tsol, Tres}(mesh, sbp, eqn, opts)
  EulerEquationMod.allocateArrays(capture, mesh, shockmesh)

  fill!(eqn.res, 0)
  # set w_el to be polynomial
  setWPoly(mesh, shockmesh, capture.w_el, degree)
  setWPolyDeriv(mesh, shockmesh, capture.grad_w, degree)

  w_face = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
  for i=1:shockmesh.numInterfaces
    iface = shockmesh.ifaces[i].iface
    idx_orig = shockmesh.ifaces[i].idx_orig
    elnumL = shockmesh.elnums_all[iface.elementL]
    elnumR = shockmesh.elnums_all[iface.elementR]

    wL = ro_sview(capture.w_el, :, :, iface.elementL)
    wR = ro_sview(capture.w_el, :, :, iface.elementR)
    nrm_face = ro_sview(mesh.nrm_face, :, :, idx_orig)
    dxidxL = ro_sview(mesh.dxidx, :, :, :, elnumL)
    dxidxR = ro_sview(mesh.dxidx, :, :, :, elnumR)
    jacL = ro_sview(mesh.jac, :, elnumL)
    jacR = ro_sview(mesh.jac, :, elnumR)
    resL = sview(eqn.res, :, :, elnumL)
    resR = sview(eqn.res, :, :, elnumR)

    # the solution is polynomial, so it doesn't matter which element we
    # interpolate from (except for the nbrperm stuff)
    boundaryFaceInterpolate!(mesh.sbpface, iface.faceL, wL, w_face)
    for j=1:mesh.numNodesPerFace
      for k=1:mesh.numDofPerNode
        w_face[k, j] *= mesh.sbpface.wface[j]
      end
    end

    EulerEquationMod.applyDgkTranspose(capture, sbp, mesh.sbpface, iface,
                  capture.diffusion, w_face, w_face, wL, wR, nrm_face, dxidxL, dxidxR,
                  jacL, jacR, resL, resR)

  end

  # now use applyE for every element, then apply Dx^T and Dy^T

  work = zeros(Tsol, mesh.numDofPerNode, mesh.numNodesPerFace)
  work2 = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerFace, mesh.dim)
  work3 = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
  E_term = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.dim)
  res2 = zeros(Tres, mesh.numDofPerNode, mesh.numNodesPerElement)
  for i=1:shockmesh.numShock
    i_full = shockmesh.elnums_all[i]
    w_i = ro_sview(capture.w_el, :, :, i)

    applyE(mesh, i_full, sview(iface_idx, :, i_full), w_i, work, work2, E_term)

    dxidx_i = ro_sview(mesh.dxidx, :, :, :, i_full)
    jac_i = ro_sview(mesh.jac, :, i_full)

    fill!(res2, 0)
    EulerEquationMod.applyDxTransposed(sbp, E_term, dxidx_i, jac_i, work3, res2)

    @test maximum(abs.(res2 - eqn.res[:, :, i_full])) < 1e-12
  end

  return nothing
end


function test_br2_ESS(mesh, sbp, eqn::EulerData{Tsol, Tres}, _opts; fullmesh=false) where {Tsol, Tres}

  opts = copy(_opts)
  # the scheme is entropy stable when the dirichlet BC = 0
  for i=1:opts["numBC"]
    opts[string("BC", i, "_name")] = "zeroBC"
  end

  # construct the shock mesh with all elements in it that are fully interior
  if fullmesh
    iface_idx, shockmesh = getEntireMesh(mesh, sbp, eqn, opts)
  else
    iface_idx, shockmesh = getInteriorMesh(mesh, sbp, eqn, opts)
  end

  capture = EulerEquationMod.SBPParabolicSC{Tsol, Tres}(mesh, sbp, eqn, opts)
  EulerEquationMod.allocateArrays(capture, mesh, shockmesh)

  # add random component to q
  # test entropy stability
  q_pert = 0.1*rand(size(eqn.q_vec))
  eqn.q_vec .+= q_pert
  array1DTo3D(mesh, sbp, eqn, opts, eqn.q_vec, eqn.q)
  fill!(eqn.res, 0)

  EulerEquationMod.calcShockCapturing(mesh, sbp, eqn, opts, capture, shockmesh)

  array3DTo1D(mesh, sbp, eqn, opts, eqn.res, eqn.res_vec)

  w_vec = zeros(Tsol, mesh.numDof)
  copy!(w_vec, eqn.q_vec)
  EulerEquationMod.convertToIR(mesh, sbp, eqn, opts, w_vec)

  val = dot(w_vec, eqn.res_vec)
  println("val = ", val)
  @test val < 0

  return nothing
end


function test_br2_serialpart(mesh, sbp, eqn::EulerData{Tsol, Tres}, _opts) where {Tsol, Tres}


  opts = copy(_opts)
  sensor = EulerEquationMod.ShockSensorEverywhere{Tsol, Tres}(mesh, sbp, opts)
  capture = EulerEquationMod.SBPParabolicSC{Tsol, Tres}(mesh, sbp, eqn, opts)

  # solve the PDE to get a solution with non-zero jump between elements
  # that can be reproduced in parallel
  opts["solve"] = true
  opts["addVolumeIntegrals"] = true
  opts["addFaceIntegrals"] = true
  opts["addBoundaryIntegrals"] = true

  solvePDE(mesh, sbp, eqn, opts)

  # the scheme is entropy stable when the dirichlet BC = 0
  for i=1:opts["numBC"]
    opts[string("BC", i, "_name")] = "zeroBC"
  end
  EulerEquationMod.getBCFunctors(mesh, sbp, eqn, opts)

  w_vec = zeros(Tsol, mesh.numDof)
  copy!(w_vec, eqn.q_vec)
  EulerEquationMod.convertToIR(mesh, sbp, eqn, opts, w_vec)

  fill!(eqn.res, 0)
  EulerEquationMod.applyShockCapturing(mesh, sbp, eqn, opts, sensor, capture)

  val = dot(w_vec, eqn.res_vec)
  println("val = ", val)
  @test val < 0

  # save this for parallel tests
  f = open("br2_entropy_serial.dat", "w")
  println(f, real(val))
  close(f)

  saveSolutionToMesh(mesh, eqn.q_vec)
  writeVisFiles(mesh, "br2_serial")

  return nothing
end


