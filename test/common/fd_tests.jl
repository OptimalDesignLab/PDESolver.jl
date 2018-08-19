# finite difference tests

"""
  Performes finite difference tests to check that coordinate/metric/normal
  calculation of the given mesh can be complex stepped

  Tests are available for the following fields of the mesh:

   * coords
   * coords_interface
   * coords_bndry
   * dxidx
   * jac
   * nrm_face
   * nrm_bndry

TODO: add parallel tests
"""
function testMeshDerivatives(mesh, sbp, opts)

  deriv_cs = computeMeshDeriv(mesh, sbp, opts, mesh.coords)
  deriv_fd = computeMeshDeriv_fd(mesh, sbp, opts, mesh.coords)
  @test maximum(abs.(deriv_cs - deriv_fd)) < 1e-5

  deriv_cs = computeMeshDeriv(mesh, sbp, opts, mesh.coords_interface)
  deriv_fd = computeMeshDeriv_fd(mesh, sbp, opts, mesh.coords_interface)
  @test maximum(abs.(deriv_cs - deriv_fd)) < 1e-5

  deriv_cs = computeMeshDeriv(mesh, sbp, opts, mesh.coords_bndry)
  deriv_fd = computeMeshDeriv_fd(mesh, sbp, opts, mesh.coords_bndry)
  @test maximum(abs.(deriv_cs - deriv_fd)) < 1e-5

  deriv_cs = computeMeshDeriv(mesh, sbp, opts, mesh.dxidx)
  deriv_fd = computeMeshDeriv_fd(mesh, sbp, opts, mesh.dxidx)
  @test maximum(abs.(deriv_cs - deriv_fd)) < 1e-5

  deriv_cs = computeMeshDeriv(mesh, sbp, opts, mesh.jac)
  deriv_fd = computeMeshDeriv_fd(mesh, sbp, opts, mesh.jac)
  @test maximum(abs.(deriv_cs - deriv_fd)) < 3e-4

  deriv_cs = computeMeshDeriv(mesh, sbp, opts, mesh.nrm_face)
  deriv_fd = computeMeshDeriv_fd(mesh, sbp, opts, mesh.nrm_face)
  @test maximum(abs.(deriv_cs - deriv_fd)) < 1e-5

  deriv_cs = computeMeshDeriv(mesh, sbp, opts, mesh.nrm_bndry)
  deriv_fd = computeMeshDeriv_fd(mesh, sbp, opts, mesh.nrm_bndry)
  @test maximum(abs.(deriv_cs - deriv_fd)) < 1e-5

  return nothing
end




"""
  Compute the derivative of the given mesh field with respect to the vertex
  coordinates (note that this is the derivative with respect to the 1D vector
  representation of the vertex coordinates, not mesh.vert_coord, which duplicates
  nodes).

  Uses the complex step method to compute the derivative.  The mesh must have
  been created with Tmsh = Complex128.

  **Inputs**

   * mesh
   * sbp
   * opts
   * field: a field of the mesh object that is in the "numerator" of the
            derivative.  Must be the array that is stored in the mesh, not
            a copy.
"""
function computeMeshDeriv(mesh, sbp, opts, field::AbstractArray)

  @assert eltype(mesh.vert_coords) == Complex128

  h = 1e-20
  pert = Complex128(0, h)
  numdof = mesh.dim*mesh.coord_numNodes

  dRdx = zeros(length(field), numdof)


  coords_vec = zeros(Complex128, numdof)
  coords3DTo1D(mesh, mesh.vert_coords, coords_vec, PdePumiInterface.AssignReduction{Complex128}())

  for i=1:numdof
    coords_vec[i] += pert
    coords1DTo3D(mesh, coords_vec, mesh.vert_coords)
    recalcCoordinatesAndMetrics(mesh, sbp, opts)
 
    for j=1:length(field)
      dRdx[j, i] = imag(field[j])/h
    end

    coords_vec[i] -= pert
    coords1DTo3D(mesh, coords_vec, mesh.vert_coords)
    recalcCoordinatesAndMetrics(mesh, sbp, opts)
  end

  return dRdx
end


"""
  Similar to [`computeMeshDeriv`](@ref), but computes the quantity using
  finite differences
"""
function computeMeshDeriv_fd(mesh::AbstractMesh{Tmsh}, sbp, opts, field::AbstractArray) where {Tmsh}

  h = 1e-6
  numdof = mesh.dim*mesh.coord_numNodes

  dRdx = zeros(length(field), numdof)

  coords_vec = zeros(Tmsh, numdof)
  coords3DTo1D(mesh, mesh.vert_coords, coords_vec, PdePumiInterface.AssignReduction{Tmsh}())

  R0 = copy(field)

  for i=1:numdof
    coords_vec[i] += h
    coords1DTo3D(mesh, coords_vec, mesh.vert_coords)
    recalcCoordinatesAndMetrics(mesh, sbp, opts)
 
    res_vec_i = sview(dRdx, :, i)
    R1 = copy(field)
    for j=1:length(R1)
      res_vec_i[j] = (R1[j] - R0[j])/h      
    end


    coords_vec[i] -= h
    coords1DTo3D(mesh, coords_vec, mesh.vert_coords)
    recalcCoordinatesAndMetrics(mesh, sbp, opts)
  end

  return dRdx
end



#------------------------------------------------------------------------------
# test derivative of residual with respect to mesh quantities

"""
  Compares finite difference and complex step derivatives of the residual
  with respect to mesh quantities.  Checks the same quantities as
  [`testMeshDerivatives`](@ref)
"""
function testResidualDerivatives(mesh, eqn, sbp, opts)

  deriv_cs = computeResidualDeriv(mesh, eqn, sbp, opts, mesh.coords)
  deriv_fd = computeResidualDeriv_fd(mesh, eqn, sbp, opts, mesh.coords)
  @test maximum(abs.(deriv_cs - deriv_fd)) < 1e-5

  deriv_cs = computeResidualDeriv(mesh, eqn, sbp, opts, mesh.coords_interface)
  deriv_fd = computeResidualDeriv_fd(mesh, eqn, sbp, opts, mesh.coords_interface)
  @test maximum(abs.(deriv_cs - deriv_fd)) < 1e-5

  deriv_cs = computeResidualDeriv(mesh, eqn, sbp, opts, mesh.coords_bndry)
  deriv_fd = computeResidualDeriv_fd(mesh, eqn, sbp, opts, mesh.coords_bndry)
  @test maximum(abs.(deriv_cs - deriv_fd)) < 1e-5

  deriv_cs = computeResidualDeriv(mesh, eqn, sbp, opts, mesh.dxidx)
  deriv_fd = computeResidualDeriv_fd(mesh, eqn, sbp, opts, mesh.dxidx)
  @test maximum(abs.(deriv_cs - deriv_fd)) < 1e-5

  deriv_cs = computeResidualDeriv(mesh, eqn, sbp, opts, mesh.jac)
  deriv_fd = computeResidualDeriv_fd(mesh, eqn, sbp, opts, mesh.jac)
  @test maximum(abs.(deriv_cs - deriv_fd)) < 1e-5

  deriv_cs = computeResidualDeriv(mesh, eqn, sbp, opts, mesh.nrm_face)
  deriv_fd = computeResidualDeriv_fd(mesh, eqn, sbp, opts, mesh.nrm_face)
  @test maximum(abs.(deriv_cs - deriv_fd)) < 1e-5

  deriv_cs = computeResidualDeriv(mesh, eqn, sbp, opts, mesh.nrm_bndry)
  deriv_fd = computeResidualDeriv_fd(mesh, eqn, sbp, opts, mesh.nrm_bndry)
  @test maximum(abs.(deriv_cs - deriv_fd)) < 1e-5

  return nothing
end


"""
  Computes the derivative of the residual with respect to a given field of the
  mesh.  See [`computeMeshDeriv`](@ref), has a similar setup, except it computes
  derivatives of mesh quantities with respect to mesh vertex coordinates.

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * field: field of the mesh object.  Cannot be a copy of the array.
"""
function computeResidualDeriv(mesh, sbp, eqn, opts, field::AbstractArray)

  h = 1e-20
  pert = Complex128(0, h)
#  numdof = mesh.dim*mesh.coord_numNodes
  numdof = length(field)
  dRdx = zeros(mesh.numDof, numdof)

  for i=1:numdof
    field[i] += pert
 
    evalResidual(mesh, sbp, eqn, opts)
    res_vec_i = sview(dRdx, :, i)
    array3DTo1D(mesh, sbp, eqn, opts, imag(eqn.res)/h, res_vec_i)


    field[i] -= pert
  end

  return dRdx
end


function computeResidualDeriv_fd(mesh, sbp, eqn, opts, field::AbstractArray)

  h = 1e-6
  numdof = length(field)
  dRdx = zeros(mesh.numDof, numdof)

  R0 = zeros(eqn.res_vec)
  R1 = zeros(eqn.res_vec)
  evalResidual(mesh, sbp, eqn, opts)
  array3DTo1D(mesh, sbp, eqn, opts, eqn.res, R0)

  for i=1:numdof
    field[i] += h
 
    evalResidual(mesh, sbp, eqn, opts)
    array3DTo1D(mesh, sbp, eqn, opts, eqn.res, R1)
    res_vec_i = sview(dRdx, :, i)
    for j=1:length(R1)
      res_vec_i[j] = (R1[j] - R0[j])/h
    end


    field[i] -= h
  end

  return dRdx
end


#------------------------------------------------------------------------------
# Test gradient computed via adjoint

"""
  Tests the total derivative of a functional with respect to the
  mesh vertex coordinates against the total derivative computed via adjoint

  This solves the PDE once for every coordinate of every node in the mesh.
  This is very expensive unless done on a tiny mesh

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
   * functional

  
"""
function testDJDx(mesh, sbp, eqn, opts, func::AbstractFunctional)

  solvePDE(mesh, sbp, eqn, opts)

  deriv_cs = computeDJDx_adjoint(mesh, sbp, eqn, opts, func)
  deriv_fd = computeDJDx_adjoint(mesh, sbp, eqn, opts, func)
  @test maximum(abs.(deriv_cs - deriv_fd)) < 1e-5

  return nothing
end


"""
  Computes total derivative of a functional J with respect ot the mesh vertex
  coordinates.
"""
function computeDJDx_fd(mesh, sbp, eqn, opts, func::AbstractFunctional)

  h = 1e-6
  J0 = evalFunctional(mesh, sbp, eqn, opts, func)

  numdof = mesh.dim*mesh.coord_numNodes
  coords_vec = zeros(Complex128, numdof)
  coords3DTo1D(mesh, mesh.vert_coords, coords_vec, PdePumiInterface.AssignReduction{Complex128}())

  DJDx = zeros(numdof)

  for i=1:numdof
    coords_vec[i] += h
    coords1DTo3D(mesh, coords_vec, mesh.vert_coords)
    recalcCoordinatesAndMetrics(mesh, sbp, opts)
    solvePDE(mesh, sbp, eqn, opts)

    J1 = evalFunctional(mesh, sbp, eqn, opts, func)
    DJDx[i] = (J1 - J0)/h

    coords_vec[i] -= h
    coords1DTo3D(mesh, coords_vec, mesh.vert_coords)
    recalcCoordinatesAndMetrics(mesh, sbp, opts)
    solvePDE(mesh, sbp, eqn, opts)
  end

  return real(DJDx)
end


"""
  Computes the total derivative of a functional J with respect to the mesh vertex
  coordinates using the (discrete) adjoint
"""
function computeDJDx_adjoint(mesh, sbp, eqn, opts, func::AbstractFunctional)

  dJdx = computedJdx_vec(mesh, sbp, eqn, opts, func)
  dRdx = computedRdx_vec(mesh, sbp, eqn, opts)
  
  psi_discrete = zeros(eqn.q_vec)
  calcAdjoint(mesh, sbp, eqn, opts, func, psi_discrete)

  DJDx = dJdx + dRdx.'*psi_discrete

  return real(DJDx)
end


"""
  Compute the vector form of dJdx.  The length is mesh.dim*mesh.coord_numNodes
"""
function computedJdx_vec(mesh, sbp, eqn, opts, func::AbstractFunctional)

  @assert eltype(mesh.vert_coords) == Complex128
  h = 1e-20
  pert = Complex128(0, h)

  numdof = mesh.dim*mesh.coord_numNodes
  coords_vec = zeros(Complex128, numdof)
  coords3DTo1D(mesh, mesh.vert_coords, coords_vec, PdePumiInterface.AssignReduction{Complex128}())

  dJdx = zeros(numdof)
  for i=1:numdof
    coords_vec[i] += pert
    coords1DTo3D(mesh, coords_vec, mesh.vert_coords)
    recalcCoordinatesAndMetrics(mesh, sbp, opts)
    
    J_val = evalFunctional(mesh, sbp, eqn, opts, func)
    dJdx[i] = imag(J_val)/h

    coords_vec[i] -= pert
    coords1DTo3D(mesh, coords_vec, mesh.vert_coords)
    recalcCoordinatesAndMetrics(mesh, sbp, opts)
  end

  return dJdx
end


"""
  Compute derivative of residual with repsect to the mesh vertex coordinates
  The result is a matrix mesh.numDof x mesh.dim*mesh.coord_numNodes
"""
function computedRdx_vec(mesh, sbp, eqn, opts)

  h = 1e-20
  pert = Complex128(0, h)
  numdof = mesh.dim*mesh.coord_numNodes
  dRdx = zeros(mesh.numDof, numdof)

  coords_vec = zeros(Complex128, numdof)
  coords3DTo1D(mesh, mesh.vert_coords, coords_vec, PdePumiInterface.AssignReduction{Complex128}())

  for i=1:numdof
    coords_vec[i] += pert
    coords1DTo3D(mesh, coords_vec, mesh.vert_coords)
    recalcCoordinatesAndMetrics(mesh, sbp, opts)
 
    evalResidual(mesh, sbp, eqn, opts)
    res_vec_i = sview(dRdx, :, i)
    array3DTo1D(mesh, sbp, eqn, opts, imag(eqn.res)/h, res_vec_i)


    coords_vec[i] -= pert
    coords1DTo3D(mesh, coords_vec, mesh.vert_coords)
    recalcCoordinatesAndMetrics(mesh, sbp, opts)
  end

  return dRdx
end


