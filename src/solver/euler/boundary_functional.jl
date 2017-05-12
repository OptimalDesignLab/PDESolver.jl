abstract FuntionalType
# Calculate the analytical force on the inner boundary of the isentropic vortex
#=
function calc_analytical_forces{Tmsh}(mesh::AbstractMesh{Tmsh}, params::ParamType{2},
	                                    coords::AbstractArray{Tmsh})

  q = zeros(mesh.numDofPerNode)
  calcIsentropicVortex(coords, params, q)  # Get analytical q ath the coordinates
  p = calcPressure(params, q) # Get the analytical pressure
  r = sqrt(coords[1]*coords[1] + coords[2]*coords[2]) # get the curcular radius
  force = 0.5*pi*r*p

  return force
end
=#
@doc """
### EulerEquationMod.calcBndryFunctional

This function calculates a functional on a geometric boundary of a the 
computational space. There is no need to call this function withing the 
nonlinear solve while computing eqn.q

**Inputs**

*  `mesh` :  Abstract mesh object
*  `sbp`  : Summation-By-Parts operator
*  `eqn`  : Euler equation object
*  `opts` : Options dictionary
*  `g_edge_number` : Geometric edge number

**Outputs**

*  `functional_val` : computed numerical force at the boundary.

"""->

function calcBndryFunctional{Tmsh, Tsol}(mesh::AbstractDGMesh{Tmsh},
                                         sbp::AbstractSBP,
                                         eqn::EulerData{Tsol}, 
                                         opts, 
                                         functor, 
                                         functional_edges)

  functional_val = zero(Tsol)


  for itr = 1:length(functional_edges)
    g_edge_number = functional_edges[itr] # Extract geometric edge number
    start_index = mesh.bndry_offsets[g_edge_number]
    end_index = mesh.bndry_offsets[g_edge_number+1]
    idx_range = start_index:(end_index-1)
    bndry_facenums = sview(mesh.bndryfaces, idx_range) # faces on geometric edge i

    nfaces = length(bndry_facenums)
    boundary_integrand = zeros(Tsol, 1, mesh.sbpface.numnodes, nfaces)

    for i = 1:nfaces
      global_facenum = idx_range[i]

      face_integrand = slice(boundary_integrand, 1, :, i)
      functor(mesh, sbp, eqn, opts, global_facenum, face_integrand)
    end    # End for i = 1:nfaces

    val_per_geom_edge = zeros(Tsol, 1)

    integratefunctional!(mesh.sbpface, mesh.bndryfaces[idx_range], 
                         boundary_integrand, val_per_geom_edge)

    functional_val += val_per_geom_edge[1]

  end  # End for itr = 1:length(functional_edges)

  # scale because of nondimentionalization
  Ma = eqn.params.Ma
  functional_val /= Ma*Ma

  return functional_val
end

@doc """
### EulerEquationMod.drag

Computes the force in the X-direction in an ADJOINT CONSISTENT way

**Inputs**

*  `params` : Parameter type
*  `q`      : Solution at a node
*  `aux_vars` : Vector of auxiliary variables
*  `nrm`    : Normal vector in the physical space

**Outputs**

*  `val`    : Momentum derivative in the X-direction

"""->

function drag_inviscid{Tsol, Tmsh}(mesh::AbstractDGMesh{Tmsh},
                                         sbp::AbstractSBP,
                                         eqn::EulerData{Tsol},
                                         opts,
                                         bndry_indx,
                                         face_integrand::AbstractArray{Tsol, 1})
  dir = 2

  # since the inviscid boundary flux is computed as
  # F_b = F(u_gamma), we should use the same boundary flux
  # in computation of functional
  num_dof = mesh.numDofPerNode
  dim = num_dof - 2
  q2 = Array(Tsol, num_dof)
  x = Array(Tmsh, 2)
  Finv = Array(Tsol, num_dof)
  sbpface = mesh.sbpface
  nrm_xi = sview(sbpface.normal, :, mesh.bndryfaces[bndry_indx].face)
  # BCType is a single node function, so we have to
  # loop over all face nodes
  for j = 1 : mesh.sbpface.numnodes
    q = sview(eqn.q_bndry, :, j, bndry_indx)
    convertToConservative(eqn.params, q, q2)
    aux_vars = sview(eqn.aux_vars_bndry, :, j, bndry_indx)
    dxidx = sview(mesh.dxidx_bndry, :, :, j, bndry_indx)

    flux_func = noPenetrationBC()
    flux_func(q, aux_vars, x, dxidx, nrm_xi, Finv, eqn.params)
    face_integrand[j] = Finv[dir]
    println(real(Finv))
  end
  return nothing
end


function drag_viscous{Tsol, Tmsh}(mesh::AbstractDGMesh{Tmsh},
                                  sbp::AbstractSBP,
                                  eqn::EulerData{Tsol},
                                  opts,
                                  bndry_indx,
                                  face_integrand::AbstractArray{Tsol, 1})

  num_dof       = mesh.numDofPerNode
  num_facenodes = mesh.numNodesPerFace
  num_elemnodes = mesh.numNodesPerElement
  Tdim  = size(eqn.q, 1) - 2
  q2    = zeros(Tsol, num_dof)
  bndry = mesh.bndryfaces[bndry_indx]
  sbpface = mesh.sbpface
  x       = Array(Tmsh, 2)
  Finv    = Array(Tsol, num_dof)
  nrm_xi  = sview(sbpface.normal, :, bndry.face)

  dir = 2

  # contribution from inviscid flux function.
  # since the inviscid boundary flux is computed as
  #   f_b = f(u_gamma), 
  # we should use the same boundary flux
  # to compute functional

  # loop over all face nodes
  for j = 1 : mesh.sbpface.numnodes
    q = sview(eqn.q_bndry, :, j, bndry_indx)
    convertToConservative(eqn.params, q, q2)
    aux_vars = sview(eqn.aux_vars_bndry, :, j, bndry_indx)
    dxidx = sview(mesh.dxidx_bndry, :, :, j, bndry_indx)

    flux_func = nonslipBC()
    flux_func(q, aux_vars, x, dxidx, nrm_xi, Finv, eqn.params)
    face_integrand[j] = Finv[dir]
  end

  # Contribution from viscous flux function.
  # Since we don't have functions for single node, first we
  # tranform data into array
  Ma = eqn.params.Ma
  Re = eqn.params.Re
  coef_nondim = Ma/Re 
  Fvis   = Array(Tsol, num_dof, num_facenodes)
  q_face = sview(eqn.q_bndry, :, :, bndry_indx)
  q_bnd  = Array(Tsol, num_dof, num_facenodes)
  dq     = Array(Tsol, num_dof, num_facenodes)
  Fvis   = Array(Tsol, Tdim, num_dof, num_facenodes)
  Gt     = zeros(Tsol, num_dof, num_dof, Tdim, Tdim, num_facenodes)
  nrm    = Array(Tmsh, Tdim, num_facenodes)
  nrm1   = Array(Tmsh, Tdim, num_facenodes)
  area   = Array(Tmsh, num_facenodes)

  for n = 1 : mesh.numNodesPerFace
    dxidx = sview(mesh.dxidx_bndry, :, :, n, bndry_indx)
    nrm[1,n] = dxidx[1, 1]*nrm_xi[1] + dxidx[2, 1]*nrm_xi[2]
    nrm[2,n] = dxidx[1, 2]*nrm_xi[1] + dxidx[2, 2]*nrm_xi[2]
    area[n] = sqrt(nrm[1,n]*nrm[1,n] + nrm[2,n]*nrm[2,n])
    nrm1[1,n] = nrm[1,n]/area[n]
    nrm1[2,n] = nrm[2,n]/area[n]
  end

  elem = bndry.element
  elem_vol = 0.0
  for n = 1:mesh.numNodesPerElement
    elem_vol +=  sbp.w[n]/mesh.jac[n, elem]
  end

  face_area = 0.0
  for n = 1:mesh.numNodesPerFace
    face_area +=  sbpface.wface[n]*area[n]
  end

  he = elem_vol/eqn.area_sum[elem]

  # for viscous flow, we need dqdx on face
  dqdx_face = Array(Tsol, Tdim, num_dof, num_facenodes)
  dqdx_elem = Array(Tsol, Tdim, num_dof, num_elemnodes)
  q_elem = sview(eqn.q, :, :, elem)
  calcGradient(mesh, sbp, elem, q_elem, dqdx_elem)
  for d = 1 : Tdim
    q_x_node = slice(dqdx_elem, d, :, :)
    q_x_face = slice(dqdx_face, d, :, :)
    boundaryinterpolate(sbpface, bndry, q_x_node, q_x_face) 
  end

  bndry_val_functor = AdiabaticWall()
  bndry_val_functor(q_face, nrm1, eqn.params, q_bnd)
  calcDiffusionTensor_adiabaticWall(q_bnd, Gt)
  calcFvis(Gt, dqdx_face, Fvis)
  
  p = opts["order"]
  Cip = opts["Cip"]
  const_tii = (p + 1.0)*(p + Tdim)/Tdim * Cip
  for n = 1 : mesh.numNodesPerFace
    for iDof = 1 : mesh.numDofPerNode
      dq[iDof, n] = q_face[iDof, n] - q_bnd[iDof, n]
    end
  end

  for n = 1 : num_facenodes
    nGtn = (nrm1[1,n]*nrm1[1,n]*Gt[dir,:,1,1,n]
          + nrm1[1,n]*nrm1[2,n]*Gt[dir,:,1,2,n]
          + nrm1[2,n]*nrm1[1,n]*Gt[dir,:,2,1,n]
          + nrm1[2,n]*nrm1[2,n]*Gt[dir,:,2,2,n])
    nGtn_dq = nGtn[1]*dq[1] + nGtn[2]*dq[2] + nGtn[3]*dq[3] + nGtn[4]*dq[4]
    integrand = const_tii /he * nGtn_dq
    integrand -= Fvis[1,dir,n] * nrm1[1,n] + Fvis[2,dir,n] * nrm1[2,n]
    face_integrand[n] += area[n] * coef_nondim * integrand
  end

  return nothing
end

@doc """
### EulerEquationMod.drag

Computes the force in the X-direction in an ADJOINT CONSISTENT way

**Inputs**

*  `params` : Parameter type
*  `q`      : Solution at a node
*  `aux_vars` : Vector of auxiliary variables
*  `nrm`    : Normal vector in the physical space

**Outputs**

*  `val`    : Momentum derivative in the X-direction

"""->

function lift_inviscid{Tsol, Tmsh}(mesh::AbstractDGMesh{Tmsh}, 
                                   sbp::AbstractSBP,
                                   eqn::EulerData{Tsol},
                                   opts,
                                   bndry_indx,
                                   face_integrand::AbstractArray{Tsol, 1})

  dir = 3
  # since the inviscid boundary flux is computed as
  # F_b = F(u_gamma), we should use the same boundary flux
  # in computation of functional
  num_dof = mesh.numDofPerNode
  dim = num_dof - 2
  q2 = zeros(Tsol, num_dof)
  x = Array(Tmsh, 2)
  Finv = Array(Tsol, num_dof)
  sbpface = mesh.sbpface
  nrm_xi = sview(sbpface.normal, :, mesh.bndryfaces[bndry_indx].face)
  # BCType is a single node function, so we have to
  # loop over all face nodes
  for j = 1 : mesh.sbpface.numnodes
    q = sview(eqn.q_bndry, :, j, bndry_indx)
    convertToConservative(eqn.params, q, q2)
    aux_vars = sview(eqn.aux_vars_bndry, :, j, bndry_indx)
    dxidx = sview(mesh.dxidx_bndry, :, :, j, bndry_indx)

    flux_func = noPenetration()
    flux_func(q, aux_vars, x, dxidx, nrm_xi, Finv, eqn.params)
    face_integrand[j] = Finv[dir]
  end
  return nothing
end

# type lift_viscous <: FunctionalType
# end

# function call{Tsol, Tres, Tmsh}(obj::lift_viscous, 

function lift_viscous{Tsol, Tmsh}(mesh::AbstractDGMesh{Tmsh},
                                        sbp::AbstractSBP,
                                        eqn::EulerData{Tsol},
                                        opts,
                                        bndry_indx,
                                        face_integrand::AbstractArray{Tsol, 1})

  dir = 3
  # since the inviscid boundary flux is computed as
  # F_b = F(u_gamma), we should use the same boundary flux
  # in computation of functional
  num_dof = mesh.numDofPerNode
  num_facenodes = mesh.numNodesPerFace
  num_elemnodes = mesh.numNodesPerElement
  Tdim = size(eqn.q, 1) - 2
  q2 = zeros(Tsol, num_dof)
  bndry = mesh.bndryfaces[bndry_indx]
  sbpface = mesh.sbpface
  x = Array(Tmsh, 2)
  Finv = Array(Tsol, num_dof)
  nrm_xi = sview(sbpface.normal, :, bndry.face)

  # contribution from inviscid flux function
  # BCType is a single node function, so we have to
  # loop over all face nodes
  for j = 1 : mesh.sbpface.numnodes
    q = sview(eqn.q_bndry, :, j, bndry_indx)
    convertToConservative(eqn.params, q, q2)
    aux_vars = sview(eqn.aux_vars_bndry, :, j, bndry_indx)
    dxidx = sview(mesh.dxidx_bndry, :, :, j, bndry_indx)

    bndry_val_functor = nonslipBC()
    bndry_val_functor(q, aux_vars, x, dxidx, nrm_xi, Finv, eqn.params)
    face_integrand[j] = Finv[dir]
    # face_integrand[j] = 0
  end

  # contribution from viscous flux function
  # since we don't have functions for single node, first we
  # tranform data into array
  Ma = eqn.params.Ma
  Re = eqn.params.Re
  coef_nondim = Ma/Re 
  Fvis = Array(Tsol, num_dof, num_facenodes)
  q_face = sview(eqn.q_bndry, :, :, bndry_indx)
  q_bnd = Array(Tsol, num_dof, num_facenodes)
  dq = Array(Tsol, num_dof, num_facenodes)
  Fvis = Array(Tsol, Tdim, num_dof, num_facenodes)
  q_bnd = Array(Tsol, num_dof, num_facenodes)
  Gt = zeros(Tsol, num_dof, num_dof, Tdim, Tdim, num_facenodes)
  nrm = Array(Tmsh, Tdim, num_facenodes)
  nrm1 = Array(Tmsh, Tdim, num_facenodes)
  area = Array(Tmsh, num_facenodes)

  for n = 1 : mesh.numNodesPerFace
    dxidx = sview(mesh.dxidx_bndry, :, :, n, bndry_indx)
    nrm[1,n] = dxidx[1, 1]*nrm_xi[1] + dxidx[2, 1]*nrm_xi[2]
    nrm[2,n] = dxidx[1, 2]*nrm_xi[1] + dxidx[2, 2]*nrm_xi[2]
    area[n] = sqrt(nrm[1,n]*nrm[1,n] + nrm[2,n]*nrm[2,n])
    nrm1[1,n] = nrm[1,n]/area[n]
    nrm1[2,n] = nrm[2,n]/area[n]
  end

  elem = bndry.element
  elem_vol = 0.0
  for n = 1:mesh.numNodesPerElement
    elem_vol +=  sbp.w[n]/mesh.jac[n, elem]
  end

  face_area = 0.0
  for n = 1:mesh.numNodesPerFace
    face_area +=  sbpface.wface[n]*area[n]
  end

  he = elem_vol/eqn.area_sum[elem]

  # for viscous flow, we need dqdx on face
  dqdx_face = Array(Tsol, Tdim, num_dof, num_facenodes)
  dqdx_elem = Array(Tsol, Tdim, num_dof, num_elemnodes)
  q_elem = sview(eqn.q, :, :, elem)
  calcGradient(mesh, sbp, elem, q_elem, dqdx_elem)
  for d = 1 : Tdim
    q_x_node = slice(dqdx_elem, d, :, :)
    q_x_face = slice(dqdx_face, d, :, :)
    boundaryinterpolate(sbpface, bndry, q_x_node, q_x_face) 
  end

  bndry_val_functor = AdiabaticWall()
  bndry_val_functor(q_face, nrm1, eqn.params, q_bnd)
  calcDiffusionTensor_adiabaticWall(q_bnd, Gt)
  calcFvis(Gt, dqdx_face, Fvis)
  p = opts["order"]
  Cip = opts["Cip"]
  const_tii = (p + 1.0)*(p + Tdim)/Tdim * Cip
  for n = 1 : mesh.numNodesPerFace
    for iDof = 1 : mesh.numDofPerNode
      dq[iDof, n] = q_face[iDof, n] - q_bnd[iDof, n]
    end
  end

  for n = 1 : num_facenodes
    nGtn = (nrm1[1,n]*nrm1[1,n]*Gt[dir,:,1,1,n]
          + nrm1[1,n]*nrm1[2,n]*Gt[dir,:,1,2,n]
          + nrm1[2,n]*nrm1[1,n]*Gt[dir,:,2,1,n]
          + nrm1[2,n]*nrm1[2,n]*Gt[dir,:,2,2,n])
    nGtn_dq = nGtn[1]*dq[1] + nGtn[2]*dq[2] + nGtn[3]*dq[3] + nGtn[4]*dq[4]
    integrand = const_tii /he * nGtn_dq
    integrand -= Fvis[1,dir,n] * nrm1[1,n] + Fvis[2,dir,n] * nrm1[2,n]
    face_integrand[n] += area[n] * coef_nondim * integrand
  end

  return nothing
end


@doc """
### EulerEquationMod.FunctionalDict

It stores the names of all possible functional options that can be computed. 
Whenever a new functional is created, it should be added to FunctionalDict.

"""->
global const FunctionalDict = Dict{ASCIIString, Function}(
# "drag" => drag(),
# "lift" => lift(),
"drag_inviscid" => drag_inviscid,
"drag_viscous" => drag_viscous,
"lift_inviscid" => lift_inviscid,
"lift_viscous" => lift_viscous,
)


@doc """
### EulerEquationMod.getFunctionalName

Gets the name of the functional that needs to be computed at a particular point

**Inputs**

*  `opts`     : Input dictionary
*  `f_number` : Number of the functional in the input dictionary

**Outputs**

*  `functional` : Returns the functional name in the dictionary

"""->
function getFunctionalName(opts, f_number)

  key = string("functional_name", f_number)
  val = opts[key]

  return functional = FunctionalDict[val]
end

