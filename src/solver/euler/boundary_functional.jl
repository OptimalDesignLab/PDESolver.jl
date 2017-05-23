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
                                         functional_edges)

  Cl = zero(Tsol)
  Cd = zero(Tsol)

  boundaryinterpolate!(mesh.sbpface, mesh.bndryfaces, eqn.q, eqn.q_bndry)

  for itr = 1:length(functional_edges)
    g_edge_number = functional_edges[itr] # Extract geometric edge number
    start_index = mesh.bndry_offsets[g_edge_number]
    end_index = mesh.bndry_offsets[g_edge_number+1]
    idx_range = start_index:(end_index-1)
    bndry_facenums = sview(mesh.bndryfaces, idx_range) # faces on geometric edge i

    nfaces = length(bndry_facenums)
    boundary_integrand = zeros(Tsol, 2, mesh.sbpface.numnodes, nfaces)

    for i = 1:nfaces
      global_facenum = idx_range[i]

      face_integrand = slice(boundary_integrand, :, :, i)
      if eqn.params.isViscous
        calcForceCoef_viscous(mesh, sbp, eqn, opts, global_facenum, face_integrand)
      else
        calcForceCoef_inviscid(mesh, sbp, eqn, opts, global_facenum, face_integrand)
      end
    end    # End for i = 1:nfaces

    val_per_geom_edge = zeros(Tsol, 2)

    integratefunctional!(mesh.sbpface, mesh.bndryfaces[idx_range], 
                         boundary_integrand, val_per_geom_edge)

    Cl += val_per_geom_edge[1]
    Cd += val_per_geom_edge[2]

  end  # End for itr = 1:length(functional_edges)

  # scale because of nondimentionalization

  return [Cl, Cd]
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

function calcForceCoef_inviscid{Tsol, Tmsh}(mesh::AbstractDGMesh{Tmsh},
                                            sbp::AbstractSBP,
                                            eqn::EulerData{Tsol},
                                            opts,
                                            bndry_indx,
                                            force::AbstractArray{Tsol, 2})
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
  Cx = zeros(Tsol, mesh.numNodesPerFace)
  Cy = zeros(Tsol, mesh.numNodesPerFace)
  for j = 1 : mesh.sbpface.numnodes
    q = sview(eqn.q_bndry, :, j, bndry_indx)
    convertToConservative(eqn.params, q, q2)

    aux_vars = sview(eqn.aux_vars_bndry, :, j, bndry_indx)
    dxidx = sview(mesh.dxidx_bndry, :, :, j, bndry_indx)

    flux_func = noPenetrationBC()
    flux_func(q, aux_vars, x, dxidx, nrm_xi, Finv, eqn.params)
    Cx[j] = Finv[2]
    Cy[j] = Finv[3]
  end
  # Ma = eqn.params.Ma
  # face_integrand /= Ma*Ma
  aoa = eqn.params.aoa
  force[1, :] = -Cx[:] * sin(aoa) + Cy[:] * cos(aoa)
  force[2, :] =  Cx[:] * cos(aoa) + Cy[:] * sin(aoa)
  return nothing
end


function calcForceCoef_viscous{Tsol, Tmsh}(mesh::AbstractDGMesh{Tmsh},
                                  sbp::AbstractSBP,
                                  eqn::EulerData{Tsol},
                                  opts,
                                  bndry_indx,
                                  force::AbstractArray{Tsol, 2})

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

  Cx = Array(Tsol, mesh.numNodesPerFace)
  Cy = Array(Tsol, mesh.numNodesPerFace)
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
    Cx[j] = Finv[2]
    Cy[j] = Finv[3]
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

  # compute boundary value
  bndry_val_functor = AdiabaticWall()
  bndry_val_functor(q_face, nrm1, eqn.params, q_bnd)
  # using boundary value to compute diffusion tensor and viscousl flux
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
    nGtn1 = (nrm1[1,n]*nrm1[1,n]*Gt[2,:,1,1,n]
           + nrm1[1,n]*nrm1[2,n]*Gt[2,:,1,2,n]
           + nrm1[2,n]*nrm1[1,n]*Gt[2,:,2,1,n]
           + nrm1[2,n]*nrm1[2,n]*Gt[2,:,2,2,n])
    nGtn2 = (nrm1[1,n]*nrm1[1,n]*Gt[3,:,1,1,n]
           + nrm1[1,n]*nrm1[2,n]*Gt[3,:,1,2,n]
           + nrm1[2,n]*nrm1[1,n]*Gt[3,:,2,1,n]
           + nrm1[2,n]*nrm1[2,n]*Gt[3,:,2,2,n])
    nGtn1_dq = nGtn1[1]*dq[1] + nGtn1[2]*dq[2] + nGtn1[3]*dq[3] + nGtn1[4]*dq[4]
    nGtn2_dq = nGtn2[1]*dq[1] + nGtn2[2]*dq[2] + nGtn2[3]*dq[3] + nGtn2[4]*dq[4]
    integrand = const_tii /he * nGtn1_dq
    integrand -= Fvis[1,2,n] * nrm1[1,n] + Fvis[2,2,n] * nrm1[2,n]
    Cx[n] += area[n] * coef_nondim * integrand
    integrand = const_tii /he * nGtn2_dq
    integrand -= Fvis[1,3,n] * nrm1[1,n] + Fvis[2,3,n] * nrm1[2,n]
    Cy[n] += area[n] * coef_nondim * integrand
  end

  aoa = eqn.params.aoa
  force[1, :] = -Cx[:] * sin(aoa) + Cy[:] * cos(aoa)
  force[2, :] =  Cx[:] * cos(aoa) + Cy[:] * sin(aoa)

  # Ma = eqn.params.Ma
  # face_integrand /= Ma*Ma
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
# "drag_inviscid" => drag_inviscid,
# "drag_viscous" => drag_viscous,
# "lift_inviscid" => lift_inviscid,
# "lift_viscous" => lift_viscous,
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

