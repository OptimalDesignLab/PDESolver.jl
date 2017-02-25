# This script holds differentiation by the mesh metrics. Write functions in here.

function getdFdm{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh},
                 sbp::AbstractSBP, eqn::EulerData{Tsol, Tres, Tdim}, opts)

  dFluxdm = zeros(Tsol, mesh.numDofPerNode, Tdim, Tdim, sbp.numnodes, mesh.numEl)
  nrm = zeros(Tmsh, Tdim)
  for i = 1:mesh.numEl
    for j = 1:mesh.numNodesPerElement
      q_vals = sview(eqn.q, :, j, i)
      aux_vars = sview(eqn.aux_vars, :, j, i)
      for k=1:Tdim  # loop over dimensions
        for p=1:Tdim
          nrm[p] = mesh.dxidx[k, p, j, i]
        end
        dfdm = sview(dFluxdm, :, :, k, j, i)
        calcdFluxdm(eqn.params, q_vals, aux_vars, nrm, dfdm)
      end
    end
  end

  return dFluxdm
end

function calcdFluxdm{Tmsh, Tsol, Tres}(params::ParamType{2, :conservative},
                      q::AbstractArray{Tsol,1},
                      aux_vars::AbstractArray{Tres, 1},
                      dir::AbstractArray{Tmsh},  dF::AbstractArray{Tsol,2})


  # Function for computing the derivative of flux w.r.tthe mesh metrics

  press = calcPressure(params, q) # Calculate pressure

  # dF = zeros(Tsol, mesh.dim+2, mesh.dim)

  fac = 1/q[1]
  u = q[2]*fac
  v = q[3]*fac

  dF[1,1] = q[1]*u
  dF[2,1] = q[2]*u + press
  dF[3,1] = q[3]*u
  dF[4,1] = (q[4] + press)*u

  dF[1,2] = q[1]*v
  dF[2,2] = q[2]*v
  dF[3,2] = q[3]*v + press
  dF[4,2] = (q[4] + press)*v

  return nothing
end

function calcdFluxdm{Tmsh, Tsol, Tres}(params::ParamType{3, :conservative},
                      q::AbstractArray{Tsol,1},
                      aux_vars::AbstractArray{Tres, 1},
                      dir::AbstractArray{Tmsh},  dF::AbstractArray{Tsol,2})

  # Function for computing the derivative of flux w.r.tthe mesh metrics

  press = calcPressure(params, q) # Calculate pressure

# dF = zeros(Tsol, mesh.dim+2, mesh.dim)

  fac = 1/q[1]
  u = q[2]*fac
  v = q[3]*fac
  w = q[4]*fac

  dF[1,1] = q[1]*u
  dF[2,1] = q[2]*u + press
  dF[3,1] = q[3]*u
  dF[4,1] = q[4]*u
  dF[5,1] = (q[5] + press)*u

  dF[1,2] = q[1]*v
  dF[2,2] = q[2]*v
  dF[3,2] = q[3]*v + press
  dF[4,2] = q[4]*v
  dF[5,2] = (q[4] + press)*v

  dF[1,3] = q[1]*w
  dF[2,3] = q[2]*w
  dF[3,3] = q[3]*w
  dF[4,3] = q[4]*w + press
  dF[5,3] = (q[5] + press)*w

  return nothing
end

@doc """
getBndryFluxdm

Get the derivative of the boundary flux w.r.t the mesh metrix dxidx matrix

"""->

function getdBndryFluxdm{Tmsh, Tsol}(mesh::AbstractMesh{Tmsh}, sbp::AbstractSBP,
                         eqn::EulerData{Tsol}, opts, deriv_bndry_funcs)

  dBndryFluxdm = zeros(Tsol, mesh.numDofPerNode, mesh.dim*mesh.dim, mesh.numNodesPerFace,
                       mesh.numBoundaryFaces)

  for i=1:mesh.numBC
    functor_i = deriv_bndry_funcs[i]
    start_index = mesh.bndry_offsets[i]
    end_index = mesh.bndry_offsets[i+1]
    idx_range = start_index:end_index  # TODO: should this be start_index:(end_index - 1) ?
    bndry_facenums_i = sview(mesh.bndryfaces, start_index:(end_index - 1))
    bndryflux_i = sview(dBndryFluxdm, :, :, :, start_index:(end_index - 1))

    # call the function that calculates the flux for this boundary condition
    # passing the functor into another function avoid type instability
    calcdBndryFluxdm(mesh, sbp, eqn, functor_i, idx_range, bndry_facenums_i, bndryflux_i)
  end

  return dBndryFluxdm
end
#=
function calcdBndryFluxdm{Tmsh,  Tsol, Tres}( mesh::AbstractCGMesh{Tmsh},
                          sbp::AbstractSBP, eqn::EulerData{Tsol},
                          functor::BCType, idx_range::UnitRange,
                          bndry_facenums::AbstractArray{Boundary,1},
                          bndryflux::AbstractArray{Tres, 4})


  nfaces = length(bndry_facenums)
  q2 = zeros(Tsol, mesh.numDofPerNode)
  for i=1:nfaces  # loop over faces with this BC
    bndry_i = bndry_facenums[i]
    global_facenum = idx_range[i]
    for j = 1:mesh.numNodesPerFace

      # get components
      q = sview(eqn.q_bndry, :, j, global_facenum)
      # convert to conservative variables if needed
      convertToConservative(eqn.params, q, q2)
      aux_vars = sview(eqn.aux_vars_bndry, :, j, global_facenum)
      x = sview(mesh.coords_bndry, :, j, global_facenum)
      dxidx = sview(mesh.dxidx_bndry, :, :, j, global_facenum)
      nrm = sview(sbp.facenormal, :, bndry_i.face)
      bndryflux_i = sview(bndryflux, :, :, j, i)

      functor(q2, aux_vars, x, dxidx, nrm, bndryflux_i, eqn.params)
    end
  end

  return nothing
end
=#

function calcdBndryFluxdm{Tmsh,  Tsol, Tres}( mesh::AbstractDGMesh{Tmsh},
                            sbp::AbstractSBP, eqn::EulerData{Tsol},
                            functor::BCType, idx_range::UnitRange,
                            bndry_facenums::AbstractArray{Boundary,1},
                            bndryflux::AbstractArray{Tres, 4})
  # calculate the boundary flux for the boundary condition evaluated by the
  # functor
  println("shape of bndryflux = ", size(bndryflux))
  nfaces = length(bndry_facenums)
  q2 = zeros(Tsol, mesh.numDofPerNode)
  for i=1:nfaces  # loop over faces with this BC
    bndry_i = bndry_facenums[i]
    global_facenum = idx_range[i]
    for j = 1:mesh.numNodesPerFace

      # get components
      q = sview(eqn.q_bndry, :, j, global_facenum)
      # convert to conservative variables if needed
      convertToConservative(eqn.params, q, q2)
      aux_vars = sview(eqn.aux_vars_bndry, :, j, global_facenum)
      x = sview(mesh.coords_bndry, :, j, global_facenum)
      dxidx = sview(mesh.dxidx_bndry, :, :, j, global_facenum)
      nrm = sview(sbp.facenormal, :, bndry_i.face)
      bndryflux_i = sview(bndryflux, :, :, j, i)

      functor(q2, aux_vars, x, dxidx, nrm, bndryflux_i, eqn.params)
    end
  end

  return nothing
end

type disentropicVortexBC_dm <: BCType
end


function call{Tmsh, Tsol, Tres}(obj::disentropicVortexBC_dm, q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, x::AbstractArray{Tmsh,1},
              dxidx::AbstractArray{Tmsh,2}, nrm::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 2}, params::ParamType{2})



  return nothing
end

type dnoPenetrationBC_dm <: BCType
end

function call{Tmsh, Tsol, Tres}(obj::dnoPenetrationBC_dm, q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, x::AbstractArray{Tmsh,1},
              dxidx::AbstractArray{Tmsh,2}, nrm::AbstractArray{Tmsh,1},
              bndryflux::AbstractArray{Tres, 2}, params::ParamType{2})

  nx = zero(Tmsh)
  ny = zero(Tmsh)
  nx = dxidx[1,1]*nrm[1] + dxidx[2,1]*nrm[2]
  ny = dxidx[1,2]*nrm[1] + dxidx[2,2]*nrm[2]
  nx = nx/(sqrt(nx*nx + ny*ny)) # Make them unit vectors
  ny = ny/(sqrt(nx*nx + ny*ny))

  Unrm = nx*q[2] + ny*q[3] # Get the normal momentum

  fac = -1.0/((nx*nx + ny*ny)^1.5)
  dnormalmomnetum_dxidx = fac*nx*nrm[1]*q[2]
  dnormalmomnetum_detadx = fac*nx*nrm[2]*q[2]
  dnormalmomnetum_dxidy = fac*ny*nrm[1]*q[3]
  dnormalmomnetum_detady = fac*ny*nrm[2]*q[3]

  qg = params.qg
  for i=1:length(q)
    qg[i] = q[i]
  end

  # DbndryFlux/dxidx
  qg[2] -= nx*dnormalmomnetum_dxidx + nrm[1]*Unrm
  qg[3] -= ny*dnormalmomnetum_dxidx
  dnx2_dxidx = nrm[1]
  dny2_dxidx = zero(Tsol)

  v_vals = params.v_vals
  convertFromNaturalToWorkingVars(params, qg, v_vals)

  calcEulerFlux(params, v_vals, aux_vars, [dnx2_dxidx, dny2_dxidx], bndryflux[:,1])

  # DbndryFlux/dxidy

  return nothing
end



global const BCDeriv_dm_Dict = Dict{ASCIIString, BCType}(
"isentropicVortexBC" => disentropicVortexBC_dm(),
"noPenetrationBC" => dnoPenetrationBC_dm(),
#"FreeStreamBC" => dFreeStreamBC_dm(),
)

function getBCDerivFunctors(mesh::AbstractMesh, sbp::AbstractSBP, eqn::EulerData, opts)

  deriv_bndry_funcs = Array(BCType, mesh.numBC)

  for i=1:mesh.numBC
    key_i = string("BC", i, "_name")
    val = opts[key_i]
    println("BCDeriv_dm_Dict[val] = ", BCDeriv_dm_Dict[val])
    deriv_bndry_funcs[i] = BCDeriv_dm_Dict[val]
  end

  return deriv_bndry_funcs
end

#------------------------------------------------------------------------------
# Routines for checking against complex step
#------------------------------------------------------------------------------


function complex_calcBoundaryFluxdm{Tmsh,  Tsol, Tres}( mesh::AbstractDGMesh{Tmsh},
                          sbp::AbstractSBP, eqn::EulerData{Tsol},
                          functor::BCType, idx_range::UnitRange,
                          bndry_facenums::AbstractArray{Boundary,1},
                          bndryflux::AbstractArray{Tres, 3})
  # calculate the boundary flux for the boundary condition evaluated by the
  # functor

  pert = complex(0, 1e-20)

  nfaces = length(bndry_facenums)
  q2 = zeros(Tsol, mesh.numDofPerNode)
  dxidx = zeros(Complex{Float64}, mesh.dim, mesh.dim)
  for i=1:nfaces  # loop over faces with this BC
    bndry_i = bndry_facenums[i]
    global_facenum = idx_range[i]
    for j = 1:mesh.numNodesPerFace

      # get components
      q = sview(eqn.q_bndry, :, j, global_facenum)
      # convert to conservative variables if needed
      convertToConservative(eqn.params, q, q2)
      aux_vars = sview(eqn.aux_vars_bndry, :, j, global_facenum)
      x = sview(mesh.coords_bndry, :, j, global_facenum)
      dxidx[:,:] = mesh.dxidx_bndry[:, :, j, global_facenum]
      dxidx[1,1] += pert
      nrm = sview(sbp.facenormal, :, bndry_i.face)
      bndryflux_i = sview(bndryflux, :, j, i)

      functor(q2, aux_vars, x, dxidx, nrm, bndryflux_i, eqn.params)
      bndryflux_i[:] = imag(bndryflux_i[:])/imag(pert)
    end
  end

  return nothing
end
