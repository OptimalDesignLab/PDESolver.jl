# Eval analytic derivative
function calcResidualDerivativeFWD(mesh::AbstractMesh, sbp::AbstractSBP, eqn::EulerData, opts)

  getDerivBCFunctors(mesh, sbp, eqn, opts)

  deriv_res = zeros(eqn.res)
  bndryflux_deriv = zeros(eqn.bndryflux)
  res_vec = zeros(eqn.res_vec)

  # Evaluate volume derivative

  # Evaluate boundary derivative
  EulerEquationMod.

  return nothing
end

function calc_dResdAlpha(mesh, sbp, eqn, opts, bndryflux_deriv)

  # Volume residual is zero. boundary faces with noPenetrationBC will also have
  # No contribution. The only contribution will come from faces with 
  # FreeStreamBC

  # Identify which faces have FreeStreamBC
  FSGeomFaces = Int[]
  for i = 1:mesh.numBC
  	key_i = string("BC", i, "_name")
  	if opts[key_i] == "FreeStreamBC"
  	  face_key = string("BC",i)
  	  BCedges_i = opts[face_key]
  	  FSGeomFaces = append!(FSGeomFaces, BCedges_i)
  	end
  end

  # Loop over those faces and fill the vector with residual
  for itr = 1:length(FSGeomFaces)
  	g_face_number = FSGeomFaces[itr]
  	# get the boundary array associated with the geometric edge
    itr2 = 0
    for itr2 = 1:mesh.numBC
      if findfirst(mesh.bndry_geo_nums[itr2], g_face_number) > 0
        break
      end
    end

    start_index = mesh.bndry_offsets[itr2]
    end_index = mesh.bndry_offsets[itr2+1]
    idx_range = start_index:(end_index-1)
    bndry_facenums = sview(mesh.bndryfaces, idx_range) # faces on geometric edge i

    nfaces = length(bndry_facenums)
    boundary_integrand = zeros(Tsol, 1, mesh.sbpface.numnodes, nfaces)
    q2 = zeros(Tsol, mesh.numDofPerNode)

    for i = 1:nfaces
      bndry_i = bndry_facenums[i]
      global_facenum = idx_range[i]
      for j = 1:mesh.sbpface.numnodes
        # println("element number = ", bndry_i.element-1, " face = ", bndry_i.face)
        q = sview(eqn.q_bndry, :, j, global_facenum)
        convertToConservative(eqn.params, q, q2)
        aux_vars = sview(eqn.aux_vars_bndry, :, j, global_facenum)
        x = sview(mesh.coords_bndry, :, j, global_facenum)
        dxidx = sview(mesh.dxidx_bndry, :, :, j, global_facenum)
        nrm = sview(sbp.facenormal, :, bndry_i.face)
        bndryflux_deriv_i = sview(bndryflux_deriv,:,j,i)
        FreeStreamBC_dAlpha(q2, aux_vars, x, dxidx, nrm, bndryflux_deriv_i, eqn.params)

      end  # End for j = 1:mesh.sbpface.numnodes
    end    # End for i = 1:nfaces
  end # End itr

  return nothing
end

@doc """
### EulerEquationMod.FreeStreamBC_dAlpha <: BCTypes

  This functor uses the Roe solver to calculate the flux for a boundary
  state corresponding to the free stream velocity, using rho_free, Ma, aoa, and E_free

  This is a low level functor
"""->

function FreeStreamBC_dAlpha(q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1},  x::AbstractArray{Tmsh,1},
              dxidx::AbstractArray{Tmsh,2}, nrm::AbstractArray{Tmsh,1}, 
              bndryflux::AbstractArray{Tres, 1}, params::ParamType{2})

  qg = params.qg

  calcFreeStream_dAlpha(x, params, qg)
  RoeSolver(params, q, qg, aux_vars, dxidx, nrm, bndryflux)
  
  return nothing
end

#=
@doc """
### EulerEquationMod.allZerosBC <: BCTypes

  This functor uses the Roe solver to calculate the flux for a boundary
  state where all the conservative variables have a value 0.0

  This is a low level functor
"""->

type allZerosBC <: BCType
end

function call{Tmsh, Tsol, Tres}(obj::allZerosBC, q::AbstractArray{Tsol,1},
              aux_vars::AbstractArray{Tres, 1}, x::AbstractArray{Tmsh,1},
              dxidx::AbstractArray{Tmsh,2}, nrm::AbstractArray{Tmsh,1}, 
              bndryflux::AbstractArray{Tres, 1}, params::ParamType{2})

  qg = zeros(Tsol, 4)
  calcZeros(x, params, qg)

  RoeSolver(params, q, qg, aux_vars, dxidx, nrm, bndryflux)

  # println("bndryflux = ", bndryflux)
  return nothing
end # end function call

global const DerivBCDict = Dict{ASCIIString, BCType}(
"FreeStreamBC_dAlpha" => FreeStreamBC_dAlpha(),
"allZerosBC" => allZerosBC(),
)

@doc """
### EulerEquationMod.getDerivBCFunctors

  This function uses the opts dictionary to populate mesh.bndry_funcs with
  the the functors

  This is a high level function.
"""->
# use this function to populate access the needed values in BCDict
function getDerivBCFunctors(mesh::AbstractMesh, sbp::AbstractSBP, eqn::EulerData, opts)
  # populate the array mesh.bndry_funcs with the functors for the boundary 
  # condition types

  for i=1:mesh.numBC
    key_i = string("DerivBC", i, "_name")
    val = opts[key_i]
    println("DerivBCDict[val] = ", DerivBCDict[val])
    mesh.bndry_funcs[i] = BCDict[val]
  end

  return nothing

end

=#