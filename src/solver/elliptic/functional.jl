abstract AbstractFunctional

type volumeAverage <: AbstractFunctional
end
function (obj::volumeAverage){Tmsh, Tsol, Tres, Tdim}(
                                      mesh::AbstractMesh{Tmsh},
                                      sbp::AbstractSBP,
                                      eqn::EllipticData{Tsol, Tres, Tdim},
                                      opts,
                                      val::Array{Tsol, 1})
  @assert(length(val) == mesh.numDofPerNode)
  val[:] = 0.0

  w = sview(sbp.w, :)
  for elem = 1 : mesh.numEl
    jac = sview(mesh.jac, :, elem)
    q = sview(eqn.q, :, :, elem)
    for n = 1 : mesh.numNodesPerElement
      val[:] += w[n]/jac[n]*q[:, n]
    end
  end

  return nothing
end


type volumeEnergy <: AbstractFunctional
end
function (obj::volumeEnergy){Tmsh, Tsol, Tres, Tdim}(
                                      mesh::AbstractMesh{Tmsh},
                                      sbp::AbstractSBP,
                                      eqn::EllipticData{Tsol, Tres, Tdim},
                                      opts,
                                      val::Array{Tsol, 1})
  @assert(length(val) == mesh.numDofPerNode)
  val[:] = 0.0

  w = sview(sbp.w, :)
  for elem = 1 : mesh.numEl
    jac = sview(mesh.jac, :, elem)
    q = sview(eqn.q, :, :, elem)
    for n = 1 : mesh.numNodesPerElement
      for dof = 1:mesh.numDofPerNode
        val[dof] += w[n]/jac[n]*q[dof, n]*q[dof,n]
      end
    end
  end

  return nothing
end
#

#
# IP Modification of target functional
#
function FunctionalModification{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh},
                                                        sbp::AbstractSBP,
                                                        eqn::EllipticData{Tsol, Tres, Tdim},
                                                        otps::Array{Tsol, 1})
  if haskey(opts, "FunctionalModification") && opts["FunctionalModification"] == "No"	
    return nothing
  end

  for iBC = 1:mesh.numBC

    bc_func = mesh.bndry_funcs[iBC]
    if isNeumann(bc_func)
      continue
    end

    indx0 = mesh.bndry_offsets[iBC]
    indx1 = mesh.bndry_offsets[iBC+1] - 1

    for f = indx0:indx1

    end

  end
  return nothing
end

global const FunctionalDict = Dict{String, AbstractFunctional}( 
  "volumeAverage" => volumeAverage(),
  "energy" => volumeEnergy()
)

function getFunctional(mesh::AbstractMesh,
                       sbp::AbstractSBP,
                       eqn::EllipticData,
                       opts)
  eqn.functional = FunctionalDict[opts["Functional"]]
end
