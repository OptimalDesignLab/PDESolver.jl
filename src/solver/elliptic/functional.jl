abstract type AbstractFunctional end

mutable struct volumeAverage <: AbstractFunctional
end
function (obj::volumeAverage)(
              mesh::AbstractMesh{Tmsh},
              sbp::AbstractOperator,
              eqn::EllipticData{Tsol, Tres, Tdim},
              opts,
              val::Array{Tsol, 1}) where {Tmsh, Tsol, Tres, Tdim}
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


mutable struct volumeEnergy <: AbstractFunctional
end
function (obj::volumeEnergy)(
              mesh::AbstractMesh{Tmsh},
              sbp::AbstractOperator,
              eqn::EllipticData{Tsol, Tres, Tdim},
              opts,
              val::Array{Tsol, 1}) where {Tmsh, Tsol, Tres, Tdim}
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
function FunctionalModification(mesh::AbstractMesh{Tmsh},
                                sbp::AbstractOperator,
                                eqn::EllipticData{Tsol, Tres, Tdim},
                                otps::Array{Tsol, 1}) where {Tmsh, Tsol, Tres, Tdim}
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
                       sbp::AbstractOperator,
                       eqn::EllipticData,
                       opts)
  eqn.functional = FunctionalDict[opts["Functional"]]
end
