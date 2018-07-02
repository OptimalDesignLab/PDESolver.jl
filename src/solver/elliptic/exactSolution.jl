
abstract type ExactSolutionType end

export ExactSolutionType, ExactTrig, calcErrorL2Norm

mutable struct ExactTrig <: ExactSolutionType
end
function (obj::ExactTrig)(xy::AbstractArray{Tmsh}, q::AbstractArray{Tsol, 1}) where {Tmsh, Tsol}  # (numDofPerNode))
  k = 2.0
  q[:] = sin(2*k*pi*xy[1])*sin(2*k*pi*xy[2])
  return nothing
end

mutable struct ExactExpTrig <: ExactSolutionType
end
function (obj::ExactExpTrig)(xy::AbstractArray{Tmsh}, q::AbstractArray{Tsol, 1}) where {Tmsh, Tsol}  # (numDofPerNode))
  k = 2.0
  ex = exp(xy[1] + xy[2])
  q[:] = ex * sin(2*k*pi*xy[1])*sin(2*k*pi*xy[2])
  return nothing
end

mutable struct ExactPoly2nd <: ExactSolutionType
end
function (obj::ExactPoly2nd)(xy::AbstractArray{Tmsh}, q::AbstractArray{Tsol, 1}) where {Tmsh, Tsol}
  a = 1.0
  b = 1.0
  c = 1.0
  q[:] =  a*xy[1]*xy[1] + b*xy[1]*xy[2] + c*xy[2]*xy[2]
  return nothing
end


global const ExactDict = Dict{String, ExactSolutionType}(
 "ExactTrig" => ExactTrig(),
 "ExactPoly2nd" => ExactPoly2nd(),
 "ExactExpTrig" => ExactExpTrig()
)


function calcErrorL2Norm(mesh::AbstractMesh{Tmsh},
                         sbp::AbstractSBP,
                         eqn::AbstractEllipticData{Tsol, Tres},
                         opts) where {Tmsh, Tsol, Tres}
  l2norm::Float64 = 0.
  lInfnorm::Float64 = 0.
  qe = Array{Tsol}(mesh.numDofPerNode)
  exactFunc = ExactDict[opts["exactSolution"]]

  lInfnorm = -1.0
  elem_with_max_dq = 0
  for e = 1 : mesh.numEl
    for n = 1 : mesh.numNodesPerElement
      xy = sview(mesh.coords, :, n, e)
      exactFunc(xy, qe)
      q = sview(eqn.q, :, n, e)
      jac = mesh.jac[n, e]
      for v = 1:mesh.numDofPerNode
        # dq = Real(q[v] - qe[v])
        dq = Float64(q[v] - qe[v])
        l2norm += dq*dq*sbp.w[n]/jac
        if lInfnorm < absvalue(dq)
          lInfnorm = absvalue(dq)
        end
      end
    end
  end

  return sqrt(l2norm), lInfnorm
end
