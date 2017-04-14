
type SRC0 <: SRCType
end
function call(obj::SRC0, 
              xy::AbstractVector, 
              f::AbstractArray)
  f[:] = 0.0
  return nothing
end

type SRC1 <: SRCType
end
function call(obj::SRC1, 
              xy::AbstractVector, 
              f::AbstractArray)
  f[:] = 1.0
  return nothing
end

#
# a scalar diffusion coefficient, 10
#
type SrcExpTrigPoly0thDiffn <: SRCType
end
function call(obj::SrcExpTrigPoly0thDiffn, 
              xy::AbstractVector, 
              src::AbstractArray)
  k = 2.0
  ss = sin(2*k*pi*xy[1])*sin(2*k*pi*xy[2])
  cs = cos(2*k*pi*xy[1])*sin(2*k*pi*xy[2])
  sc = sin(2*k*pi*xy[1])*cos(2*k*pi*xy[2])
  ex = exp(xy[1] + xy[2])
  src[:] = (2. - 8*k*k*pi*pi) * ss + 4*k*pi*(cs + sc)
  src[:] *= -10.*ex
  return nothing
end
#
# a scalar diffusion coefficient, 10
#
type SrcTrigPoly0thDiffn <: SRCType
end
function call(obj::SrcTrigPoly0thDiffn, 
              xy::AbstractVector, 
              src::AbstractArray)
  n = 2.0
  factor = 8.0*n*n*pi*pi
  ss = sin(2*n*pi*xy[1])*sin(2*n*pi*xy[2])
  src[:] = 0.0
  for i = 1:length(src)
    src[i] = factor*ss * 1.0e1
  end
  return nothing
end
#
# a 2nd order diffusion coefficient tensor
#
type SrcTrigPoly2ndDiffn <: SRCType
end
function call(obj::SrcTrigPoly2ndDiffn, 
              xy::AbstractVector, 
              src::AbstractArray)
  n = 2.0
  sc = sin(2*n*pi*xy[1])*cos(2*n*pi*xy[2])
  cs = cos(2*n*pi*xy[1])*sin(2*n*pi*xy[2])
  ss = sin(2*n*pi*xy[1])*sin(2*n*pi*xy[2])
  cc = cos(2*n*pi*xy[1])*cos(2*n*pi*xy[2])
  src[:] = 0.0
  for i = 1:length(src)
    src[i] += 6*n*pi*xy[1]*cs
    src[i] += 6*n*pi*xy[2]*sc
    src[i] += 4*n*n*pi*pi*(-2.0 - xy[1]*xy[1] - xy[2]*xy[2])*ss
    src[i] += 8*n*n*pi*pi*xy[1]*xy[2]*cc
    src[i] = -src[i]
  end
  return nothing
end



#
# a 6th order diffusion coefficient tensor
#
type SrcTrigPoly6thDiffn <: SRCType
end
function call(obj::SrcTrigPoly6thDiffn, 
              xy::AbstractVector, 
              src::AbstractArray)
  n = 2.0
  sc = sin(2*n*pi*xy[1])*cos(2*n*pi*xy[2])
  cs = cos(2*n*pi*xy[1])*sin(2*n*pi*xy[2])
  ss = sin(2*n*pi*xy[1])*sin(2*n*pi*xy[2])
  cc = cos(2*n*pi*xy[1])*cos(2*n*pi*xy[2])
  src[:] = 0.0
  for i = 1:length(src)
    src[i] += 2*n*pi*(6*xy[1]^5 + 3*xy[1]^3*xy[2]^2)*cs
    src[i] += 2*n*pi*(6*xy[2]^5 + 3*xy[1]^2*xy[2]^3)*sc
    src[i] += 4*n*n*pi*pi*(-2.0 - xy[1]^6 - xy[2]^6)*ss
    src[i] += 8*n*n*pi*pi*xy[1]^3*xy[2]^3*cc
    src[i] = -src[i]
  end
  return nothing
end

global const SRCDict = Dict{ASCIIString, SRCType}(
 "SRC0" => SRC0(),
 "SRC1" => SRC1(),
 "SrcTrigPoly0thDiffn" => SrcTrigPoly0thDiffn(),
 "SrcExpTrigPoly0thDiffn" => SrcExpTrigPoly0thDiffn(),
 "SrcTrigPoly2ndDiffn" => SrcTrigPoly2ndDiffn(),
 "SrcTrigPoly6thDiffn" => SrcTrigPoly6thDiffn(),
)


function getSrcFuntors(mesh::AbstractMesh,
                       sbp::AbstractSBP,
                       eqn::EllipticData,
                       opts)
  eqn.src_func = SRCDict[opts["SRC_name"]]
  return nothing
end

function calcSource{Tmsh, Tsol, Tres, Tdim}(mesh::AbstractMesh{Tmsh},
                                            sbp::AbstractSBP,
                                            eqn::EllipticData{Tsol, Tres, Tdim},
                                            srcFunc::SRCType,
                                            src::AbstractArray{Tres, 3})

  numElems = mesh.numEl
  numNodesPerElement = mesh.numNodesPerElement
  # xy = Array(Tmsh, Tdim)
  src_tmp = Array(eltype(src), size(src, 1))
  for elem = 1:numElems
    # @bp
    for n = 1:numNodesPerElement
      xy = sview(mesh.coords, :, n, elem)

      srcFunc(xy, sview(src, :, n, elem))
    end
  end
end
