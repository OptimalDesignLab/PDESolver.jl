# an array that returns the indices from getindex

import Base: size, getindex, setindex!

"""
  Special array type such that `A[i] = i` in 1D and `A[i, j] = i` in 2D.
  This array can be stack-allocated.  It supportes the standard Array
  constructors.
"""
struct IdentityArray{T, N} <: DenseArray{T, N}
  dims::NTuple{N, Int}
end

size(A::IdentityArray) = A.dims

function setindex!(A::IdentityArray, v, i::Int)
  throw(ErrorException("setindex! not defind for IdentityArray"))
end

# construct from individual dimensions
# These constructors can be called as IdentityArray{T}(i, j) (which is not
# documented in Julia).
function IdentityArray{T}(i::Integer) where {T}
  return IdentityArray{T, 1}( (Int(i), ))
end

function IdentityArray{T}(i::Integer, j::Integer) where {T}
  return IdentityArray{T, 2}( (Int(i), Int(j)) )
end

function IdentityArray{T}(i::Integer, j::Integer, k::Integer) where {T}
  return IdentityArray{T, 3}( (Int(i), Int(j), Int(k)) )
end

function IdentityArray{T}(i::Integer, j::Integer, k::Integer, l::Integer
                         ) where {T}
  return IdentityArray{T, 4}( (Int(i), Int(j), Int(k), Int(l)) )
end

function IdentityArray{T}(i::Integer, j::Integer, k::Integer, l::Integer,
                          m::Integer) where {T}
  return IdentityArray{T, 5}( (Int(i), Int(j), Int(k), Int(l), Int(m)) )
end

function IdentityArray{T}(i::Integer, j::Integer, k::Integer, l::Integer,
                          m::Integer, n::Integer) where {T}
  return IdentityArray{T, 6}( (Int(i), Int(j), Int(k), Int(l), Int(m), Int(n)) )
end

function IdentityArray{T}(i::Integer, j::Integer, k::Integer, l::Integer,
                          m::Integer, n::Integer, o::Integer) where {T}
  return IdentityArray{T, 7}( (Int(i), Int(j), Int(k), Int(l), Int(m), Int(n), Int(o)) )
end

# NTuple constructor
function IdentityArray{T}(d::NTuple{T, Int})
  return IdentityArray{T, N}(d)
end


#------------------------------------------------------------------------------
# Indexing

@inline function getindex(A::IdentityArray{T, 1}, i::Int) where {T}
  return T(i)
end

@inline function getindex(A::IdentityArray{T, 2}, i::Integer, j::Integer) where {T}
  return T(i)  # for 2D arrays, make all columns the same
end

# it looks like Julia automatically implements linar indexing using the 2D method

# not doing the higher dimension yet
