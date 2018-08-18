@doc """
### EulerEquationMod.absvalue

  This function computes the absolute value in a way such that the complex 
  step method works.

  There is probably a bit-twiddling way of doing this faster

"""
function absvalue(x::T) where T <: Complex

  c_part = flipsign(imag(x), real(x))
  return T(abs(real(x)), c_part)

end
#=
  x2 = complex(real(x), imag(x))

  if (real(x) < 0)
    x2 *= -1
  end

#  println("after modification, x2 = ", x2)
  return x2
end
=#
function absvalue(x)   # general fallback method
  return abs(x)
end

# array method
function absvalue(x::AbstractArray{T}) where T

  x2 = copy(x)
  for i=1:length(x)
    x2[i] = absvalue(x[i])
  end

  return x2
end

import Base.isless

@doc """
### EulerEquationMod.isless

  This function defines the isless() function (ie. <  operator) in such a way
  that algorithmic differentiation works

"""->
function isless(x::Complex, y::Complex)
  return real(x) < real(y)
end

function isless(x::Real, y::Complex)
  return x < real(y)
end

function isless(x::Complex, y::Real)
  return real(x) < y
end
# there is no isgreater, just a reuse of isless


#------------------------------------------------------------------------------
# trig functions
# most are defind in Base, but a few are needed

import Base: atan2

"""
  atan2 for complex number.  This function protects against numerical errors
  if x -> 0, but I don't know if it is as accurate as Base.atan2
"""
function atan2(y::Complex, x::Complex)

  x2 = absvalue(x)
  y2 = absvalue(y)
  # avoid numerical problems if x -> 0
  if x2 > y2
    theta = atan(y2/x2)
  else
    theta = pi/2 - atan(x2/y2)
  end

  if x < 0
    theta = pi - theta
  end

  if y < 0
    theta = -theta
  end

  return theta
end

# if only one complex number promote both to complex
function atan2(y::Number, x::T) where {T<:Complex}
  return atan2(T(y), x)
end

function atan2(y::T, x::Number) where {T<:Complex}
  return atan2(y, T(x))
end





