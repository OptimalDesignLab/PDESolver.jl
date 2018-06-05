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
