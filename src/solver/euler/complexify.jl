

function absvalue(x::Complex)
# redefine the absolute value function for complex step method
# I am pretty sure there is a bit twiddling way to doing this faster

#  println("entered absvalue")
#  println("x = ", x)

  x2 = complex(real(x), imag(x))
#  println("before modification, x2 = ", x2)

  if (real(x) < 0)
    x2 *= -1
  end

#  println("after modification, x2 = ", x2)
  return x2
end

function absvalue(x)   # general fallback method
  return abs(x)
end

import Base.isless

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
