

function absvalue(x::Complex)
# redefine the absolute value function for complex step method
# I am pretty sure there is a bit twiddling way to doing this faster
  x2 = complex(real(x), imag(x))
  if (real(x) < 0)
    x2 *= -1
  end

  return x
end

function absvalue(x)   # general fallback method
  return abs(x)
end

import Base.isless

function isless(x::Complex, y::Complex)
  return real(x) < real(y)
end

# there is no isgreater, just a reuse of isless
