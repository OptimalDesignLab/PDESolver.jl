# functions for doing complex step

import Base: flipsign

"""
  Extends Base.flipsign with method for complex numbers
"""
function flipsign(a::T, b::Number) where {T <: Complex}

  xreal = flipsign(real(a), b)
  ximag = flipsign(imag(a), b)

  return T(xreal, ximag)
end

@doc """
### EulerEquationMod.absvalue

  This function computes the absolute value in a way such that the complex 
  step method works.

  There is probably a bit-twiddling way of doing this faster

"""
function absvalue(x::T) where T <: Complex

#  c_part = flipsign(imag(x), real(x))
#  return T(abs(real(x)), c_part)
  return flipsign(x, real(x))
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


function absvalue_rev(x::Number, x_bar::Number)

#  if real(x) < 0
#    x = -x
#    x_bar = -x_bar
#  end

  # this function has the odd property that forward and reverse modes are
  # equivalent
  x2 = flipsign(x, real(x))
  x2_bar = flipsign(x_bar, real(x))


  return x2, x2_bar
end


@doc """
###Utils.absvalue_deriv

Computes the derivative of the absolute value of a variable w.r.t itself.

**Inputs**

* `val` : The variable whose derivative needs to be computed e.r.t itself

"""->

function absvalue_deriv(val::Tval) where Tval

  Tval(sign_c(val))
  #=
  if val > zero(Tval)
    return one(Tval)
  elseif val < zero(Tval)
    return -one(Tval)
  else
    return zero(Tval)
  end
  =#

end # End function absvalue_deriv

#------------------------------------------------------------------------------
# replace abs with a cubic spline fit near 0, to make it differentiable

"""
  A replacement for absvalue() that is differentiable everywhere (C^1 continuous)

  Unfortunately, it is somewhat more expensive than absvalue()
"""
function absvalue2(val::Number)

  delta = 1e-13
  if absvalue(val) > delta
    return absvalue(val)
  elseif real(val) < 0
    tmp = val + 100e-15
    return @evalpoly tmp 100e-15 -1.0 -10e12 100e24
  else  # real(val) > 0
    tmp = val
    return @evalpoly tmp 0.0 0.0 20e12 -100e24
  end
end

"""
  Derivative of abs()
"""
function absvalue2_deriv(val::T) where {T <: Number}

  delta = 1e-13
  if absvalue(val) > delta
    return sign_c(val)
  elseif real(val) < 0
    tmp = val + 100e-15
     return @evalpoly tmp -1.0 -20e12 300e24
  else  # real(val) > 0
    tmp = val
    return @evalpoly tmp 0.0 40e12 -300e24
  end
end



function max_deriv_rev(x::Tval, y::Tval, max_val_bar::Tval) where Tval

  x_bar = zero(Tval)
  y_bar = zero(Tval)
  
  if real(x) > real(y)
    x_bar += max_val_bar
  else
    y_bar += max_val_bar
  end

  return x_bar, y_bar
end # End function max_rev


"""
  Smooth abs() function based on Harten's second entropy fix
"""
function absvalue3(val::Number)

  delta = -1
  val1 = absvalue(val)
  if val1 > delta
    return val1
  else
    term1 = val*val/delta
    term2 = delta
    return ((val*val)/delta + delta)/2
  end
end

function absvalue3_deriv(val::T) where {T <: Number}

  delta = -1

  val1 = absvalue(val)
  if val1 > delta
    return absvalue_deriv(val)
  else 
    return val/delta
  end
end




"""
  Complex-step safe version of the sign function
"""
function sign_c(a::Number)

  return sign(real(a))
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



"""
  Reverse-mode differentiated version of atan2.

  **Inputs**

   * y
   * x
   * theta_bar

  **Outputs**

   * theta
   * y_bar
   * x_bar

"""
function atan2_rev(y::Number, x::Number, theta_bar::Number)

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

  # reverse sweep
  x_bar = zero(x)
  y_bar = zero(y)
  x2_bar = zero(x)
  y2_bar = zero(y)

  if y < 0
    theta_bar = -theta_bar
  end

  if x < 0
    theta_bar = -theta_bar
  end

  if x2 > y2
    t1 = theta_bar/(1 + y2*y2/(x2*x2))
    x2_bar += t1*-y2/(x2*x2)
    y2_bar += t1/x2
  else
    t2 = -theta_bar/(1 + x2*x2/(y2*y2))
    x2_bar += t2/y2
    y2_bar += t2*-x2/(y2*y2)
  end

  # there is no += for functions that return multiple values
  t3, t3_bar = absvalue_rev(x, x2_bar)
  x_bar += t3_bar

  t4, t4_bar = absvalue_rev(y, y2_bar)
  y_bar += t4_bar

  return theta, y_bar, x_bar
end







