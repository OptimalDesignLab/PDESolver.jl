function func3{T}(entity_ptr, r_::Ptr{T}, h_::Ptr{T}, m_ptr, u_::Ptr{T})
# an anisotropic function
# populates h with the desired mesh size in all three dimension

#  println("entered func3")
  # load arrays from C pointers
  r = pointer_to_array(r_, (3,3))  # REMEMBER: r is transposed when it is passed back to C
#  println("r = \n", r)
  h = pointer_to_array(h_, 3)
#  println("h = \n", h)
  u = pointer_to_array(u_, 9)
#  println("u = \n", u)


  # get vertex coords
  coords = zeros(3,1)
  getVertCoords(entity_ptr, coords, 3, 1)
  x = coords[1]
  y = coords[2]
  z = coords[3]

  # use xyz frame
  r[1,1] = 1.0
  r[2,2] = 1.0
  r[3,3] = 1.0

#  println("in julia, r = ", r)

  # calculate derivative of rho here

  drho_dx = smoothHeavisideder(x)
  h[1] = abs(0.5/drho_dx)  # make mesh size proportional to 1/drho/dx (larger gradient -> smaller mesh)
     # 2 is an emperical things with weird units


  ubound = 0.50
  lbound = 0.001
  # upper and lower bounds
  if h[1] < lbound
    h[1] = lbound
  elseif h[1] > ubound
    h[1] = ubound
  end


#  h[1] = 0.5
  h[2] = h[1]
  h[3] = 2.0


  println("x = ", x, " h = ", h, " drho_dx = ", drho_dx)
return nothing
end

# smooth heaviside function
function smoothHeavisideder(x)
# calculate the value of the smooth heaviside function at a location x
# x0 is specified within this function

  x0 = 0
  L = 5
  k = 5

#  return (L/(1 + e^(-k*(x-x0))))
  return L*(2*k*e^(-2*k*x))/(e^(-2*k*x) +1 )^2
end

