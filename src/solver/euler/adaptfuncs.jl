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



function shockRefine{T}(entity_ptr, r_::Ptr{T}, h_::Ptr{T}, m_ptr, f_ptr)
# an anisotropic function
# populates h with the desired mesh size in all three dimension
# f_ptr is a pointer to a solution field (apf::Field)

#  println("entered func3")
  # load arrays from C pointers
  r = pointer_to_array(r_, (3,3))  # REMEMBER: r is transposed when it is passed back to C
#  println("r = \n", r)
  h = pointer_to_array(h_, 3)
#  println("h = \n", h)


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

  u_node = zeros(4)
  retrieveNodeSolution(f_ptr, entity_ptr, u_node)
  # calculate derivative of rho here

  drho_dx = u_node[1]  # rho  value
#  drho_dx = smoothHeavisideder(x)
  
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

function shockRefine2{T}(entity_ptr, r_::Ptr{T}, h_::Ptr{T}, m_ptr, f_ptr, operator_ptr)
# an anisotropic function
# populates h with the desired mesh size in all three dimension
# entity_ptr is a pointer to a vertex
# r_ptr is a pointer to an array of doubles that specifiy the orthonormal
# coordinate system to use
# h_ptr is the desired mesh size in each direction of r
# m_ptr is the mesh pointer
# f_ptr is a pointer to a solution field (apf::Field)
# operator_ptr = pointer to SBP operator

#  println("entered func3")

  # load arrays from C pointers
  r = pointer_to_array(r_, (3,3))  # REMEMBER: r is transposed when it is passed back to C
#  println("r = \n", r)

  h = pointer_to_array(h_, 3)
#  println("h = \n", h)

#  println("operator_ptr = ", operator_ptr)
  sbp = unsafe_pointer_to_objref(operator_ptr)

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

  # get the elements this vertex is part of
  num_elements = countAdjacent(m_ptr, entity_ptr, 2)  # get adjacnet elements
  elements = getAdjacent(num_elements)

  # choose first element (arbitrary)
  el_ptr = elements[1]
  verts, tmp = getDownward(m_ptr, el_ptr, 0)

  # get coordinates
  coords = zeros(3,3, 1)
  sub_coords = sub(coords, :, :, 1)
  getFaceCoords(el_ptr, sub_coords, 3, 3)

  # get solution at vertices, find out the local index of the vertex
  u_vals = zeros(4,3, 1)
  vert_index = 0
  for j=1:3  # get solution value of each node, figure out which vertex is original
    subarray = sub(u_vals, :, j, 1)
    retrieveNodeSolution(f_ptr, verts[j], subarray)

    if (entity_ptr == verts[j])
      vert_index = j
    end
  end

  # get mapping jacobian
  dxi_dx = zeros(2,2,3,1)
  jac = zeros(3,1)
  mappingjacobian!(sbp, coords[1:2,:,:], dxi_dx, jac)
  

  # differentiate (mapping jacobian is constant for an element)
  u_dxidx = u_vals*dxi_dx[1,1,1,1]
  u_detadx = u_vals*dxi_dx[2,1,1,1]
  res = zeros(size(u_vals))

  differentiate!(sbp, 1, u_dxidx, res)
  differentiate!(sbp, 2, u_detadx, res)

  drho_dx = res[1, vert_index, 1]*jac[vert_index,1]
#  println("drho_dx = ", drho_dx)


  # make mesh size proportional to 1/ drhodx
  h[1] = abs(0.25/drho_dx) 


  ubound = 0.50
  lbound = 0.0001
  # upper and lower bounds
  if h[1] < lbound
    h[1] = lbound
  elseif h[1] > ubound
    h[1] = ubound
  end


#  h[1] = 0.5
  h[2] = h[1]  # make mesh square
  h[3] = 2.0


#  println("x = ", x, " h = ", h, " drho_dx = ", drho_dx)
#  print("\n")
return nothing
end


