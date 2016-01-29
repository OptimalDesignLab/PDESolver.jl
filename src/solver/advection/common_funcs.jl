# common_funcs.jl

function calc_x5plusy5{Tmsh}(coords::AbstractArray{Tmsh})

  x = coords[1]
  y = coords[2]
  u = x^5 + y^5
  
  return u
end # end function calc_x5plusy5

function calc_exp_xplusy{Tmsh}(coords::AbstractArray{Tmsh})

  x = coords[1]
  y = coords[2]
  u = exp(x+y)

  return u
end

function calc_sinwave{Tmsh}(coords::AbstractArray{Tmsh}, t)
 
  x = coords[1]
  y = coords[2]

  return sin (-x + t)
#  return sin( (-x + -y)/2 + t)
end



