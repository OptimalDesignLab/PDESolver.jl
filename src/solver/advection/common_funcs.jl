# common_funcs.jl

function calc_x5plusy5{Tmsh,Tsol}(coords::AbstractArray{Tmsh}, u::Tsol)

  x = coords[1]
  y = coords[2]
  u = x^5 + y^5
  
  return nothing
end # end function calc_x5plusy5

function calc_exp_xplusy{Tmsh,Tsol}(coords::AbstractArray{Tmsh}, u::Tsol)

  x = coords[1]
  y = coords[2]
  u = exp(x+y)

  return nothing
end