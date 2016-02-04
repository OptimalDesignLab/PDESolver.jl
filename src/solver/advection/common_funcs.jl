# common_funcs.jl

function calc_x5plusy5{Tmsh}(coords::AbstractArray{Tmsh}, alpha_x, alpha_y, t)

  x = coords[1]
  y = coords[2]
  u = x^5 + y^5
  
  return u
end # end function calc_x5plusy5

function calc_exp_xplusy{Tmsh}(coords::AbstractArray{Tmsh}, alpha_x, alpha_y, t)

  x = coords[1]
  y = coords[2]
  u = exp(x+y)

  return u
end

function calc_sinwave{Tmsh}(coords::AbstractArray{Tmsh}, alpha_x, alpha_y, t)
 
  x = coords[1]
  y = coords[2]

  return sin (-x + t)
#  return sin( (-x + -y)/2 + t)
end


function calc_sinwavey{Tmsh}(coords::AbstractArray{Tmsh}, alpha_x, alpha_y, t)
  y = coords[2]
  return sin(y)^2 + 5*sin(y) + 3/sin(y)
end

function calc_sinwavey_pert{Tmsh}(coords::AbstractArray{Tmsh}, alpha_x, alpha_y, t)
  x = coords[1]
  return calc_sinwavey(coords, t)*1000*sin(x)
end

function calc_mms1{Tmsh}(coords::AbstractArray{Tmsh}, alpha_x, alpha_y, t)
  x = coords[1]
  y = coords[2]
  px = 1.5
  py = 1.5
  return sin(px*x)*cos(py*y) + sin(px*x)*cos(py*y)
end

# x derivative of mms1
function calc_mms1dx{Tmsh}(coords::AbstractArray{Tmsh}, alpha_x, alpha_y, t)
  x = coords[1]
  y = coords[2]
  px = 1.5
  py = 1.5
  return px*cos(px*x)*cos(py*y) + px*cos(px*x)*cos(py*y)
end

