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

function calc_mms1{Tmsh}(coords::AbstractArray{Tmsh}, t)
  x = coords[1]
  y = coords[2]
  return (sin(200*pi*x)*cos(200*pi*y) + sin(50*pi*x)*cos(400*pi*y))
end

# x derivative of mms1
function calc_mms1dx{Tmsh}(coords::AbstractArray{Tmsh}, t)
  x = coords[1]
  y = coords[2]
  return 200*pi*cos(200*pi*x)*cos(200*pi*y) + 50*pi*cos(50*pi*x)*cos(400*pi*y)
end

