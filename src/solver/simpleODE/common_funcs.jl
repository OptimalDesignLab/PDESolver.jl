
@doc """
### SimpleODEMod.calc_x2_4t3

  Calculates and returns x^2 + 4*t^3
"""->
function calc_x2_4t3{Tmsh}(coords::AbstractArray{Tmsh}, params::ParamType2, t)

  x = coords[1]
  y = coords[2]

  u = x^2 + 4*(t^3)

  return u
end


@doc """
### SimpleODEMod.calc_x2_3t2

  Calculates and returns x^2 + 3*t^2
"""->
function calc_x2_3t2{Tmsh}(coords::AbstractArray{Tmsh}, params::ParamType2, t)

  x = coords[1]
  y = coords[2]

  u = x^2 + 3*(t^2)

  return u
end

@doc """
### SimpleODEMod.calc_x2_t2

  Calculates and returns x^2 + t^2
"""->
function calc_x2_t2{Tmsh}(coords::AbstractArray{Tmsh}, params::ParamType2, t)

  x = coords[1]
  y = coords[2]

  u = x^2 + t^2

  return u
end

@doc """
### SimpleODEMod.calc_4t3

  Calculates and returns 4t^3
  (du/dt of calc_x2_t4)
"""->
function calc_4t3{Tmsh}(coords::AbstractArray{Tmsh}, params::ParamType2, t)

  x = coords[1]
  y = coords[2]

  u = 4*(t^3)

  return u
end

@doc """
### SimpleODEMod.calc_3t2

  Calculates and returns 3t^2
  (du/dt of calc_x2_t3)
"""->
function calc_3t2{Tmsh}(coords::AbstractArray{Tmsh}, params::ParamType2, t)

  x = coords[1]
  y = coords[2]

  u = 3*(t^2)

  return u
end

@doc """
### SimpleODEMod.calc_2t

  Calculates and returns 2t
  (du/dt of calc_x2_t2)
"""->
function calc_2t{Tmsh}(coords::AbstractArray{Tmsh}, params::ParamType2, t)

  x = coords[1]
  y = coords[2]

  u = 2*t

  return u
end

