
@doc """
### SimpleODEMod.calc_x2_t2

  Calculates and returns x^2 and t^2
"""->
function calc_x2_t2{Tmsh}(coords::AbstractArray{Tmsh}, params::ParamType2, t)

  x = coords[1]
  y = coords[2]

  u = x^2 + t^2

  return u
end

