# macros for accessing the auxiliary variables
# note that each variable needs two macros, a get and a set

#=
macro getPressure(aux_vars, node, element)
  pressure_index = 1
  return :($aux_vars[$pressure_index, $node, $element])
end
=#

macro getPressure(aux_vars)
  pressure_index = 1
  return :($aux_vars[$pressure_index])
end

function getPressure(aux_vars)
  return aux_vars[1]
end

macro setPressure(aux_vars, node, element, val)
  pressure_index = 1
  println("aux_vars = ", aux_vars)
  println("typeof(aux_vars) = ", typeof(aux_vars))
  return :($aux_vars[$pressure_index, $node, $element] = $val)
end
