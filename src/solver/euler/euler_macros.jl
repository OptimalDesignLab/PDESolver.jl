# macros for accessing the auxiliary variables
# note that each variable needs two macros, a get and a set

#=
macro getPressure(aux_vars, node, element)
  pressure_index = 1
  return :($aux_vars[$pressure_index, $node, $element])
end
=#


@doc """
  These macros set and retrieve values of auxiliary variables.  Users
   should not acess values in the auxiliary variables directly.
"""


@doc """
### EulerEquationMod.@getPressure

  This macro returns pressure from the vector of auxiliary variables at 
  a node.  The user should use this function to access pressure rather 
  than accessing the vector of auxiliary variables directly
"""->
macro getPressure(vars)
  pressure_index = 1
  return :($vars[$pressure_index])
end


function getPressure(aux_vars)
  return aux_vars[1]
end



@doc """
### EulerEquationMod.@setPressure

  This function tests the pressure for the specified node and element,
  wheer aux_vars is the 3D array of auxiliary variables.  The users should
   use this function to set the value rather than accessing the array of 
   auxiliary variables directly
"""->
macro setPressure(aux_vars, node, element, val)
  pressure_index = 1
  println("aux_vars = ", aux_vars)
  println("typeof(aux_vars) = ", typeof(aux_vars))
  return :($aux_vars[$pressure_index, $node, $element] = $val)
end
