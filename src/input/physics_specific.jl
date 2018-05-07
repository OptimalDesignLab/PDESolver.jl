# functions that allow physics to register their own options checking

"""
  Dictionary to hold the physics-specific options checking functions.

  Do not put functions into this dictionary directly, use the accessor
  function [`registerOptionsChecker`](@ref)
"""
global const PhysicsOptionsFuncs = Dict{ASCIIString, Function}()

"""
  API function for registering an options checking function for a given
  physics.

  Note that users generally do not have to call this function directly,
  it is handled by [`register_physics`](@ref)

  This function is run after the non-physics specific default options 
  have been supplied.  The function must have a a method with signature

  ```
    func(opts::Dict)
  ```

  which modifies the dictionary in-place.

  It is illegal to register more than one options checker per physics
  or a single options checker for more than one physics.

  **Inputs**

   * physics_name: the name defined by the physics module
   * func: the options checking function
"""
function registerOptionsChecker(physics_name::ASCIIString, func::Function)

  if haskey(PhysicsOptionsFuncs, physics_name)
    error("An options checker has already been registered for physics $(physics_name)")
  end

  for (key, val) in PhysicsOptionsFuncs
    if val == func
      error("Options checker for physics $(physics_name) is already registered for physics $(key)")
    end
  end

  if length(methods(func, (Dict,))) == 0
    error("The options checker for physics $(physics_name) does not have a method matching func(opts::Dict)")
  end

  # function appears acceptable, add it to the list
  PhysicsOptionsFuncs[physics_name] = func

  return nothing
end
