# mechanism for registering physics modules, initial conditions, and boundary
# conditions

"""
  This dictionary maps physics module names to the module and the startup 
  function for the module.  The fact that this dictionary exists is an
  implementation details, it should never be accessed directly.  Instead,
  the accessor functions below should be used.
"""
global const PhysicsModDict = Dict{ASCIIString, Tuple{Module, Function}}()

"""
  This function registered a new physics module with the global list of all
  known physics modules.  Every physics module should do this as part of
  module initialization.  The name, the module, and the startup function must
  be unique (ie. they must not already exist in the list).  This function
  throws and exception if they are not.
    '
  Inputs:
    modname:  an ASCIIString name for this entry in the list.  It is used
              to retrieve the module and startup function in the 
              retrieve_physics function
    mod:  the Module itself
    startup_func: the function for running the physics.  It must have signature
                  startup_func(fname::ASCIIString), where fname is the name of
                  an input file

  Outputs: none
"""
function register_physics(modname::ASCIIString, mod::Module, startup_func::Function)

  # check if this name is already used
  if haskey(PhysicsModDict, modname)
    mod_i, func_i = PhysicsModDict[modname]
    throw(ErrorException("Physics module with name $modname is already registered to module $mod_i"))
  end

  # check if the module or function is already registered
  for (key, val) in PhysicsModDict
    mod_i = val[1]
    func_i = val[2]

    if mod ==  mod_i
      throw(ErrorException("Physics module $mod is already registered to name $key"))
    end

    if startup_func == func_i
      throw(ErrorException("Startup function $startup_func is already registered to physics module $mod_i with name $key"))
    end
  end

  # if we get here, the registration is new
  PhysicsModDict[modname] = (mod, startup_func)

  return nothing
end  # end function

function retrieve_physics(modname::ASCIIString)

  if !haskey(PhysicsModDict, modname)
    # construct the error message
    err_msg = "unknown physics module name $modname, known names are:"
    for known_modname in keys(PhysicsModDict)
      err_msg *= "\n"*known_modname
    end

    throw(ErrorException(err_msg))
  end

  # if we get here, then this is a real physics module name

  return PhysicsModDict[modname]
end



