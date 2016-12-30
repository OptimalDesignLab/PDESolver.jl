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
              retrieve_physics function. Typically the name is capitalized.
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

"""
  This function registers a new initial condition function with the specified
  physics module.  The function must have the signature:

   ICfunc(mesh::AbstractMesh, sbp::AbstractSBP, eqn:AbstractSolutionData{Tsol},
                opts:Dict, q_vec::AbstractVector{Tsol})

  where q_vec is a vector of length numDof that will be populated with the
  initial condition.  If the function is used to provide the exact solution
  for an unsteady problem (for calculating the error via the `exact_soln_func`
  key in the options dictionary), then it should use eqn.params.t as the
  current time value.

  This function does not attempt to verify that the functor has the correct
  signature, so the user should take care to ensure it is correct for the
  physics module.

  Inputs:
    mod: physics module to register the function with
    fname: name associated with the function, used as the value for the
           for any key in the options dictionary that specifies an initial
           condition function
    func: the function itself

  Outputs:
    none
"""
function registerIC(mod::Module, fname::ASCIIString, func::Function)

  # check if name is already registered
  if haskey(mod.ICDict, fname)
    throw(ErrorException("IC name $fname is already registered with physics module $mod"))
  end

  # check if function is already registered
  for (key, val) in mod.ICDict
    if val == func
      throw(ErrorException("IC function $func is already registered with physics module $mod"))
    end
  end


  # there is currently no good way of verify the signature. func.env.def.sig
  # is a tuple of argument types, but julia's type system is covarient type
  # system makes checking them difficult

  # if we get here, this name/function pair is unique
  mod.ICDict[fname] = func

  return nothing
end

"""
  This function registers a new boundary condition functor with the specified
  physics module.  The exact signature of the functor varies by physics module.
  The purpose of the functor is to calculate the flux across a boundary at
  a particular node. It usually takes in the solution values at that node, as
  well as the coordinates and mapping jacobian.  The params type is typically
  used to dispatch to different methods for 2D or 3D, if needed.  Generally,
  this functor will calculate the boundary state and then call a numerical
  flux function of some kind to compute the flux.

  This function does not attempt to verify that the functor has the correct
  signature, so the user should take care to ensure it is correct for the
  physics module.

  Inputs:
    mod: module to register the the functor with
    fname: the name associated with this function, used as the value for any
           key in the options dictionary that specifies a boundary condition,
           for example `BC1_name`
    func: the functor itself

  Outputs:
    none
"""
function registerBC(mod::Module, fname::ASCIIString, func::BCType)

  # check if name is already registered
  if haskey(mod.BCDict, fname)
    throw(ErrorException("BC name $fname is already registered with physics module $mod"))
  end

  # check if function is already registered
  for (key, val) in mod.BCDict
    if val == func
      throw(ErrorException("BC functor $func is already registered with physics module $mod"))
    end
  end

  mod.BCDict[fname] = func

  return nothing
end
