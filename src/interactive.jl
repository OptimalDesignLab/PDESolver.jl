# these functions are more useful for interactive usage of the code

"""
  Print all initial condition names for a given physics module

  Inputs:
    mod: the module to print the ICs for
    f: the IO to print to, defaults to STDOUT

  Outputs:
    none
"""
function printICNames(mod::Module, f::IO=STDOUT)

  println(f, "IC names for module ", mod, ":")
  for name in keys(mod.ICDict)
    println(f, "  ", name)
  end

  return nothing
end

"""
  Like printICNames, but for boundary conditions, see that function for details
"""
function printBCNames(mod::Module, f::IO=STDOUT)

  println(f, "BC names for module ", mod, ":")
  for name in keys(mod.BCDict)
    println(f, "  ", name)
  end

  return nothing
end



