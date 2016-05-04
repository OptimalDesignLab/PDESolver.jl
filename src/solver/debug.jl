# debug.jl
# this file contain debugging utilities

# global variable that determines debug level
# must be a literal value
global const DB_LEVEL = 1

@doc """
### PDESolver.debug1

  This macro returns the expression if the debug level is greater than or 
  equal to 1.

  Inputs/Outputs:
    expr1: an expression to be executed only if the debug level is greater
           than 1

"""->
macro debug1(expr1)
#  println("entered macro do_db")
  if DB_LEVEL >= 1
#    println("at compile time, in DB_Level >=", DB_LEVEL)
    return quote
#      println("runtime expression")
        $(esc(expr1))
    end
  else
    return nothing
  end

end  # end macro


@doc """
### PDESolver.debug2

  This macro returns the expression if the debug level is greater than or 
  equal to 2

  Inputs/Outputs:
    expr1: an expression to be executed only if the debug level is greater
           than 2

"""->
macro debug2(expr1)
#  println("entered macro do_db")
  if DB_LEVEL >= 2
#    println("at compile time, in DB_Level >=", DB_LEVEL)
    return quote
#      println("runtime expression")
        $(esc(expr1))
    end
  else
    return nothing
  end

end  # end macro



