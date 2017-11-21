# functions to create input files
using ODLCommonTools

"""
  Make an input file from an options dictionary.

  Inputs:
    dict: the options dictionary
    fname: the file name (without extension)

  Outputs:
    fname_ex: file name with extension
"""
function make_input(dict::Dict, fname::AbstractString)

  fname_ex = string(fname, ".jl")
  rmfile(fname_ex)  

  f = open(fname_ex, "w")

  println(f, "arg_dict = Dict{Any,Any}(")

  for i in keys(dict)
    show(f, i)
    print(f, " => ")
    show(f, dict[i])
    println(f, ",")
  end

  println(f, ")")
  close(f)

  return fname_ex
end
