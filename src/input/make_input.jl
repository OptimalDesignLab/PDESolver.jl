# functions to create input files
using ODLCommonTools

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

end
