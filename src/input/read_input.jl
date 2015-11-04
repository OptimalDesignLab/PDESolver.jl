#include("new_file2.jl")  # creates arg_dict
#include("../tools/misc.jl")
@doc """
### read_input

  This function reads a file which must declare a dictionary of options.
  See the documention on input variables for valid keywords.  This function
  returns the dictionary after doing some sanity checks

  Arguments:
    * fname : name of file to read
"""->
function read_input(fname::AbstractString)

println("pwd = ", pwd())

include(joinpath(pwd(), fname))  # include file in the users pwd()

# take action based on the dictionary

if haskey(arg_dict, "var1")
  global DB_LEVEL =  arg_dict["var1"]
else
  global DB_LEVEL = 0
end

# record fname in dictionary
arg_dict["fname"] = fname

# supply default values if not given 
# get() = get!(dictionary, key, default)
gamma = get!(arg_dict, "gamma", 1.4)
gamma_1 = gamma - 1
R = get!(arg_dict, "R", 287.058)
cv = R/gamma_1

Ma = get!(arg_dict, "Ma", -1.0)
Re = get!(arg_dict, "Re", -1.0)
aoa = get!(arg_dict, "aoa", -1.0)
arg_dict["aoa"] = aoa*pi/180  # convert to radians
#rho_free = get!(arg_dict, "rho_free", -1)
#E_free = get!(arg_dict, "E_free", -1)
get!(arg_dict, "vortex_x0", 0.0)
get!(arg_dict, "vortex_strength", -1.0)



get!(arg_dict, "Relfunc_name", "none")


# stabilization options
get!(arg_dict, "use_edgestab", true)
get!(arg_dict, "edgestab_gamma", -0.1)
get!(arg_dict, "use_filter", false)
get!(arg_dict, "use_res_filter", false)
get!(arg_dict, "use_dissipation", false)
get!(arg_dict, "dissipation_name", "none")
get!(arg_dict, "dissipation_const", 0.0)

# preconditioning stabilization options
# non-logical values are shared between regular, preconditioned run
get!(arg_dict, "use_edgestab_prec", false) 
get!(arg_dict, "use_filter_prec", false)
get!(arg_dict, "use_dissipation_prec", false)

if arg_dict["use_filter"]
  get!(arg_dict, "filter_name", "raisedCosineFilter")
  # the raised cosine filter has no paramters
end

# figure out coloring distances
if arg_dict["use_edgestab"]
  get!(arg_dict, "coloring_distance", 2)
else
  get!(arg_dict, "coloring_distance", 0)
end


if arg_dict["use_edgestab_prec"]
  get!(arg_dict, "coloring_distance_prec", 2)
else
  get!(arg_dict, "coloring_distance_prec", 0)
end


# misc options
get!(arg_dict, "calc_error", false)
get!(arg_dict, "calc_error_infname", "none")
get!(arg_dict, "calc_error_outfname", "error_calc.dat")
get!(arg_dict, "calc_trunc_error", false)


# Algorithmic Differentiation options
get!(arg_dict, "use_edge_res", false)

# deal with file names
smb_name = arg_dict["smb_name"]
arg_dict["smb_name"] = joinpath(Pkg.dir("PDESolver"), smb_name)

if haskey(arg_dict, "dmg_name") && arg_dict["dmg_name"] != ".null"
  dmg_name = arg_dict["dmg_name"]
  arg_dict["dmg_name"] = joinpath(Pkg.dir("PDESolver"), dmg_name)
else  # no dmg_name specified
  arg_dict["dmg_name"] = ".null"
end

get!(arg_dict, "perturb_ic", false)
get!(arg_dict, "perturb_mag", 0.0)
get!(arg_dict, "write_finalsolution", false)
get!(arg_dict, "write_finalresidual", false)

# solver debugging options
writeflux = get!(arg_dict, "writeflux", false)
writeboundary = get!(arg_dict, "writeboundary", false)
get!(arg_dict, "writeq", false)

# mesh debugging options
get!(arg_dict, "write_edge_vertnums", false)
get!(arg_dict, "write_face_vertnums", false)
get!(arg_dict, "write_boundarynums", false)
get!(arg_dict, "write_dxidx", false)
get!(arg_dict, "write_coords", false)
get!(arg_dict, "write_sparsity", false)
get!(arg_dict, "verify_coloring", true)
get!(arg_dict, "write_counts", false)

# mesh options
get!(arg_dict, "reordering_algorithm", "default")
get!(arg_dict, "reordering_start_coords", [0.0, 0.0])

# Newton's Method options
get!(arg_dict, "write_rhs", false)
get!(arg_dict, "write_jac", false)
get!(arg_dict, "print_cond", false)
get!(arg_dict, "write_sol", false)
get!(arg_dict, "write_qic", false)
get!(arg_dict, "write_vis", false)
get!(arg_dict, "write_res", false)
get!(arg_dict, "output_freq", 1)
get!(arg_dict, "recalc_prec_freq", 1)
get!(arg_dict, "jac_type", 2)
get!(arg_dict, "res_abstol", 1e-6)
get!(arg_dict, "res_reltol", 1e-6)
get!(arg_dict, "res_reltol0", -1.0)
get!(arg_dict, "print_eigs", false)
get!(arg_dict, "write_eigs", false)
get!(arg_dict, "newton_globalize_euler", true)
get!(arg_dict, "euler_tau", 1)
  # figure out Newtons method type
run_type = arg_dict["run_type"]
if run_type == 4
  arg_dict["jac_method"] = 1  # finite difference
  get!(arg_dict, "epsilon", 1e-6)
elseif run_type == 5
  arg_dict["jac_method"] = 2
  get!(arg_dict, "epsilon", 1e-20)
end

get!(arg_dict, "real_time", false)

# testing options
get!(arg_dict, "solve", true)

# write complete dictionary to file
fname = "arg_dict_output.txt"
rmfile(fname)
f = open(fname, "a+")
arg_keys = keys(arg_dict)

for key_i in arg_keys
  println(f, key_i, " => ", arg_dict[key_i])
end
close(f)


# do some sanity checks here
# deal with boundary conditions
# "numBC" must be dictionary key whose value is the number of boundary conditions
# for each boundary condition there must be keys BCi and BCi_name for i=1:numBC
numBC = arg_dict["numBC"]

#=
# sort all the BC arrays
for i=1:numBC
  key_i = string("BC", i)
  println("edge nums before sort = ", arg_dict[key_i])
  sort!(arg_dict[key_i])
  println("edge nums after sort = ", arg_dict[key_i])
  enum_key_i = string("BC", i, "_name")
  println("arg_dict[enum_key_i] = ", arg_dict[enum_key_i])
end

# check for repeated edge numbers
# this isn't cache friendly
for i=1:numBC
key_i = string("BC", i)
vals = arg_dict[key_i]

  for j=1:length(vals)
    val_j = vals[j]
    println("val_j = ", val_j)

    # check this value against all previous one
    for k=1:(i-1)
      key_k = string("BC", k)
      println("key_k = ", key_k)
      println("arg_dict[key_k] = ", arg_dict[key_k])
      index = findfirst(arg_dict[key_k], val_j)
      if index != 0
	println("Error: cannot apply more than one boundary condition to a model entity")
	println("  Model entity ", val_j, " from BC", i, " is repeated at index ", index, " of BC", k)
      end
    end
  end
end

# check fo repeated edge numbers within each array

for i=1:numBC
 key_i = string("BC", i)
 vals = arg_dict[key_i]

 if vals != unique(vals)
   println("Error: cannot apply more than one boundary condition to a model entity")
   println("BC", i, " has a repeated value")
 end
end

=#

return arg_dict

end  # end function


macro do_db(expr1)
  println("entered macro do_db")
#  println("expr1 = ", expr1)
#  println("typeof(expr1) = ", typeof(expr1))
  if DB_LEVEL < 2
    println("at compile time, in DB_Level < 2")
    return quote
            println("runtime expression")
             expr1
    end
  else
    return nothing
  end
end



