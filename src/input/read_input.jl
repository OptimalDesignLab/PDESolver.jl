#include("new_file2.jl")  # creates arg_dict
#include("../tools/misc.jl")


@doc """
### PDESolver.read_input

  This function reads a file which must  be a julia source file that declares 
  a dictionary of option keywords and values for the options named arg_dict.
  See the documention on input variables for valid keywords.  
    
  read_input() returns the dictionary after doing some sanity checks and
  supplying default values for any unspecified keys.

  After supplying default values, it prints the dictonary to arg_dict_output.jl,
  which is a valid julia source file and can be read in to re-run a simulation.

  This function checks whether the keys in arg_dict are recognized keywords 
  and prints a warning to STDERR if an unrecognized key is found.  The list of
  known keys is read from the julia source file known_keys.jl

  Inputs:
    * fname : name of file to read

  Outputs:
    arg_dict: a Dict{Any, Any} containing the option keywords and values

"""->
function read_input(fname::AbstractString)

println("pwd = ", pwd())

#include(joinpath(pwd(), fname))  # include file in the users pwd()
#include(joinpath(Pkg.dir("PDESolver"), "src/Input/known_keys.jl"))  # include the dictonary of known keys
# take action based on the dictionary

include(joinpath(pwd(), fname))  # include file in the users pwd()
include(joinpath(Pkg.dir("PDESolver"), "src/input/known_keys.jl"))  # include the dictonary of known keys

# record fname in dictionary
arg_dict["fname"] = fname

# type of variables, defaults to conservative
get!(arg_dict, "variable_type", :conservative)

# supply default values if not given 
# get() = get!(dictionary, key, default)
gamma = get!(arg_dict, "gamma", 1.4)
gamma_1 = gamma - 1
R = get!(arg_dict, "R", 287.058)
cv = R/gamma_1

get!(arg_dict, "dimensions", 2)
get!(arg_dict, "use_DG", false)

Ma = get!(arg_dict, "Ma", -1.0)
Re = get!(arg_dict, "Re", -1.0)
aoa = get!(arg_dict, "aoa", -1.0)
arg_dict["aoa"] = aoa*pi/180  # convert to radians
#rho_free = get!(arg_dict, "rho_free", -1)
#E_free = get!(arg_dict, "E_free", -1)
get!(arg_dict, "vortex_x0", 0.0)
get!(arg_dict, "vortex_strength", 1.0)

# should this really have a default value?
get!(arg_dict, "CFL", 0.4)
get!(arg_dict, "mesh_size", 1.0)  # this key should not exist
get!(arg_dict, "Relfunc_name", "none")

#SBP Options
get!(arg_dict, "Q_transpose", true)

# source term options
get!(arg_dict, "SRCname", "SRC0")
if arg_dict["SRCname"] == "SRC0"
  get!(arg_dict, "use_src_term", false)
else
  get!(arg_dict, "use_src_term", true)
end


# stabilization options
get!(arg_dict, "use_edgestab", true)
get!(arg_dict, "edgestab_gamma", -0.1)
get!(arg_dict, "use_filter", false)
get!(arg_dict, "use_res_filter", false)
get!(arg_dict, "use_dissipation", false)
get!(arg_dict, "dissipation_name", "none")
get!(arg_dict, "dissipation_const", 0.0)
get!(arg_dict, "use_GLS", false)
get!(arg_dict, "use_GLS2", false)
get!(arg_dict, "tau_type", 1)

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
get!(arg_dict, "calc_havg", false)

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

# solver options
get!(arg_dict, "write_entropy", false)
get!(arg_dict, "write_entropy_fname", "entropy.dat")
get!(arg_dict, "check_density", true)
get!(arg_dict, "check_pressure", true)

# solver debugging options
writeflux = get!(arg_dict, "writeflux", false)
writeboundary = get!(arg_dict, "writeboundary", false)
get!(arg_dict, "writeq", false)

#DEBUGGING
get!(arg_dict, "test_GLS2", false)
# mesh debugging options
get!(arg_dict, "write_edge_vertnums", false)
get!(arg_dict, "write_face_vertnums", false)
get!(arg_dict, "write_boundarynums", false)
get!(arg_dict, "write_dxidx", false)
get!(arg_dict, "write_coords", false)
get!(arg_dict, "write_sparsity", false)
get!(arg_dict, "write_sparsity_nodebnds", false)
get!(arg_dict, "write_offsets", false)
get!(arg_dict, "write_dofs", false)
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
get!(arg_dict, "write_eigdecomp", false)
get!(arg_dict, "newton_globalize_euler", false)
get!(arg_dict, "euler_tau", 1.0)
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


# Krylov options
get!(arg_dict, "krylov_reltol", 1e-2)
get!(arg_dict, "krylov_abstol", 1e-12)
get!(arg_dict, "krylov_dtol", 1e5)
get!(arg_dict, "krylov_itermax", 1000)
get!(arg_dict, "krylov_gamma", 2)

# testing options
get!(arg_dict, "solve", true)


# postprocessing options
get!(arg_dict, "do_postproc", false)
get!(arg_dict, "exact_soln_func", "nothing")

# write complete dictionary to file
fname = "arg_dict_output.jl"
rmfile(fname)
f = open(fname, "a+")

println(f, "arg_dict = Dict{Any, Any}(")
arg_keys = keys(arg_dict)



for key_i in arg_keys
  show(f, key_i)
  print(f, " => ")
  show(f, arg_dict[key_i])
  println(f, ",")
#  println(f, show(key_i), " => ", show(arg_dict[key_i]), ",")
end
println(f, ")")
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

checkKeys(arg_dict, known_keys)

return arg_dict

end  # end function

@doc """
### PDESolver.checkKeys

  This function verifies all the keys in the first argument are also keys
  of the second argument and prints a warning to STDERR if they are not.

  Inputs
    arg_dict: first dictonary
    known_keys: second dictonary

  Outputs:
    cnt: number of unrecognized keys
"""->
function checkKeys(arg_dict, known_keys)

  cnt = 0
  for key in keys(arg_dict)
    if !haskey(known_keys, key)
      println(STDERR, "Warning: Key ", key, " in input dictonary not ",
               "recognized")
       cnt += 1
    end
  end

  return cnt
end




