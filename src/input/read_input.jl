#include("new_file2.jl")  # creates arg_dict
#include("../tools/misc.jl")
"""
  This function reads an input file and turns it into a dictionary.
  Unlike [`read_input`](@ref), this function returns the contents of the
  input file exactly as a dictionary, without supplying default values.

  **Inputs**

   * fname: a string containing the file name.  If the file name does not
            start with a . or .. or /, then fname is treated as relative to
            the pwd, otherwise fname is interpreted as an absolute path

  **Outputs**

   * arg_dict: the dictionary containing all the options.

  Exceptions: if the input file does not contain a dictionary, an exception
              is thrown.

  Notes: currently, the input file format is a literal dictionary that is
         processed using eval().  This is likely to change to change in the
         future to support static compilation.
"""
function read_input_file(fname::AbstractString)

  if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    println("pwd = ", pwd())
    println("fname = ", fname)
  end

  fpath = joinpath(pwd(), fname)

  # this uses eval, which is evil (and not statically compilable)
  arg_dict = evalfile(fpath)  # include file in the users pwd()

  if !( typeof(arg_dict) <: Dict )
    throw(ErrorException("Input file does not contain a dictionary"))
  end

  return arg_dict
end


@doc """
### PDESolver.read_input

  This function read an input file, supplies default values whenever possible,
  and does sanity checks on the options.
  A dictionary (referred to as the options dictionary) is returned.
  See [`read_input_file`](@ref) for the description of the input file format.

  After default values are supplied, the dictionary is printed to
  arg_dict_output.jl (by MPI rank 0) in the format of an input file.
  This is useful for rerunning a simulation.

  This function prints warnings if keys in the dictionary do not correspond
  to known options.  The file known_keys.jl lists all known keys.
  See input_vals.txt for the list of keys, possible values, and their meanings.

  This function is idempotent; this is essential for using arg_dict_output.jl
  to rerun a simulation.

  Inputs:
    * fname : name of file to read, can be relative or absolute path.

  Outputs:
    arg_dict: a Dict{Any, Any} containing the option keywords and values

"""->
function read_input(fname::AbstractString)

  arg_dict = read_input_file(fname)
  arg_dict = read_input(arg_dict)

  return arg_dict
end

"""
  This method supplies default values to a given dictionary.  See the other
  method for details.

  **Inputs**

   * arg_dict: a dictionary

  **Outputs**

   * arg_dict: the same dictionary, with default values supplied
"""
function read_input(arg_dict::Dict)


#include(joinpath(pwd(), fname))  # include file in the users pwd()
#include(joinpath(Pkg.dir("PDESolver"), "src/Input/known_keys.jl"))  # include the dictonary of known keys
# take action based on the dictionary

#arg_dict = evalfile(fpath)  # include file in the users pwd()a

#arg_dict = read_input_file(fname)
# TODO: make this a global const
#known_keys = evalfile(joinpath(Pkg.dir("PDESolver"), "src/input/known_keys.jl"))  # include the dictonary of known keys


# new (201612) options checking function
checkForIllegalOptions_pre(arg_dict)

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
if arg_dict["use_DG"]
  get!(arg_dict, "operator_type", "SBPOmega")
else
  get!(arg_dict, "operator_type", "SBPGamma")
end

get!(arg_dict, "operator_type2", "SBPNone")
get!(arg_dict, "use_staggered_grid", arg_dict["operator_type2"] != "SBPNone")


Ma = get!(arg_dict, "Ma", -1.0)
Re = get!(arg_dict, "Re", -1.0)
aoa = get!(arg_dict, "aoa", -1.0)
arg_dict["aoa"] = aoa*pi/180  # convert to radians
#rho_free = get!(arg_dict, "rho_free", -1)
#E_free = get!(arg_dict, "E_free", -1)
get!(arg_dict, "vortex_x0", 0.0)
get!(arg_dict, "vortex_strength", 1.0)

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

# volume integral options
get!(arg_dict, "volume_integral_type", 1)
get!(arg_dict, "Volume_flux_name", "StandardFlux")
get!(arg_dict, "Viscous_flux_name", "ErrorFlux")
get!(arg_dict, "face_integral_type", 1)
get!(arg_dict, "FaceElementIntegral_name", "ESLFFaceIntegral")

# viscous specific options
get!(arg_dict, "isViscous", false)
get!(arg_dict, "SAT_type", "Hartman")

# timestepping options
get!(arg_dict, "t_max", 0.0)

if !haskey(arg_dict, "delta_t") && (arg_dict["run_type"] == 1 || arg_dict["run_type"] == 20 || arg_dict["run_type"] == 30)
  arg_dict["calc_dt"] = true
else
  arg_dict["calc_dt"] = false
end
# should this really have a default value?
get!(arg_dict, "CFL", 0.4)
get!(arg_dict, "use_itermax", haskey(arg_dict, "itermax"))

# stabilization options
get!(arg_dict, "use_edgestab", false)
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
  # the raised cosine filter has no parameters
end

# figure out coloring distances
if arg_dict["use_edgestab"] || arg_dict["use_DG"]
  get!(arg_dict, "coloring_distance", 2)
else
  get!(arg_dict, "coloring_distance", 0)
end


if arg_dict["use_edgestab_prec"]
  get!(arg_dict, "coloring_distance_prec", 2)
else
  get!(arg_dict, "coloring_distance_prec", 0)
end

# parallel options
if arg_dict["run_type"] == 1 || arg_dict["run_type"] == 30
  get!(arg_dict, "parallel_type", 1)
else
  get!(arg_dict, "parallel_type", 2)
end

if arg_dict["run_type"] == 1 || arg_dict["run_type"] == 30
  if arg_dict["face_integral_type"] == 2  # entropy stable
    get!(arg_dict, "parallel_data", "element")
  else
    get!(arg_dict, "parallel_data", "face")
  end
else
  get!(arg_dict, "parallel_data", "element")
end

#-----------------------------------------------
# physics module options
get!(arg_dict, "use_Minv", false)       # apply inverse mass matrix to residual calc in physics module. needed for CN

if arg_dict["use_Minv"] == false && arg_dict["run_type"] == 20
  println("INPUT: User did not specify use_Minv but selected run_type is CN. Setting use_Minv = true.")
  arg_dict["use_Minv"] = true
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
arg_dict["smb_name"] = update_path(smb_name)

if haskey(arg_dict, "dmg_name") && arg_dict["dmg_name"] != ".null"
  dmg_name = arg_dict["dmg_name"]
  arg_dict["dmg_name"] = update_path(dmg_name)
else  # no dmg_name specified
  arg_dict["dmg_name"] = ".null"
end

get!(arg_dict, "perturb_ic", false)
get!(arg_dict, "perturb_mag", 0.0)
get!(arg_dict, "write_finalsolution", false)
get!(arg_dict, "write_finalresidual", false)

# solver option
get!(arg_dict, "addVolumeIntegrals", true)
get!(arg_dict, "addBoundaryIntegrals", true)
get!(arg_dict, "addFaceIntegrals", true)
get!(arg_dict, "addStabilization", true)

# logging options

# entropy
get!(arg_dict, "write_entropy", false)
get!(arg_dict, "write_entropy_freq", 1)
get!(arg_dict, "write_entropy_fname", "entropy.dat")

# integral q
get!(arg_dict, "write_integralq", false)
get!(arg_dict, "write_integralq_fname", "integralq.dat")

# enstrophy
get!(arg_dict, "write_enstrophy", false)
get!(arg_dict, "write_enstrophy_fname", "enstrophy.dat")
get!(arg_dict, "write_enstrophy_freq", 1)

# kinetic energy
get!(arg_dict, "write_kinetic_energy", false)
get!(arg_dict, "write_kinetic_energy_fname", "kinetic_energy.dat")
get!(arg_dict, "write_kinetic_energy_freq", 1)

# kinetic energy dt
get!(arg_dict, "write_kinetic_energydt", false)
get!(arg_dict, "write_kinetic_energydt_fname", "kinetic_energydt.dat")
get!(arg_dict, "write_kinetic_energydt_freq", 1)


get!(arg_dict, "check_density", true)
get!(arg_dict, "check_pressure", true)


# DG Flux options
get!(arg_dict, "LFalpha", 0.0)


get!(arg_dict, "precompute_volume_flux", true)
if arg_dict["volume_integral_type"] != 2
  get!(arg_dict, "precompute_face_flux", true)
else
  get!(arg_dict, "precompute_face_flux", false)
  get!(arg_dict, "precompute_q_face", false)
end
get!(arg_dict, "precompute_boundary_flux", true)

if !arg_dict["use_DG"]
  get!(arg_dict, "precompute_q_face", true)
  get!(arg_dict, "precompute_q_bndry", true)
end

# if not already specified, set to true just in case
get!(arg_dict, "precompute_q_face", true)
get!(arg_dict, "precompute_q_bndry", true)

# solver debugging options
writeflux = get!(arg_dict, "writeflux", false)
writeboundary = get!(arg_dict, "writeboundary", false)
get!(arg_dict, "writeq", false)
get!(arg_dict, "writeqface", false)
get!(arg_dict, "write_fluxface", false)

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
get!(arg_dict, "write_interfaces", false)
get!(arg_dict, "write_boundaries", false)
get!(arg_dict, "write_sharedboundaries", false)
get!(arg_dict, "use_linear_metrics", false)

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
get!(arg_dict, "write_vorticity_vis", false)
get!(arg_dict, "exact_visualization", false)
get!(arg_dict, "write_res", false)
get!(arg_dict, "output_freq", 1)
get!(arg_dict, "recalc_prec_freq", 1)
get!(arg_dict, "jac_type", 2)
get!(arg_dict, "use_jac_precond", false)
get!(arg_dict, "res_abstol", 1e-6)
get!(arg_dict, "res_reltol", 1e-6)
get!(arg_dict, "res_reltol0", -1.0)
get!(arg_dict, "step_tol", 0.0)
get!(arg_dict, "print_eigs", false)
get!(arg_dict, "write_eigs", false)
get!(arg_dict, "write_eigdecomp", false)
get!(arg_dict, "newton_globalize_euler", false)
get!(arg_dict, "euler_tau", 1.0)
get!(arg_dict, "use_volume_preconditioner", false)

if arg_dict["run_type"] == 5  # steady newton
  get!(arg_dict, "newton_verbosity", 5)
else
  get!(arg_dict, "newton_verbosity", 4)
end
  # figure out Newtons method type
run_type = arg_dict["run_type"]

if haskey(arg_dict, "jac_method")
  if arg_dict["jac_method"] == 1
    get!(arg_dict, "epsilon", 1e-6)
  elseif arg_dict["jac_method"] == 2
    get!(arg_dict, "epsilon", 1e-20)
  end
end

# clean-sheet Newton's method (internal to CN) option - only for debugging
get!(arg_dict, "cleansheet_CN_newton", false)

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
get!(arg_dict, "write_timing", false)

myrank = MPI.Comm_rank(MPI.COMM_WORLD)

# Functional computational options
get!(arg_dict, "calc_functional", false)
get!(arg_dict, "objective_function", "none")
get!(arg_dict, "num_functionals", 0)
get!(arg_dict, "functional_error", false)
get!(arg_dict, "functional_error_outfname", "functional_error")
get!(arg_dict, "analytical_functional_val", 0.0)

# Adjoint computation options
get!(arg_dict, "need_adjoint", false)
get!(arg_dict, "write_adjoint", false)
get!(arg_dict, "write_adjoint_vis", false)

# Unsteady adjoint (CN) computation options --- EXPERIMENTAL, NONWORKING CODE
get!(arg_dict, "adjoint_revolve", false)
get!(arg_dict, "adjoint_straight", false)
get!(arg_dict, "uadj_global", false)
get!(arg_dict, "use_Minv_override_for_uadj", false)


# checkpointing/restart options
get!(arg_dict, "most_recent_checkpoint", -1)
get!(arg_dict, "most_recent_checkpoint_path", "")
get!(arg_dict, "writing_checkpoint", -1)
get!(arg_dict, "writing_checkpoint_path", "")
get!(arg_dict, "is_restart", false)
get!(arg_dict, "ncheckpoints", 2)
get!(arg_dict, "checkpoint_freq", 200)
get!(arg_dict, "use_checkpointing", false)

# Options passed directly to Petsc
petsc_opts = Dict{AbstractString, AbstractString}(
  "-malloc" => "",
  "-malloc_debug" => "",
  "-ksp_monitor" => "",
  "-pc_type" => "bjacobi",
  "-sub_pc_type" => "ilu",
  "-sub_pc_factor_levels" => "4",
  "-ksp_gmres_modifiedgramschmidt" => "",
  "-ksp_pc_side" => "right",
  "-ksp_gmres_restart" => "30",
)

if arg_dict["use_volume_preconditioner"]
  get!(petsc_opts, "-pc_type", "shell")
end


petsc_opts = get!(arg_dict, "petsc_options", petsc_opts)

# set the options in Petsc, if Petsc will be used
if arg_dict["jac_type"] == 3 || arg_dict["jac_type"] == 4
  PetscSetOptions(arg_dict["petsc_options"])
end

# Advection specific options
# TODO; move these into physics module
get!(arg_dict, "advection_velocity", [1.0, 1.0, 1.0])

checkForIllegalOptions_post(arg_dict)

# write complete dictionary to file
myrank = MPI.Comm_rank(MPI.COMM_WORLD)
commsize = MPI.Comm_size(MPI.COMM_WORLD)
if myrank == 0
  fname = "arg_dict_output"
  make_input(arg_dict, fname)
  #=
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
  =#
end
# do some sanity checks here


checkKeys(arg_dict, KNOWN_KEYS)

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
      println(STDERR, "Warning: Key ", key, " in input dictionary not ",
               "recognized")
       cnt += 1
    end
  end

  return cnt
end


function update_path(path)

  if startswith(path, "SRCMESHES")
    path = joinpath(Pkg.dir("PDESolver"), "src", "mesh_files", lstrip(path[10:end], '/'))
  end

  return path
end

"""
  Check the user supplied options for errors before supplying default values
"""
function checkForIllegalOptions_pre(arg_dict)

  # Ensure that jac-method is not specified
  if haskey(arg_dict, "jac_method")
    if arg_dict["run_type"] == 1
      warn("jac_method specified, but run_type is RK4.")
    end
  end

  if haskey(arg_dict, "jac_type")
    if arg_dict["jac_type"] == 3 || arg_dict["jac_type"] == 4
      if arg_dict["jac_method"] != 2
        warn("PETSc jac_type specified, but jac_method is not 2 (complex step)")
      end
    end
  end

  if get(arg_dict, "use_volume_preconditioner", false)
    petsc_opts = get(arg_dict, "petsc_options", Dict{Any, Any}())
    val = get(petsc_opts, "-pc_type", "shell")
    if val != "shell"
      error("when use_volume_preconditioner, the petsc_opts -pc_type must be either unspecified or \"shell\"")
    end
  end


  return nothing

end


"""
  Check the user supplied options for errors after supplying default options.
"""
function checkForIllegalOptions_post(arg_dict)

  myrank = MPI.Comm_rank(MPI.COMM_WORLD)
  commsize = MPI.Comm_size(MPI.COMM_WORLD)

  jac_type = arg_dict["jac_type"]
  if commsize > 1 && !( jac_type == 3 || jac_type == 4) && (arg_dict["run_type"] != 1 && arg_dict["run_type"] != 30)
  error("Invalid jacobian type for parallel run")
end


  if !arg_dict["use_DG"]
    keys = ["precompute_volume_flux", "precompute_face_flux", "precompute_boundary_flux"]
    for key in keys
      !arg_dict[key] && error("$key not supported for CG")
    end
  end

  if !arg_dict["precompute_volume_flux"] && !arg_dict["Q_transpose"]
    error("cannot combine non-transposed Q and not precomputing the volume flux")
  end

  # error if checkpointing not supported
  checkpointing_run_types = [1, 30, 20]
  if arg_dict["use_checkpointing"] && !(arg_dict["run_type"] in checkpointing_run_types)
    error("checkpointing only supported with RK4 and LSERK and CN")
  end

  if arg_dict["use_checkpointing"]
    if arg_dict["ncheckpoints"] <= 0
      error("checkpointing requires ncheckpoints > 0")
    end

    if arg_dict["checkpoint_freq"] <= 0
      error("checkpointing requires checkpoint_freq > 0")
    end
  end

  jac_type = arg_dict["jac_type"]
  if arg_dict["use_volume_preconditioner"] && (jac_type != 3 && jac_type != 4)
    error("cannot precondition non-iterative method")
  end

  if arg_dict["use_volume_preconditioner"] && arg_dict["run_type"] != 20
    error("cannot use volume preconditioner with any method except CN")
  end


  checkBCOptions(arg_dict)

  return nothing
end

"""
  Check that a model entity does not have more than one boundary condition
  applied to it
"""
function checkBCOptions(arg_dict)

  numBC = arg_dict["numBC"]

  # sort all the BC arrays
  for i=1:numBC
    key_i = string("BC", i)
    sort!(arg_dict[key_i])
  end

  # check for repeated edge numbers
  # this isn't cache friendly
  for i=1:numBC
  key_i = string("BC", i)
  vals = arg_dict[key_i]

    for j=1:length(vals)
      val_j = vals[j]

      # check this value against all previous one
      for k=1:(i-1)
        key_k = string("BC", k)
        index = findfirst(arg_dict[key_k], val_j)
        if index != 0
          throw(ErrorException("cannot apply more than one boundary condition to a model entity:\n Model entity $val_j from BC$i is repated in BC$k"))
        end
      end
    end
  end

  # check for repeated edge numbers within each array

  for i=1:numBC
    key_i = string("BC", i)
    vals = arg_dict[key_i]

    if vals != unique(vals)
      throw(ErrorException("cannot apply a boundary condition to a model entity more than once: BC$i has repeated model entities"))
    end
  end

  return nothing
end
