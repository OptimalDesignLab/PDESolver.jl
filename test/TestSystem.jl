#module TestSystem
# this module implements a mechanism for running tests, including 
# subsets of tests defined by tags
module TestSystem
using Input

global const TAG_DEFAULT = "tag_default"  # tag that every function must have
include("tags.jl")  # list of tags

global const EmptyDict = Dict{Any, Any}()

export TAG_DEFAULT, TestList, add_func1!, add_func2!, add_func3!, run_testlist,
       runTestSystem, not_isapprox

"""
  This type stores all the data that describes a set of tests and their
  associated tags.

  ** Fields**

   * funcs: list of functions that contain the tests.  The order in which
           the tests run is the same as the order of insertion

   * func_tags: Vector, same length as funcs, where each element is a vector
               of tags associated with each func.

   * func_type:  Vector, same length as funcs, describing how to run each
                function. 1 = no arguments, 2 = args = (mesh, sbp, eqn, opts),
                3 = (mesh, sbp, eqn, opts), using a modified version of an
                existing input file
   * input_name: Vector, same length as funcs, containing the name of the input
                file to use with each function.  For func_type == 1 this value
                is unused, for func type == 2 it is used directly, for 
                func_type == 3 it is modified and written to a new name
                according to mod_input

   * mod_input: Vector, same length as funcs, where each element is a dictionary.
               For func_type == 3 the input file specified by input_name is
               modified with these keys.  There must also be a key called
               "new_fname" that specifies name to write the modified file
               to (not including file extension.  
               This modified file will be loaded before the function 
               is called.

   * tag_list: collection of all known tags
"""
mutable struct TestList
  funcs::Vector{Function}  # functions containing the tests
  func_tags::Vector{Vector{String}}  # tags associated with each function
  func_type::Vector{Int}  # enum defining how to run the functino
  input_name::Vector{String} # the name of the input file for each function
  input_mod::Vector{Dict}  # dictionary of modifications to input_name
  tag_list::Set{String}  # all known tags

  function TestList()
    sizehint = 0  # initial length of vectors
    funcs = Array{Function}(sizehint)
    func_tags = Array{Vector{String}}(sizehint)
    func_type = Array{Int}(sizehint)
    input_name = Array{String}(sizehint)
    input_mod = Array{Dict}(sizehint)
    tag_list = Set{String}()
    push!(tag_list, TAG_DEFAULT)

    return new(funcs, func_tags, func_type, input_name, input_mod, tag_list)
  end  # end constructor
end  # end type definition

#------------------------------------------------------------------------------
# TestList functions

"""
  This function adds a new test function of func_type == 1 to the list
  (a function that takes no arguments).

  **Inputs**

   * test_list: the list of tests to append the function to
   * func: the function
   * tags: an array of tags associated with this function.  The user is free
          to modify this array afterwards, but not to mutate the strings within
          the array. (optional)

"""
function add_func1!(testlist::TestList, func::Function, tags::Array{String}=String[])

  check_tags(func, tags)

  # function
  push!(testlist.funcs, func)

  # tags
  tag_list_full = Array{String}(length(tags)+1)
  tag_list_full[1] = TAG_DEFAULT
  tag_list_full[2:end] = tags
  push!(testlist.func_tags, tag_list_full)

  # function type
  push!(testlist.func_type, 1)

  # input name
  push!(testlist.input_name,  "")

  # input_mod
  push!(testlist.input_mod, EmptyDict)

  for i=1:length(tags)
    push!(testlist.tag_list, tags[i])
  end

end


"""
  This function adds a new test function of `func_type == 2` to the list
  (a function that takes the arguments (mesh, sbp, eqn, opts), as
  obtained from the input file specified).

  **Inputs**

   * test_list: the list of tests to append the function to
   * func: the function
   * input_name: name of input file to be used with this function
   * tags: an array of tags associated with this function.  The user is free
          to modify this array afterwards, but not to mutate the strings within
          the array. (optional)

"""
function add_func2!(testlist::TestList, func::Function, input_name::String,
                    tags::Array{String}=String[])

  check_tags(func, tags)

  # function
  push!(testlist.funcs, func)

  # tags
  tag_list_full = Array{String}(length(tags)+1)
  tag_list_full[1] = TAG_DEFAULT
  tag_list_full[2:end] = tags
  push!(testlist.func_tags, tag_list_full)

  # function type
  push!(testlist.func_type, 2)

  # input name
  push!(testlist.input_name,  input_name)

  # input_mod
  push!(testlist.input_mod, EmptyDict)

  for i=1:length(tags)
    push!(testlist.tag_list, tags[i])
  end

end

"""
  This function adds a new test function of `func_type == 3` to the list
  (a function that takes the arguments `(mesh, sbp, eqn, opts)`, as
  obtained from modifying the input file specified).

  **Inputs**

   * test_list: the list of tests to append the function to
   * func: the function
   * input_name: name of input file to be modified for use with this function
   * mod_dict: dictionary of keys to be added (or replaced) in the input file
   * tags: an array of tags associated with this function.  The user is free
           to modify this array afterwards, but not to mutate the strings within
           the array. (optional)

"""
function add_func3!(testlist::TestList, func::Function, input_name::String,
                    mod_dict::Dict, tags::Array{String}=String[])

  check_tags(func, tags)

  # function
  push!(testlist.funcs, func)

  # tags
  tag_list_full = Array{String}(length(tags)+1)
  tag_list_full[1] = TAG_DEFAULT
  tag_list_full[2:end] = tags
  push!(testlist.func_tags, tag_list_full)

  # function type
  push!(testlist.func_type, 3)

  # input name
  push!(testlist.input_name,  input_name)

  # input_mod
  push!(testlist.input_mod, mod_dict)

  for i=1:length(tags)
    push!(testlist.tag_list, tags[i])
  end

end

"""
  This function checks that all necessary tags are present and throws
  an exception is they are not.

  Currently this only verifies that one of LengthTags is present

  **Inputs**

   * func: the function (only needed to produce error message)
   * tags: the tags for this function

  Outputs: none
"""
function check_tags(func::Function, tags::Array{String})

  n = length(intersect(tags, LengthTags))
  if n > 1
    throw(ErrorException("function $func has $n > 1 LengthTags"))
  elseif n < 1
    throw(ErrorException("function $func has $n < 1 LengthTags"))
  end

  return nothing
end

"""
  This function runs a test list.  Tests are run in the order they were
  loaded into the TestList object.  This implementation handles the case of
  several tests sets sharing an input file efficiencly, ie. if several test 
  functions use the same input file and are placed in the test list 
  consecutively, the input file will be loaded only once.

  A list of tags can be optionally supplied.  In this case, only the tests
  that have the specified tags will be run.  If no list is supplied, all tags
  are run.

  **Inputs**

   * testlist: a TestList loaded with functions
   * prep_func = function used with test functions of type 2 or 3.  It must
                 have signature prep_func(fname::String). fname is the
                 name of hte input file associated with the test function
   * tags: an array of tags (optional)

  **Outputs**

   * ntests: number of tests that were run.
"""
function run_testlist(testlist::TestList, prep_func::Function, tags::Vector{String}=String[TAG_DEFAULT])

  println("Running tests with tags matching = ", tags)
  tags_set = Set(tags)
  ntests = length(testlist.funcs)

  # need to declare these here so they have the right scope
  # these variables are type unstable anways, the default the type of the
  # default value doesn't matter
  mesh = 0
  sbp = 0
  eqn = 0
  opts = 0

  ftiming = open("timing.dat", "w")
  input_name_prev = ""  # holds the previous input name
  n_run_tests = 0
  for i=1:ntests
    func_i = testlist.funcs[i]
    func_tags_i = testlist.func_tags[i]
    functype_i = testlist.func_type[i]
    input_name_i = testlist.input_name[i]
    input_mod_i = testlist.input_mod[i]

    for j=1:length(func_tags_i)
      tag_j = func_tags_i[j]
      if tag_j in tags_set
        n_run_tests += 1
        println("running function ", func_i)
        println("function type = ", functype_i)
        # run this test
        t_test = @elapsed if functype_i == 1  # function with no arguments
          func_i()

        elseif functype_i == 2  # function with all 4 arguments
#          if input_name_i != input_name_prev
            println("about to run prep function ", prep_func, " on input ", input_name_i)
            mesh, sbp, eqn, opts = prep_func(input_name_i)
#          end
          func_i(mesh, sbp, eqn, opts)
          input_name_prev = input_name_i

        elseif functype_i == 3  # modify input before running
          new_fname = input_mod_i["new_fname"]*".jl"

          if input_name_prev != new_fname  # don't create a file if it was done already
            arg_dict = evalfile(joinpath(pwd(), input_name_i))

            for (key, val) in input_mod_i
              if key != "new_fname"
                arg_dict[key] = val
              end
            end

            make_input(arg_dict, new_fname)

            mesh, sbp, eqn, opts = prep_func(new_fname)
            input_name_prev = new_fname
          end  # end if file already created

          func_i(mesh, sbp, eqn, opts)
        else
          throw(ErrorException("unsupported function type"))
        end

        println(ftiming, func_i, ": ", t_test, " seconds")
        println("time for test function ", func_i, " = ", t_test)

        break  # don't run this test more than once even if it matches
                  # multiple tags

      end  # end if tag matches
    end  # end loop over tags

  end  # end loop over test functions

  close(ftiming)

  return n_run_tests
end

"""
  This function runs a test list, according to the tags supplied.  Unlike 
  [`run_tests`](@ref), this function supplies default behavior for the `tags`
  list.  If the list is empty, then all tests are run, otherwise only the
  specified tag are run.  See [`run_tests`](@ref) for argument descriptions.

  **Inputs**

   * testlist
   * prep_func
   * tags
"""
function runTestSystem(testlist::TestList, prep_func::Function, tags::Vector{String})

  nargs = length(tags)
  if nargs == 0
    tags2 = String[TAG_DEFAULT]
  else
    tags2 = Array{String}(nargs)
    copy!(tags2, tags)
  end

  n = run_testlist(testlist, prep_func, tags2)

  if n == 0
    error("No tests run matching tags: ", tags2)
  end

  return nothing
end


"""
  Helper function to compute !isapprox(args..., kwargs) in @test macros because

  @test !isapprox(args..., kwargs...)

  doesn't work.
"""
function not_isapprox(args...; kwargs...)
  return !isapprox(args...; kwargs...)
end



end  # end module
