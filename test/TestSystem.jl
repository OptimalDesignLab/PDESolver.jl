#module TestSystem
# this module implements a mechanism for running tests, including 
# subsets of tests defined by tags

include( joinpath(Pkg.dir("PDESolver"), "src/input/make_input.jl"))

global const TAG_DEFAULT = "default_tag"  # tag that every function must have
global const EmptyDict = Dict{Any, Any}()
"""
  This type stores all the data that describes a set of tests and their
  associated tags.

  Fields:
    funcs: list of functions that contain the tests.  The order in which
           the tests run is the same as the order of insertion

    func_tags: Vector, same length as funcs, where each element is a vector
               of tags associated with each func.

    func_type:  Vector, same length as funcs, describing how to run each
                function. 1 = no arguments, 2 = args = (mesh, sbp, eqn, opts),
                3 = (mesh, sbp, eqn, opts), using a modified version of an
                existing input file
    input_name: Vector, same length as funcs, containing the name of the input
                file to use with each function.  For func_type == 1 this value
                is unused, for func type == 2 it is used directly, for 
                func_type == 3 it is modified and written to a new name
                according to mod_input

    mod_input: Vector, same length as funcs, where each element is a dictionary.
               For func_type == 3 the input file specified by input_name is
               modified with these keys.  There must also be a key called
               "new_fname" that specifies name to write the modified file
               to (not including file extension.  
               This modified file will be loaded before the function 
               is called.

    tag_list: collection of all known tags
"""
type TestList
  funcs::Vector{Function}  # functions containing the tests
  func_tags::Vector{Vector{ASCIIString}}  # tags associated with each function
  func_type::Vector{Int}  # enum defining how to run the functino
  input_name::Vector{ASCIIString} # the name of the input file for each function
  input_mod::Vector{Dict}  # dictionary of modifications to input_name
  tag_list::Set{ASCIIString}  # all known tags

  function TestList()
    sizehint = 0  # initial length of vectors
    funcs = Array(Function, sizehint)
    func_tags = Array(Vector{ASCIIString}, sizehint)
    func_type = Array(Int, sizehint)
    input_name = Array(ASCIIString, sizehint)
    input_mod = Array(Dict, sizehint)
    tag_list = Set{ASCIIString}()
    push!(tag_list, TAG_DEFAULT)

    return new(funcs, func_tags, func_type, input_name, input_mod, tag_list)
  end  # end constructor
end  # end type definition

#------------------------------------------------------------------------------
# TestList functions

"""
  This function adds a new test function of func_type == 1 to the list
  (a function that takes no arguments).

  Inputs
    test_list: the list of tests to append the function to
    func: the function
    tags: an array of tags associated with this function.  The user is free
          to modify this array afterwards, but not to mutate the strings within
          the array. (optional)

"""
function add_func1!(testlist::TestList, func::Function, tags::Array{ASCIIString}=ASCIIString[])

  # function
  push!(testlist.funcs, func)

  # tags
  tag_list_full = Array(ASCIIString, length(tags)+1)
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
  This function adds a new test function of func_type == 2 to the list
  (a function that takes the arguments (mesh, sbp, eqn, opts), as
  obtained from the input file specified).

  Inputs
    test_list: the list of tests to append the function to
    func: the function
    input_name: name of input file to be used with this function
    tags: an array of tags associated with this function.  The user is free
          to modify this array afterwards, but not to mutate the strings within
          the array. (optional)

"""
function add_func2!(testlist::TestList, func::Function, input_name::ASCIIString,
                    tags::Array{ASCIIString}=ASCIIString[])

  # function
  push!(testlist.funcs, func)

  # tags
  tag_list_full = Array(ASCIIString, length(tags)+1)
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
  This function adds a new test function of func_type == 2 to the list
  (a function that takes the arguments (mesh, sbp, eqn, opts), as
  obtained from modifying the input file specified).

  Inputs
    test_list: the list of tests to append the function to
    func: the function
    input_name: name of input file to be modified for use with this function
    mod_dict: dictionary of keys to be added (or replaced) in the input file
    tags: an array of tags associated with this function.  The user is free
          to modify this array afterwards, but not to mutate the strings within
          the array. (optional)

"""
function add_func3!(testlist::TestList, func::Function, input_name::ASCIIString,
                    mod_dict::Dict, tags::Array{ASCIIString}=ASCIIString[])

  # function
  push!(testlist.funcs, func)

  # tags
  tag_list_full = Array(ASCIIString, length(tags)+1)
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
  This function runs a test list.  Tests are run in the order they were
  loaded into the TestList object.  This implementation handles the case of
  several tests sets sharing an input file efficiencly, ie. if several test 
  functions use the same input file and are placed in the test list 
  consecutively, the input file will be loaded only once.

  A list of tags can be optionally supplied.  In this case, only the tests
  that have the specified tags will be run.  If no list is supplied, all tags
  are run.

  Inputs:
    testlist: a TestList loaded with functions
    tags: an array of tags (optional)

"""
function run_testlist(testlist::TestList, tags::Vector{ASCIIString}=ASCIIString[TAG_DEFAULT])

  println("Running tests with tags matching = ", tags)
  tags_set = Set(tags)
  ntests = length(testlist.funcs)
  for i=1:ntests
    func_i = testlist.funcs[i]
    func_tags_i = testlist.func_tags[i]
    functype_i = testlist.func_type[i]
    input_name_i = testlist.input_name[i]
    input_mod_i = testlist.input_mod[i]

    for j=1:length(func_tags_i)
      tag_j = func_tags_i[j]
      if tag_j in tags_set
        println("running function ", func_i)
        # run this test
        if functype_i == 1  # function with no arguments
          func_i()
        elseif functype_i == 2  # function with all 4 arguments
          if ARGS[1] != input_name_i
            ARGS[1] = input_name_i
            include(STARTUP_PATH)
          end
          func_i(mesh, sbp, eqn, opts)
        elseif functype_i == 3  # modify input before running
          new_fname = input_mod_i["new_fname"]*".jl"

          if ARGS[1] != new_fname  # don't create a file if it was done already
            include(joinpath(pwd(), input_name_i))

            for (key, val) in input_mod_i
              if key != "new_fname"
                arg_dict[key] = val
              end
            end

            make_input(arg_dict, input_mod_i["new_fname"])

            ARGS[1] = new_fname
            include(STARTUP_PATH)
          end  # end if file already created

          func_i(mesh, sbp, eqn, opts)
        else
          throw(ErrorException("unsupported function type"))
        end

        continue  # don't run this test more than once even if it matches
                  # multiple tags

      end  # end if tag matches
    end  # end loop over tags

  end  # end loop over test functions

  return nothing
end


#end  # end module
