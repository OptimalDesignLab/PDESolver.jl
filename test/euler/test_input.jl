# verify processing of input files works
"""
  This function tests whether the known keys list works correctly
"""
function test_input()
  facts("---- Testing input processing ------") do

    extract_path = joinpath(Pkg.dir("PDESolver"), "src/input/extract_keys.jl")
    extract_path2 = joinpath(Pkg.dir("PDESolver"), "src/input/read_input.jl")

    include(extract_path)
    include(extract_path2)
    include("known_keys.jl")
    @fact haskey(known_keys, "key1") --> true
    @fact haskey(known_keys, "smb_name") --> true
    @fact haskey(known_keys, "var1") --> false

    #TODO: make this statically compilable
    include("input_test.jl")
    @fact checkKeys(arg_dict, known_keys) --> 1

  end
end

#test_input()
add_func1!(EulerTests, test_input, [TAG_INPUT])
