# verify processing of input files works
"""
  This function tests whether the known keys list works correctly
"""
function test_input()
  facts("---- Testing input processing ------") do

    extract_path = joinpath(Pkg.dir("PDESolver"), "src/input/extract_keys.jl")

    include(extract_path)
#    include("known_keys.jl")
    @fact haskey(Input.KNOWN_KEYS, "key1") --> false
    @fact haskey(Input.KNOWN_KEYS, "smb_name") --> true
    @fact haskey(Input.KNOWN_KEYS, "var1") --> false

    arg_dict = Input.read_input_file("input_test.jl")
    @fact Input.checkKeys(arg_dict, Input.KNOWN_KEYS) --> 1

    arg_dict["numBC"] = 1
    arg_dict["BC1"] = [0, 0]

    @fact_throws Inputs.checkBCOptions(arg_dict)

    # pick another random input file
    opts = Input.read_input("input_vals_channel.jl")

    # verify that supplying default values is idempotent
    opts_orig = copy(opts)
    Input.read_input(opts)

    for (key, val) in opts
      @fact val --> opts_orig[key]
    end
    


  end
end

#test_input()
add_func1!(EulerTests, test_input, [TAG_INPUT, TAG_SHORTTEST])
