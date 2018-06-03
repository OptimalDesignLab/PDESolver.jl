# test complexify.jl

"""
  This function tests whether compexify.jl functions work correctly
"""
function test_complexify()

  @testset "----- Testing Complexify.jl ------" begin

    a = 1.0
    b = -1.0
    @test isapprox( absvalue(a), abs(a)) atol=1e-15
    @test isapprox( absvalue(b), abs(b)) atol=1e-15

    c = complex(1.0, 1.0)
    d = complex(1.0, -1.0)
    f = complex(-1.0, 1.0)
    g = complex(-1.0, -1.0)

    @test isapprox( absvalue(c), complex(1.0, 1.0)) atol=1e-15
    @test isapprox( absvalue(d), complex(1.0, -1.0)) atol=1e-15
    @test isapprox( absvalue(f), complex(1.0, -1.0)) atol=1e-15
    @test isapprox( absvalue(g), complex(1.0, 1.0)) atol=1e-15

  end
end

#test_complexify()
add_func1!(EulerTests, test_complexify, [TAG_COMPLEX, TAG_SHORTTEST])
