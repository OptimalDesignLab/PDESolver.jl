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

    x = 1.0; y = 1.0
    checkdiff(x, y)
    checkderiv(x, y)


    x = 0.1; y = 1.0
    checkdiff(x, y)
    checkderiv(x, y)


    x = 0.01; y = 1.0
    checkdiff(x, y)
    checkderiv(x, y)


    x = 0.1; y = -1.0
    checkdiff(x, y)
    checkderiv(x, y)



    x = -0.001; y = -1.0
    checkdiff(x, y)
    checkderiv(x, y)


  end
end

"""
  Helper function for testing atan2 value
"""
function checkdiff(x::Real, y::Real)

  x2 = Complex(x, 0)
  y2 = Complex(y, 0)

  theta1 = atan2(y, x)
  theta2 = atan2(y2, x2)

  @test isapprox(theta1, theta2) atol=2*eps()

  return theta1 - theta2
end

"""
  Helper function for testing atan2 derivative (finite difference)
"""
function checkderiv(x::Real, y::Real)

  h = 1e-20
  hfd = 1e-6
  x2 = Complex128(x, h)
  
  d1 = imag(atan2(y, x2))/h
  d1fd = (atan2(y, x + hfd) - atan2(y, x))/hfd

  @test isapprox(d1, d1fd) atol=1e-6

  y2 = Complex128(y, h)
  d1 = imag(atan2(y2, x))/h
  d1fd = (atan2(y + hfd, x) - atan2(y, x))/hfd

  @test isapprox(d1, d1fd) atol=1e-6

  return nothing
end



#test_complexify()
add_func1!(EulerTests, test_complexify, [TAG_COMPLEX, TAG_SHORTTEST])
