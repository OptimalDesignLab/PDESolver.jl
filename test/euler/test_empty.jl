"""
  A test of the test system itself
"""
function test_empty()
  @testset "Testing empty..." begin

    @testset "Testing if 1 == 1" begin
      @test ( 1 )== 1
    end

  end
end

#test_empty()
add_func1!(EulerTests, test_empty, [TAG_SHORTTEST])

