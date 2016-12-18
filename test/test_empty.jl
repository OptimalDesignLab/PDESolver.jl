
function test_empty()
  facts("Testing empty...") do

    context("Testing if 1 == 1") do
      @fact 1 --> 1
    end

  end
end

#test_empty()
add_func1!(EulerTests, test_empty)

