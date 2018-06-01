function test_elliptic_conv_rate()
    facts("---Test convergence rates of Elliptic Solver---") do
        input_file = "./input_conv_study1.jl"
        mesh, sbp, eqn, opts = solvePDE(input_file)
        l2nrm0 = 0.07403342662560318
        func0  = 6.435906901469668e-5
        l2nrm = readdlm("l2norm.dat") 
        func  = readdlm("functional.dat") 
        @fact l2nrm[1] --> roughly(l2nrm0, atol=1e-10)
        @fact func[1]  --> roughly(func0,  atol=1e-10)

        input_file = "./input_conv_study2.jl"
        mesh, sbp, eqn, opts = solvePDE(input_file)
        l2nrm0 = 0.007932898520294955
        func0  = 1.6477731494431852e-5
        l2nrm = readdlm("l2norm.dat") 
        func  = readdlm("functional.dat") 
        @fact l2nrm[1] --> roughly(l2nrm0, atol=1e-10)
        @fact func[1]  --> roughly(func0,  atol=1e-10)

        input_file = "./input_conv_study3.jl"
        mesh, sbp, eqn, opts = solvePDE(input_file)
        l2nrm0 = 0.0009120826559896547
        func0  = 1.4640116780054169e-6
        l2nrm = readdlm("l2norm.dat") 
        func  = readdlm("functional.dat") 
        println(l2nrm)
        println(func)
        @fact l2nrm[1] --> roughly(l2nrm0, atol=1e-10)
        @fact func[1]  --> roughly(func0,  atol=1e-10)
    end
end

add_func1!(EllipticTests, test_elliptic_conv_rate, [TAG_SHORTTEST])
